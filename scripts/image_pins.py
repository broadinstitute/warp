#!/usr/bin/env python3
"""Release and pin warp-tools images in WARP's WDLs.

A WARP pin names a release tag and the digest it points at, so the digest fixes what
runs and the tag says what it is:

    us.gcr.io/broad-gotc-prod/<image>:<release tag>@sha256:<digest>

Commands (the script finds the WARP repository from its own location):

  release IMAGE SOURCE VERSION [PATH ...]
      Tag the build at SOURCE (the warp-tools branch name, any existing tag, or a
      sha256 digest) as VERSION-<unix timestamp>, then pin every WDL reference to it.
      VERSION is <image-version>[-<tool-version>], as in the image's docker_build.sh.
      If that build already carries a VERSION-<timestamp> tag, it is reused.
      Adding a tag needs registry write access (gcloud auth).

  pin IMAGE TAG [PATH ...]
      Pin every WDL reference to an existing tag, without tagging anything.

  Both rewrite every reference to IMAGE, or only those under the given PATHs
  (files or directories, e.g. pipelines/wdl/scanvi).

  check [--base REF | --all]
      Report pins that are not <image>:<tag>@sha256:<digest> with the tag resolving to
      that digest. By default only lines added since REF (default origin/develop).

Procedure: website/docs/contribution/contribute_to_warp/update_warp_tools_image.md
Standard library only.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request

REGISTRY = "us.gcr.io/broad-gotc-prod"
ROOTS = ("pipelines", "tasks", "verification", "all_of_us", "beta-pipelines", "projects")
# Shared by ~27 CI path filters: moving its pin starts all of those test suites.
SHARED = "verification/VerifyTasks.wdl"
ACCEPT = ", ".join((
    "application/vnd.oci.image.index.v1+json",
    "application/vnd.docker.distribution.manifest.list.v2+json",
    "application/vnd.oci.image.manifest.v1+json",
    "application/vnd.docker.distribution.manifest.v2+json",
))

NAME = r"[a-z0-9]+(?:[._-][a-z0-9]+)*"
TAG = r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}"
DIGEST = r"sha256:[0-9a-f]{64}"
REF_TAIL = r"(?::(?P<tag>%s))?(?:@(?P<digest>%s))?\s*" % (TAG, DIGEST)
VERSION_RE = re.compile(r"[0-9]+\.[0-9]+\.[0-9]+(?:-[A-Za-z0-9_.]+)*")
LITERAL_RE = re.compile(r"\"([^\"\n]*)\"|'([^'\n]*)'")


class Error(Exception):
    pass


# --- registry (anonymous reads; the registry is public) ---------------------------

def _get(url, accept=ACCEPT):
    request = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "warp-image-pins"})
    return urllib.request.urlopen(request, timeout=30)


def _api(registry, path):
    host, _, project = registry.partition("/")
    return "https://%s/v2/%s/%s" % (host, project, path)


def digest_of(registry, image, ref):
    """Digest that <registry>/<image>:<ref> (or @<ref>) points at; None if it doesn't exist."""
    try:
        with _get(_api(registry, "%s/manifests/%s" % (image, ref))) as response:
            return response.headers["Docker-Content-Digest"]
    except urllib.error.HTTPError as err:
        if err.code == 404:
            return None
        raise Error("registry answered HTTP %d for %s/%s:%s" % (err.code, registry, image, ref))
    except urllib.error.URLError as err:
        raise Error("cannot reach the registry: %s" % err.reason)


def tags_by_digest(registry, image):
    """{digest: [tags]} for an image, from the registry's tag list (empty if no such image)."""
    try:
        with _get(_api(registry, "%s/tags/list" % image), "application/json") as response:
            listing = json.load(response)
    except urllib.error.HTTPError as err:
        if err.code == 404:
            return {}
        raise Error("registry answered HTTP %d listing %s/%s" % (err.code, registry, image))
    except urllib.error.URLError as err:
        raise Error("cannot reach the registry: %s" % err.reason)
    return {d: m.get("tag", []) for d, m in listing.get("manifest", {}).items()}


def add_tag(registry, image, digest, tag):
    """Add <tag> to <digest> with the caller's gcloud credentials."""
    command = ["gcloud", "container", "images", "add-tag", "--quiet",
               "%s/%s@%s" % (registry, image, digest), "%s/%s:%s" % (registry, image, tag)]
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError:
        raise Error("gcloud is not installed; tagging needs it (and registry write access)")
    except subprocess.CalledProcessError:
        raise Error("gcloud could not add the tag; check 'gcloud auth list' and your registry write access")


# --- finding references in WDLs ---------------------------------------------------

def wdl_files(root):
    for top in ROOTS:
        for directory, _, names in sorted(os.walk(os.path.join(root, top))):
            for name in sorted(names):
                if name.endswith(".wdl"):
                    yield os.path.relpath(os.path.join(directory, name), root)


def references(text, registry, image=None):
    """Yield (line, start, end, name, tag, digest, is_full) for each image reference.

    A reference is a whole string literal: "<registry>/<name>[:tag][@digest]", or a bare
    "<name>[:tag][@digest]" in a file that holds the literal "<registry>/" (pipelines that
    prepend a docker prefix). Comment lines are skipped. With image set, only that image.
    """
    name = re.escape(image) if image else NAME
    full = re.compile(r"\s*%s/(?P<name>%s)%s" % (re.escape(registry), name, REF_TAIL))
    bare = re.compile(r"\s*(?P<name>%s)%s" % (name, REF_TAIL)) if '"%s/"' % registry in text else None
    for number, line in enumerate(text.splitlines(), 1):
        if line.lstrip().startswith("#"):
            continue
        for literal in LITERAL_RE.finditer(line):
            group = 1 if literal.group(1) is not None else 2
            content = literal.group(group)
            match = full.fullmatch(content) or (bare.fullmatch(content) if bare else None)
            if match:
                yield (number, literal.start(group), literal.end(group), match.group("name"),
                       match.group("tag"), match.group("digest"), match.re is full)


# --- pin --------------------------------------------------------------------------

def pin(root, registry, image, tag, digest, include_shared=False, dry_run=False, only=()):
    """Point references to <image> at <tag>@<digest>, all or only under <only>.

    Returns (changed, skipped) lines."""
    changed, skipped = [], []
    prefixes = [os.path.normpath(p) for p in only]
    for path in wdl_files(root):
        if prefixes and not any(path == p or path.startswith(p + os.sep) for p in prefixes):
            continue
        with open(os.path.join(root, path), encoding="utf-8") as handle:
            text = handle.read()
        refs = list(references(text, registry, image))
        if not refs:
            continue
        if path == SHARED and not include_shared:
            skipped += ["%s:%d" % (path, r[0]) for r in refs]
            continue
        lines = text.splitlines(keepends=True)
        for number, start, end, name, _, _, is_full in reversed(refs):
            new = "%s%s:%s@%s" % (registry + "/" if is_full else "", name, tag, digest)
            line = lines[number - 1]
            if line[start:end] != new:
                lines[number - 1] = line[:start] + new + line[end:]
                changed.append("%s:%d: %s -> %s" % (path, number, line[start:end].strip(), new))
        if not dry_run:
            with open(os.path.join(root, path), "w", encoding="utf-8") as handle:
                handle.write("".join(lines))
    return sorted(changed), skipped


def report_pin(registry, image, tag, digest, changed, skipped, dry_run):
    print("%s %s/%s:%s@%s" % ("Would pin" if dry_run else "Pinned", registry, image, tag, digest))
    for line in changed:
        print("  " + line)
    if not changed:
        print("  (no reference needed a change)")
    for line in skipped:
        print("  skipped %s (shared verification pin; --include-shared-verification moves it)" % line)
    if changed and not dry_run:
        print("Next: bump every affected pipeline (minor by default) and its changelog, commit, then run "
              "scripts/validate_release.sh -g origin/staging.")


# --- check ------------------------------------------------------------------------

def added_lines(root, base):
    """{path: {line numbers}} for WDL lines added between <base> and HEAD."""
    diff = subprocess.run(["git", "-C", root, "diff", "-U0", "%s...HEAD" % base, "--", "*.wdl"],
                          capture_output=True, text=True)
    if diff.returncode:
        raise Error("git diff against %s failed: %s" % (base, diff.stderr.strip()))
    added, path, number = {}, None, 0
    for line in diff.stdout.splitlines():
        if line.startswith("+++ "):
            path = line[6:] if line.startswith("+++ b/") else None
        elif line.startswith("@@"):
            number = int(re.search(r"\+(\d+)", line).group(1))
        elif line.startswith("+") and path:
            added.setdefault(path, set()).add(number)
            number += 1
    return added


def check(root, registry, base=None):
    """Return (path, line, message) for each checked pin that breaks the rule."""
    wanted = None if base is None else added_lines(root, base)
    ours, problems = {}, []
    for path in (sorted(wanted) if wanted is not None else wdl_files(root)):
        if not os.path.exists(os.path.join(root, path)):
            continue
        with open(os.path.join(root, path), encoding="utf-8") as handle:
            text = handle.read()
        for number, _, _, name, tag, digest, is_full in references(text, registry):
            if wanted is not None and number not in wanted[path]:
                continue
            if not is_full and not digest:
                # A bare "name:tag" may come from another registry; count it only if the image is ours.
                if name not in ours:
                    ours[name] = bool(tags_by_digest(registry, name))
                if not ours[name]:
                    continue
            if not tag:
                problems.append((path, number, "%s has a digest but no release tag" % name))
            elif not digest:
                problems.append((path, number, "%s:%s has no digest; pin <tag>@sha256:<digest>" % (name, tag)))
            else:
                actual = digest_of(registry, name, tag)
                if actual is None:
                    problems.append((path, number, "tag %s:%s not found in %s" % (name, tag, registry)))
                elif actual != digest:
                    problems.append((path, number, "tag %s:%s points at %s, not the pinned digest"
                                     % (name, tag, actual)))
    return problems


# --- command line -----------------------------------------------------------------

def repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def release_tag(registry, image, source, version, dry_run):
    """Tag the build at <source> as <version>-<timestamp> (or reuse one). Returns (tag, digest)."""
    if not VERSION_RE.fullmatch(version):
        raise Error("VERSION must be <image-version>[-<tool-version>], e.g. 1.7.3-1.5.1")
    digest = source if re.fullmatch(DIGEST, source) else digest_of(registry, image, source)
    if digest is None:
        raise Error("%s/%s:%s not found; is the warp-tools build finished, and is %s the branch name?"
                    % (registry, image, source, source))
    same_version = re.compile(r"%s-[0-9]{10}" % re.escape(version))
    reuse = sorted(t for t in tags_by_digest(registry, image).get(digest, []) if same_version.fullmatch(t))
    if reuse:
        print("Reusing release tag %s, already on %s" % (reuse[0], digest))
        return reuse[0], digest
    tag = "%s-%d" % (version, int(time.time()))
    if dry_run:
        print("Would tag %s/%s@%s as %s" % (registry, image, digest, tag))
        return tag, digest
    add_tag(registry, image, digest, tag)
    if digest_of(registry, image, tag) != digest:
        raise Error("tag %s did not resolve to %s after tagging" % (tag, digest))
    print("Tagged %s/%s@%s as %s" % (registry, image, digest, tag))
    return tag, digest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--registry", default=REGISTRY, help="default: %(default)s")
    commands = parser.add_subparsers(dest="command", required=True)
    release = commands.add_parser("release", help="tag a build and pin it")
    release.add_argument("image", help="registry repository, e.g. warp-tools")
    release.add_argument("source", help="warp-tools branch name, existing tag, or sha256 digest")
    release.add_argument("version", help="<image-version>[-<tool-version>], e.g. 1.7.3-1.5.1")
    pin_cmd = commands.add_parser("pin", help="pin to an existing tag")
    pin_cmd.add_argument("image")
    pin_cmd.add_argument("tag")
    for sub in (release, pin_cmd):
        sub.add_argument("paths", nargs="*", metavar="PATH", help="only rewrite references under these paths")
        sub.add_argument("--dry-run", action="store_true", help="show the changes; write and tag nothing")
        sub.add_argument("--include-shared-verification", action="store_true", help="also move pins in " + SHARED)
    check_cmd = commands.add_parser("check", help="report pins that break the rule")
    scope = check_cmd.add_mutually_exclusive_group()
    scope.add_argument("--base", default="origin/develop", help="check lines added since REF (default %(default)s)")
    scope.add_argument("--all", action="store_true", help="check every WDL")
    args = parser.parse_args(argv)
    root, registry = repo_root(), args.registry.rstrip("/")

    try:
        if args.command == "check":
            problems = check(root, registry, None if args.all else args.base)
            for path, number, message in problems:
                print("%s:%d: %s" % (path, number, message))
                if os.environ.get("GITHUB_ACTIONS") == "true":
                    print("::error file=%s,line=%d::%s" % (path, number, message))
            print("%d pin problem(s)" % len(problems) if problems else "All checked pins are tag@digest and resolve.")
            return 1 if problems else 0
        if args.command == "pin":
            tag, digest = args.tag, digest_of(registry, args.image, args.tag)
            if digest is None:
                raise Error("%s/%s:%s not found" % (registry, args.image, args.tag))
        else:
            tag, digest = release_tag(registry, args.image, args.source, args.version, args.dry_run)
        changed, skipped = pin(root, registry, args.image, tag, digest,
                               args.include_shared_verification, args.dry_run, args.paths)
        report_pin(registry, args.image, tag, digest, changed, skipped, args.dry_run)
        return 0
    except Error as err:
        print("error: %s" % err, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
