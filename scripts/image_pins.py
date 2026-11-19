#!/usr/bin/env python3
"""Tag and pin warp-tools images in WARP's WDLs.

A WARP pin names a version tag and the digest it points at, so the digest fixes what
runs and the tag says which image version it is:

    us.gcr.io/broad-gotc-prod/<image>:<version tag>@sha256:<digest>

Commands (the script finds the WARP repository from its own location):

  pin PR [PATH ...]
      PR is a warp-tools pull request: its number, URL or branch name. For each image
      that the PR's latest commit built, give that build the version tag
      <version>-<unix timestamp> and pin the WDLs to it. <version> is the tag the image's
      docker_build.sh builds, less the timestamp (e.g. 1.7.2-1.5.1). First it checks that
      the build finished, that the branch tag is that build, and that DOCKER_IMAGE_VERSION
      was bumped. A build that already has a <version>-<timestamp> tag keeps it. Tagging
      needs registry write access (gcloud auth).
        --digest D   pin this earlier build of the image instead (the one you tested)
        --version V  for the images whose docker_build.sh computes part of the tag at build
                     time (glimpse2, sv-imputation-rust-tools, zcall, illumina-iaap-autocall)

  pin IMAGE:TAG [PATH ...]
      Pin every WDL reference to an existing tag, without tagging anything.

  Both rewrite every reference to the image, or only those under the given PATHs
  (files or directories, e.g. pipelines/wdl/scanvi).

  check [--base REF | --all]
      Report pins that are not <image>:<tag>@sha256:<digest> with the tag resolving to
      that digest. By default only lines added since REF (default origin/develop).

Procedure: website/docs/contribution/contribute_to_warp/update_warp_tools_image.md
Standard library only.
"""

import argparse
import calendar
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

REGISTRY = "us.gcr.io/broad-gotc-prod"
WARP_TOOLS = "broadinstitute/warp-tools"
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


def manifests(registry, image):
    """{digest: {"tag": [...], "timeUploadedMs": ...}} from the registry's tag list ({} if no such image)."""
    try:
        with _get(_api(registry, "%s/tags/list" % image), "application/json") as response:
            return json.load(response).get("manifest", {})
    except urllib.error.HTTPError as err:
        if err.code == 404:
            return {}
        raise Error("registry answered HTTP %d listing %s/%s" % (err.code, registry, image))
    except urllib.error.URLError as err:
        raise Error("cannot reach the registry: %s" % err.reason)


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


# --- warp-tools pull requests (public GitHub reads) ---------------------------------

def _github(url, missing_ok=False):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "warp-image-pins"}),
                                    timeout=30) as response:
            return response.read().decode("utf-8")
    except urllib.error.HTTPError as err:
        if err.code == 404 and missing_ok:
            return None
        raise Error("GitHub answered HTTP %d for %s" % (err.code, url))
    except urllib.error.URLError as err:
        raise Error("cannot reach GitHub: %s" % err.reason)


def api(path):
    return json.loads(_github("https://api.github.com/repos/%s/%s" % (WARP_TOOLS, path)))


def raw(sha, path, missing_ok=False):
    return _github("https://raw.githubusercontent.com/%s/%s/%s" % (WARP_TOOLS, sha, path), missing_ok)


def find_pr(ref):
    """A warp-tools PR from its number, URL or branch name."""
    number = re.fullmatch(r"(?:https://github\.com/%s/pull/)?#?([0-9]+)(?:/.*)?" % WARP_TOOLS, ref)
    if number:
        return api("pulls/" + number.group(1))
    prs = api("pulls?state=all&head=%s:%s" % (WARP_TOOLS.split("/")[0], urllib.parse.quote(ref, safe="")))
    if not prs:
        raise Error("no warp-tools pull request from a branch named %s" % ref)
    return prs[0]


def builds(pr, registry):
    """[(image, directory, run)] for each image in <registry> that the PR's latest commit built."""
    sha, found, seen = pr["head"]["sha"], [], set()
    for run in api("actions/runs?event=pull_request&per_page=100&head_sha=" + sha)["workflow_runs"]:
        if run["path"] in seen or not re.fullmatch(r"\.github/workflows/build-[^/]+\.yml", run["path"]):
            continue
        seen.add(run["path"])  # runs come newest first
        workflow = raw(sha, run["path"])
        host, path, directory = (re.search(r"^\s*%s:\s*(\S+)" % key, workflow, re.M)
                                 for key in ("DOCKER_REGISTRY", "GCR_PATH", "working-directory"))
        if not (host and path and directory):
            continue  # pushes nothing (e.g. build-warp-tools-inverse.yml)
        repository = "%s/%s" % (host.group(1), path.group(1))
        if not repository.startswith(registry + "/"):
            raise Error("%s pushes to %s, outside %s" % (run["path"], repository, registry))
        image = repository[len(registry) + 1:]
        if run["status"] != "completed":
            raise Error("the %s build of %s is still running; rerun when it finishes" % (image, sha[:9]))
        if run["conclusion"] != "success":
            raise Error("the %s build of %s ended in %s; fix it in the warp-tools PR" % (image, sha[:9], run["conclusion"]))
        found.append((image, directory.group(1), run))
    if not found:
        raise Error("warp-tools PR #%d's latest commit (%s) built no image" % (pr["number"], sha[:9]))
    return found


def assignments(script):
    """Top-level NAME=value assignments with literal values in a shell script."""
    values = {}
    for name, value in re.findall(r'^([A-Z0-9_]+)="?([^"$\s]*)"?[ \t]*(?:#.*)?$', script, re.M):
        values.setdefault(name, value)
    return values


def image_version(script, where):
    """<image-version>[-<tool-version>]: the tag docker_build.sh builds, less its timestamp."""
    values = assignments(script)
    template = re.search(r'^\s*IMAGE_TAG="([^"]*TIMESTAMP[^"]*)"', script, re.M)
    if template:
        tag = re.sub(r"\$\{?([A-Z0-9_]+)\}?", lambda m: values.get(m.group(1), m.group(0)), template.group(1))
        tag = re.sub(r"-\$\{?TIMESTAMP\}?$", "", tag.rpartition(":")[2])
        if VERSION_RE.fullmatch(tag):
            return tag
    raise Error("cannot work out the version from %s; pass --version" % where)


def pushed_by(run, uploaded_ms):
    """Whether a manifest uploaded at <uploaded_ms> was pushed during <run>."""
    start, end = (calendar.timegm(time.strptime(run[key], "%Y-%m-%dT%H:%M:%SZ"))
                  for key in ("created_at", "updated_at"))
    return start <= int(uploaded_ms) / 1000 <= end + 60


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
                    ours[name] = bool(manifests(registry, name))
                if not ours[name]:
                    continue
            if not tag:
                problems.append((path, number, "%s has a digest but no tag" % name))
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


def version_tag(registry, image, digest, version, dry_run):
    """Tag <digest> as <version>-<timestamp>, or reuse such a tag already on it. Returns the tag."""
    if not VERSION_RE.fullmatch(version):
        raise Error("the version must be <image-version>[-<tool-version>], e.g. 1.7.3-1.5.1, not %s" % version)
    same_version = re.compile(r"%s-[0-9]{10}" % re.escape(version))
    reuse = sorted(t for t in manifests(registry, image).get(digest, {}).get("tag", []) if same_version.fullmatch(t))
    if reuse:
        print("Reusing version tag %s, already on %s/%s@%s" % (reuse[0], registry, image, digest))
        return reuse[0]
    tag = "%s-%d" % (version, int(time.time()))
    if dry_run:
        print("Would tag %s/%s@%s as %s" % (registry, image, digest, tag))
        return tag
    add_tag(registry, image, digest, tag)
    if digest_of(registry, image, tag) != digest:
        raise Error("tag %s did not resolve to %s after tagging" % (tag, digest))
    print("Tagged %s/%s@%s as %s" % (registry, image, digest, tag))
    return tag


def pin_pr(root, registry, pr_ref, only=(), digest=None, version=None, include_shared=False, dry_run=False):
    """Tag the builds of a warp-tools PR and pin the WDLs to them."""
    pr = find_pr(pr_ref)
    found = builds(pr, registry)
    if (digest or version) and len(found) > 1:
        raise Error("--digest and --version need a PR that builds one image; #%d built %s"
                    % (pr["number"], ", ".join(f[0] for f in found)))
    if digest and not re.fullmatch(DIGEST, digest):
        raise Error("--digest must be sha256:<64 hex digits>")
    plan = []
    for image, directory, run in found:
        script = raw(pr["head"]["sha"], directory + "/docker_build.sh")
        base = raw(pr["base"]["sha"], directory + "/docker_build.sh", missing_ok=True)
        bumped = assignments(script).get("DOCKER_IMAGE_VERSION")
        if base is not None and assignments(base).get("DOCKER_IMAGE_VERSION") == bumped:
            raise Error("%s/docker_build.sh still has DOCKER_IMAGE_VERSION=%s, as on %s; bump it in the PR"
                        % (directory, bumped, pr["base"]["ref"]))
        listing = manifests(registry, image)
        if digest:
            build = digest
            if build not in listing:
                raise Error("%s/%s@%s not found" % (registry, image, build))
        else:
            build = digest_of(registry, image, pr["head"]["ref"])
            if build is None or not pushed_by(run, listing.get(build, {}).get("timeUploadedMs", 0)):
                raise Error("%s/%s:%s is not the build of %s yet; rerun when its push finishes, or pass --digest"
                            % (registry, image, pr["head"]["ref"], pr["head"]["sha"][:9]))
        plan.append((image, build, version or image_version(script, directory + "/docker_build.sh")))
    for image, build, image_version_ in plan:
        tag = version_tag(registry, image, build, image_version_, dry_run)
        changed, skipped = pin(root, registry, image, tag, build, include_shared, dry_run, only)
        report_pin(registry, image, tag, build, changed, skipped, dry_run)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--registry", default=REGISTRY, help="default: %(default)s")
    commands = parser.add_subparsers(dest="command", required=True)
    pin_cmd = commands.add_parser("pin", help="tag a warp-tools PR's build and pin it, or pin an existing tag")
    pin_cmd.add_argument("target", metavar="PR|IMAGE:TAG",
                         help="a warp-tools PR's number, URL or branch name; or IMAGE:TAG, an existing tag")
    pin_cmd.add_argument("paths", nargs="*", metavar="PATH", help="only rewrite references under these paths")
    pin_cmd.add_argument("--digest", help="with a PR: pin this earlier build of the image instead")
    pin_cmd.add_argument("--version", help="with a PR: <image-version>[-<tool-version>], if docker_build.sh "
                                           "computes part of its tag")
    pin_cmd.add_argument("--dry-run", action="store_true", help="show the changes; write and tag nothing")
    pin_cmd.add_argument("--include-shared-verification", action="store_true", help="also move pins in " + SHARED)
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
        # IMAGE:TAG never looks like a PR: git branch names cannot contain ":".
        existing = re.fullmatch(r"(?:%s/)?(%s):(%s)" % (re.escape(registry), NAME, TAG), args.target)
        if not existing:
            pin_pr(root, registry, args.target, args.paths, args.digest, args.version,
                   args.include_shared_verification, args.dry_run)
            return 0
        if args.digest or args.version:
            raise Error("--digest and --version go with a warp-tools PR, not an existing tag")
        image, tag = existing.groups()
        digest = digest_of(registry, image, tag)
        if digest is None:
            raise Error("%s/%s:%s not found" % (registry, image, tag))
        changed, skipped = pin(root, registry, image, tag, digest,
                               args.include_shared_verification, args.dry_run, args.paths)
        report_pin(registry, image, tag, digest, changed, skipped, args.dry_run)
        return 0
    except Error as err:
        print("error: %s" % err, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
