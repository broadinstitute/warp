# FISS-MCP runbook — debugging a failed WARP Terra test

**Scope:** debugging WARP's **Terra-based WDL pipeline tests** (the Plumbing/Scientific CI
tests that run a `Test<Pipeline>` workflow on Terra and compare outputs to truth). Nothing
else. It assumes the `fiss-mcp` tools are available. Goal: go from a Terra submission URL to
the exact line that failed, fast, without fumbling.

## 0. Parse the URL

```
https://app.terra.bio/#workspaces/<NAMESPACE>/<WORKSPACE>/submission_history/<SUBMISSION_ID>
```
e.g. `warp-pipelines / WARP Tests / 22618e75-...`. Note `%20` = space in the workspace name ("WARP Tests").

## 1. Submission → which workflow failed

`get_submission_status(namespace, workspace, submission_id)`
- `status_summary` tells you Failed/Succeeded counts.
- Grab each `workflows[].workflowId`. A submission usually has ONE workflow (the `Test<Pipeline>` wrapper).
- If `status: Succeeded` you're done — that URL is green, stop.

## 2. Workflow → which task failed (summary mode ALWAYS first)

`get_job_metadata(..., workflow_id, mode="summary")`
- Cheap (~1-2K tokens). Returns `failed_tasks[]` and `workflow_failures[].causedBy` chain.
- The `causedBy` message names the real culprit, e.g.
  `Job VerifyMultiome.CompareAtacLibraryMetrics:NA:1 exited with return code 1`.
- **Gotcha:** the top-level failed task is often `Test<Pipeline>.Verify` — that is a *sub-workflow*, not the real failure. Its `stderr_url`/`stdout_url` are empty at the parent level. The real task lives one level down (`VerifyMultiome.<Task>`). Go to step 3.

## 3. Find the real task's logs (subworkflow drill-down via GCS)

`get_job_metadata` extract mode CANNOT path into subworkflows here — call keys contain dots
(`TestMultiome.Verify`) and the dot-path parser chokes (`Key 'TestMultiome' not found`). **Don't fight it.** Walk the GCS tree instead with `list_gcs_objects(recursive=false)`:

```
gs://<bucket>/submissions/intermediates/<SUB_ID>/<Wrapper>/<WF_ID>/call-Verify/
  └─ <SubWorkflowName>/            (e.g. VerifyMultiome/)
       └─ <SUBWF_ID>/              (a fresh UUID, different from WF_ID)
            └─ call-<TaskName>/    (e.g. call-CompareAtacLibraryMetrics/)
                 ├─ stdout   ← comparison output usually goes HERE
                 ├─ stderr   ← often 0 bytes; don't assume the error is here
                 ├─ rc       ← 2-byte file: the return code
                 └─ script   ← the exact bash that ran, if you need it
```
Get the bucket from any `*_url` in step 2 (`gs://fc-...`). List each level with `recursive=false` and read the `prefixes[]` to descend one hop at a time.

## 4. Read the log

`read_gcs_object(gcs_uri, offset, max_bytes)`.

Prefer `stdout` over `stderr` for WARP verification tasks — the python comparison scripts
print PASS/FAIL to stdout; stderr is frequently empty. For a comparison task, the lines you
want are `Error: Metric ... exceeds threshold`.

**Sandboxed-agent caveat (does NOT apply to a human with gsutil):** in the agent container,
`read_gcs_object` returns an opaque `<<ccr:...>>` reference for large reads (a context-
compression layer you can't actually read), and it *dedups identical re-reads* to the same
ref. So: keep `max_bytes` small (≤ ~1400 → real UTF-8 text), and walk the file in
**contiguous non-overlapping windows** (`offset` 0, 400, 1800, 2500, …) — never re-read the
same range; if a window comes back compressed, shift the offset instead of retrying it.

## 5. Other gotchas

- **No `gsutil`/`gcloud` in the agent container.** Use the MCP GCS tools (`list_gcs_objects`, `read_gcs_object`, `get_gcs_object_metadata`), not shell. (A human on a laptop can just `gsutil cat` the stdout path from step 3 and skip step 4's windowing entirely.)
- `get_workflow_logs(fetch_content=false)` gives you all task stderr/stdout URLs in one shot — handy for the *parent* tasks, but it won't reach into subworkflows (their entries have empty URLs). Use step 3 for those.
- For infra-looking failures (RC 137 / "stopped before command finished" / 0-second tasks) use `get_batch_job_status` — those errors are NOT in the GCS stderr.
- `cost` in the status is real money spent; a green submission still cost a few cents.
- A **stale** failure is common: confirm the submission's date against the branch HEAD. A run from before a fix landed will show the old failure.

## 6. One-screen cheat sheet

```
get_submission_status            → workflowId, Failed?
get_job_metadata mode=summary    → failed task + causedBy (spot the .Verify subworkflow)
list_gcs_objects recursive=false → descend call-Verify/<SubWF>/<uuid>/call-<Task>/
read_gcs_object                  → read stdout (human: gsutil cat; agent: small windows, see §4)
```
