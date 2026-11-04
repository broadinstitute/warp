# End-of-Session Reflection Prompt

Paste the prompt below to an AI coding agent at the end of a substantial,
AI-assisted session — typically just before opening a PR for review. Its purpose
is to turn what the agent learned into durable guidance in [AGENTS.md](AGENTS.md)
for future agents, while keeping that file lean.

This is a **user-invoked** prompt, not a standing instruction: run it deliberately,
review the proposed diff, and only then apply it.

---

```
Before I open this PR, reflect on everything you learned this session about this
repository and how to work in it effectively.

1. Identify durable, generalizable lessons — conventions, gotchas, dependency
   chains, validation steps, or mistakes you made and corrected — that would save
   a FUTURE agent time. Ignore anything specific to this one task.

2. Review AGENTS.md. For each lesson, decide:
   - Already covered (in AGENTS.md or a doc it links to) → skip it.
   - Agent-operational and missing → propose a concise addition to the right
     section of AGENTS.md. Keep it short and link out to the canonical doc rather
     than restating rules.
   - Really a human-facing rule/rationale → don't put it in AGENTS.md; flag it as
     a suggested change to the relevant website doc instead.

3. While you're in AGENTS.md, eliminate any redundancy you spot — duplicated
   guidance, stale paths, or content that now belongs in a linked doc.

4. Verify every claim against the actual repo (file paths, task names, commands)
   before writing it — do not add aspirational or unverified notes.

Show me the proposed AGENTS.md diff and a short rationale for each change before
applying anything.
```
