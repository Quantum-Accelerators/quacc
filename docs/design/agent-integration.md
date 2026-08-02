# Design: Phased Agent Integration for quacc

**Status:** Draft / proposal
**Tracking issue:** [#3154](https://github.com/Quantum-Accelerators/quacc/issues/3154)
**Audience:** quacc maintainers

quacc adopts agent tooling — context files, skills, CI automation, machine-readable docs, and (maybe) an MCP server — to help two different populations:

| Population | Needs help with |
| --- | --- |
| **Contributors** | Recipe/flow conventions, docstrings, `ruff`, test topology, per-engine test suites |
| **Users** | Finding the right recipe among ~78 jobs across 75 modules, composing flows, choosing among 6 workflow engines, deciphering a failed run |

An integration that only helps contributors does nothing for users, and vice versa. Each phase below names which population it serves and ships independently — no phase depends on a later one being built.

---

## Guiding principles

1. **No hard dependency on any agent vendor.** Nothing in quacc's runtime, tests, or docs build may require an API key.
2. **Provider-neutral by default.** Root context file is `AGENTS.md`, not `CLAUDE.md`; skills live in a top-level `skills/` folder using the open [agentskills.io](https://agentskills.io/what-are-skills) format, not `.claude/skills/`. Vendor-specific files (`.claude/settings.json`, `.github/workflows/claude.yml`) stay, but only where the vendor-specific ergonomics are the point — a CI job has to run some concrete agent, and quacc's already runs Claude.
3. **Reuse the source of truth.** Recipe metadata, settings docs, and schemas already exist in the code. Agent-facing artifacts must be generated from or point at them — never hand-copied, since a hand-copied fact rots the first time the code changes.
4. **A skill ships only if it closes a real gap** (see below) — otherwise it's ceremony.
5. **Every skill needs a verification check in CI.** No skill ships without one (Phase 6).

**Non-goals:** replacing human code review or auto-merging; letting an agent submit/manage jobs on shared HPC resources; committing agent-generated recipes without a human author who understands the science.

---

## The three-gap test

A skill (or index, or context file) is worth building only if it closes one of these — this is the filter applied to every deliverable below:

1. **The knowledge isn't in the code.** A judgment call, not a fact — e.g. "prefer EMT for solids, LJ for molecules as the cheap-calculator-first choice" (from `contributing.md`), or "use Parsl on a Slurm cluster, Dask on a single multicore node."
2. **The knowledge is in the code, but pattern-matching it is error-prone.** quacc has two functions named `relax_job` (`recipes/vasp/core.py` and `recipes/vasp/slabs.py`) — an agent can often pick the right one by inferring from context, but that's a guess, not a guarantee, across ~78 jobs.
3. **The knowledge changes faster than an agent's memory, or the agent has no way to look it up at all.** Exact CLI flags, current schema shape — stale in a general model's training data, and unavailable entirely to an agent with no filesystem access to the repo (e.g. a hosted assistant, no local clone).

**A finding that reshapes Phase 5, from a live test:** asking a local Claude Code session (with the quacc repo cloned) "how do I relax a slab with VASP" produced a fully correct answer — it grepped `def relax_job`, read the matching file, and explained the slab-specific settings (`isif=2`, `auto_dipole=True`) — with **no skill, index, or `AGENTS.md` in place at all**. For a user who already has the repo cloned and a coding agent with Bash/grep, gap #3 is often already closed for free, and even gap #2 can resolve correctly by lucky inference (it picked `slabs.py` over `core.py` from the filename matching "slab"). This means the highest-value remaining work for users is **not** "look up a recipe by name" — grep already does that — it's the cases grep doesn't cover: vague queries with no obvious search term, genuine disambiguation between near-duplicates, and users with no local checkout to grep at all.

---

## Phases

| Phase | Serves | Ships | Portable? | Effort |
| --- | --- | --- | --- | --- |
| 0 | Maintainers | Policy, spend cap, CI trigger policy | n/a | S |
| 1 | Contributors | `AGENTS.md` | Yes | S |
| 2 | Contributors | `skills/{new-recipe,new-calculator,recipe-tests,docs-sync}/` | Yes | M |
| 3 | Maintainers | Tuned `claude.yml` + advisory PR-review workflow | No | M |
| 4 | Users (any assistant) | `llms.txt`, generated recipe index, generated settings reference | Yes | M |
| 5 | Users | `skills/quacc-workflow/`, scoped to discovery/disambiguation/engine choice; MCP only if that proves insufficient | Yes → maybe not | M → L |
| 6 | Everyone | Eval suite, per-skill CI verification, ownership, sunset rule | n/a | M |

### Phase 0 — Decide, and set guardrails
Maintainer sign-off before anything lands: an `AGENTS_POLICY` note in `contributing.md` (agent-assisted PRs welcome; human author owns correctness; visibly unreviewed agent output gets closed); confirm `ANTHROPIC_API_KEY` ownership and spend cap; artifacts live in-repo; CI trigger restricted to collaborators.

### Phase 1 — `AGENTS.md`
Root file, <120 lines: orientation to `src/quacc/{recipes,calculators,runners,schemas,wflow_tools,settings.py}`; the recipe contract (`@job`/`@flow`, pickle-serializable I/O, `Atoms` first arg, `quacc.schemas` return, no `.` in dict keys); hard rules an agent violates by default (absolute paths only, never `os.chdir`/`os.getcwd()`, NumPy docstrings, type hints, `ruff check --fix && ruff format`, test required, gzip large fixtures); test topology (`tests/core`, per-engine dirs, monkeypatched `conftest.py`, `--noconftest`, don't run licensed codes locally); pointers (not copies) to `docs/dev/contributing.md` and `docs/dev/recipes/{jobs,flows}.md`. Plus `.claude/settings.json` (allowlist `pytest`, `ruff`, `git status/diff/log`, `pip install -e`) and a `.gitignore` entry for `.claude/settings.local.json`.

**Acceptance:** a fresh agent session asked to "add an `md_job` to `recipes/lj`" produces `@job`, NumPy docstrings, a schema return, no `os.getcwd()`, and a test — verified by running it, not by reading the file.

### Phase 2 — Contributor skills
Top-level `skills/<name>/SKILL.md` (agentskills.io format), installed per-agent via `npx skills add ./skills -a <agent>`:

| Skill | Gap closed | Does |
| --- | --- | --- |
| `new-recipe` | #2 | Scaffold a `@job`/`@flow` from the closest existing recipe, correct schema/runner, matching test; encodes the cheap-calculator-first rule (#1). |
| `new-calculator` | #2 | Calculator dir + presets + `pyproject.toml` entry + optional extra + `requirements-<code>.txt` + monkeypatched `conftest.py` if needed. |
| `recipe-tests` | #2 | Write/repair tests for an existing recipe: small system, low runtime, regression coverage for bug fixes. |
| `docs-sync` | #2 | Update `docs/user/`/`docs/dev/` and nav config after a code change — the most common review comment on recipe PRs. |

Cap at ~4 skills; delete any unused for 3 months. Each needs a Phase 6 verification check before it ships.

### Phase 3 — CI automation
Tune the already-live `.github/workflows/claude.yml` (add a prompt pointing at `AGENTS.md`, explicit `claude_args` allowlist, collaborator gate); add a non-blocking `claude-review.yml` against the Phase 2 checklist (never a required check, never auto-approves). This phase is inherently vendor-specific — CI has to run some concrete agent.

### Phase 4 — Machine-readable docs
`llms.txt`/`llms-full.txt` generated at docs-build time; a generated recipe index (import path, one-line summary, calculator, calc type — introspected from `quacc.recipes`, can't drift from code); a generated settings reference from `QuaccSettings`. This is what serves users with no local checkout to grep, and is the input Phase 5's skill consumes.

**Acceptance:** given only `llms.txt`, an assistant with no quacc training data writes a runnable EMT relax + slab flow and a correct `quacc set WORKFLOW_ENGINE` call.

### Phase 5 — User-facing skill (MCP only if needed)
`skills/quacc-workflow/SKILL.md`, scoped to what grep-on-a-local-checkout genuinely can't do well (see the three-gap test finding above):
- **Discovery without an obvious search term** — "what's the best way to compute an adsorption energy" has no guessable `grep` string; the skill searches the Phase 4 index instead.
- **Disambiguation** — quacc has two `relax_job`s (`core.py`, `slabs.py`) and similar near-duplicates elsewhere; state the rule explicitly instead of relying on an agent inferring correctly from a filename.
- **Workflow-engine choice** — decision logic for `dask`/`parsl`/`prefect`/`ray`/`redun`/`jobflow` based on the user's infra, plus the config gotchas (e.g. `WORKFLOW_ENGINE` can't change inside a context manager — must use `export QUACC_WORKFLOW_ENGINE=...`).
- **Flow composition** — the `@job`/`@flow` composition pattern from `docs/user/basics/` for multi-step workflows (relax → static → phonon).

An MCP server (`quacc mcp`, read-only, no job submission, no credential handling) is built later, only for what a static index can't cover — live settings resolution, reading an actual failed run directory — and only if a named maintainer will own it long-term.

**Acceptance:** a user with no docs read can get from a structure file to a correctly-composed (not necessarily submitted) workflow with correct imports and no hallucinated kwargs.

### Phase 6 — Evaluation and upkeep
A per-skill CI check that referenced paths/commands still resolve (deepmd-kit's own `skills/` has zero such wiring — worth not repeating); a ~15-task eval set (5 contributor, 5 user, 5 diagnostic) run before/after each phase; a named owner per artifact in `docs/maintainers/internal.md`; sunset rule — unused for two release cycles, delete it; quarterly check that `AGENTS.md`/skills still match `docs/dev/`.

---

## Open questions

1. Restrict `AGENTS.md`'s scope to what's true today, or also to what Claude Code specifically reads (it reads `AGENTS.md` natively; other agents' support varies)?
2. In-repo `skills/` vs. a separate repo? (Recommend in-repo — must version with the conventions it describes.)
3. CI trigger policy: collaborator-gated or open to all?
4. Does anyone want to own an MCP server long-term? If not, stop at Phase 5.
5. Spend ceiling for the GitHub Action?
6. Should `docs/design/` be in the published nav, or excluded entirely?

---

## Decision log

| Date | Decision | Rationale |
| --- | --- | --- |
| 2026-08-01 | `AGENTS.md`, not `CLAUDE.md` | Provider-neutral; Claude Code reads it natively anyway |
| 2026-08-01 | Skills in top-level `skills/` (agentskills.io format), not `.claude/skills/` | Portable via `npx skills add -a <agent>` |
| 2026-08-01 | Phase 5 skill scoped to discovery/disambiguation/engine-choice, not general recipe lookup | Live test showed grep+read on a local checkout already solves named-recipe lookup for free |
| 2026-08-01 | Every skill ships with a CI verification check | Avoid the staleness gap deepmd-kit's own `skills/` has (no CI wiring at all) |
| 2026-08-01 | Any MCP server is read-only, no job submission, no credential handling | Safety/support burden on shared HPC resources |
| 2026-08-01 | Nothing in quacc may require an API key | Core installs and CI must work unchanged |
