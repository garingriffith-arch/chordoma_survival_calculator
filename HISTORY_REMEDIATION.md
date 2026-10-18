# Historical restricted-data remediation

**Status: unresolved.** The current revision working tree excludes NCDB-derived fitted RDS objects, but historical commits on the public repository previously contained such objects. An ordinary file deletion does not remove an object from Git history.

## Verified current exposure (18 October 2026)

Repository visibility was **public** during this audit. The current `main` branch itself includes two fitted RDS blobs in `data/processed/` (about 224 KB and 212 KB); the current 2024 revision branch has no such RDS blobs, but its ancestors include the old `main` lineage. The repository has two ordinary branches (`main`, `2024-tripod-reproducibility`), one open draft pull request (#1), and no reported forks. **Temporarily set repository visibility to private** as the first action. This document cannot guarantee the removal of any copies already downloaded.

A two-phase Windows PowerShell cleanup package was prepared separately for the repository owner: a fresh mirror clone followed by `git-filter-repo --sensitive-data-removal`, comprehensive read-only verification of the rewritten mirror, and a distinct administrator-approved force-push step. None of those history-changing or visibility-changing steps has been executed by the assistant. This repository must not be advertised as clean before the new remote refs, PR caches and affected local clones are addressed.

## Scope

Before publication, audit **all refs and historical objects** for restricted NCDB records, not only the current branch files. Known paths in the old repository include `data/processed/chordoma_model_objects.rds` and `data/processed/chordoma_treatment_model_objects.rds`. The old analysis saved complete-case datasets inside fitted RDS objects; these cannot be treated as aggregate coefficients or safely redistributed.

The repository’s latest 2024 pipeline uses CSV coefficients and reference survival and does not need patient-level fitted objects for the Shiny calculator.

## Remediation procedure (administrator action; not yet executed)

1. Back up any needed source code to a **private, restricted location**. Do not publish a mirror archive or historical Git bundle containing NCDB-derived records.
2. Wait until the new manuscript-exhibit stage has been run and reviewed locally; check every present and historical file for patient-level records, restricted model objects, and any additional embedded values.
3. Use a new, fresh mirror clone and `git-filter-repo` to remove every restricted object from *all branches and tags*. Review the filtered refs and object inventory **before** any forced push.
4. Coordinate with collaborators because forced updates to public refs invalidate existing local commit histories and can disrupt pull requests.
5. Force-update the intended branches and tags only after explicit administrator approval. Re-clone from scratch and audit the resulting public refs.
6. Contact GitHub Support about cached views, PR refs, and server-side object retention. A force-push cannot guarantee deletion of copies already downloaded, cached, or referenced from pull requests.
7. If complete removal cannot be established or DUA requirements demand it, consider a freshly initialized repository with no connection to the old Git history, together with repository takedown/support escalation as appropriate.

**Do not rewrite or force-push history casually.** A successful syntax check or new clean branch is not evidence that old restricted records are no longer accessible. The repository should not be described as fully remediated until these steps have been verified.

## Proposed local verification commands (read-only)

From a **fresh mirror clone**, not from the active Windows analysis directory:

~~~powershell
git rev-list --objects --all | Select-String '(?i)\.(rds|rda|RData)$'
git rev-list --objects --all | Select-String '(?i)(data/raw/|data/processed/|internal_do_not_share/|outputs/)'
~~~

A zero-result scan is necessary but not sufficient: encrypted archives, mislabeled files, GitHub caches and external clones require separate consideration. Publication does not require depositing any restricted NCDB participant-level objects.
