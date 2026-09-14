# Public release preflight for v2.0.0

Date: 2026-09-14

## Scope

This candidate aligns the public repository with the corrected Scientific Reports submission candidate. The v1.0.0 tag and Zenodo record remain unchanged. The public copy adds final v3.1 correction records, bounded v4 repair/closeout artifacts and the compact v5 public clinical-extension archive.

## Checks completed

- Required final strict membership, association, module-deletion, spatial and clinical-evidence records are present.
- The v5 compact archive manifest verifies after portable-path conversion.
- The repository-wide SHA256 manifest covers all release-candidate files except itself.
- All included Python sources compile.
- No release-candidate file is 95 MB or larger.
- No user-specific local project path or temporary container path remains in public candidate files.
- Internal rendering files, author-review manuscripts, contact sheets and execution logs are excluded through `.gitignore`.
- No FASTQ, IDAT, H5, H5AD, CLOUPE, controlled-access file or complete large public matrix is added.
- The clinical-extension semantic checks remain 10/10 passing.

## Publication state

This is a local release candidate. It has not been pushed, tagged, archived by Zenodo or assigned a new DOI. `CITATION.cff` intentionally omits a v2 DOI until Zenodo creates the immutable version record.
