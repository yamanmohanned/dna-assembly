# tools/

Third-party programs that are built or unpacked locally and are **not committed**.

## LKH-3 (reference ATSP solver)

- Source: the official LKH-3 page of Keld Helsgaun. Last checked 2026-10-22: `webhotel4.ruc.dk` and `akira.ruc.dk` were unreachable from this machine (see `docs/phases/PHASE_0.md`).
- Place the downloaded tarball here, for example `tools/LKH-3.<version>.tgz`. The next step is to record its SHA-256 in `docs/DECISIONS.md` and build it (WSL once `build-essential` is installed, or native MinGW if it compiles; both are TODO-VERIFY).
- Licence: academic use. TODO-VERIFY the exact wording on the official page before redistributing anything.
- Rule R9: nothing other than this program is ever called "LKH-3" in results, figures or slides.
