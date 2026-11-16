# Phase 9: molecular-biology pages, 3D viewer, Presentation Mode (D-055)

Requested by the user on 2026-11-14 as an enhancement after the frozen release. Nothing about the
experiments changed.

## What was added
- `genomeasm/molbio/`: `sequence` (complement, transcription, genetic code, open reading frames,
  mutation classifier), `structure` (PDB reader, covalent bonds, Watson-Crick pairs, helix
  measurement), `model` (scenes for the viewer), `processes` (replication, transcription,
  translation scenarios), `diagrams` (SVG), `library` (references, facts, genes), `quiz`.
- `genomeasm/appkit/`: `biology` (what the pages show), `lecture` (28 slides), `studio`
  (self-contained HTML), `web/` (3D engine, molecule view, animation player, deck).
- Pages: Presentation, DNA structure, DNA to protein, Genes and mutations, References.
- Data: `data/raw/1BNA.pdb` (+ metadata, checksum), `configs/phix174_genes.yaml`,
  `docs/references.yaml`, `docs/molbio_facts.yaml`.
- Command: `genomeasm app lecture [--author NAME]` writes `reports/lecture/dna_lecture.html`.

## Checked, and how
| Item | Check | Result |
|---|---|---|
| Genetic code | compared with Biopython tables 1 and 11 (test) | equal |
| 11 genes of phiX174 | start codon, whole codons, one final stop (test) | all clean |
| 1BNA | 12 pairs, 2 or 3 hydrogen bonds per pair, lengths 2.5 to 3.5 angstrom (test) | 32 bonds, 2.61 to 3.25 |
| Helix | rise 3.4 +/- 0.15, 9.6 to 10.6 bp per turn, positive twist (test) | 3.37, 10.1, right-handed |
| 39 journal references | Crossref record of each DOI, 2026-11-14 | all match |
| Quoted numbers | abstract of the PubMed record where it states them | see `docs/molbio_facts.yaml` |
| Deck in the browser | every slide measured through the DOM at 1600 x 900: no overflow, no clashing labels | 28 of 28 |
| Viewer in the browser | WebGL 2 drawing, picking, selection card, 2D fallback, freeing on slide change | worked |
| Animations in the browser | frames of every step looked at for the three processes | correct against the sequences |

## Not done or not verified
- Textbook values (fork speed, fidelity, Okazaki fragment length, 46 chromosomes, the conclusions
  of four classic papers) are TODO-VERIFY: NCBI Bookshelf answered with a CAPTCHA.
- Full screen (key F) was not tried: it needs a real key press in a visible window.
- Screenshots of the new Streamlit pages were taken only for DNA structure and Presentation; the
  Browser pane was hidden afterwards, so the other pages were checked by the page tests only.
- The deck was looked at in one browser engine (Chromium) at one size (1600 x 900).
- `genomeasm docs screenshots` and `docs verify` were not rerun; the report, the pptx and
  `docs/GUIDE.md` were not changed and do not mention the new pages.
