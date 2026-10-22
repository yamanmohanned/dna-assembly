# Data loading: from an NCBI accession to a checked genome in memory

Code: [genomeasm/io.py](../genomeasm/io.py), [genomeasm/datacard.py](../genomeasm/datacard.py),
[genomeasm/cli.py](../genomeasm/cli.py). Config: [configs/genomes.yaml](../configs/genomes.yaml).
Tests: `tests/test_fasta_parser.py`, `test_genome.py`, `test_cache.py`, `test_datacard_cli.py`,
`test_real_genomes.py`. Numbers about the genomes are **not** repeated here: they live in the data
cards `data/processed/<accession>.card.md`, generated from the saved files (rule R3).

## The story in plain words

A genome is one very long word written with four letters: A, C, G, T. NCBI, a big public library
for DNA, keeps the official copies. Each copy has a catalogue number (an *accession*) and an edition
number (the *version*), for example `NC_001422.1`.

We go to the library **once**, copy the word, and let the library tell us a few facts about it
(name, kind of molecule, length). We check that our copy and the library's facts agree. Then we save
the copy with a fingerprint (a *SHA-256 checksum*). From then on the program only reads the saved
copy, after checking the fingerprint. So the program works with no internet, and if even one letter
of the saved file changes, the program notices.

## The pipeline, step by step

```
configs/genomes.yaml            list of versioned accessions
        |  load_accessions()
        v
  validate_accession()          format check, e.g. NC_001422.1 (no version = rejected)
        |
        |  ---- only for `genomeasm data fetch` (needs internet) ----
        v
  fetch_from_ncbi()             Entrez esummary (facts) + efetch (FASTA text)
        v
  parse_single_record()         FASTA rules (below); exactly one record
        v
  Genome.from_record()          header id == accession; only A,C,G,T (rule D-016)
        v
  cross-checks with NCBI        accession version, molecule type = dna, length = NCBI's length
        v
  save (atomic)                 data/raw/<acc>.fasta  +  <acc>.meta.json  +  SHA256SUMS
        |
        |  ---- every normal run (offline) ----
        v
  load_genome()                 SHA-256 of the file == SHA256SUMS ?  else ChecksumMismatchError
        v                       parse -> Genome.from_record -> compare with <acc>.meta.json
  LoadedGenome                  genome + NCBI metadata + file fingerprint
        v
  build_card()                  data/processed/<acc>.card.json and <acc>.card.md
```

If **any** check in the download part fails, **nothing is written**. If NCBI cannot be reached,
the program stops and says so. It never makes up sequence data (rule R1).

## Files

| File | What it is | Written by |
|---|---|---|
| `data/raw/<acc>.fasta` | The exact bytes NCBI sent (header line plus 70-letter lines) | `fetch` |
| `data/raw/<acc>.meta.json` | NCBI's own summary of the record (organism, molecule type, length, topology, strand, dates), plus when and how we downloaded it | `fetch` |
| `data/raw/SHA256SUMS` | One line per FASTA file: `<64 hex digits>  <file name>`. Same format as the `sha256sum -c` tool | `fetch` |
| `data/processed/<acc>.card.json` | The data card as data (read by the report and the app) | `cards` |
| `data/processed/<acc>.card.md` | The data card as a table for people | `cards` |

Two fingerprints are recorded. The **file SHA-256** covers the exact bytes of the FASTA file, line
breaks included. The **sequence SHA-256** covers only the uppercase letters, so it stays the same if
a file is re-wrapped to a different line width.

## FASTA parsing rules (`parse_fasta`)

| Input | What happens |
|---|---|
| Unix, Windows (`\r\n`) or old Mac (`\r`) line endings; blank lines; spaces around a line | Accepted |
| Lowercase letters (`acgt`) | Accepted, turned into uppercase; the number of lowercase letters is recorded |
| IUPAC ambiguity codes `R Y S W K M B D H V N` | Accepted **by the parser** (so they can be counted), rejected **by `Genome`** (see below) |
| Empty text or only blank lines | `FastaFormatError` |
| Text before the first `>`; empty header (`>` alone); header with no sequence | `FastaFormatError` |
| Digits, `-`, `*`, `U`, `;` comments, non-ASCII, spaces or tabs inside a sequence line | `FastaFormatError` with the line number |
| More than one record where one genome is expected | `FastaFormatError` |

## The alphabet rule for ambiguity codes (decision D-016)

A letter such as `N` means "the machine could not tell which base this is". The project's rule is:
**a genome must contain only A, C, G, T.** A genome with any other letter raises
`AmbiguousBaseError`, which lists how many there are and where the first ones are. Nothing is
repaired, dropped (that would shift all positions and break the ground truth of simulated reads) or
guessed. If a future genome needs ambiguity handling, the safe extension is to split it into
ACGT-only pieces; that is deliberately not built yet.

## What each check protects against

| Check | Protects against |
|---|---|
| Versioned accession required | Silently getting a newer edition of the sequence later |
| FASTA identifier equals the requested accession | Downloading the wrong record |
| Pure A/C/G/T | Ambiguity codes corrupting overlap computations and ground truth |
| NCBI summary agrees (accession version, `dna`, length) | A truncated or damaged download |
| SHA-256 verified on every load | A changed or corrupted file on disk |
| Atomic writes, nothing saved after a failed check | Half-written caches |
| Loading never touches the network | Surprises during the offline presentation |

## Things this pipeline does not do (declared simplifications)

- **Circular and double-stranded genomes are treated as one linear string.** NCBI lists phiX174 as
  circular single-stranded DNA and lambda as linear double-stranded DNA (see the data cards). A read
  that would cross the origin of a circular genome, and reads from the reverse-complement strand, are
  not simulated in the baseline (decision D-017).
- Only single-record DNA FASTA files are loaded.
- The Entrez contact email is not sent (decision D-010); this is allowed but NCBI would not be able
  to contact us.

## Run it yourself

```bash
genomeasm data fetch      # needs internet; skips genomes that are already cached and verified
genomeasm data verify     # offline: re-checks every checksum and every cross-check
genomeasm data cards      # offline: rewrites data/processed/*.card.json and *.card.md
```
