# genome-annotation-pipeline

**Version 1.0.0** · Python ≥ 3.10 · **No third-party dependencies**
(pure standard library) · MIT license

An end-to-end genome annotation pipeline for FASTA files, in four
subcommands:

| Subcommand | Input | What it does |
| --- | --- | --- |
| `qc` | nucleotide FASTA | Assembly quality control: record counts, lengths, GC content, N50/L50 |
| `orf` | nucleotide FASTA | Six-frame open reading frame (ORF) finding with the standard genetic code |
| `props` | protein FASTA | Protein properties: molecular weight, amino-acid composition, GRAVY, theoretical pI |
| `annotate` | nucleotide FASTA | The whole pipeline at once, rendered as Markdown **and** HTML annotation reports |

Under the subcommands sit a strict multi-record FASTA parser, the
standard genetic code built programmatically from NCBI's compact
representation, and property calculations whose constants are
transcribed from the primary sources (ExPASy mass and hydropathy
tables, the Bjellqvist pKa set) — see `docs/RESEARCH.md` for every
source and design note. The `annotate` report is verified, not
merely produced: its example genome is synthetic and *constructed*,
with three genes planted at known coordinates, and the test suite
re-derives that gene set with an independent scanning method before
comparing it against the report. An end-to-end integration test
additionally carries one genome through all four subcommands via the
public CLI and requires their outputs to agree with each other.

The project was built in five stages — parser + QC, ORF finding,
protein properties, annotation reports, and this 1.0.0 integration
and polish release; the roadmap at the bottom lists them.

## Background and theory

### FASTA

FASTA is the de-facto plain-text format for nucleotide and protein
sequences. Each record starts with a single-line *definition line* (header)
beginning with a `>` character in the first column; the first word after
`>` is the record identifier and the rest of the line is a free-text
description. The sequence itself follows on one or more lines until the
next header or end of file, so one file can hold many records
(multi-FASTA). Letter case carries no meaning, and sequences are widely
wrapped at 60–80 characters — a parser must therefore join wrapped lines
and normalise case before doing anything else.

This parser accepts the **IUPAC nucleotide alphabet**: the canonical bases
`A C G T`, `U` (uracil, normalised to `T`), and the ambiguity codes
`R Y S W K M B D H V N`. Anything else — for example `X`, digits, or gap
symbols — is rejected with an error that names the record, the line, and
the offending character. Failing loudly at parse time is cheaper than
annotating a corrupted record later.

### GC content

For a sequence `s` of length `L`:

```
GC% = 100 * (count of G + count of C) / L
```

The assembly-wide figure uses the same formula over all records combined
(total G+C over total bases), matching the QUAST definition. Ambiguity
symbols are counted in the length but not as G or C; the per-record report
shows how many ambiguity symbols each record contains so the effect is
never hidden.

### N50 and L50

N50/L50 summarise how contiguous an assembly is. Sort the contig lengths
from longest to shortest and walk down the list, accumulating lengths.
The first point where the running total reaches at least half of the
total assembly length gives:

- **N50** — the length of the contig at that point (equivalently: the
  length for which all contigs of that length or longer cover at least
  half of the assembly), and
- **L50** — how many contigs were needed to reach that point.

Worked example, also used as a unit test: lengths `[10, 20, 30, 40]`,
total 100, half is 50. Sorted descending: 40 (cumulative 40 < 50), then
30 (cumulative 70 ≥ 50). So **N50 = 30** and **L50 = 2**.

Note that N50 is computed against the *assembly* size. The related NG50
uses a known/estimated genome size instead; this tool reports N50 only.

### Open reading frames and translation

An **open reading frame (ORF)** is a sequence of successive nucleotide
triplets, read as codons, that begins with an initiator codon and ends
with a stop codon (the NCBI MeSH definition). Because a gene can sit on
either strand, a double-stranded sequence has **six reading frames**:
three on the given strand, three on its reverse complement. The ORF
finder searches all six for ORFs that begin with the start codon `ATG`
and end at the first in-frame stop codon — `TAA`, `TAG`, or `TGA` in the
standard genetic code (NCBI translation table 1).

Translation uses that same standard code. The codon table is *built
programmatically* from NCBI's compact representation of the code (one
64-character amino-acid string whose positions are enumerated by three
base strings) instead of 64 hand-typed dictionary entries, and unit
tests pin well-known codons so a construction error cannot hide. Two
edge cases are explicit: a trailing partial codon (sequence length not
a multiple of three) is ignored, and a codon containing an IUPAC
ambiguity symbol cannot be resolved to one amino acid, so it translates
as `X`.

How ORFs are reported:

- **Coordinates** are 1-based and inclusive on the *original* sequence,
  always with start < end, plus the strand (`+`/`-`) and the frame
  (1–3, counted on the strand the ORF was found on).
- **Minimum length** is measured in amino acids *excluding* the stop
  codon (`--min-aa`, default 30). Random sequence produces a stop codon
  roughly every 21 codons, so short ATG-to-stop stretches are mostly
  chance; the threshold keeps them out of the default report.
- **Nested starts:** within one frame, an ORF runs from the *first*
  `ATG` after the previous stop to the *first* stop that follows it —
  one ORF per stop per frame, the longest one. (The alternative,
  reporting every ATG-to-stop pair including all nested sub-ORFs, was
  rejected: gene callers annotate the longest ORF per stop.)

### Protein properties

Stage 3 of the build added protein properties for a sequence itself —
the natural input is a translation produced by the ORF finder, read
back from a protein FASTA file
with a dedicated parser mode that accepts the 20 standard amino acids
plus `X` (unknown residue) and rejects everything else (including the
ambiguity codes `B`/`Z` and the stop symbol `*`). Four properties are
computed:

- **Molecular weight (MW).** Each amino acid inside a chain contributes
  its *residue mass* — the free amino acid minus the water lost when its
  peptide bonds formed — so a chain of *n* residues is short exactly one
  water overall:

  ```
  MW = sum of residue masses + mass of one H2O
  ```

  The masses are average isotopic values from ExPASy's amino-acid mass
  table (the table ExPASy Compute pI/Mw uses). Example, worked by hand:
  `ACD` = 71.0788 + 103.1388 + 115.0886 + 18.01524 = **307.32 Da**.
- **Amino-acid composition.** Counts of each of the 20 standard
  residues, with percentages of the total sequence length.
- **GRAVY** (grand average of hydropathy): the arithmetic mean of the
  Kyte & Doolittle (1982) hydropathy values of the residues.
  Positive means hydrophobic overall, negative means hydrophilic;
  poly-alanine is exactly 1.800, poly-lysine exactly −3.900.
- **Isoelectric point (pI).** The pI is the pH at which the protein's
  net charge is zero. The charge model sums Henderson–Hasselbalch
  fractional charges over the N- and C-termini and the ionisable side
  chains (K, R, H positive; D, E, C, Y negative) using the **Bjellqvist
  pKa set** — the values ExPASy Compute pI/Mw uses, including
  residue-specific terminal pKa corrections. Net charge falls strictly
  monotonically as pH rises, so the zero crossing is found by bisection
  on [0, 14]. This is an *estimate*: the Bjellqvist values were
  calibrated on polypeptides focusing in immobilised pH gradients, and
  no sequence-only method can see the local environments of a folded
  protein. As an external check, the implementation reproduces the
  published reference values pI(`INGAR`) = 9.75 and pI(`PETER`) = 4.53
  for this pKa set exactly.

Unknown residues (`X`) are handled explicitly, never averaged away:
they are counted and reported, contribute no mass and no charge, and
are excluded from both the GRAVY sum and its count (so a reported MW
for an `X`-containing sequence is the mass of its known part, and the
report says so). A record consisting only of `X` has no defined
properties and is an error.

### Annotation reports

The `annotate` subcommand runs the pipeline end to end over one
file: QC the assembly,
find ORFs in every record, compute each ORF protein's properties,
and present the result as a report. Each annotated gene gets a
sequential identifier (`gene_001`, `gene_002`, ...) and one table
row carrying its location and properties. Locations follow the
conventions of the GFF3 specification: coordinates are 1-based and
inclusive on the record, start is always less than end, and the
strand is `+` or `-`. The frame column is the reading frame (1–3)
on the strand the ORF was found on — the ORFfinder convention of
stage 2 — and is deliberately *not* labelled "phase": in GFF3,
phase is the number of bases to skip to reach the next codon (0–2),
a different quantity.

The same in-memory result is rendered twice — Markdown for reading
next to the code and data, and a single self-contained HTML file
(inline CSS, no external resources) for sharing or archiving — so
the two formats cannot disagree. Both end with a methods section
that states the definitions and sources behind the numbers, in the
report itself rather than only in this README.

## Project structure

```
genome-annotation-pipeline/
├── main.py                        # run the CLI without installing
├── pyproject.toml                 # packaging; installs `genome-annotate`
├── requirements.txt               # no dependencies — see the note inside
├── LICENSE                        # MIT
├── src/
│   └── genome_annotation_pipeline/
│       ├── __init__.py
│       ├── fasta.py               # strict multi-record FASTA parser
│       ├── qc.py                  # GC, length, N50/L50 statistics
│       ├── genetics.py            # standard genetic code, translation
│       ├── orf.py                 # six-frame ORF finder
│       ├── protein.py             # MW, composition, GRAVY, pI estimate
│       ├── report.py              # annotation reports (Markdown + HTML)
│       └── cli.py                 # `genome-annotate` CLI (qc, orf, props, annotate)
├── examples/
│   ├── tiny.fasta                 # 4 records, hand-checkable statistics
│   ├── orfs.fasta                 # 2 records with planted ORFs (34 / 30 aa)
│   ├── proteins.fasta             # 3 short proteins, hand-checkable properties
│   └── synthetic_genome.fasta     # 2 records, 3 planted genes, known coordinates
├── tests/
│   └── run_tests.py               # plain-Python test runner (no pytest)
└── docs/
    └── RESEARCH.md                # sources and design notes
```

## Usage

From a source checkout (no installation needed):

```bash
python3 main.py qc examples/tiny.fasta
python3 main.py qc examples/tiny.fasta --json
python3 main.py orf examples/orfs.fasta
python3 main.py orf examples/orfs.fasta --min-aa 5 --json
python3 main.py props examples/proteins.fasta
python3 main.py props examples/proteins.fasta --json
python3 main.py annotate examples/synthetic_genome.fasta --out-prefix synthetic_genome
```

Or install the package and use the console script:

```bash
pip install .
genome-annotate qc examples/tiny.fasta
genome-annotate orf examples/orfs.fasta
genome-annotate props examples/proteins.fasta
genome-annotate annotate examples/synthetic_genome.fasta --out-prefix synthetic_genome
```

Run the tests:

```bash
python3 tests/run_tests.py
# All 48 tests passed.
```

### One genome, end to end

The four subcommands compose into a single workflow over one file —
QC it, find its ORFs, inspect the ORF proteins, then render the full
annotation (this exact flow, run through the public CLI, is what the
end-to-end integration test verifies):

```bash
python3 main.py qc examples/synthetic_genome.fasta
python3 main.py orf examples/synthetic_genome.fasta
# (write the reported proteins to a protein FASTA file, then:)
python3 main.py props orf-proteins.fasta
python3 main.py annotate examples/synthetic_genome.fasta --out-prefix synthetic_genome
```

## Example

`examples/tiny.fasta` contains four short records whose statistics can be
verified by hand (lengths 40, 30, 20, 10; hand-counted GC counts 20, 19,
9, 6):

```fasta
>contig_040 synthetic longest fragment
ACGTACGTAC
ACGTACGTAC
ACGTACGTAC
ACGTACGTAC
>contig_030 synthetic middle fragment
GGGGCCCCAA
ACGTACGTAC
AAAACCCCGG
>contig_020 lowercase and wrapped fragment
acgtacgtac
ttttccccaa
>contig_010 synthetic short fragment
GCGCGCATAT
```

Actual output of `python3 main.py qc examples/tiny.fasta`:

```
Assembly QC
-----------
Records        : 4
Total bases    : 100
GC (%)         : 54.00
Length min/max : 10 / 40
Mean length    : 25.00
N50            : 30
L50            : 2

Per-record
----------
ID          Length    GC   GC (%)  Ambig
----------------------------------------
contig_040      40    20    50.00      0
contig_030      30    19    63.33      0
contig_020      20     9    45.00      0
contig_010      10     6    60.00      0
```

With `--json` the same run prints the statistics as structured JSON
(record array plus assembly summary with `record_count`, `total_bases`,
`gc_count`, `gc_percent`, `min_length`, `max_length`, `mean_length`,
`n50`, `l50`), suitable for piping into other tools.

Malformed input fails cleanly with a non-zero exit code and a precise
message, for example an `X` in a sequence reports the record, line, and
character instead of producing statistics.

### ORF finding

`examples/orfs.fasta` contains two synthetic records with one planted
ORF each: a 34-aa ORF on the plus strand of `demo_plus` and a 30-aa ORF
on the minus strand of `demo_minus` (planted as the reverse complement
of a known cassette, so the reported coordinates exercise the
strand-mirroring code path). Actual output of
`python3 main.py orf examples/orfs.fasta`:

```
Open reading frames
-------------------
ORFs found     : 2
Min length     : 30 aa (excluding the stop codon)

Record       Start     End  Strand  Frame  Length (aa)  Protein
---------------------------------------------------------------
demo_plus       13     117       +      1           34  MAKFPGDEHQRSTVWYCNILMAKFPGDEHQRSTV
demo_minus      11     103       -      1           30  MAKFPGDEHQRSTVWYCNILMAKFPGDEHQ
```

Lowering the threshold shows why it exists: `--min-aa 5` reports 7
ORFs, the 5 extra ones being short chance ORFs (9–13 aa) that random
sequence produces in the other frames. With `--json` the same run
prints one object per ORF with `record`, `start`, `end`, `strand`,
`frame`, `length_nt`, `length_aa`, and `protein`. In the text table a
protein longer than 60 characters is truncated with `...`; its full
length always remains visible in the Length (aa) column.

### Protein properties

`examples/proteins.fasta` contains three short synthetic proteins whose
properties can be checked by hand against the mass and hydropathy
tables (poly-alanine, the tripeptide `ACD`, poly-lysine):

```fasta
>ala8 synthetic poly-alanine demo
AAAA
AAAA
>acd synthetic tripeptide demo
ACD
>lys4 synthetic poly-lysine demo
KKKK
```

Actual output of `python3 main.py props examples/proteins.fasta`:

```
Protein properties
------------------
Records        : 3

ala8 — synthetic poly-alanine demo
  Length (aa) : 8
  MW (Da)     : 586.65
  pI          : 5.57
  GRAVY       : 1.800
  Composition : A 8 (100.0%)

acd — synthetic tripeptide demo
  Length (aa) : 3
  MW (Da)     : 307.32
  pI          : 4.30
  GRAVY       : 0.267
  Composition : A 1 (33.3%), C 1 (33.3%), D 1 (33.3%)

lys4 — synthetic poly-lysine demo
  Length (aa) : 4
  MW (Da)     : 530.71
  pI          : 10.48
  GRAVY       : -3.900
  Composition : K 4 (100.0%)
```

The three records bracket the scales on purpose: GRAVY runs from
hydrophilic (−3.900) to hydrophobic (1.800), and the pI from acidic
(4.30 — the Asp side chain and C-terminus dominate `ACD`) to basic
(10.48 — four Lys side chains). With `--json` the same run prints, per
protein, `length`, `standard_residues`, `unknown_residues`,
`molecular_weight`, `isoelectric_point`, `gravy`, and the full
20-residue `composition` with counts and percentages.

### Annotation reports

`examples/synthetic_genome.fasta` is a two-record synthetic genome
(486 nt in total) that was *constructed*, not sampled, so its genes
are known exactly. Three gene cassettes were planted at computed
positions in filler made only of the bases A and C (the repeating
pattern `ACCA`): no codon in any frame of an A/C-only stretch, on
either strand, can be a start or a stop codon, so the filler is
ORF-free by construction and nothing but the planted genes can be
found. The layout:

| Record | Length | Planted genes |
| --- | --- | --- |
| `synth_chr1` | 318 nt | 34-aa gene at 38–142 (+ strand); 30-aa gene at 185–277 (− strand) |
| `synth_chr2` | 168 nt | 31-aa gene at 26–121 (+ strand) |

Actual output of `python3 main.py annotate examples/synthetic_genome.fasta --out-prefix synthetic_genome`:

```
Annotation report
-----------------
Records        : 2
Genes annotated: 3
Min length     : 30 aa (excluding the stop codon)
Markdown report: synthetic_genome.md
HTML report    : synthetic_genome.html
```

The Markdown report it writes (`synthetic_genome.md`, shown in
full; the HTML file carries the same tables and methods text):

```markdown
# Genome annotation report

Source file: examples/synthetic_genome.fasta
Pipeline: genome-annotation-pipeline 1.0.0
Minimum ORF length: 30 aa (excluding the stop codon)

## Assembly summary

| Metric | Value |
| --- | --- |
| Records | 2 |
| Total bases | 486 |
| GC (%) | 39.09 |
| Length min/max (nt) | 168 / 318 |
| Mean length (nt) | 243.00 |
| N50 (nt) | 318 |
| L50 | 1 |
| Annotated genes | 3 |

## Annotated genes

| Gene | Record | Start | End | Strand | Frame | Length (aa) | MW (Da) | pI | GRAVY |
| --- | --- | ---: | ---: | :---: | ---: | ---: | ---: | ---: | ---: |
| gene_001 | synth_chr1 | 38 | 142 | + | 2 | 34 | 3980.51 | 6.70 | -0.726 |
| gene_002 | synth_chr1 | 185 | 277 | - | 3 | 30 | 3537.00 | 5.98 | -0.763 |
| gene_003 | synth_chr2 | 26 | 121 | + | 2 | 31 | 3622.20 | 9.19 | -0.535 |

## Methods

- Input was parsed with the pipeline's strict FASTA parser (IUPAC nucleotide alphabet, upper-case normalised, U converted to T); malformed records are errors, not guesses.
- Assembly statistics follow the QUAST definitions: GC (%) is the total number of G and C nucleotides divided by the total length, and N50/L50 are computed against the assembly size.
- Genes are open reading frames found by scanning all six reading frames for the first ATG after the previous stop codon through the first in-frame stop codon, using the standard genetic code (NCBI translation table 1). One ORF is reported per stop codon per frame, and ORFs shorter than the report's minimum length are excluded.
- Coordinates are 1-based and inclusive on the original record, always with start < end, together with the strand (+/-) and the reading frame (1-3) of the strand the ORF was found on — the location and strand conventions of the GFF3 specification. The frame is the reading frame in the ORFfinder sense, not the GFF3 phase field.
- Molecular weight is the sum of the average isotopic residue masses (ExPASy amino-acid mass table) plus one water molecule; GRAVY is the mean of the Kyte & Doolittle (1982) hydropathy values; the pI is an estimate from sequence alone, computed with the Bjellqvist pKa set.
- Gene identifiers (gene_001, gene_002, ...) are assigned sequentially in file order of the records and, within a record, by start coordinate, strand, and frame.
```

The recovered table matches the construction exactly — every
coordinate, strand, frame, and protein — and the test suite proves
that two ways: an independent scan (translating each frame and
scanning the protein strings for M…stop runs, a different method
from the finder's codon scan) derives the same three genes from the
fixture, and the report is compared against both. Note `gene_002`
is exactly 30 aa long: `--min-aa 31` removes precisely that gene
and renumbers the rest, which the tests also pin. With `--json` the
console prints the whole report as structured JSON instead of the
text summary (the two report files are still written).

## Design decisions

- **Strict parser over silent repair.** Invalid characters, sequence data
  before the first header, empty files, and records without sequence data
  are all errors with line/record context. A QC tool that guesses hides
  exactly the problems QC exists to find.
- **`U` is normalised to `T`** so every later stage can assume one DNA
  alphabet; the conversion is documented rather than silent.
- **GC% denominator is the full record length**, including ambiguity
  symbols (QUAST convention), with the ambiguity count reported alongside
  so low-complexity or uncertain records stay visible.
- **Standard library only.** Parsing and these statistics need nothing
  else; zero dependencies keep the foundation trivially reproducible.
- **Statistics are pure functions of parsed records** (`qc.py` never
  touches the filesystem), which makes every number unit-testable against
  hand-computed fixtures.
- **The codon table is generated, not transcribed.** Building it from
  NCBI's compact `AAs`/`Base1`/`Base2`/`Base3` strings removes 64
  hand-typed entries as a failure mode; spot-check tests (ATG→M, TGG→W,
  the stop set) pin the construction.
- **One ORF per stop per frame, from the first ATG.** Nested ATGs do
  not spawn extra rows; the reported ORF is the longest one ending at
  each stop, matching how gene callers annotate.
- **Minimum ORF length is in amino acids, excluding the stop** — the
  unit biologists reason in — rather than nucleotides including it.
- **One parser core, two alphabets.** Protein parsing reuses the strict
  FASTA parser with an amino-acid alphabet and the `U`→`T` conversion
  disabled, so error behaviour (record/line/character context) is
  identical across subcommands instead of drifting between two parsers.
- **`X` is data, not an error — and never silently averaged.** Unknown
  residues are counted and reported, excluded from GRAVY's sum *and*
  its denominator, and given no mass or charge; the report states the
  exclusion. Only an all-`X` record, whose properties are genuinely
  undefined, is an error.
- **pI is presented as an estimate, with its pKa set named.** The
  Bjellqvist set was chosen because it is the set behind ExPASy Compute
  pI/Mw and its calibration is published; the bisection searches the
  full [0, 14] range rather than a narrower window, so very acidic or
  basic peptides are not silently clamped.
- **One result, two renderings.** The Markdown and HTML reports are
  both formatted from a single `AnnotationReport` value, and neither
  renderer computes anything — every number comes from the stage 1–3
  modules — so the two formats cannot drift apart, and a report
  cannot disagree with the `qc`/`orf`/`props` subcommands.
- **Gene identifiers are sequential, not meaningful.** `gene_001`,
  `gene_002`, ... are assigned in file order within one run and
  renumber when `--min-aa` changes the gene set. They identify rows
  of one report; no stable locus identity across runs is promised.
- **The HTML report is a single file on purpose.** Inline CSS and
  zero external references mean it renders offline and archives
  byte-for-byte. Every input-derived string (source name, record
  identifiers) is HTML-escaped before it lands in the document.
- **The stage-4 fixture is a constructive proof.** The synthetic
  genome was assembled from parts whose ORF content is derivable by
  hand — A/C-only filler is ORF-free in all six frames — and the
  tests re-derive its gene set with an independent scanning method
  before comparing it with the report, so "the pipeline finds the
  planted genes" is checked, not assumed.

## Limitations and future work

- Gene calling is ORF finding: there are no gene models — no
  promoters, splice sites, or scoring of candidates — so a "gene" in
  a report is a candidate coding region exactly as stage 2 defines
  it, and on real multi-megabase genomes the unfiltered output would
  need downstream triage this pipeline does not do.
- Reports are point-in-time snapshots: re-running `annotate`
  overwrites the same output paths, and there is no diffing between
  runs or tracking of a gene across them (identifiers renumber, see
  above).
- The pI is a sequence-only estimate from one pKa set (Bjellqvist):
  it ignores folded structure, post-translational modifications, and
  buffer effects, and the underlying calibration (pH 4.5–7.3, in urea)
  makes very short peptides and very basic proteins the least reliable
  cases. Two decimals are reported because the arithmetic has them,
  not because the prediction does.
- Molecular weight is the *average isotopic* mass; no monoisotopic
  masses are computed.
- The protein parser accepts only the 20 standard amino acids plus
  `X`. ExPASy Compute pI/Mw additionally accepts `B`/`Z` with averaged
  masses and pKa values of 0; this pipeline instead rejects them at
  parse time, consistent with its strict-parser policy. A stop symbol
  (`*`) at the end of a translation must be stripped before `props`.
- ORF finding uses the standard code only, requires the start codon
  `ATG` (alternative initiation codons such as GTG/TTG are not
  considered), and requires the stop codon to lie inside the record:
  an ORF that runs off the end of a record is not reported.
- No compressed input (`.gz`/`.zip`) and no FASTQ support.
- GC counts canonical `G`/`C` only; the IUPAC code `S` (G or C) is
  reported as an ambiguity symbol, not split between the two.
- N50 is reported, not NG50 (no reference genome size is accepted yet).

Natural next steps, none of them in this release: GFF3 output for
the annotated genes; alternative genetic codes and initiation
codons (GTG/TTG) for the ORF finder; NG50 given a genome size;
monoisotopic masses alongside the average isotopic ones; FASTQ and
gzip-compressed input; and gene models or candidate scoring so real
genomes need less downstream triage.

## Development stages

1. Stage 1 — foundation: FASTA parser + QC statistics.
2. Stage 2 — ORF finding: six-frame ORF finder + translation
   (standard genetic code).
3. Stage 3 — protein properties: molecular weight, amino-acid
   composition, GRAVY, Bjellqvist-pKa pI estimate.
4. Stage 4 — annotation reports: Markdown + HTML reports over a
   synthetic test genome with known embedded genes.
5. **Stage 5 — release 1.0.0 (this release):** end-to-end
   integration tests over the public CLI, final documentation, and
   release polish. The five stages are complete; further work is
   the future work above.

## References

Every definition and constant in this pipeline was checked against
a primary or authoritative source during development — QUAST's
metric definitions, NCBI's genetic code and MeSH ORF definition,
ExPASy's mass/hydropathy tables and Compute pI/Mw documentation,
Kyte & Doolittle (1982), Bjellqvist et al. (1993, 1994), and the
GFF3 column conventions. The full annotated source list, with the
design notes for each stage, is in
[`docs/RESEARCH.md`](docs/RESEARCH.md).

## License

MIT — see `LICENSE`.
