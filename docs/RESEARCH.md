# Research notes — genome-annotation-pipeline

Stage 1 covers two things: what a FASTA file is, and how the QC numbers
(GC%, N50, L50) are defined. Each definition below was checked against the
listed source during development; the implementation follows the source,
not the other way round.

Stage 2 (sections 6–8) covers the standard genetic code and the ORF
definition the finder implements, checked the same way against NCBI
sources during the stage-2 pass (2026-11-08).

Stage 3 (sections 9–12) covers the protein-property constants —
residue masses, the Kyte & Doolittle hydropathy scale, and the
Bjellqvist pKa set behind the pI estimate — checked against ExPASy
documentation and the primary papers during the stage-3 pass
(2026-11-08). Every numeric table in `protein.py` is transcribed from
one of these sources, not from memory.

Stage 4 (sections 13–14) covers the annotation report's location
conventions, checked against the GFF3 specification's column
definitions (via the GMOD tutorial that reproduces them, and NCBI's
GFF3 format page) during the stage-4 pass (2026-11-08), and the
construction of the synthetic-genome fixture the report is verified
against.

Stage 5 (section 15) is the 1.0.0 integration and polish pass. It
introduces no new science and therefore no new references: every
definition and constant still traces to the sources of stages 1–4
listed below, and the section records only how the assembled
pipeline is verified end to end.

## 1. The FASTA format

A FASTA record consists of a single-line description (the *definition
line*), distinguished from the sequence data by a `>` symbol in the first
column, followed by raw sequence data on the following lines. The first
word of the definition line is the sequence identifier (SeqID); the rest
is an optional description. Sequence data for one record continues until
the next definition line, which is what makes multi-record (multi-FASTA)
files work. NCBI's submission guidance adds that the definition line must
be a single line, the SeqID must not contain spaces, sequence lines
should be 80 characters or shorter, and letter case is ignored.

Consequences baked into `fasta.py`:

- The parser is line-oriented and joins wrapped sequence lines.
- The identifier is the first whitespace-delimited token of the header;
  everything after it is kept as the description.
- Case is normalised to upper case (case carries no meaning).
- A file with sequence data before the first `>`, an empty header, or a
  record with no sequence lines is malformed and is rejected.

## 2. Sequence alphabet: IUPAC nucleotide codes

The accepted alphabet is the IUPAC nucleotide nomenclature: canonical
bases `A`, `C`, `G`, `T`, uracil `U`, and the ambiguity codes
`R` (G or A), `Y` (T or C), `M` (A or C), `K` (G or T), `S` (G or C),
`W` (A or T), `H` (A/C/T), `B` (G/T/C), `V` (G/C/A), `D` (G/A/T) and
`N` (any base), as tabulated by the UCSC Genome Browser's IUPAC page.
`U` is normalised to `T` at parse time so downstream stages see a single
DNA alphabet.

Deliberate exclusions: `X` and gap characters (`.`/`-`) are **not** in
this alphabet, so the parser reports them as invalid characters with the
record id and line number instead of guessing at their meaning.

## 3. GC content

Following the QUAST manual: GC (%) is the total number of G and C
nucleotides divided by the total length, expressed as a percentage.
Applied per record and to the whole assembly (totals over totals, not an
average of percentages — the two differ whenever record lengths differ).

## 4. N50 and L50

Following the QUAST manual's metric definitions:

- **N50** is the length for which the collection of all contigs of that
  length or longer covers at least half of the assembly.
- **L50** is the number of contigs equal to or longer than N50 — in other
  words, the minimal number of contigs that cover half of the assembly.

Algorithm used in `qc.py`: sort lengths descending, accumulate until the
running total first reaches ≥ 50% of the total length; the current
contig's length is N50 and its rank is L50. Hand-verified fixture:
lengths [10, 20, 30, 40] → total 100 → cumulative 40, then 70 → N50 = 30,
L50 = 2 (asserted in `tests/run_tests.py`).

N50 is defined against the assembly size. NG50 substitutes a known or
estimated genome size; stage 1 does not accept a genome size, so only
N50/L50 are reported.

## 5. Initial design decisions

- **Strictness first.** Every malformed-input path raises `FastaError`
  with source, line, and record context, and the CLI maps that to exit
  code 2 with the message on stderr. Silent repair (dropping unknown
  characters, skipping empty records) was rejected: in a QC stage it
  would manufacture false confidence in exactly the files that need
  scrutiny.
- **Ambiguity symbols are data, not errors.** They are valid IUPAC, so
  they parse; they are excluded from the GC numerator, included in the
  length denominator, and counted per record (`ambiguous_count`) so
  their influence on GC% is auditable.
- **No dependencies.** The Python standard library covers line parsing,
  `dataclasses`, and `argparse`; keeping stage 1 dependency-free makes
  the hand-verified fixtures the only source of truth for correctness.

## 6. The standard genetic code (NCBI transl_table=1)

NCBI's Genetic Codes page publishes every code in a compact four-string
form. For transl_table=1 (The Standard Code), read directly from the
page during this pass:

```
AAs    = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M------**--*----M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
```

Column *i* of the four strings jointly defines codon *i*:
`Base1[i] + Base2[i] + Base3[i]` maps to `AAs[i]`, with `*` marking the
stop codons (TAA, TAG, TGA). The page also states that the initiation
codon is AUG (ATG in the DNA alphabet the page uses by GenBank
convention).

`genetics.py` builds `CODON_TABLE` by zipping exactly these strings —
the `AAs` string verbatim, and the Base strings regenerated from their
obvious patterns (16× each base for Base1; `TTTTCCCCAAAAGGGG` ×4 for
Base2; `TCAG` ×16 for Base3) — instead of a 64-entry hand-typed
dictionary. Because a pattern bug would corrupt the table silently, the
tests pin hand-known assignments (ATG→M, TGG→W, TTT→F, CAT→H, AAA→K,
CCC→P, GGG→G, and the stop set {TAA, TAG, TGA}), and the derived
`STOP_CODONS` set is taken from the table itself rather than restated.

Translation policy in `translate()`: stop codons translate to `*` and
remain visible in the output; a trailing partial codon (length not a
multiple of three) is ignored; a codon containing an IUPAC ambiguity
symbol is unresolvable and translates to `X`.

## 7. What counts as an ORF

The finder implements the NCBI MeSH definition of *Open Reading
Frames* (MeSH Unique ID D016366), read directly from the MeSH page:
"A sequence of successive nucleotide triplets that are read as CODONS
specifying AMINO ACIDS and begin with an INITIATOR CODON and end with
a stop codon (CODON, TERMINATOR)." The NCBI Handbook glossary states
the same definition ("Sequence of DNA that begins with an initiation
codon and ends with a termination codon"); its Bookshelf page was
behind a bot check during this pass, so it was read via NCBI's indexed
text and cross-checked against the MeSH wording above.

Tool behaviour follows NCBI's ORFfinder, whose home page describes a
search for ORFs in the entered DNA sequence that returns each ORF's
range and protein translation, with parameters of the same shape this
CLI exposes: a choice of start codons (its examples include genetic
code 1 with 'ATG only') and a minimal ORF length. Concretely, this
stage searches all six reading frames (three per strand, the second
strand searched as the reverse complement), requires ATG as the only
initiator, requires an in-frame stop inside the record, and reports
coordinates on the original sequence.

## 8. Stage-2 design decisions

- **Minimum length in amino acids, excluding the stop** (`--min-aa`,
  default 30). ORFfinder expresses its minimum in nucleotides; amino
  acids were chosen because protein length is the unit the downstream
  stages (protein properties, reports) reason in. The default of 30
  reflects the chance background: with 3 stop codons among 64, a random
  frame hits a stop every ~21 codons on average, so sub-30-aa ORFs are
  mostly noise — visible on demand, e.g. `--min-aa 5` on
  `examples/orfs.fasta` reports 7 ORFs instead of the 2 planted ones.
- **One ORF per stop per frame, from the first ATG.** When several
  ATGs precede one stop, the reported ORF starts at the first of them
  (the longest). Rejected alternative: one ORF per ATG (every nested
  sub-ORF), which multiplies rows without adding information — all
  nested variants share the same stop and differ only by an N-terminal
  truncation of the reported protein.
- **Partial ORFs are not reported.** An ATG with no in-frame stop
  before the record ends yields nothing: the MeSH definition requires
  a terminating stop codon, and a truncated record edge is a data
  problem to surface at QC time, not a gene to invent.
- **Ambiguous codons translate to `X` but cannot start or end an
  ORF.** A codon such as ATH contains no literal ATG, so it never
  opens an ORF; inside an ORF it contributes an `X` residue, so a
  low-quality base does not silently split an otherwise solid ORF.

## 9. Residue masses and molecular weight

ExPASy Compute pI/Mw defines the molecular weight it reports as
obtained "by the addition of average isotopic masses of amino acids in
the protein and the average isotopic mass of one water molecule"
(documentation, section "Molecular Weight"). The masses themselves
come from ExPASy's amino-acid mass table, which lists for each amino
acid the *residue* mass — the free amino acid minus one water, i.e.
what the residue contributes inside a chain — together with the
average isotopic mass of water, 18.01524 Da.

`protein.py` therefore computes

```
MW = sum of residue masses + 18.01524
```

with `RESIDUE_MASSES` transcribed from that table (e.g. A 71.0788,
G 57.0519, W 186.2132; Ile and Leu share 113.1594). The "+ one water"
term is not a fudge: each residue mass already had its peptide-bond
water subtracted, and a chain of *n* residues forms only *n* − 1
peptide bonds, leaving exactly one water unaccounted for. The tests
verify the tripeptide `ACD` two independent ways — directly
(71.0788 + 103.1388 + 115.0886 + 18.01524 = 307.32144 Da) and via free
amino-acid masses minus two waters — and both routes agree.

Note the convention is *average isotopic* mass (the quantity Compute
pI/Mw and ProtParam report); monoisotopic masses are a different table
and are not used anywhere in this stage.

## 10. GRAVY and the Kyte & Doolittle hydropathy scale

The ProtParam documentation defines the property exactly as
implemented: "The GRAVY value for a peptide or protein is calculated
as the sum of hydropathy values of all the amino acids, divided by the
number of residues in the sequence." The values are the hydropathy
scale of Kyte & Doolittle (1982), transcribed from ExPASy ProtScale's
table for that scale (and identical to Biopython's `ProtParamData`
table): Ile 4.5 is the most hydrophobic residue, Arg −4.5 the most
hydrophilic, Ala 1.8. Anchors pinned in the tests: poly-Ala gives
exactly 1.800, poly-Lys exactly −3.900, and `ACD` gives
(1.8 + 2.5 − 3.5)/3 = 0.2667.

The paper: Kyte & Doolittle, *J. Mol. Biol.* 157(1):105–132 (1982),
DOI 10.1016/0022-2836(82)90515-0, PMID 7108955.

## 11. The isoelectric point and the Bjellqvist pKa set

Compute pI/Mw's documentation states that the tool "uses the pK
values of amino acids described in Bjellqvist et al." Bjellqvist and
co-workers calibrated those values against the measured focusing
positions of polypeptides in immobilised pH gradients (1993 paper,
29 polypeptides of known sequence; 1994 paper refining the pH scale).
The same documentation carries the method's own caveats, which this
project repeats rather than hides: the prediction for very basic
proteins may be inadequate, and buffer capacity affects accuracy, so
pI values of small proteins/peptides are the most problematic.

The pKa set as implemented in `protein.py` was read from Biopython's
`Bio/SeqUtils/IsoelectricPoint.py`, whose header attributes the method
to Bjellqvist (1993, 1994); a second, independent calculator preset
for "Bjellqvist" agreed on every value. The set:

| Group | pKa |
|---|---|
| N-terminus (generic) | 7.5 |
| N-terminus when first residue is A / M / S / P / T / V / E | 7.59 / 7.00 / 6.93 / 8.36 / 6.82 / 7.44 / 7.70 |
| Lys (K) | 10.0 |
| Arg (R) | 12.0 |
| His (H) | 5.98 |
| C-terminus (generic) | 3.55 |
| C-terminus when last residue is D / E | 4.55 / 4.75 |
| Asp (D) | 4.05 |
| Glu (E) | 4.45 |
| Cys (C) | 9.0 |
| Tyr (Y) | 10.0 |

Each group contributes a Henderson–Hasselbalch fractional charge —
1/(1 + 10^(pH−pKa)) for a basic group, 1/(1 + 10^(pKa−pH)) for an
acidic one — and the pI is the pH where the fractions sum to zero.
Because every term moves monotonically with pH, the net charge is
strictly decreasing in pH, so bisection is guaranteed to find the
unique crossing. `protein.py` bisects the full [0, 14] interval;
Biopython's implementation bisects a narrower window ([4.05, 12])
matched to its calibration range — the full interval was chosen here
so that very acidic or basic peptides report their modelled pI
instead of a clamped boundary value, at the documented cost of
extrapolating the calibration.

Validation is external, not self-referential: Biopython's module
publishes doctest values pI(`INGAR`) = 9.75 and pI(`PETER`) = 4.53
for this pKa set, and this implementation reproduces both exactly
(asserted in the tests). The tests additionally check, with a
separately written per-residue charge function, that the modelled
charge at each reported pI is ~0, and that an Asp-rich peptide's pI
(3.32 for `DDDDDD`) sits far below a Lys-rich one's (10.48 for
`KKKK`).

## 12. Stage-3 design decisions

- **One parser core, two alphabets.** Protein FASTA parsing calls the
  stage-1 parser with the protein alphabet (the 20 standard amino
  acids plus `X`) and the `U`→`T` conversion disabled, so protein
  input gets the same strictness and the same record/line/character
  error context as nucleotide input. Deliberate contrast with ExPASy
  Compute pI/Mw, which accepts `B`/`Z`/`X` using averaged masses and
  pKa values of 0: averaging an unknown into a number contradicts
  this pipeline's strictness policy, so `B`, `Z`, `U`, `O`, `J`, the
  stop symbol `*`, digits and gaps are rejected at parse time.
- **The `X` policy.** An unknown residue is counted and reported
  (`unknown_count`, and an explicit note in the text report). It
  contributes no mass (MW is then the mass of the known part), no
  side-chain charge, and is excluded from GRAVY's sum *and* count —
  GRAVY is "the average hydropathy of the residues whose hydropathy
  is known". An `X` at a terminus uses the generic terminal pKa. A
  record with no standard residues at all has no defined properties
  and is an error, mirroring the parser's refusal to guess.
- **Composition percentages use the total length** (including `X`)
  as denominator, as ProtParam-style composition does; with `X`
  present the standard percentages therefore sum to less than 100,
  and the shortfall *is* the `X` share — visible rather than
  normalised away.
- **Two decimals for pI, labelled an estimate.** The bisection
  converges to 0.0001 pH units, but the prediction's accuracy is
  bounded by the pKa calibration (immobilised gradients, pH 4.5–7.3,
  in urea), so the README and the module docstring both call the pI
  an estimate and name the pKa set. Reporting more decimals would
  claim precision the method does not have.

## 13. Report location conventions: GFF3

An annotation report is only useful if its coordinates mean the same
thing to every reader, so the gene table adopts the field semantics
of GFF3, the standard feature format, as published by the Sequence
Ontology and reproduced in the GMOD GFF3 tutorial (reference 17) and
NCBI's GFF3 format page (reference 18), both read during this pass:

- Columns 4 and 5 ("start" and "end") are "1-based integer
  coordinates" relative to the column-1 landmark (the *seqid*), and
  "start is always less than or equal to end" — the convention
  `orf.py` already reports (1-based, inclusive, start < end on the
  original record), with the record ID playing the seqid role.
- Column 7 ("strand") is "+ for positive strand (relative to the
  landmark), − for minus strand" — the same two symbols the report
  uses.
- Column 8 ("phase"), required for CDS features, is "one of the
  integers 0, 1, or 2, indicating the number of bases that should be
  removed from the beginning of this feature to reach the first base
  of the next codon". That is a *skip count for multi-line CDS
  features*, not a reading-frame number: for a complete CDS that
  begins exactly at a codon boundary the GFF3 phase is 0, whatever
  frame of the landmark it sits in. The report's Frame column
  instead keeps the stage-2/ORFfinder convention — the reading frame
  (1–3) counted on the strand the ORF was found on — and is labelled
  "Frame", never "phase", so the two quantities cannot be confused.

The reports are not GFF3 files (no GFF3 output is produced; the
Markdown and HTML tables are for reading), but every location in
them can be transcribed into GFF3 fields without re-interpretation,
and the methods section of each report says so.

## 14. Stage-4 design decisions

- **The fixture is designed so its truth is derivable, not just
  asserted.** `examples/synthetic_genome.fasta` (486 nt, 2 records)
  plants three gene cassettes in filler consisting only of the bases
  A and C (the repeat `ACCA`). The filler argument is closed-form:
  any codon drawn from {A, C} is one of AAA/AAC/ACA/ACC/CAA/CAC/CCA/
  CCC — never ATG and never a stop — and on the reverse strand every
  codon is drawn from {G, T}, while ATG and all three stop codons
  (TAA, TAG, TGA) contain an A; so filler alone yields no ORF in any
  of the six frames, regardless of frame alignment. Only
  junction-spanning codons and cassette interiors can produce ORFs,
  and the test suite scans for exactly those with a second,
  independent method (translate each frame, scan the protein string
  for M…stop runs) and requires the result to equal the planted set
  — coordinates, strands, frames, and proteins — before comparing
  the report against it. The cassette layout was adjusted once under
  this check: with the first candidate filler length, a shifted
  frame starting at an internal ATG of cassette A ran 35 codons
  through the downstream cassette before hitting a stop; the
  middle-filler length (42 nt) in the shipped fixture was selected
  because the independent scan proves that run — and every other
  unplanted one — stays below the 30-aa threshold.
- **One in-memory result, two renderers.** `report.py` first builds
  an `AnnotationReport` (assembly statistics from `qc.py`, ORFs from
  `orf.py`, properties from `protein.py`); `format_markdown` and
  `format_html` only format it. The methods text is a single shared
  constant rendered by both, so the formats cannot diverge in
  numbers *or* in wording.
- **Gene identifiers are per-report.** `gene_001`, ... are assigned
  in file order (and the finder's within-record order) and renumber
  when `--min-aa` changes the set — the 30-aa boundary case is
  pinned by a test: at `--min-aa 31` exactly the 30-aa gene drops
  out. Stable cross-run locus identity is explicitly out of scope.
- **HTML is one file, escaped.** Inline CSS, no external references
  (a test asserts no `http(s)://`, `src=`, or `href=` occurs in the
  output), and `html.escape` on every input-derived string, so the
  file is archivable as-is and a hostile record ID cannot inject
  markup into it.

## 15. Stage-5 verification: end-to-end, through the public CLI

Stages 1–4 were each verified at module level and through the
in-process CLI entry point. Stage 5 adds the missing angle: the
whole user workflow in a subprocess, exactly as a new user runs it
(`python3 main.py ...` from the repository root), carrying one file
— `examples/synthetic_genome.fasta` — through all four subcommands
and requiring the outputs to agree *with each other*, not only with
hand-computed fixtures:

1. `qc --json` reports the assembly (2 records, 486 bases, GC count
   190, N50 318, L50 1).
2. `orf --json` recovers exactly the three planted genes (the same
   set section 14 derives independently).
3. The reported proteins, written to a protein FASTA — the hand-off
   a user performs between the `orf` and `props` subcommands — go
   through `props --json`.
4. `annotate --json` must then reproduce all of it: its assembly
   block is byte-identical to the `qc` output, its gene coordinates
   are the `orf` output, and each gene's molecular weight, pI, and
   GRAVY equal the `props` output for the same protein (both sides
   round to four decimals in their `to_dict` forms, so the
   comparison is exact). Both written reports (Markdown and HTML)
   are finally checked to carry the same numbers and the 1.0.0
   version string.

## References

All URLs below were opened or returned by search and read during the
stage-1 research pass (2026-11-07).

1. **QUAST 5.0.2 manual** — metric definitions used for GC (%), N50 and
   L50 (section "Understanding the output / Metrics description").
   Institut Pasteur Gensoft mirror.
   https://gensoft.pasteur.fr/docs/quast/5.0.2/
2. **Using BLAT to Find Sequence Similarity in Closely Related Genomes**
   — *Current Protocols in Bioinformatics* (2012), hosted on PubMed
   Central. States the FASTA sequence format as a header line beginning
   with the `>` symbol and a unique name, followed by the sequence.
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4101998
3. **The GenBank Submissions Handbook — "Formatting your Submission"**
   — NCBI Bookshelf. Definition line construction, SeqID rules, 80-char
   sequence lines, case-insensitivity, IUPAC symbols for nucleotide
   submissions. (Read via NCBI's indexed text; NCBI's direct page view
   was behind a bot check during this pass, and the definition was
   cross-checked against reference 2.)
   https://www.ncbi.nlm.nih.gov/books/NBK566986/
4. **Genome Browser IUPAC Codes** — UCSC Genome Browser help. The
   single-character IUPAC representation of single and ambiguous bases
   used as this parser's accepted alphabet.
   https://genome.soe.ucsc.edu/goldenPath/help/iupac.html

Stage-2 references (opened and read during the stage-2 research pass,
2026-11-08):

5. **The Genetic Codes** — NCBI (compiled by Elzanowski & Ostell;
   page last updated 2024-09-23). Section 1, "The Standard Code
   (transl_table=1)", provided the `AAs`/`Starts`/`Base1`/`Base2`/
   `Base3` strings the codon table is built from, verbatim.
   https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
6. **Open Reading Frames — MeSH** — NCBI, MeSH Unique ID D016366. The
   definition of an ORF quoted in section 7 (initiator codon to stop
   codon), read directly from the MeSH page.
   https://www.ncbi.nlm.nih.gov/mesh?term=Open+Reading+Frames
7. **ORFfinder** — NCBI. Tool description and example parameter sets
   (genetic code, 'ATG only' start codons, minimal ORF length) that
   the `orf` subcommand's behaviour mirrors.
   https://www.ncbi.nlm.nih.gov/orffinder/
8. **The NCBI Handbook — glossary, "open reading frame"** — NCBI
   Bookshelf NBK470040. Same definition as reference 6, cross-checked
   via NCBI's indexed text (the direct page view was behind a bot
   check during this pass, as in stage 1).
   https://www.ncbi.nlm.nih.gov/books/NBK470040/def-item/open-reading-frame/

Stage-3 references (opened and read during the stage-3 research pass,
2026-11-08):

9. **ProtParam documentation** — ExPASy (SIB). Defines GRAVY as the
   sum of hydropathy values divided by the number of residues, and
   states that the molecular weight and pI are computed as in Compute
   pI/Mw (references 10–11).
   https://web.expasy.org/protparam/protparam-doc.html
10. **Amino acid masses** — ExPASy (FindMod mass tables; the table
    Compute pI/Mw links to for its average isotopic masses). Source
    of the 20 residue masses and of the average isotopic water mass
    (18.01524 Da) transcribed into `protein.py`.
    https://web.expasy.org/findmod/findmod_masses.html
11. **Compute pI/Mw documentation** — ExPASy (SIB). Molecular weight
    defined as the sum of average isotopic amino-acid masses plus one
    water; pI stated to use "the pK values of amino acids described
    in Bjellqvist et al."; also the source of the caveats quoted in
    section 11 (very basic proteins, buffer capacity/small proteins).
    https://web.expasy.org/compute_pi/pi_tool-doc.html
12. **Kyte, J. & Doolittle, R. F. (1982).** A simple method for
    displaying the hydropathic character of a protein. *Journal of
    Molecular Biology* 157(1), 105–132.
    DOI: 10.1016/0022-2836(82)90515-0. PMID: 7108955.
13. **ProtScale — hydropathicity scale (Kyte & Doolittle)** — ExPASy
    (SIB). The per-residue hydropathy table transcribed into
    `protein.py` (Ile 4.5 … Arg −4.5).
    https://web.expasy.org/protscale/pscale/Hphob.Doolittle.html
14. **Bjellqvist, B., Hughes, G. J., Pasquali, C., Paquet, N.,
    Ravier, F., Sanchez, J.-C., Frutiger, S. & Hochstrasser, D.
    (1993).** The focusing positions of polypeptides in immobilized
    pH gradients can be predicted from their amino acid sequences.
    *Electrophoresis* 14(10), 1023–1031.
    DOI: 10.1002/elps.11501401163. PMID: 8125050.
15. **Bjellqvist, B., Basse, B., Olsen, E. & Celis, J. E. (1994).**
    Companion paper refining the pH scale so that isoelectric points
    correlate with polypeptide compositions. *Electrophoresis* 15,
    529–539. DOI: 10.1002/elps.1150150171. PMID: 8055880.
16. **IsoelectricPoint module source** — Biopython
    (`Bio/SeqUtils/IsoelectricPoint.py`). The Bjellqvist pKa set as
    implemented was read from this file (its header attributes the
    method to references 14–15); its doctests publish the reference
    values pI(`INGAR`) = 9.75 and pI(`PETER`) = 4.53 that the stage-3
    tests reproduce.
    https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/IsoelectricPoint.py

Stage-4 references (returned by search and read during the stage-4
research pass, 2026-11-08):

17. **GFF3 tutorial** — GMOD. Reproduces the Sequence Ontology GFF3
    specification's column definitions used in section 13: columns
    4–5 as 1-based coordinates with start ≤ end, column 7 strand as
    +/−, and column 8 phase as the 0/1/2 bases-to-skip count.
    https://gmod.org/wiki/GFF3_2011_tutorial
18. **GFF3 format** — NCBI Datasets documentation. NCBI's statement
    of the same column specifications (start/end in 1-based
    coordinates; strand; phase related to /codon_start), used as the
    cross-check for reference 17.
    https://www.ncbi.nlm.nih.gov:443/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/
