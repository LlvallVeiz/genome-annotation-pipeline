"""Annotation report generation (stage 4).

Stage 4 runs the earlier stages as one pipeline over a FASTA file —
assembly QC (stage 1), six-frame ORF finding (stage 2), and protein
properties (stage 3) — and renders the combined result as a report in
two formats that carry identical numbers:

- **Markdown**: a title, an assembly-summary table, a table of the
  annotated genes, and a methods section.
- **HTML**: a single self-contained file (inline CSS, no external
  resources) with the same content, so it can be opened anywhere or
  archived next to the data it describes.

A *gene* here is exactly what stage 2 reports: an ORF from the first
ATG after the previous stop codon to the first in-frame stop codon.
Genes are numbered sequentially (``gene_001``, ``gene_002``, ...) in
file order of the records and, within a record, in the ORF finder's
sort order (start coordinate, then strand, then frame). Gene rows use
the location conventions of the GFF3 specification — coordinates are
1-based and inclusive on the landmark (here: the record), start is
always less than end, and the strand is ``+`` or ``-`` — while the
frame column reports the reading frame (1–3) on the strand the ORF
was found on, the ORFfinder convention; that frame is *not* the GFF3
``phase`` field, which counts bases to skip rather than naming a
frame.

Both renderers are pure functions of an :class:`AnnotationReport`,
and every number in them comes from the stage 1–3 modules, so a
report cannot drift from the numbers the other subcommands print.
"""

from __future__ import annotations

import html
from dataclasses import dataclass

from . import __version__
from .fasta import FastaRecord
from .orf import DEFAULT_MIN_AA, Orf, find_orfs_in_records
from .protein import ProteinProps, protein_props
from .qc import AssemblyStats, assembly_stats

#: Method notes shared verbatim by the Markdown and HTML reports.
#: They describe what the numbers mean and where their definitions
#: come from; both renderers must stay this honest, so the text lives
#: in exactly one place.
METHOD_NOTES = (
    "Input was parsed with the pipeline's strict FASTA parser (IUPAC "
    "nucleotide alphabet, upper-case normalised, U converted to T); "
    "malformed records are errors, not guesses.",
    "Assembly statistics follow the QUAST definitions: GC (%) is the "
    "total number of G and C nucleotides divided by the total length, "
    "and N50/L50 are computed against the assembly size.",
    "Genes are open reading frames found by scanning all six reading "
    "frames for the first ATG after the previous stop codon through "
    "the first in-frame stop codon, using the standard genetic code "
    "(NCBI translation table 1). One ORF is reported per stop codon "
    "per frame, and ORFs shorter than the report's minimum length "
    "are excluded.",
    "Coordinates are 1-based and inclusive on the original record, "
    "always with start < end, together with the strand (+/-) and the "
    "reading frame (1-3) of the strand the ORF was found on — the "
    "location and strand conventions of the GFF3 specification. The "
    "frame is the reading frame in the ORFfinder sense, not the GFF3 "
    "phase field.",
    "Molecular weight is the sum of the average isotopic residue "
    "masses (ExPASy amino-acid mass table) plus one water molecule; "
    "GRAVY is the mean of the Kyte & Doolittle (1982) hydropathy "
    "values; the pI is an estimate from sequence alone, computed with "
    "the Bjellqvist pKa set.",
    "Gene identifiers (gene_001, gene_002, ...) are assigned "
    "sequentially in file order of the records and, within a record, "
    "by start coordinate, strand, and frame.",
)


@dataclass(frozen=True)
class AnnotatedGene:
    """One annotated gene: an ORF plus its protein's properties."""

    gene_id: str
    orf: Orf
    props: ProteinProps

    @property
    def record_id(self) -> str:
        """Identifier of the record the gene was found on."""
        return self.orf.record_id

    @property
    def start(self) -> int:
        """1-based inclusive start coordinate on the record."""
        return self.orf.start

    @property
    def end(self) -> int:
        """1-based inclusive end coordinate on the record."""
        return self.orf.end

    @property
    def strand(self) -> str:
        """Strand the ORF was found on (``+`` or ``-``)."""
        return self.orf.strand

    @property
    def frame(self) -> int:
        """Reading frame (1–3) on the strand the ORF was found on."""
        return self.orf.frame

    @property
    def length_aa(self) -> int:
        """Protein length in amino acids, excluding the stop codon."""
        return self.orf.length_aa

    @property
    def protein(self) -> str:
        """The ORF's translation, without the terminal stop codon."""
        return self.orf.protein

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serialisable representation of the gene."""
        payload = dict(self.orf.to_dict())
        payload["gene_id"] = self.gene_id
        payload["molecular_weight"] = round(self.props.molecular_weight, 4)
        payload["isoelectric_point"] = round(self.props.isoelectric_point, 4)
        payload["gravy"] = round(self.props.gravy, 4)
        payload["unknown_residues"] = self.props.unknown_count
        return payload


@dataclass(frozen=True)
class AnnotationReport:
    """The full result of annotating one FASTA file."""

    source: str
    min_aa: int
    stats: AssemblyStats
    genes: tuple[AnnotatedGene, ...]

    @property
    def gene_count(self) -> int:
        """Number of annotated genes in the report."""
        return len(self.genes)

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serialisable representation of the report."""
        return {
            "source": self.source,
            "min_aa": self.min_aa,
            "assembly": self.stats.to_dict(),
            "gene_count": self.gene_count,
            "genes": [gene.to_dict() for gene in self.genes],
        }


def annotate_records(
    records: list[FastaRecord],
    min_aa: int = DEFAULT_MIN_AA,
    source: str = "",
) -> AnnotationReport:
    """Run the full pipeline over parsed records.

    Computes the assembly statistics, finds ORFs in every record (in
    file order), and computes the stage-3 properties of each ORF's
    protein. The statistics are pure functions of the records, so
    they are computed even when no ORF is found.

    Args:
        records: Parsed FASTA records, as returned by
            :func:`~genome_annotation_pipeline.fasta.parse_fasta`.
        min_aa: Minimum ORF length in amino acids, excluding the
            stop codon.
        source: Name of the input the records came from, shown in
            the rendered reports (usually the file path).

    Raises:
        ValueError: If ``records`` is empty or ``min_aa`` is
            smaller than 1.
    """
    stats = assembly_stats(records)
    orfs = find_orfs_in_records(records, min_aa=min_aa)
    genes = []
    for index, orf in enumerate(orfs, start=1):
        gene_id = f"gene_{index:03d}"
        protein_record = FastaRecord(
            record_id=gene_id,
            description="",
            sequence=orf.protein,
        )
        genes.append(
            AnnotatedGene(
                gene_id=gene_id,
                orf=orf,
                props=protein_props(protein_record),
            )
        )
    return AnnotationReport(
        source=source,
        min_aa=min_aa,
        stats=stats,
        genes=tuple(genes),
    )


def _summary_rows(report: AnnotationReport) -> list[tuple[str, str]]:
    """(metric, value) rows of the assembly summary, shared by renderers."""
    stats = report.stats
    return [
        ("Records", str(stats.record_count)),
        ("Total bases", str(stats.total_bases)),
        ("GC (%)", f"{stats.gc_percent:.2f}"),
        ("Length min/max (nt)", f"{stats.min_length} / {stats.max_length}"),
        ("Mean length (nt)", f"{stats.mean_length:.2f}"),
        ("N50 (nt)", str(stats.n50)),
        ("L50", str(stats.l50)),
        ("Annotated genes", str(report.gene_count)),
    ]


def _meta_lines(report: AnnotationReport) -> list[str]:
    """Report header lines (source, version, threshold), shared."""
    lines = []
    if report.source:
        lines.append(f"Source file: {report.source}")
    lines.append(f"Pipeline: genome-annotation-pipeline {__version__}")
    lines.append(
        f"Minimum ORF length: {report.min_aa} aa "
        "(excluding the stop codon)"
    )
    return lines


def format_markdown(report: AnnotationReport) -> str:
    """Render the report as a Markdown document."""
    lines = ["# Genome annotation report", ""]
    lines.extend(_meta_lines(report))
    lines.extend(["", "## Assembly summary", ""])
    lines.append("| Metric | Value |")
    lines.append("| --- | --- |")
    for metric, value in _summary_rows(report):
        lines.append(f"| {metric} | {value} |")
    lines.extend(["", "## Annotated genes", ""])
    if report.genes:
        lines.append(
            "| Gene | Record | Start | End | Strand | Frame | "
            "Length (aa) | MW (Da) | pI | GRAVY |"
        )
        lines.append(
            "| --- | --- | ---: | ---: | :---: | ---: | ---: | "
            "---: | ---: | ---: |"
        )
        for gene in report.genes:
            lines.append(
                f"| {gene.gene_id} | {gene.record_id} | {gene.start} | "
                f"{gene.end} | {gene.strand} | {gene.frame} | "
                f"{gene.length_aa} | {gene.props.molecular_weight:.2f} | "
                f"{gene.props.isoelectric_point:.2f} | "
                f"{gene.props.gravy:.3f} |"
            )
    else:
        lines.append(
            f"_No ORFs met the minimum length of {report.min_aa} aa, "
            "so no genes are annotated._"
        )
    lines.extend(["", "## Methods", ""])
    lines.extend(f"- {note}" for note in METHOD_NOTES)
    lines.append("")
    return "\n".join(lines)


_HTML_STYLE = """\
body { font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
       margin: 2em auto; max-width: 60em; padding: 0 1em;
       color: #1a1a1a; }
table { border-collapse: collapse; margin: 1em 0; }
th, td { border: 1px solid #bbb; padding: 0.3em 0.7em;
         text-align: left; }
th { background: #f0f0f0; }
td.num, th.num { text-align: right; }
.meta { color: #555; }
"""


def format_html(report: AnnotationReport) -> str:
    """Render the report as one self-contained HTML file.

    The document carries its styling in an inline ``<style>`` block
    and references no external resources, so it renders identically
    offline and can be archived as-is. All values taken from the
    input (source name, record identifiers) are HTML-escaped.
    """
    lines = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        "<title>Genome annotation report</title>",
        "<style>",
        _HTML_STYLE.rstrip(),
        "</style>",
        "</head>",
        "<body>",
        "<h1>Genome annotation report</h1>",
        '<p class="meta">',
    ]
    meta = _meta_lines(report)
    lines.extend(
        html.escape(line) + ("<br>" if index < len(meta) - 1 else "")
        for index, line in enumerate(meta)
    )
    lines.extend(["</p>", "<h2>Assembly summary</h2>", "<table>"])
    lines.append("<tr><th>Metric</th><th>Value</th></tr>")
    for metric, value in _summary_rows(report):
        lines.append(
            f"<tr><td>{html.escape(metric)}</td>"
            f"<td>{html.escape(value)}</td></tr>"
        )
    lines.extend(["</table>", "<h2>Annotated genes</h2>"])
    if report.genes:
        lines.append("<table>")
        lines.append(
            "<thead><tr>"
            "<th>Gene</th><th>Record</th>"
            '<th class="num">Start</th><th class="num">End</th>'
            "<th>Strand</th>"
            '<th class="num">Frame</th><th class="num">Length (aa)</th>'
            '<th class="num">MW (Da)</th><th class="num">pI</th>'
            '<th class="num">GRAVY</th>'
            "</tr></thead>"
        )
        lines.append("<tbody>")
        for gene in report.genes:
            lines.append(
                "<tr>"
                f"<td>{html.escape(gene.gene_id)}</td>"
                f"<td>{html.escape(gene.record_id)}</td>"
                f'<td class="num">{gene.start}</td>'
                f'<td class="num">{gene.end}</td>'
                f"<td>{html.escape(gene.strand)}</td>"
                f'<td class="num">{gene.frame}</td>'
                f'<td class="num">{gene.length_aa}</td>'
                f'<td class="num">{gene.props.molecular_weight:.2f}</td>'
                f'<td class="num">{gene.props.isoelectric_point:.2f}</td>'
                f'<td class="num">{gene.props.gravy:.3f}</td>'
                "</tr>"
            )
        lines.append("</tbody>")
        lines.append("</table>")
    else:
        lines.append(
            f"<p>No ORFs met the minimum length of {report.min_aa} aa, "
            "so no genes are annotated.</p>"
        )
    lines.extend(["<h2>Methods</h2>", "<ul>"])
    lines.extend(f"<li>{html.escape(note)}</li>" for note in METHOD_NOTES)
    lines.extend(["</ul>", "</body>", "</html>", ""])
    return "\n".join(lines)
