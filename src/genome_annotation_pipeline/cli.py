"""Command-line interface for the genome annotation pipeline.

Four subcommands over one strict FASTA parser:

- ``qc`` — per-record and assembly-level statistics (lengths, GC,
  N50/L50) for a nucleotide FASTA file;
- ``orf`` — six-frame ORF finding over the same parser, with the
  standard genetic code;
- ``props`` — physicochemical properties (molecular weight,
  amino-acid composition, GRAVY, theoretical pI) for protein FASTA
  files;
- ``annotate`` — the whole pipeline over a nucleotide FASTA file,
  writing Markdown and HTML annotation reports.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

from .fasta import FastaError, parse_fasta
from .orf import DEFAULT_MIN_AA, Orf, find_orfs_in_records
from .protein import ProteinProps, parse_protein_fasta, protein_props_for_records
from .qc import AssemblyStats, assembly_stats
from .report import AnnotationReport, annotate_records, format_html, format_markdown

#: Proteins longer than this are truncated in the text table (the full
#: length is always shown in the Length (aa) column).
_PROTEIN_DISPLAY_LIMIT = 60


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="genome-annotate",
        description=(
            "Genome annotation pipeline: FASTA parsing + assembly QC, "
            "six-frame ORF finding, protein properties, and Markdown/HTML "
            "annotation reports. Standard library only."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    qc_parser = subparsers.add_parser(
        "qc",
        help="Report length, GC and N50/L50 statistics for a FASTA file.",
    )
    qc_parser.add_argument("input", help="Path to the input FASTA file.")
    qc_parser.add_argument(
        "--json",
        action="store_true",
        help="Print the statistics as JSON instead of a text table.",
    )

    orf_parser = subparsers.add_parser(
        "orf",
        help="Find open reading frames in all six reading frames.",
    )
    orf_parser.add_argument("input", help="Path to the input FASTA file.")
    orf_parser.add_argument(
        "--min-aa",
        type=int,
        default=DEFAULT_MIN_AA,
        help=(
            "Minimum protein length in amino acids, excluding the stop "
            f"codon (default: {DEFAULT_MIN_AA})."
        ),
    )
    orf_parser.add_argument(
        "--json",
        action="store_true",
        help="Print the ORFs as JSON instead of a text table.",
    )

    props_parser = subparsers.add_parser(
        "props",
        help="Report molecular weight, composition, GRAVY and pI for proteins.",
    )
    props_parser.add_argument("input", help="Path to the input protein FASTA file.")
    props_parser.add_argument(
        "--json",
        action="store_true",
        help="Print the properties as JSON instead of a text report.",
    )

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate a FASTA file and write Markdown + HTML reports.",
    )
    annotate_parser.add_argument("input", help="Path to the input FASTA file.")
    annotate_parser.add_argument(
        "--min-aa",
        type=int,
        default=DEFAULT_MIN_AA,
        help=(
            "Minimum protein length in amino acids, excluding the stop "
            f"codon (default: {DEFAULT_MIN_AA})."
        ),
    )
    annotate_parser.add_argument(
        "--out-prefix",
        required=True,
        help=(
            "Output path prefix: the reports are written to "
            "<prefix>.md and <prefix>.html."
        ),
    )
    annotate_parser.add_argument(
        "--json",
        action="store_true",
        help="Print the report summary as JSON instead of text.",
    )
    return parser


def format_qc_table(stats: AssemblyStats) -> str:
    """Format assembly statistics as an aligned plain-text report."""
    lines = [
        "Assembly QC",
        "-----------",
        f"Records        : {stats.record_count}",
        f"Total bases    : {stats.total_bases}",
        f"GC (%)         : {stats.gc_percent:.2f}",
        f"Length min/max : {stats.min_length} / {stats.max_length}",
        f"Mean length    : {stats.mean_length:.2f}",
        f"N50            : {stats.n50}",
        f"L50            : {stats.l50}",
        "",
        "Per-record",
        "----------",
    ]
    id_width = max(len(record.record_id) for record in stats.records)
    header = (
        f"{'ID':<{id_width}}  {'Length':>6}  {'GC':>4}  "
        f"{'GC (%)':>7}  {'Ambig':>5}"
    )
    lines.append(header)
    lines.append("-" * len(header))
    for record in stats.records:
        lines.append(
            f"{record.record_id:<{id_width}}  {record.length:>6}  "
            f"{record.gc_count:>4}  {record.gc_percent:>7.2f}  "
            f"{record.ambiguous_count:>5}"
        )
    return "\n".join(lines)


def run_qc(input_path: str, as_json: bool) -> int:
    """Run the ``qc`` subcommand and return a process exit code."""
    records = parse_fasta(input_path)
    stats = assembly_stats(records)
    if as_json:
        print(json.dumps(stats.to_dict(), indent=2))
    else:
        print(format_qc_table(stats))
    return 0


def _display_protein(protein: str) -> str:
    """Truncate a protein for the text table, keeping the display tidy."""
    if len(protein) <= _PROTEIN_DISPLAY_LIMIT:
        return protein
    keep = _PROTEIN_DISPLAY_LIMIT - 3
    return protein[:keep] + "..."


def format_orf_table(orfs: list[Orf], min_aa: int) -> str:
    """Format ORF results as an aligned plain-text report."""
    lines = [
        "Open reading frames",
        "-------------------",
        f"ORFs found     : {len(orfs)}",
        f"Min length     : {min_aa} aa (excluding the stop codon)",
        "",
    ]
    if not orfs:
        lines.append("No ORFs found.")
        return "\n".join(lines)
    id_width = max(len(orf.record_id) for orf in orfs)
    header = (
        f"{'Record':<{id_width}}  {'Start':>6}  {'End':>6}  "
        f"{'Strand':>6}  {'Frame':>5}  {'Length (aa)':>11}  Protein"
    )
    lines.append(header)
    lines.append("-" * len(header))
    for orf in orfs:
        lines.append(
            f"{orf.record_id:<{id_width}}  {orf.start:>6}  {orf.end:>6}  "
            f"{orf.strand:>6}  {orf.frame:>5}  {orf.length_aa:>11}  "
            f"{_display_protein(orf.protein)}"
        )
    return "\n".join(lines)


def run_orfs(input_path: str, min_aa: int, as_json: bool) -> int:
    """Run the ``orf`` subcommand and return a process exit code."""
    if min_aa < 1:
        raise ValueError(f"--min-aa must be at least 1, got {min_aa}")
    records = parse_fasta(input_path)
    orfs = find_orfs_in_records(records, min_aa=min_aa)
    if as_json:
        payload = {
            "record_count": len(records),
            "min_aa": min_aa,
            "orf_count": len(orfs),
            "orfs": [orf.to_dict() for orf in orfs],
        }
        print(json.dumps(payload, indent=2))
    else:
        print(format_orf_table(orfs, min_aa))
    return 0


def _compact_composition(props: ProteinProps) -> str:
    """One-line composition: nonzero residues as ``A 8 (100.0%)``."""
    parts = [
        f"{residue} {count} ({percent:.1f}%)"
        for residue, count, percent in props.composition
        if count
    ]
    if props.unknown_count:
        percent = 100.0 * props.unknown_count / props.length
        parts.append(f"X {props.unknown_count} ({percent:.1f}%)")
    return ", ".join(parts)


def format_props_report(props_list: list[ProteinProps]) -> str:
    """Format protein properties as a plain-text report, one block each."""
    lines = [
        "Protein properties",
        "------------------",
        f"Records        : {len(props_list)}",
        "",
    ]
    for index, props in enumerate(props_list):
        if index:
            lines.append("")
        title = props.record_id
        if props.description:
            title += f" — {props.description}"
        lines.append(title)
        lines.append(f"  Length (aa) : {props.length}")
        lines.append(f"  MW (Da)     : {props.molecular_weight:.2f}")
        lines.append(f"  pI          : {props.isoelectric_point:.2f}")
        lines.append(f"  GRAVY       : {props.gravy:.3f}")
        lines.append(f"  Composition : {_compact_composition(props)}")
        if props.unknown_count:
            lines.append(
                f"  Note        : {props.unknown_count} unknown (X) "
                "residue(s) excluded from MW, pI and GRAVY"
            )
    return "\n".join(lines)


def run_props(input_path: str, as_json: bool) -> int:
    """Run the ``props`` subcommand and return a process exit code."""
    records = parse_protein_fasta(input_path)
    props_list = protein_props_for_records(records)
    if as_json:
        payload = {
            "record_count": len(records),
            "proteins": [props.to_dict() for props in props_list],
        }
        print(json.dumps(payload, indent=2))
    else:
        print(format_props_report(props_list))
    return 0


def format_annotate_summary(
    report: AnnotationReport, markdown_path: str, html_path: str
) -> str:
    """Format the ``annotate`` run's outcome as a plain-text summary."""
    lines = [
        "Annotation report",
        "-----------------",
        f"Records        : {report.stats.record_count}",
        f"Genes annotated: {report.gene_count}",
        f"Min length     : {report.min_aa} aa (excluding the stop codon)",
        f"Markdown report: {markdown_path}",
        f"HTML report    : {html_path}",
    ]
    return "\n".join(lines)


def run_annotate(
    input_path: str, min_aa: int, out_prefix: str, as_json: bool
) -> int:
    """Run the ``annotate`` subcommand and return a process exit code.

    Writes ``<out_prefix>.md`` and ``<out_prefix>.html``; the console
    shows a short text summary of the run, or the report itself as
    JSON with ``--json``.
    """
    if min_aa < 1:
        raise ValueError(f"--min-aa must be at least 1, got {min_aa}")
    records = parse_fasta(input_path)
    report = annotate_records(records, min_aa=min_aa, source=input_path)
    markdown_path = f"{out_prefix}.md"
    html_path = f"{out_prefix}.html"
    Path(markdown_path).write_text(format_markdown(report), encoding="utf-8")
    Path(html_path).write_text(format_html(report), encoding="utf-8")
    if as_json:
        print(json.dumps(report.to_dict(), indent=2))
    else:
        print(format_annotate_summary(report, markdown_path, html_path))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point for the ``genome-annotate`` command."""
    args = build_parser().parse_args(argv)
    try:
        if args.command == "qc":
            return run_qc(args.input, args.json)
        if args.command == "orf":
            return run_orfs(args.input, args.min_aa, args.json)
        if args.command == "props":
            return run_props(args.input, args.json)
        if args.command == "annotate":
            return run_annotate(args.input, args.min_aa, args.out_prefix, args.json)
    except (FastaError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 2  # Unreachable while subcommands are required.


if __name__ == "__main__":
    raise SystemExit(main())
