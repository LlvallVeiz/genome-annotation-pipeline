#!/usr/bin/env python3
"""Plain-Python test runner for genome-annotation-pipeline (no pytest).

Run from the repository root:

    python3 tests/run_tests.py

Every expected number in this file was computed by hand from the fixture
sequences first, and only then written down as an assertion.
"""

from __future__ import annotations

import io
import json
import subprocess
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from genome_annotation_pipeline import __version__

from genome_annotation_pipeline.cli import format_orf_table, main as cli_main
from genome_annotation_pipeline.fasta import (
    FastaError,
    FastaRecord,
    parse_fasta,
    parse_fasta_text,
)
from genome_annotation_pipeline.genetics import (
    CODON_TABLE,
    START_CODON,
    STOP_CODONS,
    reverse_complement,
    translate,
)
from genome_annotation_pipeline.orf import Orf, find_orfs, find_orfs_in_records
from genome_annotation_pipeline.protein import (
    HYDROPATHY,
    PKA_CTERM,
    PKA_CTERM_DEFAULT,
    PKA_NEGATIVE,
    PKA_NTERM,
    PKA_NTERM_DEFAULT,
    PKA_POSITIVE,
    RESIDUE_MASSES,
    WATER_MASS,
    gravy,
    isoelectric_point,
    molecular_weight,
    parse_protein_fasta_text,
    protein_props,
)
from genome_annotation_pipeline.qc import assembly_stats, n50_l50, record_stats
from genome_annotation_pipeline.report import (
    annotate_records,
    format_html,
    format_markdown,
)

TINY = ROOT / "examples" / "tiny.fasta"
ORFS_EXAMPLE = ROOT / "examples" / "orfs.fasta"
PROTEINS_EXAMPLE = ROOT / "examples" / "proteins.fasta"
SYNTHETIC_GENOME = ROOT / "examples" / "synthetic_genome.fasta"


def expect_fasta_error(text: str, *fragments: str) -> None:
    """Assert that parsing ``text`` raises FastaError mentioning fragments."""
    try:
        parse_fasta_text(text, source="test-input")
    except FastaError as exc:
        message = str(exc)
        for fragment in fragments:
            assert fragment in message, f"{fragment!r} not in error: {message}"
        return
    raise AssertionError(f"expected FastaError for input: {text!r}")


def test_parser_multiline_lowercase_and_description() -> None:
    records = parse_fasta(TINY)
    assert [record.record_id for record in records] == [
        "contig_040",
        "contig_030",
        "contig_020",
        "contig_010",
    ]
    assert records[0].description == "synthetic longest fragment"
    assert records[0].length == 40
    # Lower-case input is normalised; the two wrapped lines are joined.
    assert records[2].sequence == "ACGTACGTACTTTTCCCCAA"
    assert records[2].length == 20


def test_parser_normalises_u_to_t() -> None:
    (record,) = parse_fasta_text(">rna_fragment\nACGU\n")
    assert record.sequence == "ACGT"


def test_parser_invalid_character_is_reported() -> None:
    expect_fasta_error(">bad\nACGTXAC\n", "invalid character 'X'", "'bad'", "line 2")


def test_parser_empty_input_is_an_error() -> None:
    expect_fasta_error("  \n\n", "empty FASTA input")


def test_parser_sequence_before_header_is_an_error() -> None:
    expect_fasta_error("ACGT\n>late\nACGT\n", "before the first FASTA header")


def test_parser_record_without_sequence_is_an_error() -> None:
    expect_fasta_error(">empty_record\n>other\nACGT\n", "'empty_record'", "no sequence data")


def test_gc_counts_match_hand_computation() -> None:
    records = {record.record_id: record for record in parse_fasta(TINY)}
    stats_040 = record_stats(records["contig_040"])  # "ACGT" x 10
    assert stats_040.gc_count == 20
    assert stats_040.gc_percent == 50.0
    stats_030 = record_stats(records["contig_030"])  # 8 + 5 + 6 G/C per line
    assert stats_030.gc_count == 19
    assert abs(stats_030.gc_percent - (100.0 * 19 / 30)) < 1e-9
    stats_020 = record_stats(records["contig_020"])  # 5 + 4 G/C per line
    assert stats_020.gc_count == 9
    assert stats_020.gc_percent == 45.0
    stats_010 = record_stats(records["contig_010"])  # GCGCGC + ATAT
    assert stats_010.gc_count == 6
    assert stats_010.gc_percent == 60.0


def test_ambiguous_bases_are_counted_not_rejected() -> None:
    (record,) = parse_fasta_text(">with_ns\nACGTNN\n")
    stats = record_stats(record)
    assert stats.length == 6
    assert stats.ambiguous_count == 2
    assert abs(stats.gc_percent - (100.0 * 2 / 6)) < 1e-9


def test_n50_l50_known_example() -> None:
    # Lengths 10, 20, 30, 40 -> total 100, half is 50.
    # Sorted descending: 40 (cum 40 < 50), 30 (cum 70 >= 50) -> N50=30, L50=2.
    assert n50_l50([10, 20, 30, 40]) == (30, 2)
    assert n50_l50([7]) == (7, 1)
    assert n50_l50([]) == (0, 0)


def test_assembly_stats_for_tiny_example() -> None:
    stats = assembly_stats(parse_fasta(TINY))
    assert stats.record_count == 4
    assert stats.total_bases == 100
    assert stats.gc_count == 54  # 20 + 19 + 9 + 6, counted by hand above.
    assert stats.gc_percent == 54.0
    assert stats.min_length == 10
    assert stats.max_length == 40
    assert stats.mean_length == 25.0
    assert stats.n50 == 30
    assert stats.l50 == 2


def test_cli_qc_json_output() -> None:
    stdout = io.StringIO()
    with redirect_stdout(stdout):
        exit_code = cli_main(["qc", str(TINY), "--json"])
    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["record_count"] == 4
    assert payload["total_bases"] == 100
    assert payload["n50"] == 30
    assert payload["l50"] == 2
    assert payload["records"][1]["gc_count"] == 19


def test_cli_missing_file_returns_error_code() -> None:
    stderr = io.StringIO()
    with redirect_stderr(stderr):
        exit_code = cli_main(["qc", str(ROOT / "examples" / "does-not-exist.fasta")])
    assert exit_code == 2
    assert "error:" in stderr.getvalue()


def test_codon_table_spot_checks() -> None:
    # Well-known assignments of the standard code, checked by hand
    # against NCBI transl_table=1 — these pin the programmatic build.
    assert len(CODON_TABLE) == 64
    assert CODON_TABLE["ATG"] == "M"  # the start codon
    assert CODON_TABLE["TGG"] == "W"  # the only Trp codon
    assert CODON_TABLE["TTT"] == "F"
    assert CODON_TABLE["CAT"] == "H"
    assert CODON_TABLE["AAA"] == "K"
    assert CODON_TABLE["CCC"] == "P"
    assert CODON_TABLE["GGG"] == "G"
    assert STOP_CODONS == {"TAA", "TAG", "TGA"}
    assert START_CODON == "ATG"


def test_translate_hand_computed_codons() -> None:
    assert translate("ATGGCTTAA") == "MA*"  # ATG=M, GCT=A, TAA=stop
    assert translate("TTTGGGCAT") == "FGH"  # TTT=F, GGG=G, CAT=H


def test_translate_trailing_partial_codon_is_ignored() -> None:
    assert translate("ATGGCTT") == "MA"  # trailing T cannot form a codon
    assert translate("ATGGC") == "M"  # trailing GC likewise
    assert translate("") == ""


def test_translate_ambiguous_codon_becomes_x() -> None:
    assert translate("ATGNNNTAA") == "MX*"  # NNN is unresolvable -> X
    assert translate("ATHGCT") == "XA"  # ATH (H = A/C/T) -> X, GCT -> A


def test_reverse_complement_hand_computed() -> None:
    # ATGGCT: complement TACCGA, reversed AGCCAT.
    assert reverse_complement("ATGGCT") == "AGCCAT"
    # Ambiguity codes complement pairwise (R<->Y): ACGTRY -> RYACGT.
    assert reverse_complement("ACGTRY") == "RYACGT"


def test_orf_finder_plus_strand_nested_start() -> None:
    # Hand-built fixture: GGGG + ATG GCT ATG AAA TTT CCC GGG TAA + CCCC.
    # The ORF starts at the first ATG (position 5, frame 2), the second
    # ATG (position 11) is nested inside it, and the stop TAA sits at
    # 26-28. An independent scan (translating each frame and scanning
    # the protein string) confirms this is the only ORF at min_aa=1.
    (record,) = parse_fasta_text(">plus\nGGGGATGGCTATGAAATTTCCCGGGTAACCCC\n")
    orfs = find_orfs(record, min_aa=1)
    assert len(orfs) == 1
    orf = orfs[0]
    assert (orf.start, orf.end, orf.strand, orf.frame) == (5, 28, "+", 2)
    assert orf.protein == "MAMKFPG"  # includes the nested ATG's M
    assert orf.length_nt == 24
    assert orf.length_aa == 7


def test_orf_finder_minus_strand_coordinates() -> None:
    # The cassette ATG GCT TTT GGG TAA (protein MAFG) was reverse-
    # complemented and embedded as AAAA + TTACCCAAAAGCCAT + GGGG, so
    # the ORF lives on the minus strand. On the reverse complement it
    # occupies 5-19; mirrored back to the original 23-nt record the
    # coordinates are again 5-19 (start < end).
    (record,) = parse_fasta_text(">minus\nAAAATTACCCAAAAGCCATGGGG\n")
    orfs = find_orfs(record, min_aa=4)
    assert len(orfs) == 1
    orf = orfs[0]
    assert (orf.start, orf.end, orf.strand, orf.frame) == (5, 19, "-", 2)
    assert orf.protein == "MAFG"
    assert orf.length_nt == 15


def test_orf_finder_no_orf_cases() -> None:
    (no_start,) = parse_fasta_text(">no_start\nACGTACGTACGT\n")
    assert find_orfs(no_start, min_aa=1) == []  # no ATG anywhere
    (no_stop,) = parse_fasta_text(">no_stop\nATGAAAGGGCCC\n")
    assert find_orfs(no_stop, min_aa=1) == []  # ATG but no in-frame stop
    (plus,) = parse_fasta_text(">plus\nGGGGATGGCTATGAAATTTCCCGGGTAACCCC\n")
    assert find_orfs(plus, min_aa=8) == []  # the 7-aa ORF is filtered out


def test_cli_orf_json_on_example() -> None:
    stdout = io.StringIO()
    with redirect_stdout(stdout):
        exit_code = cli_main(["orf", str(ORFS_EXAMPLE), "--json"])
    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["record_count"] == 2
    assert payload["orf_count"] == 2
    by_record = {orf["record"]: orf for orf in payload["orfs"]}
    plus = by_record["demo_plus"]
    assert (plus["start"], plus["end"], plus["strand"], plus["frame"]) == (
        13,
        117,
        "+",
        1,
    )
    assert plus["length_nt"] == 105
    assert plus["length_aa"] == 34
    assert plus["protein"] == "MAKFPGDEHQRSTVWYCNILMAKFPGDEHQRSTV"
    minus = by_record["demo_minus"]
    assert (minus["start"], minus["end"], minus["strand"], minus["frame"]) == (
        11,
        103,
        "-",
        1,
    )
    assert minus["length_nt"] == 93
    assert minus["length_aa"] == 30
    assert minus["protein"] == "MAKFPGDEHQRSTVWYCNILMAKFPGDEHQ"


def test_orf_table_truncates_long_protein() -> None:
    orf = Orf(
        record_id="long",
        start=1,
        end=213,
        strand="+",
        frame=1,
        protein="M" + "A" * 69,  # 70 aa
    )
    table = format_orf_table([orf], min_aa=30)
    assert "70" in table  # the full length stays visible
    assert ("M" + "A" * 56 + "...") in table  # 57 chars + ellipsis = 60
    assert orf.protein not in table  # the untruncated protein is not shown


def test_cli_orf_no_orfs_found_message() -> None:
    # tiny.fasta records are at most 40 nt, so nothing can reach the
    # default minimum of 30 aa; the run succeeds and says so.
    stdout = io.StringIO()
    with redirect_stdout(stdout):
        exit_code = cli_main(["orf", str(TINY)])
    assert exit_code == 0
    assert "No ORFs found." in stdout.getvalue()


def test_cli_orf_min_aa_must_be_positive() -> None:
    stderr = io.StringIO()
    with redirect_stderr(stderr):
        exit_code = cli_main(["orf", str(TINY), "--min-aa", "0"])
    assert exit_code == 2
    assert "error:" in stderr.getvalue()


def expect_protein_error(text: str, *fragments: str) -> None:
    """Assert that protein-parsing ``text`` raises FastaError with fragments."""
    try:
        parse_protein_fasta_text(text, source="test-input")
    except FastaError as exc:
        message = str(exc)
        for fragment in fragments:
            assert fragment in message, f"{fragment!r} not in error: {message}"
        return
    raise AssertionError(f"expected FastaError for input: {text!r}")


def test_protein_parser_accepts_amino_acids_and_x() -> None:
    # E, F and X are impossible in the nucleotide alphabet; the protein
    # parser accepts them, normalises case, and joins wrapped lines.
    (record,) = parse_protein_fasta_text(">enzyme demo\nacdx\nEF\n")
    assert record.sequence == "ACDXEF"
    assert record.description == "demo"
    # The same text is invalid as nucleotide input (X is not IUPAC).
    expect_fasta_error(">enzyme\nACDXEF\n", "invalid character 'X'")


def test_protein_parser_rejects_non_standard_letters() -> None:
    # B (Asx) and Z (Glx) are ambiguity codes, U is selenocysteine, and
    # * is a stop symbol: none of them is a standard residue or X.
    expect_protein_error(">p\nACDB\n", "invalid character 'B'", "protein alphabet")
    expect_protein_error(">p\nACDZ\n", "invalid character 'Z'")
    expect_protein_error(">p\nACDU\n", "invalid character 'U'")
    expect_protein_error(">p\nACD*\n", "invalid character '*'")


def test_property_tables_spot_checks() -> None:
    # Values pinned against the checked sources (ExPASy mass and
    # ProtScale hydropathy tables; Bjellqvist pKa set), so an accidental
    # edit of a constant fails loudly.
    assert len(RESIDUE_MASSES) == 20
    assert RESIDUE_MASSES["A"] == 71.0788
    assert RESIDUE_MASSES["G"] == 57.0519
    assert RESIDUE_MASSES["W"] == 186.2132
    assert RESIDUE_MASSES["I"] == RESIDUE_MASSES["L"] == 113.1594
    assert WATER_MASS == 18.01524
    assert len(HYDROPATHY) == 20
    assert HYDROPATHY["I"] == 4.5  # most hydrophobic on the KD scale
    assert HYDROPATHY["R"] == -4.5  # most hydrophilic
    assert HYDROPATHY["A"] == 1.8
    assert PKA_POSITIVE == {"K": 10.0, "R": 12.0, "H": 5.98}
    assert PKA_NEGATIVE == {"D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
    assert PKA_NTERM["P"] == 8.36 and PKA_NTERM["A"] == 7.59
    assert PKA_CTERM == {"D": 4.55, "E": 4.75}
    assert PKA_NTERM_DEFAULT == 7.5 and PKA_CTERM_DEFAULT == 3.55


def test_molecular_weight_hand_computed() -> None:
    # ACD from the ExPASy residue masses plus one water:
    # 71.0788 + 103.1388 + 115.0886 + 18.01524 = 307.32144.
    expected = 71.0788 + 103.1388 + 115.0886 + 18.01524
    assert abs(molecular_weight("ACD") - expected) < 1e-9
    assert round(molecular_weight("ACD"), 2) == 307.32
    # Independent route: free amino acids (residue + water each) minus
    # one water per peptide bond (two bonds in a tripeptide).
    free = (71.0788 + 18.01524) + (103.1388 + 18.01524) + (115.0886 + 18.01524)
    assert abs(molecular_weight("ACD") - (free - 2 * 18.01524)) < 1e-9
    # A single residue is just the free amino acid: Ala = 89.09 Da.
    assert abs(molecular_weight("A") - 89.09404) < 1e-9
    # The example protein: 8 x 71.0788 + 18.01524 = 586.64564.
    assert abs(molecular_weight("AAAAAAAA") - 586.64564) < 1e-9


def test_gravy_known_values() -> None:
    assert gravy("AAAAAAAA") == 1.8  # poly-Ala: the Ala value itself
    assert gravy("III") == 4.5  # poly-Ile: the top of the KD scale
    assert gravy("KKKK") == -3.9
    # (1.8 + 2.5 - 3.5) / 3 = 0.2666...
    assert abs(gravy("ACD") - (1.8 + 2.5 - 3.5) / 3) < 1e-12


def test_composition_hand_computed() -> None:
    (record,) = parse_protein_fasta_text(">mix\nAACD\n")
    props = protein_props(record)
    assert props.length == 4
    assert props.standard_count == 4
    assert props.unknown_count == 0
    composition = {aa: (count, pct) for aa, count, pct in props.composition}
    assert len(composition) == 20  # all residues reported, zeros included
    assert composition["A"] == (2, 50.0)
    assert composition["C"] == (1, 25.0)
    assert composition["D"] == (1, 25.0)
    assert composition["W"] == (0, 0.0)
    assert sum(pct for _, pct in composition.values()) == 100.0


def test_unknown_residue_x_policy() -> None:
    (record,) = parse_protein_fasta_text(">with_x\nAXA\n")
    props = protein_props(record)
    assert props.length == 3
    assert props.standard_count == 2
    assert props.unknown_count == 1
    # X has no mass: the weight is the two alanines plus water.
    assert abs(props.molecular_weight - (2 * 71.0788 + 18.01524)) < 1e-9
    # X is excluded from GRAVY's numerator and denominator alike.
    assert props.gravy == 1.8
    # Composition percentages are shares of the total length (3), so
    # the Ala share is 2/3 and the missing third is the X.
    composition = {aa: (count, pct) for aa, count, pct in props.composition}
    assert composition["A"][0] == 2
    assert abs(composition["A"][1] - 200.0 / 3) < 1e-9


def test_record_with_only_x_is_an_error() -> None:
    (record,) = parse_protein_fasta_text(">only_x\nXXX\n")
    try:
        protein_props(record)
    except ValueError as exc:
        assert "only_x" in str(exc)
        assert "no standard amino" in str(exc)
    else:
        raise AssertionError("expected ValueError for an all-X record")


def test_pi_matches_published_bjellqvist_examples() -> None:
    # Reference values published in the doctests of Biopython's
    # IsoelectricPoint module, which implements the same Bjellqvist
    # pKa set (its values are rounded to two decimals, hence 0.02).
    assert abs(isoelectric_point("INGAR") - 9.75) < 0.02
    assert abs(isoelectric_point("PETER") - 4.53) < 0.02


def _independent_charge(sequence: str, ph: float) -> float:
    """Net charge recomputed per residue, independent of protein.net_charge.

    Uses the module's pKa *constants* (pinned separately above) but a
    separate aggregation: a plain loop over residues instead of the
    module's count-based sums.
    """
    charge = 0.0
    nterm_pka = PKA_NTERM.get(sequence[0], PKA_NTERM_DEFAULT)
    charge += 1.0 / (1.0 + 10.0 ** (ph - nterm_pka))
    cterm_pka = PKA_CTERM.get(sequence[-1], PKA_CTERM_DEFAULT)
    charge -= 1.0 / (1.0 + 10.0 ** (cterm_pka - ph))
    for residue in sequence:
        if residue in PKA_POSITIVE:
            charge += 1.0 / (1.0 + 10.0 ** (ph - PKA_POSITIVE[residue]))
        if residue in PKA_NEGATIVE:
            charge -= 1.0 / (1.0 + 10.0 ** (PKA_NEGATIVE[residue] - ph))
    return charge


def test_pi_charge_is_zero_at_pi() -> None:
    for sequence in ["ACD", "AAAAAAAA", "KKKK", "DDDDDD", "INGAR", "PETER"]:
        pi = isoelectric_point(sequence)
        assert abs(_independent_charge(sequence, pi)) < 0.01, sequence


def test_pi_acidic_vs_basic_sanity() -> None:
    # Hand estimates from the pKa values: poly-Asp's ionisable groups
    # are all acidic (side chain 4.05, C-term 4.55), so its pI sits
    # near 3.3; poly-Lys carries four Lys side chains (10.0), so its
    # pI sits near 10.5.
    assert isoelectric_point("DDDDDD") < 4.0
    assert isoelectric_point("KKKK") > 10.0
    assert isoelectric_point("DDDDDD") < isoelectric_point("ACD")
    assert isoelectric_point("ACD") < isoelectric_point("KKKK")


def test_cli_props_json_on_example() -> None:
    stdout = io.StringIO()
    with redirect_stdout(stdout):
        exit_code = cli_main(["props", str(PROTEINS_EXAMPLE), "--json"])
    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["record_count"] == 3
    by_id = {entry["id"]: entry for entry in payload["proteins"]}
    ala8 = by_id["ala8"]
    assert ala8["molecular_weight"] == 586.6456
    assert ala8["isoelectric_point"] == 5.57
    assert ala8["gravy"] == 1.8
    assert ala8["composition"]["A"] == {"count": 8, "percent": 100.0}
    acd = by_id["acd"]
    assert acd["molecular_weight"] == 307.3214
    assert acd["gravy"] == 0.2667
    lys4 = by_id["lys4"]
    assert lys4["gravy"] == -3.9
    assert lys4["isoelectric_point"] == 10.4778


def _write_temp_fasta(content: str, directory: str) -> str:
    path = Path(directory) / "input.fasta"
    path.write_text(content, encoding="utf-8")
    return str(path)


def test_cli_props_text_report_and_x_note() -> None:
    with tempfile.TemporaryDirectory() as directory:
        path = _write_temp_fasta(">axa demo\nAXA\n", directory)
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            exit_code = cli_main(["props", path])
    assert exit_code == 0
    report = stdout.getvalue()
    assert "MW (Da)     : 160.17" in report  # 2 x Ala + water
    assert "A 2 (66.7%), X 1 (33.3%)" in report
    assert "unknown (X)" in report  # the exclusion is stated, not silent


def test_cli_props_rejects_rna_and_all_x_inputs() -> None:
    with tempfile.TemporaryDirectory() as directory:
        rna_path = _write_temp_fasta(">rna\nACGU\n", directory)
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            exit_code = cli_main(["props", rna_path])
        assert exit_code == 2
        assert "invalid character 'U'" in stderr.getvalue()
        x_path = _write_temp_fasta(">only_x\nXXX\n", directory)
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            exit_code = cli_main(["props", x_path])
        assert exit_code == 2
        assert "error:" in stderr.getvalue()


# --- Stage 4: the synthetic genome and its planted genes -----------------
#
# examples/synthetic_genome.fasta was *constructed*, not sampled:
#
#   synth_chr1 (318 nt) = filler(37) + cassette A (105 nt)
#                         + filler(42) + revcomp(cassette B) (93 nt)
#                         + filler(41)
#   synth_chr2 (168 nt) = filler(25) + cassette C (96 nt) + filler(47)
#
# The filler is the repeat "ACCA" (bases A and C only). Every codon of
# an A/C-only stretch, in any frame, is one of AAA/AAC/ACA/ACC/CAA/
# CAC/CCA/CCC on the plus strand — never ATG, never a stop — and one
# of GGG/GGT/GTG/GTT/TGG/TGT/TTG/TTT on the minus strand (all three
# stop codons contain an A, and ATG needs one too), so filler alone is
# ORF-free in all six frames; only codons spanning the junctions and
# the cassette interiors can form ORFs, and the independent scan in
# test_annotate_recovers_exactly_the_planted_genes proves none do.
# The planted coordinates below are the arithmetic of the layout:
# gene A starts right after 37 nt of filler (frame (38-1) mod 3 + 1
# = 2) and gene B's cassette occupies 185-277, mirrored onto the minus
# strand (frame 3).

_PLANTED_PROTEIN_A = "MAKFPGDEHQRSTVWYCNILMAKFPGDEHQRSTV"  # 34 aa
_PLANTED_PROTEIN_B = "MAKFPGDEHQRSTVWYCNILMAKFPGDEHQ"  # 30 aa
_PLANTED_PROTEIN_C = "MADEKLMNPQRSTVWYACFGHIKLNPQRSTV"  # 31 aa

#: (gene_id, record, start, end, strand, frame, protein), in report order.
PLANTED_GENES = [
    ("gene_001", "synth_chr1", 38, 142, "+", 2, _PLANTED_PROTEIN_A),
    ("gene_002", "synth_chr1", 185, 277, "-", 3, _PLANTED_PROTEIN_B),
    ("gene_003", "synth_chr2", 26, 121, "+", 2, _PLANTED_PROTEIN_C),
]


def _independent_orf_scan(record: FastaRecord, min_aa: int) -> set:
    """ORFs found by a second method: scan translated frames for M..* runs.

    orf.py scans codon by codon; this scans the *protein string* of
    each translated frame for a run from the first M after a stop to
    the first stop that follows it, then converts the run's position
    back to record coordinates. Agreement between the two methods —
    and with the planted set — is what pins the fixture.
    """
    found = set()
    for strand, sequence in (
        ("+", record.sequence),
        ("-", reverse_complement(record.sequence)),
    ):
        for offset in range(3):
            protein = translate(sequence[offset:])
            run_start = None
            for index, residue in enumerate(protein):
                if residue == "*":
                    if run_start is not None:
                        orf_protein = protein[run_start:index]
                        if len(orf_protein) >= min_aa:
                            nt_start = offset + 3 * run_start
                            nt_end = offset + 3 * (index + 1)  # exclusive
                            if strand == "+":
                                start, end = nt_start + 1, nt_end
                            else:
                                start = record.length - nt_end + 1
                                end = record.length - nt_start
                            found.add(
                                (
                                    record.record_id,
                                    start,
                                    end,
                                    strand,
                                    offset + 1,
                                    orf_protein,
                                )
                            )
                    run_start = None
                elif residue == "M" and run_start is None:
                    run_start = index
    return found


def test_annotate_recovers_exactly_the_planted_genes() -> None:
    records = parse_fasta(SYNTHETIC_GENOME)
    assert [record.length for record in records] == [318, 168]
    planted = {
        (record, start, end, strand, frame, protein)
        for _, record, start, end, strand, frame, protein in PLANTED_GENES
    }
    # The independent protein-string scan finds exactly the planted
    # genes in the fixture — no chance ORF sneaks in anywhere.
    scanned = set()
    for record in records:
        scanned |= _independent_orf_scan(record, min_aa=30)
    assert scanned == planted
    # And the annotation reports exactly those genes, in order, with
    # sequential identifiers.
    report = annotate_records(records, min_aa=30)
    assert report.gene_count == 3
    actual = [
        (
            gene.gene_id,
            gene.record_id,
            gene.start,
            gene.end,
            gene.strand,
            gene.frame,
            gene.protein,
        )
        for gene in report.genes
    ]
    assert actual == PLANTED_GENES


def test_annotation_numbers_match_stage_modules() -> None:
    records = parse_fasta(SYNTHETIC_GENOME)
    report = annotate_records(records, min_aa=30, source="check")
    # Assembly numbers: the report's statistics are the qc module's,
    # and the hand-derivable ones match the fixture's construction
    # (lengths 318 + 168 = 486; GC counted character by character).
    assert report.stats == assembly_stats(records)
    assert report.stats.record_count == 2
    assert report.stats.total_bases == 486
    assert report.stats.min_length == 168
    assert report.stats.max_length == 318
    assert report.stats.mean_length == 243.0
    assert (report.stats.n50, report.stats.l50) == (318, 1)
    hand_gc = sum(
        record.sequence.count("G") + record.sequence.count("C")
        for record in records
    )
    assert report.stats.gc_count == hand_gc == 190
    # Gene rows: same ORFs as the finder, properties recomputed by
    # the stage-3 functions on each gene's protein.
    orfs = find_orfs_in_records(records, min_aa=30)
    assert len(report.genes) == len(orfs) == 3
    for gene, orf in zip(report.genes, orfs):
        assert (gene.record_id, gene.start, gene.end, gene.strand) == (
            orf.record_id,
            orf.start,
            orf.end,
            orf.strand,
        )
        assert gene.protein == orf.protein
        assert gene.props.molecular_weight == molecular_weight(orf.protein)
        assert gene.props.isoelectric_point == isoelectric_point(orf.protein)
        assert gene.props.gravy == gravy(orf.protein)
    assert [gene.gene_id for gene in report.genes] == [
        "gene_001",
        "gene_002",
        "gene_003",
    ]


def _expected_gene_row(gene_row: tuple) -> str:
    """The Markdown table row a planted gene must produce."""
    gene_id, record, start, end, strand, frame, protein = gene_row
    return (
        f"| {gene_id} | {record} | {start} | {end} | {strand} | {frame} | "
        f"{len(protein)} | {molecular_weight(protein):.2f} | "
        f"{isoelectric_point(protein):.2f} | {gravy(protein):.3f} |"
    )


def test_markdown_report_carries_the_same_numbers() -> None:
    records = parse_fasta(SYNTHETIC_GENOME)
    report = annotate_records(records, min_aa=30, source="examples/synthetic_genome.fasta")
    markdown = format_markdown(report)
    assert markdown.startswith("# Genome annotation report")
    assert "Source file: examples/synthetic_genome.fasta" in markdown
    assert "| Records | 2 |" in markdown
    assert "| Total bases | 486 |" in markdown
    assert "| GC (%) | 39.09 |" in markdown
    assert "| N50 (nt) | 318 |" in markdown
    assert "| Annotated genes | 3 |" in markdown
    gene_rows = [
        line for line in markdown.splitlines() if line.startswith("| gene_")
    ]
    assert gene_rows == [_expected_gene_row(row) for row in PLANTED_GENES]
    assert "## Methods" in markdown
    assert "Bjellqvist" in markdown


def test_html_report_is_wellformed_and_carries_the_same_numbers() -> None:
    records = parse_fasta(SYNTHETIC_GENOME)
    report = annotate_records(records, min_aa=30, source="examples/synthetic_genome.fasta")
    document = format_html(report)
    assert document.startswith("<!DOCTYPE html>")
    assert document.rstrip().endswith("</html>")
    # Basic well-formedness: every opened table/row/cell is closed.
    assert document.count("<table") == document.count("</table>") == 2
    assert document.count("<tr") == document.count("</tr>")
    assert document.count("<td") == document.count("</td>")
    assert document.count("<thead>") == document.count("</thead>") == 1
    assert document.count("<tbody>") == document.count("</tbody>") == 1
    # Self-contained: inline style only, no external references.
    assert "<style>" in document
    assert "http://" not in document and "https://" not in document
    assert "src=" not in document and "href=" not in document
    # No unfilled placeholder survived rendering.
    assert "TODO" not in document and "%s" not in document
    assert "None" not in document
    # The same numbers as the modules produce, in the same cells.
    for _, _, _, _, _, _, protein in PLANTED_GENES:
        assert f">{molecular_weight(protein):.2f}<" in document
        assert f">{isoelectric_point(protein):.2f}<" in document
        assert f">{gravy(protein):.3f}<" in document
    for gene_id, record, start, end, *_rest in PLANTED_GENES:
        assert f"<td>{gene_id}</td>" in document
        assert f"<td>{record}</td>" in document
        assert f'<td class="num">{start}</td>' in document
        assert f'<td class="num">{end}</td>' in document
    assert "<td>486</td>" in document  # total bases in the summary table


def test_cli_annotate_writes_both_reports() -> None:
    with tempfile.TemporaryDirectory() as directory:
        prefix = str(Path(directory) / "synth-report")
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            exit_code = cli_main(
                ["annotate", str(SYNTHETIC_GENOME), "--out-prefix", prefix]
            )
        assert exit_code == 0
        summary = stdout.getvalue()
        assert "Genes annotated: 3" in summary
        assert f"Markdown report: {prefix}.md" in summary
        assert f"HTML report    : {prefix}.html" in summary
        markdown = Path(prefix + ".md").read_text(encoding="utf-8")
        document = Path(prefix + ".html").read_text(encoding="utf-8")
    assert "| gene_002 | synth_chr1 | 185 | 277 | - | 3 | 30 |" in markdown
    assert "<td>gene_003</td>" in document


def test_cli_annotate_json_summary() -> None:
    with tempfile.TemporaryDirectory() as directory:
        prefix = str(Path(directory) / "synth-report")
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            exit_code = cli_main(
                ["annotate", str(SYNTHETIC_GENOME), "--out-prefix", prefix, "--json"]
            )
    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["gene_count"] == 3
    assert payload["min_aa"] == 30
    assert payload["assembly"]["total_bases"] == 486
    assert payload["assembly"]["n50"] == 318
    gene_002 = payload["genes"][1]
    assert gene_002["gene_id"] == "gene_002"
    assert (gene_002["start"], gene_002["end"], gene_002["strand"]) == (185, 277, "-")
    assert gene_002["molecular_weight"] == round(
        molecular_weight(_PLANTED_PROTEIN_B), 4
    )
    assert gene_002["protein"] == _PLANTED_PROTEIN_B


def test_cli_annotate_min_aa_filters_boundary_gene() -> None:
    # Gene B is exactly 30 aa, so raising the threshold to 31 must
    # remove exactly it; genes A (34 aa) and C (31 aa) remain and the
    # identifiers renumber.
    with tempfile.TemporaryDirectory() as directory:
        prefix = str(Path(directory) / "filtered")
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            exit_code = cli_main(
                [
                    "annotate",
                    str(SYNTHETIC_GENOME),
                    "--min-aa",
                    "31",
                    "--out-prefix",
                    prefix,
                    "--json",
                ]
            )
    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["gene_count"] == 2
    remaining = {
        (gene["record"], gene["start"], gene["end"]) for gene in payload["genes"]
    }
    assert remaining == {("synth_chr1", 38, 142), ("synth_chr2", 26, 121)}
    assert [gene["gene_id"] for gene in payload["genes"]] == [
        "gene_001",
        "gene_002",
    ]


def test_annotate_with_no_genes_found() -> None:
    # tiny.fasta's records are far too short for a 30-aa ORF: both
    # report formats must still render, with the empty result stated.
    with tempfile.TemporaryDirectory() as directory:
        prefix = str(Path(directory) / "empty-report")
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            exit_code = cli_main(["annotate", str(TINY), "--out-prefix", prefix])
        assert exit_code == 0
        assert "Genes annotated: 0" in stdout.getvalue()
        markdown = Path(prefix + ".md").read_text(encoding="utf-8")
        document = Path(prefix + ".html").read_text(encoding="utf-8")
    assert "| Annotated genes | 0 |" in markdown
    assert "No ORFs met the minimum length of 30 aa" in markdown
    assert "| gene_" not in markdown
    assert "No ORFs met the minimum length of 30 aa" in document
    # No gene table rows (the methods note may still *mention* the
    # gene_001 identifier scheme, so check for table cells, not text).
    assert "<td>gene_001</td>" not in document

# --- Stage 5: release integration ----------------------------------------
#
# The tests above exercise modules and the in-process CLI entry point.
# The two below exercise the pipeline exactly as a new user meets it:
# the public CLI in a subprocess (`python3 main.py ...` from the repo
# root), one genome carried through all four subcommands, with the
# numbers cross-checked between the subcommands' own outputs — no
# module internals involved after the command line.


def test_package_version_is_release_1_0_0() -> None:
    assert __version__ == "1.0.0"
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'version = "1.0.0"' in pyproject
    # The version reaches the user: every annotation report names it.
    report = annotate_records(parse_fasta(SYNTHETIC_GENOME), min_aa=30)
    assert "Pipeline: genome-annotation-pipeline 1.0.0" in format_markdown(report)


def _run_public_cli(*arguments: str) -> subprocess.CompletedProcess:
    """Run ``python3 main.py ...`` from the repo root, as a user would."""
    return subprocess.run(
        [sys.executable, str(ROOT / "main.py"), *arguments],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )


def test_end_to_end_public_cli_workflow() -> None:
    genome = str(SYNTHETIC_GENOME)
    # 1. Assembly QC of the raw genome file.
    qc_run = _run_public_cli("qc", genome, "--json")
    assert qc_run.returncode == 0, qc_run.stderr
    qc_payload = json.loads(qc_run.stdout)
    assert qc_payload["record_count"] == 2
    assert qc_payload["total_bases"] == 486
    assert qc_payload["gc_count"] == 190
    assert (qc_payload["n50"], qc_payload["l50"]) == (318, 1)

    # 2. ORF finding on the same file recovers the planted genes.
    orf_run = _run_public_cli("orf", genome, "--json")
    assert orf_run.returncode == 0, orf_run.stderr
    orf_payload = json.loads(orf_run.stdout)
    assert orf_payload["orf_count"] == 3
    found = {
        (orf["record"], orf["start"], orf["end"], orf["strand"], orf["frame"])
        for orf in orf_payload["orfs"]
    }
    planted = {
        (record, start, end, strand, frame)
        for _, record, start, end, strand, frame, _ in PLANTED_GENES
    }
    assert found == planted

    with tempfile.TemporaryDirectory() as directory:
        # 3. The hand-off a user performs: the ORF proteins, written
        # to a protein FASTA, go through `props`.
        proteins_path = Path(directory) / "orf-proteins.fasta"
        proteins_path.write_text(
            "".join(
                f">orf_{index} from {orf['record']}\n{orf['protein']}\n"
                for index, orf in enumerate(orf_payload["orfs"], start=1)
            ),
            encoding="utf-8",
        )
        props_run = _run_public_cli("props", str(proteins_path), "--json")
        assert props_run.returncode == 0, props_run.stderr
        props_payload = json.loads(props_run.stdout)
        assert props_payload["record_count"] == 3
        props_by_protein = {
            orf["protein"]: entry
            for orf, entry in zip(
                orf_payload["orfs"], props_payload["proteins"]
            )
        }

        # 4. The full annotation, and three-way consistency: the
        # report's assembly block is the qc output, its gene rows are
        # the orf output, and its property numbers are the props
        # output for the same proteins (all rounded identically).
        prefix = str(Path(directory) / "genome-report")
        annotate_run = _run_public_cli(
            "annotate", genome, "--out-prefix", prefix, "--json"
        )
        assert annotate_run.returncode == 0, annotate_run.stderr
        report_payload = json.loads(annotate_run.stdout)
        assert report_payload["assembly"] == qc_payload
        assert report_payload["gene_count"] == 3
        for gene in report_payload["genes"]:
            entry = props_by_protein[gene["protein"]]
            assert gene["molecular_weight"] == entry["molecular_weight"]
            assert gene["isoelectric_point"] == entry["isoelectric_point"]
            assert gene["gravy"] == entry["gravy"]
        gene_keys = {
            (gene["record"], gene["start"], gene["end"])
            for gene in report_payload["genes"]
        }
        assert gene_keys == {
            (record, start, end)
            for _, record, start, end, *_ in PLANTED_GENES
        }

        # 5. Both written reports carry the same verified numbers.
        markdown = Path(prefix + ".md").read_text(encoding="utf-8")
        document = Path(prefix + ".html").read_text(encoding="utf-8")
    assert "Pipeline: genome-annotation-pipeline 1.0.0" in markdown
    assert "| Total bases | 486 |" in markdown
    assert "| GC (%) | 39.09 |" in markdown
    for gene_row in PLANTED_GENES:
        assert _expected_gene_row(gene_row) in markdown
    assert "<td>gene_001</td>" in document
    assert "<td>486</td>" in document


TESTS = [
    test_parser_multiline_lowercase_and_description,
    test_parser_normalises_u_to_t,
    test_parser_invalid_character_is_reported,
    test_parser_empty_input_is_an_error,
    test_parser_sequence_before_header_is_an_error,
    test_parser_record_without_sequence_is_an_error,
    test_gc_counts_match_hand_computation,
    test_ambiguous_bases_are_counted_not_rejected,
    test_n50_l50_known_example,
    test_assembly_stats_for_tiny_example,
    test_cli_qc_json_output,
    test_cli_missing_file_returns_error_code,
    test_codon_table_spot_checks,
    test_translate_hand_computed_codons,
    test_translate_trailing_partial_codon_is_ignored,
    test_translate_ambiguous_codon_becomes_x,
    test_reverse_complement_hand_computed,
    test_orf_finder_plus_strand_nested_start,
    test_orf_finder_minus_strand_coordinates,
    test_orf_finder_no_orf_cases,
    test_cli_orf_json_on_example,
    test_orf_table_truncates_long_protein,
    test_cli_orf_no_orfs_found_message,
    test_cli_orf_min_aa_must_be_positive,
    test_protein_parser_accepts_amino_acids_and_x,
    test_protein_parser_rejects_non_standard_letters,
    test_property_tables_spot_checks,
    test_molecular_weight_hand_computed,
    test_gravy_known_values,
    test_composition_hand_computed,
    test_unknown_residue_x_policy,
    test_record_with_only_x_is_an_error,
    test_pi_matches_published_bjellqvist_examples,
    test_pi_charge_is_zero_at_pi,
    test_pi_acidic_vs_basic_sanity,
    test_cli_props_json_on_example,
    test_cli_props_text_report_and_x_note,
    test_cli_props_rejects_rna_and_all_x_inputs,
    test_annotate_recovers_exactly_the_planted_genes,
    test_annotation_numbers_match_stage_modules,
    test_markdown_report_carries_the_same_numbers,
    test_html_report_is_wellformed_and_carries_the_same_numbers,
    test_cli_annotate_writes_both_reports,
    test_cli_annotate_json_summary,
    test_cli_annotate_min_aa_filters_boundary_gene,
    test_annotate_with_no_genes_found,
    test_package_version_is_release_1_0_0,
    test_end_to_end_public_cli_workflow,
]


def main() -> int:
    """Run every test function and report a plain summary."""
    failures = 0
    for test in TESTS:
        try:
            test()
        except Exception as exc:  # noqa: BLE001 - report and continue.
            failures += 1
            print(f"FAIL {test.__name__}: {exc!r}")
        else:
            print(f"PASS {test.__name__}")
    if failures:
        print(f"{failures} of {len(TESTS)} tests failed.")
        return 1
    print(f"All {len(TESTS)} tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
