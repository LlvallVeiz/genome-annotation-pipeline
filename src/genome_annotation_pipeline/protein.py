"""Physicochemical properties of protein sequences.

Stage 3 consumes amino-acid sequences — typically the ORF translations
reported by stage 2 — and computes the properties an annotation report
needs: molecular weight, amino-acid composition, the GRAVY hydropathy
average, and a theoretical isoelectric point (pI) estimate.

Every constant below is transcribed from a checked source (see
``docs/RESEARCH.md`` for the full notes):

- Residue masses and the water mass: the average isotopic values in
  ExPASy's amino-acid mass table (the table ExPASy Compute pI/Mw links
  to for its masses).
- Hydropathy values: Kyte & Doolittle (1982), as tabulated on ExPASy's
  ProtScale page for that scale.
- pKa values: the Bjellqvist set that ExPASy Compute pI/Mw uses,
  including its residue-specific terminal corrections.

Unknown residues (``X``) are handled explicitly, never silently: they
are counted and reported, they contribute nothing to the molecular
weight or the charge, and they are excluded from both the numerator and
the denominator of the GRAVY average. A record with no standard amino
acids at all has no defined properties and is an error.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .fasta import FastaRecord, parse_fasta, parse_fasta_text

#: The 20 standard amino acids, in alphabetical (display) order.
STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

#: Symbol for a residue whose identity is unknown. It is also what the
#: stage-2 translator emits for a codon it cannot resolve.
UNKNOWN_RESIDUE = "X"

#: Alphabet accepted by the protein parser: the 20 standard amino acids
#: plus ``X``. Everything else — the ambiguity codes ``B``/``Z``, the
#: non-standard residues ``U``/``O``/``J``, stop symbols, digits, gaps —
#: is rejected at parse time with record and line context.
PROTEIN_ALPHABET = frozenset(STANDARD_AMINO_ACIDS + UNKNOWN_RESIDUE)

_ALPHABET_DESCRIPTION = (
    "the protein alphabet (the 20 standard amino acids plus X "
    "for an unknown residue)"
)

#: Average isotopic masses of the amino-acid *residues*, in daltons
#: (ExPASy amino-acid mass table). A residue mass is what an amino acid
#: contributes inside a chain: the free amino acid minus the water lost
#: when its peptide bonds formed.
RESIDUE_MASSES = {
    "A": 71.0788,
    "C": 103.1388,
    "D": 115.0886,
    "E": 129.1155,
    "F": 147.1766,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "K": 128.1741,
    "L": 113.1594,
    "M": 131.1926,
    "N": 114.1038,
    "P": 97.1167,
    "Q": 128.1307,
    "R": 156.1875,
    "S": 87.0782,
    "T": 101.1051,
    "V": 99.1326,
    "W": 186.2132,
    "Y": 163.1760,
}

#: Average isotopic mass of one water molecule, in daltons (same ExPASy
#: table). A chain of residue masses is short exactly one water: the
#: free N-terminal H and C-terminal OH that no peptide bond consumed.
WATER_MASS = 18.01524

#: Hydropathy values of Kyte & Doolittle (1982): positive is
#: hydrophobic, negative is hydrophilic. Ile (4.5) is the most
#: hydrophobic residue on the scale, Arg (-4.5) the most hydrophilic.
HYDROPATHY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

#: Bjellqvist pKa values for the basic (positively charging) groups:
#: the side chains of Lys, Arg and His. The N-terminus is handled
#: separately because Bjellqvist gives it residue-specific values.
PKA_POSITIVE = {"K": 10.0, "R": 12.0, "H": 5.98}

#: Bjellqvist pKa values for the acidic (negatively charging) groups:
#: the side chains of Asp, Glu, Cys and Tyr.
PKA_NEGATIVE = {"D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}

#: Generic Bjellqvist pKa of the N-terminal amino group, used when the
#: first residue has no entry in :data:`PKA_NTERM`.
PKA_NTERM_DEFAULT = 7.5

#: Residue-specific Bjellqvist pKa values for the N-terminus, keyed by
#: the first residue of the chain.
PKA_NTERM = {
    "A": 7.59,
    "E": 7.70,
    "M": 7.00,
    "P": 8.36,
    "S": 6.93,
    "T": 6.82,
    "V": 7.44,
}

#: Generic Bjellqvist pKa of the C-terminal carboxyl group, used when
#: the last residue has no entry in :data:`PKA_CTERM`.
PKA_CTERM_DEFAULT = 3.55

#: Residue-specific Bjellqvist pKa values for the C-terminus, keyed by
#: the last residue of the chain (only Asp and Glu differ).
PKA_CTERM = {"D": 4.55, "E": 4.75}


def parse_protein_fasta(path: str | Path) -> list[FastaRecord]:
    """Parse a protein FASTA file from disk.

    Uses the shared parser core with the protein alphabet: the 20
    standard amino acids plus ``X``. Anything else raises
    :class:`~genome_annotation_pipeline.fasta.FastaError` with the
    record, line, and character — so a nucleotide file fed to ``props``
    fails on its first ambiguity code instead of being misread.
    """
    return parse_fasta(
        path,
        alphabet=PROTEIN_ALPHABET,
        alphabet_description=_ALPHABET_DESCRIPTION,
        convert_u=False,
    )


def parse_protein_fasta_text(text: str, source: str = "<string>") -> list[FastaRecord]:
    """Parse protein FASTA content already loaded as text."""
    return parse_fasta_text(
        text,
        source=source,
        alphabet=PROTEIN_ALPHABET,
        alphabet_description=_ALPHABET_DESCRIPTION,
        convert_u=False,
    )


def amino_acid_counts(sequence: str) -> dict[str, int]:
    """Count the 20 standard amino acids in a protein sequence.

    The result always has all 20 keys (zero-filled), in
    :data:`STANDARD_AMINO_ACIDS` order. ``X`` residues are not counted
    here; see :func:`unknown_residue_count`.
    """
    counts = dict.fromkeys(STANDARD_AMINO_ACIDS, 0)
    for residue in sequence:
        if residue in counts:
            counts[residue] += 1
    return counts


def unknown_residue_count(sequence: str) -> int:
    """Count unknown (``X``) residues in a protein sequence."""
    return sequence.count(UNKNOWN_RESIDUE)


def _validate_sequence(sequence: str) -> None:
    """Reject empty sequences and characters outside the protein alphabet."""
    if not sequence:
        raise ValueError("protein sequence is empty")
    invalid = sorted(set(sequence) - PROTEIN_ALPHABET)
    if invalid:
        raise ValueError(
            f"invalid residue(s) {', '.join(invalid)} in protein "
            "sequence; expected the 20 standard amino acids plus X"
        )


def _require_standard_residues(sequence: str, quantity: str) -> int:
    """Return the number of standard residues, or raise if there are none."""
    standard = len(sequence) - unknown_residue_count(sequence)
    if standard == 0:
        raise ValueError(
            f"{quantity} is undefined for a sequence with no standard "
            "amino acids (only X)"
        )
    return standard


def molecular_weight(sequence: str) -> float:
    """Return the average isotopic molecular weight, in daltons.

    The weight is the sum of the residue masses plus one water
    molecule (the ExPASy Compute pI/Mw definition): each residue mass
    already has the peptide-bond water subtracted, so a chain of *n*
    residues is short exactly one water overall. ``X`` residues have no
    known mass and contribute nothing, so for a sequence containing
    ``X`` the result is the mass of its known part only.

    Raises:
        ValueError: If the sequence is empty, contains characters
            outside the protein alphabet, or has no standard amino
            acids.
    """
    _validate_sequence(sequence)
    _require_standard_residues(sequence, "molecular weight")
    return (
        sum(RESIDUE_MASSES[residue] for residue in sequence if residue != UNKNOWN_RESIDUE)
        + WATER_MASS
    )


def gravy(sequence: str) -> float:
    """Return the GRAVY (grand average of hydropathy) of a sequence.

    The mean of the Kyte & Doolittle hydropathy values over the
    sequence's standard residues. (ProtParam divides by the number of
    residues; an ``X`` has no hydropathy value, so here both the sum
    and the count cover standard residues only.)

    Raises:
        ValueError: On an empty or invalid sequence, or one with no
            standard amino acids.
    """
    _validate_sequence(sequence)
    standard = _require_standard_residues(sequence, "GRAVY")
    total = sum(HYDROPATHY[residue] for residue in sequence if residue != UNKNOWN_RESIDUE)
    return total / standard


def _fraction_protonated(ph: float, pka: float) -> float:
    """Fraction of a basic group that is protonated (charge +1) at ``ph``."""
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def _fraction_deprotonated(ph: float, pka: float) -> float:
    """Fraction of an acidic group that is deprotonated (charge -1)."""
    return 1.0 / (1.0 + 10.0 ** (pka - ph))


def net_charge(sequence: str, ph: float) -> float:
    """Return the modelled net charge of a protein at a given pH.

    Henderson–Hasselbalch fractional charges are summed over the two
    termini and the ionisable side chains (K, R, H positive; D, E, C,
    Y negative), using the Bjellqvist pKa set. The terminal pKa values
    depend on the terminal residues, per Bjellqvist; an ``X`` at a
    terminus falls back to the generic terminal pKa, and internal
    ``X`` residues contribute no charge.

    Raises:
        ValueError: On an empty sequence or one containing characters
            outside the protein alphabet.
    """
    _validate_sequence(sequence)
    counts = amino_acid_counts(sequence)
    nterm_pka = PKA_NTERM.get(sequence[0], PKA_NTERM_DEFAULT)
    cterm_pka = PKA_CTERM.get(sequence[-1], PKA_CTERM_DEFAULT)
    positive = _fraction_protonated(ph, nterm_pka) + sum(
        counts[residue] * _fraction_protonated(ph, pka)
        for residue, pka in PKA_POSITIVE.items()
    )
    negative = _fraction_deprotonated(ph, cterm_pka) + sum(
        counts[residue] * _fraction_deprotonated(ph, pka)
        for residue, pka in PKA_NEGATIVE.items()
    )
    return positive - negative


def isoelectric_point(sequence: str) -> float:
    """Estimate the isoelectric point (pI) of a protein sequence.

    The pI is the pH at which :func:`net_charge` is zero. Net charge
    decreases strictly monotonically with pH (basic groups only lose
    protons, acidic groups only gain negative charge), and the charge
    is positive at pH 0 and negative at pH 14, so bisection on
    [0, 14] is guaranteed to bracket the root; the search stops once
    the interval is narrower than 0.0001 pH units.

    This is an *estimate* from sequence alone: the Bjellqvist pKa
    values were calibrated on polypeptides focusing in immobilised pH
    gradients (pH 4.5–7.3, in urea), and no pKa set can see a folded
    protein's local environments, so treat the second decimal as
    precision of the arithmetic, not of the prediction.

    Raises:
        ValueError: On an empty sequence or one containing characters
            outside the protein alphabet.
    """
    _validate_sequence(sequence)
    low, high = 0.0, 14.0
    while high - low > 0.0001:
        mid = (low + high) / 2.0
        if net_charge(sequence, mid) > 0.0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


@dataclass(frozen=True)
class ProteinProps:
    """Computed properties of one protein record."""

    record_id: str
    description: str
    length: int
    standard_count: int
    unknown_count: int
    molecular_weight: float
    isoelectric_point: float
    gravy: float
    #: (amino acid, count, percent of total length) for all 20 standard
    #: amino acids, in :data:`STANDARD_AMINO_ACIDS` order.
    composition: tuple[tuple[str, int, float], ...]

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serialisable representation of the properties."""
        return {
            "id": self.record_id,
            "description": self.description,
            "length": self.length,
            "standard_residues": self.standard_count,
            "unknown_residues": self.unknown_count,
            "molecular_weight": round(self.molecular_weight, 4),
            "isoelectric_point": round(self.isoelectric_point, 4),
            "gravy": round(self.gravy, 4),
            "composition": {
                residue: {"count": count, "percent": round(percent, 4)}
                for residue, count, percent in self.composition
            },
        }


def protein_props(record: FastaRecord) -> ProteinProps:
    """Compute the stage-3 properties of one parsed protein record.

    Composition percentages are shares of the *total* length,
    including any ``X`` residues (so with ``X`` present the standard
    percentages sum to less than 100, and the shortfall is the ``X``
    share, reported via ``unknown_count``).

    Raises:
        ValueError: If the record contains no standard amino acids.
    """
    sequence = record.sequence
    counts = amino_acid_counts(sequence)
    standard = sum(counts.values())
    unknown = unknown_residue_count(sequence)
    if standard == 0:
        raise ValueError(
            f"record '{record.record_id}' contains no standard amino "
            "acids (only X), so its properties are undefined"
        )
    composition = tuple(
        (residue, counts[residue], 100.0 * counts[residue] / len(sequence))
        for residue in STANDARD_AMINO_ACIDS
    )
    return ProteinProps(
        record_id=record.record_id,
        description=record.description,
        length=len(sequence),
        standard_count=standard,
        unknown_count=unknown,
        molecular_weight=molecular_weight(sequence),
        isoelectric_point=isoelectric_point(sequence),
        gravy=gravy(sequence),
        composition=composition,
    )


def protein_props_for_records(records: list[FastaRecord]) -> list[ProteinProps]:
    """Compute properties for every record, in file order."""
    return [protein_props(record) for record in records]
