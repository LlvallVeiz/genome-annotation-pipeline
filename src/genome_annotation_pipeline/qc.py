"""Quality-control statistics for parsed FASTA records.

The definitions here follow the QUAST manual: GC (%) is the number of G and
C nucleotides divided by the total length, N50 is the length for which the
collection of all contigs of that length or longer covers at least half of
the assembly, and L50 is the minimal number of contigs that cover half of
the assembly.
"""

from __future__ import annotations

from dataclasses import dataclass

from .fasta import AMBIGUOUS_BASES, FastaRecord


@dataclass(frozen=True)
class RecordStats:
    """QC statistics for a single record."""

    record_id: str
    description: str
    length: int
    gc_count: int
    gc_percent: float
    ambiguous_count: int

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serialisable representation of the statistics."""
        return {
            "id": self.record_id,
            "description": self.description,
            "length": self.length,
            "gc_count": self.gc_count,
            "gc_percent": round(self.gc_percent, 4),
            "ambiguous_count": self.ambiguous_count,
        }


@dataclass(frozen=True)
class AssemblyStats:
    """QC statistics for a whole set of records (an assembly)."""

    record_count: int
    total_bases: int
    gc_count: int
    gc_percent: float
    min_length: int
    max_length: int
    mean_length: float
    n50: int
    l50: int
    records: tuple[RecordStats, ...]

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serialisable representation of the statistics."""
        return {
            "record_count": self.record_count,
            "total_bases": self.total_bases,
            "gc_count": self.gc_count,
            "gc_percent": round(self.gc_percent, 4),
            "min_length": self.min_length,
            "max_length": self.max_length,
            "mean_length": round(self.mean_length, 4),
            "n50": self.n50,
            "l50": self.l50,
            "records": [record.to_dict() for record in self.records],
        }


def gc_count(sequence: str) -> int:
    """Count canonical G and C bases in a normalised sequence."""
    return sequence.count("G") + sequence.count("C")


def gc_percent(sequence: str) -> float:
    """Return the percentage of canonical G/C bases in a sequence."""
    if not sequence:
        return 0.0
    return 100.0 * gc_count(sequence) / len(sequence)


def ambiguous_count(sequence: str) -> int:
    """Count IUPAC ambiguity symbols (anything beyond A, C, G and T)."""
    return sum(1 for base in sequence if base in AMBIGUOUS_BASES)


def record_stats(record: FastaRecord) -> RecordStats:
    """Compute QC statistics for one parsed record."""
    return RecordStats(
        record_id=record.record_id,
        description=record.description,
        length=record.length,
        gc_count=gc_count(record.sequence),
        gc_percent=gc_percent(record.sequence),
        ambiguous_count=ambiguous_count(record.sequence),
    )


def n50_l50(lengths: list[int]) -> tuple[int, int]:
    """Compute the (N50, L50) pair for a collection of contig lengths.

    Lengths are sorted from longest to shortest and accumulated until the
    running total first reaches at least half of the total length. N50 is
    the length of the contig at that point; L50 is how many contigs were
    needed to reach it.
    """
    if not lengths:
        return (0, 0)
    total = sum(lengths)
    half = total / 2
    cumulative = 0
    for rank, length in enumerate(sorted(lengths, reverse=True), start=1):
        cumulative += length
        if cumulative >= half:
            return (length, rank)
    return (0, 0)  # Unreachable: the full sum always reaches the half.


def assembly_stats(records: list[FastaRecord]) -> AssemblyStats:
    """Compute per-record and assembly-wide QC statistics.

    Raises:
        ValueError: If ``records`` is empty.
    """
    if not records:
        raise ValueError("at least one FASTA record is required")
    per_record = tuple(record_stats(record) for record in records)
    lengths = [record.length for record in records]
    total_bases = sum(lengths)
    total_gc = sum(stats.gc_count for stats in per_record)
    n50, l50 = n50_l50(lengths)
    return AssemblyStats(
        record_count=len(records),
        total_bases=total_bases,
        gc_count=total_gc,
        gc_percent=100.0 * total_gc / total_bases,
        min_length=min(lengths),
        max_length=max(lengths),
        mean_length=total_bases / len(records),
        n50=n50,
        l50=l50,
        records=per_record,
    )
