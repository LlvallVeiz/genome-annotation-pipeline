"""Six-frame open reading frame (ORF) finding.

An open reading frame is, following the NCBI MeSH definition, a sequence
of successive nucleotide triplets read as codons that begins with an
initiator codon and ends with a stop codon. This module searches all
six reading frames of a record — three on the given strand, three on
its reverse complement — for ORFs that begin with the start codon ATG
and end at the first in-frame stop codon (TAA, TAG or TGA) of the
standard genetic code.

Reporting rule for nested starts: within one frame, an ORF runs from
the *first* ATG seen since the previous stop codon (or the start of the
frame) to the *first* stop codon that follows it, so exactly one ORF is
reported per stop codon per frame — the longest ORF ending at that
stop. The alternative (reporting every ATG-to-stop pair, which lists
every nested sub-ORF) was rejected: gene callers annotate the longest
ORF for a stop, and the nested variants share its stop and protein
suffix, so nothing is hidden by the choice.
"""

from __future__ import annotations

from dataclasses import dataclass

from .fasta import FastaRecord
from .genetics import START_CODON, STOP_CODONS, reverse_complement, translate

#: Default minimum ORF size, in amino acids excluding the stop codon.
#: Random sequence hits a stop codon roughly every 21 codons, so ORFs of
#: a few amino acids occur by chance everywhere; 30 keeps the report
#: focused while staying short enough for small synthetic examples.
DEFAULT_MIN_AA = 30


@dataclass(frozen=True)
class Orf:
    """One open reading frame, in coordinates of the original record.

    ``start`` and ``end`` are 1-based and inclusive, always with
    ``start < end``, regardless of strand. ``frame`` is the reading
    frame (1, 2 or 3) on the strand the ORF was found on — the offset
    of the ORF's first base within that strand, modulo 3, plus 1, the
    convention used by NCBI's ORFfinder. ``protein`` is the translation
    of the ORF *without* the terminal stop codon.
    """

    record_id: str
    start: int
    end: int
    strand: str
    frame: int
    protein: str

    @property
    def length_nt(self) -> int:
        """Length of the ORF in nucleotides, including the stop codon."""
        return self.end - self.start + 1

    @property
    def length_aa(self) -> int:
        """Length of the protein in amino acids, excluding the stop."""
        return len(self.protein)

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serialisable representation of the ORF."""
        return {
            "record": self.record_id,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "frame": self.frame,
            "length_nt": self.length_nt,
            "length_aa": self.length_aa,
            "protein": self.protein,
        }


def find_orfs(record: FastaRecord, min_aa: int = DEFAULT_MIN_AA) -> list[Orf]:
    """Find ORFs in all six reading frames of one record.

    Args:
        record: The parsed FASTA record to search.
        min_aa: Minimum protein length, in amino acids *excluding* the
            terminal stop codon. ORFs shorter than this are not
            reported.

    Returns:
        The ORFs found, sorted by start coordinate, then strand
        (``+`` before ``-``), then frame.

    Raises:
        ValueError: If ``min_aa`` is smaller than 1.
    """
    if min_aa < 1:
        raise ValueError(f"min_aa must be at least 1, got {min_aa}")
    found: list[Orf] = []
    strands = (
        ("+", record.sequence),
        ("-", reverse_complement(record.sequence)),
    )
    for strand, sequence in strands:
        for frame_offset in range(3):
            found.extend(
                _scan_frame(
                    record.record_id,
                    sequence,
                    record.length,
                    strand,
                    frame_offset,
                    min_aa,
                )
            )
    found.sort(key=lambda orf: (orf.start, orf.strand, orf.frame))
    return found


def find_orfs_in_records(
    records: list[FastaRecord], min_aa: int = DEFAULT_MIN_AA
) -> list[Orf]:
    """Find ORFs in every record of a parsed FASTA file, in file order."""
    found: list[Orf] = []
    for record in records:
        found.extend(find_orfs(record, min_aa=min_aa))
    return found


def _scan_frame(
    record_id: str,
    sequence: str,
    record_length: int,
    strand: str,
    frame_offset: int,
    min_aa: int,
) -> list[Orf]:
    """Scan one reading frame of one strand for ATG-to-stop ORFs.

    ``sequence`` is the strand being searched (the record itself for
    ``+``, its reverse complement for ``-``); coordinates are converted
    back to the original record before the ORF is built.
    """
    orfs: list[Orf] = []
    orf_start: int | None = None  # 0-based offset of the current ATG.
    for pos in range(frame_offset, len(sequence) - 2, 3):
        codon = sequence[pos : pos + 3]
        if codon in STOP_CODONS:
            if orf_start is not None:
                protein = translate(sequence[orf_start:pos])
                if len(protein) >= min_aa:
                    orfs.append(
                        _build_orf(
                            record_id,
                            strand,
                            frame_offset,
                            orf_start,
                            pos + 3,
                            record_length,
                            protein,
                        )
                    )
                orf_start = None
        elif orf_start is None and codon == START_CODON:
            orf_start = pos
    return orfs


def _build_orf(
    record_id: str,
    strand: str,
    frame_offset: int,
    orf_start: int,
    orf_end: int,
    record_length: int,
    protein: str,
) -> Orf:
    """Convert searched-strand offsets to original-record coordinates.

    ``orf_start``/``orf_end`` are 0-based half-open offsets on the
    searched strand, running from the A of the start codon to just past
    the stop codon. On the ``-`` strand they are mirrored back so the
    reported coordinates always satisfy start < end on the original
    record.
    """
    if strand == "+":
        start = orf_start + 1
        end = orf_end
    else:
        start = record_length - orf_end + 1
        end = record_length - orf_start
    return Orf(
        record_id=record_id,
        start=start,
        end=end,
        strand=strand,
        frame=frame_offset + 1,
        protein=protein,
    )
