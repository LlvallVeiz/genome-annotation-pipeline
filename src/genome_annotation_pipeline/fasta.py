"""FASTA parsing for nucleotide sequences.

The parser is deliberately strict: a downstream annotation pipeline should
fail early with a precise message instead of silently annotating a malformed
record. It accepts the IUPAC nucleotide alphabet, normalises sequences to
upper case, and converts uracil (``U``) to thymine (``T``) so the rest of the
pipeline can assume a DNA alphabet.

Stage 3 reuses this parser core for *protein* records by passing an
amino-acid alphabet and disabling the ``U`` -> ``T`` conversion (see
``protein.py``); the default behaviour described above is unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class FastaError(ValueError):
    """Raised when FASTA input is malformed."""


#: Canonical DNA bases plus the IUPAC ambiguity codes (Cornish-Bowden, 1985).
#: ``U`` is accepted and normalised to ``T``. ``X`` and gap symbols are not
#: part of that nucleotide alphabet, so they are reported as invalid.
VALID_BASES = frozenset("ACGTURYSWKMBDHVN")

#: IUPAC ambiguity codes: every accepted symbol that is not A, C, G, T or U.
AMBIGUOUS_BASES = frozenset("RYSWKMBDHVN")

#: Alphabet description used in invalid-character errors by default.
_DEFAULT_ALPHABET_DESCRIPTION = "the IUPAC nucleotide alphabet (ACGTU plus RYSWKMBDHVN)"


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record: an identifier, an optional description, a sequence."""

    record_id: str
    description: str
    sequence: str

    @property
    def header(self) -> str:
        """Return the full definition line, without the leading ``>``."""
        if self.description:
            return f"{self.record_id} {self.description}"
        return self.record_id

    @property
    def length(self) -> int:
        """Return the number of bases in the sequence."""
        return len(self.sequence)


def parse_fasta(
    path: str | Path,
    *,
    alphabet: frozenset[str] = VALID_BASES,
    alphabet_description: str = _DEFAULT_ALPHABET_DESCRIPTION,
    convert_u: bool = True,
) -> list[FastaRecord]:
    """Parse a FASTA file from disk.

    Args:
        path: Path to a ``.fasta``/``.fa`` text file.
        alphabet: Accepted sequence characters (after upper-casing).
            Defaults to the IUPAC nucleotide alphabet; ``protein.py``
            passes the amino-acid alphabet instead.
        alphabet_description: Human-readable alphabet name used in
            invalid-character error messages.
        convert_u: Normalise ``U`` to ``T`` (nucleotide semantics).
            Protein parsing disables this: in a protein, ``U`` is
            selenocysteine, not a spelling of another residue, and it is
            outside the accepted protein alphabet anyway.

    Returns:
        The records in file order, with sequences normalised to upper case
        and ``U`` converted to ``T``.

    Raises:
        FastaError: If the file is empty or malformed.
        OSError: If the file cannot be read.
    """
    source = str(path)
    try:
        text = Path(path).read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise FastaError(f"{source}: file is not valid UTF-8 text ({exc})") from exc
    return parse_fasta_text(
        text,
        source=source,
        alphabet=alphabet,
        alphabet_description=alphabet_description,
        convert_u=convert_u,
    )


def parse_fasta_text(
    text: str,
    source: str = "<string>",
    *,
    alphabet: frozenset[str] = VALID_BASES,
    alphabet_description: str = _DEFAULT_ALPHABET_DESCRIPTION,
    convert_u: bool = True,
) -> list[FastaRecord]:
    """Parse FASTA content already loaded as text.

    Args:
        text: Raw FASTA text; sequence lines may be wrapped at any width
            and may use lower case.
        source: Name used in error messages (usually the file path).
        alphabet: Accepted sequence characters (after upper-casing).
        alphabet_description: Human-readable alphabet name used in
            invalid-character error messages.
        convert_u: Normalise ``U`` to ``T`` after validation.

    Returns:
        The parsed records in file order.

    Raises:
        FastaError: With the line number and record context when the input
            is empty, has sequence data before the first header, contains
            an invalid character, or has a record without sequence data.
    """
    if not text.strip():
        raise FastaError(f"{source}: empty FASTA input")

    records: list[FastaRecord] = []
    current_id: str | None = None
    current_description = ""
    current_header_line = 0
    chunks: list[str] = []

    def flush() -> None:
        """Append the record being accumulated, validating it has sequence."""
        if current_id is None:
            return
        if not chunks:
            raise FastaError(
                f"{source}: record '{current_id}' (header at line "
                f"{current_header_line}) has no sequence data"
            )
        records.append(
            FastaRecord(
                record_id=current_id,
                description=current_description,
                sequence="".join(chunks),
            )
        )

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            if not header:
                raise FastaError(f"{source}: empty FASTA header at line {line_number}")
            parts = header.split(None, 1)
            current_id = parts[0]
            current_description = parts[1] if len(parts) == 2 else ""
            current_header_line = line_number
            chunks = []
            continue
        if current_id is None:
            raise FastaError(
                f"{source}: sequence data before the first FASTA header "
                f"at line {line_number}"
            )
        chunk = "".join(line.split()).upper()
        for char in chunk:
            if char not in alphabet:
                raise FastaError(
                    f"{source}: invalid character {char!r} in record "
                    f"'{current_id}' at line {line_number}; expected "
                    f"{alphabet_description}"
                )
        chunks.append(chunk.replace("U", "T") if convert_u else chunk)

    flush()
    if not records:
        raise FastaError(f"{source}: no FASTA records found")
    return records
