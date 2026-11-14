"""Core package for the genome annotation pipeline.

An end-to-end, dependency-free pipeline over FASTA files:

- strict multi-record FASTA parsing and assembly quality control
  (record counts, lengths, GC content, N50/L50);
- the standard genetic code, translation, and six-frame open reading
  frame (ORF) finding;
- protein properties for the ORF translations (molecular weight,
  amino-acid composition, GRAVY, theoretical pI);
- annotation reports in Markdown and HTML built from the same
  in-memory result, so the two formats cannot disagree.

The command-line interface lives in :mod:`genome_annotation_pipeline.cli`
and exposes the four subcommands ``qc``, ``orf``, ``props`` and
``annotate``.
"""

__version__ = "1.0.0"

__all__ = ["__version__"]
