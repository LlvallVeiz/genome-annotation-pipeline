#!/usr/bin/env python3
"""Run the genome-annotate CLI from a source checkout, without installing.

Usage:
    python3 main.py qc examples/tiny.fasta
    python3 main.py qc examples/tiny.fasta --json
    python3 main.py orf examples/orfs.fasta --min-aa 5
    python3 main.py props examples/proteins.fasta
    python3 main.py annotate examples/synthetic_genome.fasta --out-prefix report
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

from genome_annotation_pipeline.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
