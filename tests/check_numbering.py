#!/usr/bin/env python3
"""The results/ step directories must be numbered 01..N with no gaps.

Removing a stage is exactly when this breaks, and a gap in the numbering reads
to anyone opening results/ as though a stage had failed or been skipped.
"""
import re
import sys
from pathlib import Path

home = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
cfg = (home / "config" / "config.sh").read_text()
nums = [int(m) for m in re.findall(r"RESULTS_DIR\}/(\d\d)_", cfg)]
if nums != list(range(1, len(nums) + 1)):
    raise SystemExit(f"gap in results/ numbering: {nums}")
if len(nums) != 6:
    raise SystemExit(f"expected 6 step directories, found {len(nums)}: {nums}")
print(f"results/ numbering contiguous: 01..{nums[-1]:02d}")
