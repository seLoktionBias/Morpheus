#!/usr/bin/env python3
"""Catch shell constructs that work on one bash and not another.

The pipeline is developed on macOS (bash 3.2) and run on RHEL-family clusters
(bash 4.x). Both directions have bitten:

  ${#arr[@]-0}  bash 3.2 tolerated it; bash 4.x rejects it as a bad
                substitution, and the job died on line 415 with no other clue.
                It is not valid anywhere: the # length operator cannot take a
                -default. Declare the array and use ${#arr[@]}.
  mapfile       bash 4+ only; on macOS a local run dies immediately with
                "mapfile: command not found".

`bash -n` catches neither: the first is a runtime expansion error, the second a
missing builtin. So grep for them.
"""
import re
import sys
from pathlib import Path

BANNED = [
    (re.compile(r"\$\{#\w+\[[@*]\]\s*[-:+?]"),
     "${#arr[@]-default}: the # length operator cannot take a modifier. "
     "Declare the array up front and use ${#arr[@]}."),
    (re.compile(r"(?<![\w-])mapfile\b"),
     "mapfile is bash 4+; macOS ships bash 3.2. Use a while-read loop."),
    (re.compile(r"(?<![\w-])readarray\b"),
     "readarray is bash 4+; macOS ships bash 3.2. Use a while-read loop."),
    (re.compile(r"\$\{\w+\^\^|\$\{\w+,,"),
     "${var^^} / ${var,,} are bash 4+. Use tr."),
    (re.compile(r"(?<![\w-])declare\s+-A\b"),
     "associative arrays are bash 4+."),
]

home = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
targets = sorted(
    p for d in ("scripts", "slurm", "config", "bin", "tests")
    for p in (home / d).rglob("*")
    if p.is_file() and (p.suffix in {".sh", ".sbatch"} or p.name == "morpheus")
)
if not targets:
    raise SystemExit("no shell files found to check")

problems = []
for p in targets:
    for n, line in enumerate(p.read_text().splitlines(), 1):
        if line.lstrip().startswith("#"):
            continue          # the explanations above live in comments
        for pattern, why in BANNED:
            if pattern.search(line):
                problems.append(f"{p.relative_to(home)}:{n}: {why}\n    {line.strip()}")

if problems:
    raise SystemExit("non-portable shell:\n" + "\n".join(problems))
print(f"portable: {len(targets)} shell files clean")
