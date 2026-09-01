#!/usr/bin/env python3
"""Consistency checks for the Slurm scripts.

These are the files least likely to be exercised on a laptop and most likely to
rot when the pipeline changes shape, so check them statically:

  * numbering contiguous from 00, matching the order they are submitted in
  * every script uses the one shared preamble
  * every script names a job and logs to logs/
  * every step a job asks the pipeline to run is a step the pipeline has
  * the submit chain in run_pipeline.sh --mode slurm submits exactly these
"""
import re
import sys
from pathlib import Path

home = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
pipeline = (home / "scripts" / "run_pipeline.sh").read_text()
steps = set(re.search(r"ALL_STEPS=\(([^)]*)\)", pipeline, re.S).group(1).split())

scripts = sorted(p for p in (home / "slurm").glob("*.sbatch"))
if not scripts:
    raise SystemExit("no .sbatch scripts found")
names = [p.name for p in scripts]

prefixes = [int(n[:2]) for n in names]
if prefixes != list(range(len(prefixes))):
    raise SystemExit(f"gap in slurm numbering: {names}")

for p in scripts:
    text = p.read_text()
    if "_common.sh" not in text:
        raise SystemExit(f"{p.name} does not source the shared preamble")
    if not re.search(r"^#SBATCH --job-name", text, re.M):
        raise SystemExit(f"{p.name} has no --job-name")
    if not re.search(r"^#SBATCH --output=logs/", text, re.M):
        raise SystemExit(f"{p.name} does not write its log under logs/")

    asked = set(re.findall(r"morpheus_run --only \"?([\w-]+)", text))
    asked |= {s for m in re.findall(r"for step in ([^;]+); do", text)
              for s in m.split()}
    unknown = asked - steps - {"${step}"}
    if unknown:
        raise SystemExit(f"{p.name} runs unknown step(s): {sorted(unknown)}")

# --mode slurm must submit every script, in order, and no others.
submitted = re.findall(r'"\$\{S\}/(\S+?\.sbatch)"', pipeline)
if submitted != names:
    raise SystemExit(
        f"--mode slurm submits {submitted}, but slurm/ holds {names}")

print("slurm consistent:", ", ".join(names))
