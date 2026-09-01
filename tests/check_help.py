#!/usr/bin/env python3
"""The help text must describe the program that exists.

Help drifts: a flag gets renamed, an example keeps the old spelling, and the
first thing a new user copies is the one command that cannot work. So pull every
flag out of both help screens and require the parser to know it, and require
every flag the parser knows to be documented.

Also guards against the failure this replaced: help built with `sed -n 'M,Np'`
over the script's own comment block, whose line range drifts until it starts
printing shell source at the user.
"""
import re
import subprocess
import sys
from pathlib import Path

home = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
morpheus = home / "bin" / "morpheus"
pipeline = home / "scripts" / "run_pipeline.sh"

def run(*args):
    r = subprocess.run([str(morpheus), *args], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"`morpheus {' '.join(args)}` exited {r.returncode}")
    return r.stdout

top = run("--help")
for alias in ("-h", "help"):
    if run(alias) != top:
        raise SystemExit(f"`morpheus {alias}` differs from `morpheus --help`")

run_help = run("run", "--help")

# No shell source may leak into either screen.
for name, text in (("morpheus --help", top), ("morpheus run --help", run_help)):
    for leak in ("set -euo", "#!/usr/bin/env", "####", 'BASH_SOURCE'):
        if leak in text:
            raise SystemExit(f"{name} leaks shell source: {leak!r}")
    if not text.strip():
        raise SystemExit(f"{name} is empty")

# Every flag the parser accepts, from its own case statement.
parser = re.search(r"while \[\[ \$# -gt 0 \]\]; do(.*?)\ndone",
                   pipeline.read_text(), re.S).group(1)
known = set(re.findall(r"(--[\w-]+)\)", parser))
known |= set(re.findall(r"(--[\w-]+)\|", parser))
known.add("-h")
# Flags spelled --_something are internal plumbing (the per-gene driver
# re-enters this script for one gene) and are deliberately not user-facing.
known = {f for f in known if not f.startswith("--_")}
if len(known) < 15:
    raise SystemExit(f"only found {len(known)} flags in the parser; parsing broke")

# Every flag mentioned in either help screen must be one of them.
mentioned = set(re.findall(r"(?<![\w-])(--[a-z][\w-]*)", top + run_help))
# sbatch's own flags, named in the text as things to hand to --slurm_extra or
# as the per-job defaults living in the sbatch files. Not Morpheus flags.
SBATCH_FLAGS = {"--mem", "--qos", "--partition", "--account", "--export"}
unknown = mentioned - known - SBATCH_FLAGS
if unknown:
    raise SystemExit(f"help mentions flags the parser rejects: {sorted(unknown)}")

# ...and every flag the parser accepts must appear in `run --help`, so nothing
# is reachable but undiscoverable. Underscore spellings are canonical; the
# hyphenated aliases are deliberately undocumented.
undocumented = {f for f in known
                if f not in run_help and f.replace("-", "_", 1) not in run_help}
undocumented -= {f for f in known if "-" in f[2:]}   # hyphenated aliases
if undocumented:
    raise SystemExit(f"parser accepts undocumented flags: {sorted(undocumented)}")

print(f"help consistent: {len(mentioned)} flags shown, {len(known)} accepted")
