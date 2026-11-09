# UCSF Chimera executes this script with its embedded Python 2.7 interpreter
# (`chimera --nogui`), never ANTIOCH-VS's own Python. Python 2.7 syntax only: no
# f-strings, no walrus, no `X | None` annotations. Do not modernize this file,
# by hand or by tool.
# ASCII ONLY: any non-ASCII byte here is a SyntaxError under Python 2 unless a
# PEP 263 encoding line is declared. Constitution Principle XIII.

import os

# =====================================================================
# CONFIGURATION BLOCK (Edit these values using any text editor)
# =====================================================================
OUTPUT_PDB         = os.environ.get("CHIMERA_OUT_PDB", "receptor_wH.pdb")

OPTIMIZE_HBONDS    = True
PRESERVE_PKA_NAMES = True
# =====================================================================

from chimera import runCommand

# Step 1: Dynamically build the optimized addh arguments list
addh_args = []

if OPTIMIZE_HBONDS:
    addh_args.append("hbond true inIsolation false")
else:
    addh_args.append("hbond false inIsolation true")

if PRESERVE_PKA_NAMES:
    addh_args.append("useHisName true useGluName true useAspName true useLysName true useCysName true")
else:
    addh_args.append("useHisName false useGluName false useAspName false useLysName false useCysName false")

# Step 2: Combine and execute the hydrogen addition command
addh_cmd = "addh {}".format(" ".join(addh_args))
runCommand(addh_cmd)

# Step 3: Save the protonated structure to a new PDB file
# '#0' targets the active model in the session
write_cmd = "write format pdb #0 %s" % OUTPUT_PDB
runCommand(write_cmd)
