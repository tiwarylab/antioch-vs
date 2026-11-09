# UCSF Chimera executes this script with its embedded Python 2.7 interpreter
# (`chimera --nogui`), never ANTIOCH-VS's own Python. Python 2.7 syntax only: no
# f-strings, no walrus, no `X | None` annotations. Do not modernize this file,
# by hand or by tool.
# ASCII ONLY: any non-ASCII byte here is a SyntaxError under Python 2 unless a
# PEP 263 encoding line is declared. Constitution Principle XIII.

import os

from chimera import MSMSModel, openModels, runCommand
from WriteDMS import writeDMS

# =====================================================================
# CONFIGURATION BLOCK (Edit these values using any text editor)
# =====================================================================
DENSITY_VALUE = int(os.environ.get('SURF_DENSITY', 5))

# Look for the pipeline's requested filename, or default to the standard
default_dms = "receptor_surf_d%s.dms" % DENSITY_VALUE
OUTPUT_FILENAME = os.environ.get("CHIMERA_OUT_DMS", default_dms)
# =====================================================================

# Run surface command
runCommand('surf')

# UCSF Chimera executes this helper with Python 2.
runCommand('setattr s density %s' % DENSITY_VALUE)

# Extract model and write out DMS
surf = openModels.list(modelTypes=[MSMSModel])[0]
writeDMS(surf, OUTPUT_FILENAME)
