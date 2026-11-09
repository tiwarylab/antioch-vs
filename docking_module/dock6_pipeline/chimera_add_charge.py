# UCSF Chimera executes this script with its embedded Python 2.7 interpreter
# (`chimera --nogui`), never ANTIOCH-VS's own Python. Python 2.7 syntax only: no
# f-strings, no walrus, no `X | None` annotations. Do not modernize this file,
# by hand or by tool.
# ASCII ONLY: any non-ASCII byte here is a SyntaxError under Python 2 unless a
# PEP 263 encoding line is declared. Constitution Principle XIII.

import os

# =====================================================================
# CONFIGURATION BLOCK
# =====================================================================
OUTPUT_MOL2     = os.environ.get("CHIMERA_OUT_MOL2", "receptor_wH_wcH.mol2")
CHARGE_SCOPE    = "all"
CHARGE_MODEL    = "ff14SB"
NONSTD_METHOD   = "am1"
# =====================================================================

import chimera
from chimera import replyobj, runCommand


# ---------------------------------------------------------------------
# CLEVER TOPOLOGICAL REPAIR FUNCTION
# ---------------------------------------------------------------------
def fix_amber_termini():
    """
    Scans the molecular graph for N-termini. 
    If it finds the MDAnalysis naming quirk (H, H2, H3), it renames 'H' to 'H1'.
    If it finds a physically missing hydrogen, it deletes them so addh can rebuild.
    """
    # Use Chimera's native replyobj instead of print()
    replyobj.info("[*] Scanning for N-termini naming quirks from MD tools...\n")
    
    for m in chimera.openModels.list(modelTypes=[chimera.Molecule]):
        for res in m.residues:
            if 'N' in res.atomsMap:
                n_atom = res.atomsMap['N'][0]
                
                is_n_term = True
                for neighbor in n_atom.neighbors:
                    if neighbor.name == 'C' and neighbor.residue != res:
                        is_n_term = False
                        break
                
                if is_n_term:
                    hydrogens = [neighbor for neighbor in n_atom.neighbors if neighbor.element.name == 'H']
                    
                    if len(hydrogens) == 3:
                        for h in hydrogens:
                            if h.name == 'H':
                                target = h.oslIdent()
                                # UCSF Chimera executes this helper with Python 2.
                                replyobj.info("[+] Renaming N-terminal hydrogen %s to H1\n" % target)
                                runCommand("setattr a name H1 " + target)
                                
                    elif 0 < len(hydrogens) < 3:
                        replyobj.info(
                            "[!] Found invalid number of hydrogens at N-terminus: %s%s\n"
                            % (res.type, res.id)
                        )
                        for h in hydrogens:
                            target = h.oslIdent()
                            replyobj.info("    -> Deleting %s for rebuild.\n" % target)
                            runCommand("delete " + target)

# 1. Execute the repair 
fix_amber_termini()

# 2. Fill in missing hydrogens (Only triggers if atoms were deleted)
runCommand("addh")
# ---------------------------------------------------------------------

# 3. Construct and execute the charge assignment command 
charge_cmd = "addcharge %s chargeModel %s method %s" % (
    CHARGE_SCOPE,
    CHARGE_MODEL,
    NONSTD_METHOD,
)
runCommand(charge_cmd)

# 4. Write out the parameterized structural data to Mol2 format
write_cmd = "write format mol2 #0 %s" % OUTPUT_MOL2
runCommand(write_cmd)
