import os

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
PDB_LIBRARY = os.path.join(PROJECT_ROOT, "compounds")
FF_LIBRARY = os.path.join(PROJECT_ROOT, "forcefields")
FOYER_FF_FORMATS = ["xml"]
EXTERNAL_FF_FORMATS = ["eam.fs", "eam.alloy"]
# Atom masses obtained from NIST
ATOM_MASSES = {
    "H": 1.0078400,
    "He": 4.002600,
    "C": 12.009600,
    "N": 14.033074,
    "O": 15.999030,
    "Al": 26.98154,
    "Si": 28.08400,
    "S": 31.972071,
    "V": 50.941500,
    "Mn": 54.93804,
    "Ni": 58.69344,
    "Cu": 63.54600,
    "Ga": 69.72310,
    "Nb": 92.90638,
    "Mo": 95.96000,
    "Ag": 107.8682,
    "Te": 127.6000,
    "Au": 196.9666,
}
ATOMIC_NUMBERS = {
    "H": 1,
    "He": 2,
    "C": 6,
    "N": 7,
    "O": 8,
    "Al": 13,
    "Si": 14,
    "S": 16,
    "V": 23,
    "Mn": 25,
    "Ni": 28,
    "Cu": 29,
    "Ga": 31,
    "Nb": 41,
    "Mo": 42,
    "Ag": 47,
    "Te": 52,
    "Au": 79,
}

