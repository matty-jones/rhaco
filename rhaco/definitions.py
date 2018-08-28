import os

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
PDB_LIBRARY = os.path.join(PROJECT_ROOT, 'compounds')
FF_LIBRARY = os.path.join(PROJECT_ROOT, 'forcefields')
FOYER_FF_FORMATS = ['xml']
EXTERNAL_FF_FORMATS = ['eam.fs', 'eam.alloy']
# Atom masses obtained from NIST
ATOM_MASSES = {'H': 1.0078400, 'He': 4.002600, 'C': 12.009600, 'O': 15.999030,
               'Ni': 58.69344, 'Mn': 54.93804, 'Ga': 69.72310, 'Si': 28.08400,
               'S': 31.972071, 'N': 14.033074, 'Mo': 95.96000, 'Nb': 92.90638,
               'Te': 127.6000, 'V': 50.941500, 'Cu': 63.54600, 'Ag': 107.8682,
               'Au': 196.9666,}
