import os

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
PDB_LIBRARY = os.path.join(PROJECT_ROOT, 'compounds')
FF_LIBRARY = os.path.join(PROJECT_ROOT, 'forcefields')
# Atom masses obtained from NIST
ATOM_MASSES = {'H': 1.00784, 'He': 4.00260, 'C': 12.00960, 'O': 15.99903,
               'Ni': 58.69344, 'Mn': 54.93804, 'Ga': 69.72310, 'Si': 28.08400}
