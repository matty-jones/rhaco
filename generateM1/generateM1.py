import mbuild as mb
import numpy as np

class m1_unit_cell(mb.Compound):
    # This class will contain the unit cell for manipulation and replication
    def __init__(self, stoichiometry_dict):
        # Call the mb.Compound initialisation
        super().__init__()
        # Load the unit cell
        mb.load('templateM1.pdb', compound=self)
        # TODO: Assign the atoms based on the input stoichiometry
        # Replacable atoms in the matrix are assigned as type 'X'
        # Note: In both Py2 and Py3, subsequent calls to keys() and values() with no
        # intervening modifications will directly correspond \cite{Documentation}
        atom_types = list(stoichiometry_dict.keys())
        atom_ratios = np.array(list(stoichiometry_dict.values()))
        probabilities = list(atom_ratios / np.sum(atom_ratios))
        for particle in self.particles():
            if particle.name == 'X':
                # `Randomly' select an atom type based on the biases given in stoichiometry_dict
                particle.name = np.random.choice(atom_types, p=probabilities)
        # Check all the 'X' atom_types got updated
        assert('X' not in [particle.name for particle in self.particles()])
        # TODO: Create ports to allow the unit cell to be repeated along each dimension (2D to start)


class m1_surface(mb.Compound):
    # This class will describe the surface and consist of several m1_unit_cell instances in a specified dimension
    # Default stoichiometry found in: Nanostructured Catalysts: Selective Oxidations (Hess and Schl\"ogl, 2011, RSC)
    def __init__(self, surface_dimensions=[1, 1], stoichiometry_dict={'Mo': 1, 'V': 0.2, 'Nb': 0.17}):
        # Call the mb.Compound initialisation
        super().__init__()
        current_cell = m1_unit_cell(stoichiometry_dict)
        self.add(current_cell)


if __name__ == "__main__":
    m1Crystal = m1_surface()
    m1Crystal.save('test.hoomdxml', overwrite=True)
