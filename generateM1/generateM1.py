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
        # Specify some dimensions for the templateM1.pdb file that we're using (these are esimated from the atom positions
        # in the template)
        x_extent = 2.1189  # mbuild assumes nm
        y_extent = 2.6738  # mbuild assumes nm
        # OUTER LOOP: Multiply up each x_row to create as many y repeats as specified
        previous_row = None
        complete_cell_matrix = []  # This is required for new bonds across diagonal elements
        for y_repeat in range(surface_dimensions[1]):
            current_row = []
            # INNER LOOP: First, create as many x repeats as specified
            # Note: Each cell has 159 atoms in it
            previous_cell = None
            for x_repeat in range(surface_dimensions[0]):
                current_cell = m1_unit_cell(stoichiometry_dict)
                current_row.append(current_cell)
                current_cell.translate([x_repeat * x_extent, y_repeat * y_extent, 0])
                self.add(current_cell)
                if previous_cell is not None:
                    self.add_bond([previous_cell[60], current_cell[21]])
                    self.add_bond([previous_cell[137], current_cell[13]])
                    self.add_bond([previous_cell[134], current_cell[16]])
                    self.add_bond([previous_cell[122], current_cell[16]])
                    self.add_bond([previous_cell[19], current_cell[119]])
                    self.add_bond([previous_cell[66], current_cell[33]])
                    self.add_bond([previous_cell[18], current_cell[65]])
                previous_cell = current_cell
            if previous_row is not None:
                for cell_ID, current_cell_y in enumerate(current_row):
                    previous_cell_y = previous_row[cell_ID]
                    self.add_bond([previous_cell_y[72], current_cell_y[27]])
                    self.add_bond([previous_cell_y[1], current_cell_y[58]])
                    self.add_bond([previous_cell_y[1], current_cell_y[73]])
                    self.add_bond([previous_cell_y[4], current_cell_y[123]])
                    self.add_bond([previous_cell_y[4], current_cell_y[141]])
                    self.add_bond([previous_cell_y[6], current_cell_y[141]])
                    self.add_bond([previous_cell_y[6], current_cell_y[156]])
                    self.add_bond([previous_cell_y[114], current_cell_y[12]])
                    self.add_bond([previous_cell_y[159], current_cell_y[12]])
            complete_cell_matrix.append(current_row)
            previous_row = current_row
        # Now that the cell_matrix is complete, there might be a few more bonds to add in
        # Across the diagonals (i.e. [0, 0] bonded to [1, 1]; [0, 1] bonded to [1, 2] etc.)
        # Go across all rows first
        for y_coord in range(surface_dimensions[1]):
            for x_coord in range(surface_dimensions[0]):
                if (x_coord + 1 < surface_dimensions[0]) and (y_coord + 1 < surface_dimensions[1]):
                    first_cell = complete_cell_matrix[x_coord][y_coord]
                    second_cell = complete_cell_matrix[x_coord + 1][y_coord + 1]
                    self.add_bond([first_cell[61], second_cell[21]])


if __name__ == "__main__":
    m1Crystal = m1_surface(surface_dimensions=[2, 2])
    m1Crystal.save('test.hoomdxml', overwrite=True)
