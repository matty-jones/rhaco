import mbuild as mb
import numpy as np
import argparse

class m1_unit_cell(mb.Compound):
    # This class will contain the unit cell for manipulation and replication
    def __init__(self, template, stoichiometry_dict):
        # Call the mb.Compound initialisation
        super().__init__()
        # Load the unit cell
        mb.load(template, compound=self)
        # Replacable atoms in the matrix are assigned as type `X'
        # Note: In both Py2 and Py3, subsequent calls to keys() and values() with no
        # intervening modifications will directly correspond \cite{PyDocumentation}
        atom_types = list(stoichiometry_dict.keys())
        atom_ratios = np.array(list(stoichiometry_dict.values()))
        probabilities = list(atom_ratios / np.sum(atom_ratios))
        for particle in self.particles():
            if particle.name == 'X':
                # `Randomly' select an atom type based on the biases given in stoichiometry_dict
                particle.name = np.random.choice(atom_types, p=probabilities)
        # Check all the 'X' atom_types got updated
        #assert('X' not in [particle.name for particle in self.particles()])


class m1_surface(mb.Compound):
    # This class will describe the surface and consist of several m1_unit_cell instances in a specified dimension
    # Default stoichiometry found in: Nanostructured Catalysts: Selective Oxidations (Hess and Schl\"ogl, 2011, RSC)
    def __init__(self, surface_dimensions, template, stoichiometry_dict):
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
                print("Adding " + [x_repeat, y_repeat] + " to system...\r", end=' ')
                current_cell = m1_unit_cell(template, stoichiometry_dict)
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
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--stoichiometry", type=lambda s: {str(key[1:-1]): float(val) for [key, val] in [splitChar for splitChar in [cell.split(':') for cell in [_ for _ in s[1:-1].split(',') if len(_) > 0]] if len(splitChar) > 0]}, default={'Mo': 1, 'V': 0.2, 'Nb': 0.17}, required=False,
                        help='''Specify a stoichiometry for the surface.\n
                        Atoms marked as type 'X' in the template will be replaced with atoms of the type given by the keys in the input dictionary, with probabilities determined from the ratios given by the values of the input dictionary.\n
                        For example: -s "{'Mo': 1, 'V': 0.2, 'Nb': 0.17}" will create a surface where there are 5 Mo atoms for every V, and 0.85 Nb.\n
                        If not specified, the default stoichiometry is set to {'Mo': 1, 'V': 0.2, 'Nb': 0.17}''')
    parser.add_argument("-d", "--dimensions", type=lambda d: [int(_) for _ in d.split(',') if len(_) > 0], default=[1, 1], required=False,
                       help='''Specify the number of cells to stitch into the surface (integers), along the x- and y-directions.
                        For example: -d "6,5" will create a surface containing 30 unit cells, with 6 along the x-axis and 5 up the y-axis.
                        If not specified, the default dimensions are a single unit cell producing a 1x1 surface.
                       ''')
    parser.add_argument("-t", "--template", type=str, default='templateM1.pdb', required=False,
                       help='''Identify the unit cell file to be used to create the surface.
                        For example: -t "templateM1.pdb".
                        If not specified, the default ./templateM1.pdb is used.
                       ''')
    parser.add_argument("-o", "--output", type=str, default='output.hoomdxml', required=False,
                       help='''Identify the location of the output file containing the final surface structure.
                        For example: -o "output.hoomdxml".
                        If not specified, the default ./output.hoomdxml is used.
                       ''')
    args = parser.parse_args()
    m1Crystal = m1_surface(args.dimensions, args.template, args.stoichiometry)
    print("Crystal generated. Saving as", args.output + "...")
    m1Crystal.save(args.output, overwrite=True)
    print("Output generated. Exitting...")
