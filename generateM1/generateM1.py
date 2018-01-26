import mbuild as mb
import numpy as np
import argparse


# Specify some dimensions for the templateM1.pdb file that we're using (these are esimated from the atom positions
# in the template)
x_extent = 2.1189  # mbuild assumes nm
y_extent = 2.6738  # mbuild assumes nm


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
    def __init__(self, surface_dimensions, template, stoichiometry_dict, bonds_periodic):
        # Call the mb.Compound initialisation
        super().__init__()
        # OUTER LOOP: Multiply up each x_row to create as many y repeats as specified
        previous_row = None
        complete_cell_matrix = []  # This is required for new bonds across diagonal elements
        for y_repeat in range(surface_dimensions[1]):
            current_row = []
            # INNER LOOP: First, create as many x repeats as specified
            # Note: Each cell has 159 atoms in it
            previous_cell = None
            for x_repeat in range(surface_dimensions[0]):
                print("\rAdding " + repr([x_repeat, y_repeat]) + " to system...", end=' ')
                current_cell = m1_unit_cell(template, stoichiometry_dict)
                current_row.append(current_cell)
                current_cell.translate([x_repeat * x_extent, y_repeat * y_extent, 0])
                self.add(current_cell)
                if previous_cell is not None:
                    self.add_x_connecting_bonds(previous_cell, current_cell)
                previous_cell = current_cell
            if previous_row is not None:
                for cell_ID, current_cell_y in enumerate(current_row):
                    previous_cell_y = previous_row[cell_ID]
                    self.add_y_connecting_bonds(previous_cell_y, current_cell_y)
            complete_cell_matrix.append(current_row)
            previous_row = current_row
        # Now that the cell_matrix is complete, there might be a few more bonds to add in
        # Go across all rows first
        for y_coord in range(surface_dimensions[1]):
            # Bottom bonding to top over periodic boundary conditions
            if bonds_periodic and y_coord == 0:
                # Iterate over every cell in the row
                for x_coord in range(surface_dimensions[0]):
                    first_cell = complete_cell_matrix[surface_dimensions[1] - 1][x_coord]
                    second_cell = complete_cell_matrix[0][x_coord]
                    self.add_y_connecting_bonds(first_cell, second_cell)
            # Now each column
            for x_coord in range(surface_dimensions[0]):
                # Left hand side bonding to right hand side over periodic boundary conditions
                if bonds_periodic and x_coord == 0:
                    print(complete_cell_matrix)
                    for rowNo, row in enumerate(complete_cell_matrix):
                        for colNo, column in enumerate(row):
                            print(rowNo, colNo, column.pos)
                    first_cell = complete_cell_matrix[y_coord][surface_dimensions[0] - 1]
                    second_cell = complete_cell_matrix[y_coord][0]
                    self.add_x_connecting_bonds(first_cell, second_cell)
                # Bonds located across the diagonals (i.e. [0, 0] bonded to [1, 1]; [0, 1] bonded to [1, 2] etc.)
                if (x_coord + 1 < surface_dimensions[0]) and (y_coord + 1 < surface_dimensions[1]):
                    first_cell = complete_cell_matrix[x_coord][y_coord]
                    second_cell = complete_cell_matrix[x_coord + 1][y_coord + 1]
                    self.add_diagonal_connecting_bonds(first_cell, second_cell)
        print()

    def add_x_connecting_bonds(self, cell1, cell2):
        self.add_bond([cell1[60], cell2[21]])
        self.add_bond([cell1[137], cell2[13]])
        self.add_bond([cell1[134], cell2[16]])
        self.add_bond([cell1[122], cell2[16]])
        self.add_bond([cell1[19], cell2[119]])
        self.add_bond([cell1[66], cell2[33]])
        self.add_bond([cell1[18], cell2[65]])

    def add_y_connecting_bonds(self, cell1, cell2):
        self.add_bond([cell1[72], cell2[27]])
        self.add_bond([cell1[1], cell2[58]])
        self.add_bond([cell1[1], cell2[73]])
        self.add_bond([cell1[4], cell2[123]])
        self.add_bond([cell1[4], cell2[141]])
        self.add_bond([cell1[6], cell2[141]])
        self.add_bond([cell1[6], cell2[156]])
        self.add_bond([cell1[114], cell2[12]])
        self.add_bond([cell1[159], cell2[12]])

    def add_diagonal_connecting_bonds(self, cell1, cell2):
        self.add_bond([cell1[61], cell2[21]])


class m1_system(mb.Compound):
    # This class will describe the surface and consist of several m1_unit_cell instances in a specified dimension
    # Default stoichiometry found in: Nanostructured Catalysts: Selective Oxidations (Hess and Schl\"ogl, 2011, RSC)
    def __init__(self, bottom_crystal, top_crystal, plane_separation):
        # Call the mb.Compound initialisation
        super().__init__()
        # First, center both crystals according to their center of geometry (we don't care about masses here)
        top_crystal.translate(-top_crystal.pos)
        bottom_crystal.translate(-bottom_crystal.pos)
        # Flip the top crystal around
        top_crystal.rotate(np.pi, [1, 0, 0])
        # Now shift both crystals away from each other by the plane_separation
        bottom_crystal.translate([0, 0, plane_separation/2.0])
        top_crystal.translate([0, 0, -plane_separation/2.0])
        # Add both crystal planes to the system
        self.add(bottom_crystal)
        self.add(top_crystal)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--stoichiometry", type=lambda s: {str(key[1:-1]): float(val) for [key, val] in [splitChar for splitChar in [cell.split(':') for cell in [_ for _ in s[1:-1].split(',') if len(_) > 0]] if len(splitChar) > 0]}, default={'Mo': 1, 'V': 0.2, 'Nb': 0.17}, required=False,
                        help='''Specify a stoichiometry for the surface.\n
                        Atoms marked as type 'X' in the template will be replaced with atoms of the type given by the keys in the input dictionary, with probabilities determined from the ratios given by the values of the input dictionary.\n
                        For example: -s "{'Mo': 1, 'V': 0.2, 'Nb': 0.17}" will create a surface where there are 5 Mo atoms for every V, and 0.85 Nb.\n
                        If not specified, the default stoichiometry is set to {'Mo': 1, 'V': 0.2, 'Nb': 0.17}''')
    parser.add_argument("-d", "--dimensions", type=lambda d: [int(_) for _ in d.split('x') if len(_) > 0], default=[2, 2], required=False,
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
    parser.add_argument("-p", "--plane_separation", type=float, default=1.0, required=False,
                       help='''Assign a pysical separation (in nm) to the two planes corresponding to the top and bottom of the simulation volume within the periodic box.
                        Note that this is not the same as the z_extent, which describes the region available to ethane molecules in the simulation.
                        This value should be larger than the interaction cut-off specified in the forcefield (pair or Coulombic) to prevent the self-interaction of the surface.
                        For example: -p 1.0.
                        If not specified, the default value of one nanometre is used.
                       ''')
    parser.add_argument("-z", "--z_extent", type=float, default=20.0, required=False,
                       help='''Assign the z-axis extent of the simulation (in nm). This defines the region available for ethane molecules to move around in, between the two catalyst plates (region depth = z_extent - plane_separation).
                        Note that this is not the same as the plane_separation, which describes the physical separation between the two M1 crystal planes.
                        For example: -z 20.0.
                        If not specified, the default value of 20 nanometres is used.
                       ''')
    parser.add_argument("-b", "--bonds_periodic", action='store_false', default='output.hoomdxml', required=False,
                       help='''A boolean that determines whether periodic bonds are inserted across the x-y plane that connect the left side of the system with the right side, and the top to the bottom across the simulation's periodic boundary conditions.
                        For example: -b.
                        The default for this parameter is 'True', but passing this flag will change it to 'False' (which makes it look prettier in VMD for outputs.
                       ''')
    args = parser.parse_args()
    print("Generating first surface (bottom)...")
    surface1 = m1_surface(args.dimensions, args.template, args.stoichiometry, args.bonds_periodic)
    print("Generating second surface (top)...")
    surface2 = m1_surface(args.dimensions, args.template, args.stoichiometry, args.bonds_periodic)
    print("Surface generated. Saving as", args.output + "...")
    system = m1_system(surface1, surface2, args.plane_separation)
    system_box = mb.Box(mins = [-(x_extent * args.dimensions[0])/2.0, -(y_extent * args.dimensions[1])/2.0, -args.z_extent/2.0],
                        maxs = [(x_extent * args.dimensions[0])/2.0, (y_extent * args.dimensions[1])/2.0, args.z_extent/2.0])
    system.save(args.output, overwrite=True, box = system_box)
    print("Output generated. Exitting...")
