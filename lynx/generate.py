import mbuild as mb
import numpy as np
import argparse
import copy
import os
from lynx.definitions import PDB_LIBRARY, FF_LIBRARY

# The crystallographic unit cell parameters
# Taken from Desanto2006 (10.1007/s11244-006-0068-8)
x_extent = 2.148490
y_extent = 2.664721
z_extent = 0.400321

# Set the defaults for all the required arguments
defaults_dict = {'stoichiometry': {'Mo': 1, 'V': 0.3, 'Nb': 0.15, 'Te': 0.15},
                 'dimensions': [1, 1, 1],
                 'template': 'templateM1.pdb',
                 'organic': 'ethane.pdb',
                 'crystal_separation': 25.0,
                 'z_box_size': 20.0,
                 'bonds_periodic': True,
                 'number_of_organic_mols': 200,
                 'forcefield': 'FF_opls_uff.xml'
                }


class m1_unit_cell(mb.Compound):
    # This class will contain the unit cell for manipulation and replication
    def __init__(self, template, stoichiometry_dict):
        # Call the mb.Compound initialisation
        super().__init__()
        # Load the unit cell
        mb.load(os.path.join(PDB_LIBRARY, template), compound=self)
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
    # This class will describe the surface and consist of several m1_unit_cell instances in a
    # specified dimension
    # Default stoichiometry found in: Nanostructured Catalysts: Selective Oxidations
    # (Hess and Schl\"ogl, 2011, RSC)
    def __init__(self, surface_dimensions, template, stoichiometry_dict):
        # Call the mb.Compound initialisation
        super().__init__()
        # OUTER LOOP: Create multiple layers based on the input dimensions
        for z_repeat in range(surface_dimensions[2]):
            # MIDDLE LOOP: Multiply up each x_row to create as many y repeats as specified
            complete_cell_matrix = []  # This is required for new bonds across diagonal elements
            previous_row = None
            for y_repeat in range(surface_dimensions[1]):
                current_row = []
                # INNER LOOP: First, create as many x repeats as specified
                # Note: Each cell has 159 atoms in it
                previous_cell = None
                for x_repeat in range(surface_dimensions[0]):
                    print("\rAdding " + repr([x_repeat, y_repeat, z_repeat])\
                          + " to system...", end=' ')
                    current_cell = m1_unit_cell(template, stoichiometry_dict)
                    current_row.append(current_cell)
                    current_cell.translate([x_repeat * x_extent,
                                            y_repeat * y_extent,
                                            z_repeat * z_extent])
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
            # Now that the cell_matrix is complete for this layer, there might be a few
            # more bonds to add in
            # Go across all rows first
            for y_coord in range(surface_dimensions[1]):
                ## Skip creating periodic bonds because I'm pretty sure this isn't correct
                ## Bottom bonding to top over periodic boundary conditions
                #if bonds_periodic and y_coord == 0:
                #    # Iterate over every cell in the row
                #    for x_coord in range(surface_dimensions[0]):
                #        first_cell = complete_cell_matrix[surface_dimensions[1] - 1][x_coord]
                #        second_cell = complete_cell_matrix[0][x_coord]
                #        print('Y_connecting_bonds for', [x_coord, y_coord, z_repeat])
                #        self.add_y_connecting_bonds(first_cell, second_cell)
                # Now each column
                for x_coord in range(surface_dimensions[0]):
                    ## Skip creating periodic bonds because I'm pretty sure this isn't correct
                    ## Left hand side bonding to right hand side over periodic boundary conditions
                    #if bonds_periodic and x_coord == 0:
                    #    first_cell = complete_cell_matrix[y_coord][surface_dimensions[0] - 1]
                    #    second_cell = complete_cell_matrix[y_coord][0]
                    #    print('X_connecting_bonds for', [x_coord, y_coord, z_repeat])
                    #    self.add_x_connecting_bonds(first_cell, second_cell)
                    # Bonds located across the diagonals (i.e. [0, 0] bonded to [1, 1];
                    # [0, 1] bonded to [1, 2] etc.)
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
    # This class will describe the surface and consist of several m1_unit_cell
    # instances in a specified dimension
    # Default stoichiometry found in: Nanostructured Catalysts: Selective Oxidations
    # (Hess and Schl\"ogl, 2011, RSC)
    def __init__(self, bottom_crystal, top_crystal, crystal_separation, solvent):
        # Call the mb.Compound initialisation
        super().__init__()
        # Firstly, get the current COM positions for each plane. This will be
        # important later
        top_COM = copy.deepcopy(top_crystal.pos)
        bottom_COM = copy.deepcopy(bottom_crystal.pos)
        # Then, center both crystals according to their center of geometry
        # (we don't care about masses here)
        top_crystal.translate(-top_crystal.pos)
        bottom_crystal.translate(-bottom_crystal.pos)
        # Now the top crystal is centered at the origin, we have no issues
        # flipping it around
        top_crystal.rotate(np.pi, [1, 0, 0])
        # Now shift both crystals in the z-direction away from each other by the
        # (crystal_separation/2.0)
        # Note that crystal_separation is given in Angstroems but currently in nm
        bottom_crystal.translate([0, 0, crystal_separation/20.0])
        top_crystal.translate([0, 0, -crystal_separation/20.0])
        # Add both crystal planes to the system
        self.add(bottom_crystal)
        self.add(top_crystal)
        # Finally, add the solvent
        if solvent is not None:
            self.add(solvent)


class mbuild_template(mb.Compound):
    # This class will contain the mb compound for hydrocarbon
    def __init__(self, template):
        # Call the mb.Compound initialisation
        super().__init__()
        # Load the unit cell
        mb.load(os.path.join(PDB_LIBRARY, template), compound=self)


def create_morphology(args):
    output_file = create_output_file_name(args)
    print("Generating first surface (bottom)...")
    surface1 = m1_surface(args.dimensions, args.template, args.stoichiometry)
    print("Generating second surface (top)...")
    surface2 = m1_surface(args.dimensions, args.template, args.stoichiometry)
    if args.number_of_organic_mols > 0:
        # Now we can populate the box with hydrocarbons
        print("Surfaces generated. Generating hydrocarbons...")
        hydrocarbon = mbuild_template(args.organic)
        # Define the regions that the hydrocarbons can go in, so we don't end up with
        # them between layers
        box_top = mb.Box(mins = [-(x_extent * args.dimensions[0])/2.0,
                                 -(y_extent * args.dimensions[1])/2.0,
                                 args.crystal_separation/20.0 + (z_extent * args.dimensions[2])],
                            maxs = [(x_extent * args.dimensions[0])/2.0,
                                    (y_extent * args.dimensions[1])/2.0,
                                    args.z_box_size/2.0])
        box_bottom = mb.Box(mins = [-(x_extent * args.dimensions[0])/2.0,
                                    -(y_extent * args.dimensions[1])/2.0,
                                    -args.z_box_size/2.0],
                            maxs = [(x_extent * args.dimensions[0])/2.0,
                                    (y_extent * args.dimensions[1])/2.0,
                                    -args.crystal_separation/20.0 - (z_extent * args.dimensions[2])])
        solvent = mb.packing.fill_region([hydrocarbon] * 2, [args.number_of_organic_mols // 2] * 2,
                                         [box_bottom, box_top])
    else:
        solvent = None
    # Now create the system by combining the two surfaces and the solvent
    system = m1_system(surface1, surface2, args.crystal_separation, solvent)
    # Generate the morphology box based on the input parameters
    system_box = mb.Box(mins = [-(x_extent * args.dimensions[0])/2.0,
                                -(y_extent * args.dimensions[1])/2.0,
                                -args.z_box_size/2.0],
                        maxs = [(x_extent * args.dimensions[0])/2.0,
                                (y_extent * args.dimensions[1])/2.0,
                                args.z_box_size/2.0])
    print("Morphology generated. Saving as", output_file + "...")
    if args.forcefield is not None:
        try:
            # Check the FF library first
            forcefield_file = os.path.join(FF_LIBRARY, args.forcefield)
            with open(forcefield_file, 'r') as file_handle:
                pass
        except FileNotFoundError:
            # Otherwise use the cwd
            forcefield_file = args.forcefield
        mbuild.formats.gsdwriter.write_gsd(system, output_file, overwrite=True, box=system_box,
                                           forcefield_name=forcefield_file)
    else:
        system.save(output_file, overwrite=True, box=system_box)
    print("Output generated. Exitting...")


def create_output_file_name(args, file_type='hoomdxml'):
    output_file = "out"
    for (arg_name, arg_val) in sorted(args._get_kwargs()):
        print(arg_name, arg_val)
        if (arg_val == defaults_dict[arg_name]):
            continue
        output_file += "_"
        if arg_name == 'stoichiometry':
            output_file += "S"
            for key, val in arg_val.items():
                output_file += str(key) + str(val)
        elif arg_name == 'dimensions':
            output_file += "D" + "x".join(list(map(str, arg_val)))
        elif arg_name == 'template':
            output_file += "T" + args.template.split('/')[-1].split('.')[0]
        elif arg_val is False:
            output_file += arg_name[0].upper() + "Off"
        elif arg_val is True:
            output_file += arg_name[0].upper() + "On"
        else:
            output_file += arg_name[0].upper() + str(arg_val)
    return output_file + '.' + file_type


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--stoichiometry",
                        type=lambda s: {str(key[1:-1]): float(val) for [key, val] in [
                            splitChar for splitChar in [cell.split(':') for cell in [
                                _ for _ in s[1:-1].split(',') if len(_) > 0]]\
                            if len(splitChar) > 0]},
                        default={'Mo': 1, 'V': 0.3, 'Nb': 0.15, 'Te': 0.15},
                        required=False,
                        help='''Specify a stoichiometry for the surface.\n
                        Atoms marked as type 'X' in the template will be replaced with atoms of the
                        type given by the keys in the input dictionary, with probabilities determined
                        from the ratios given by the values of the input dictionary.\n
                        For example: -s "{'Mo': 1, 'V': 0.3, 'Nb': 0.15, 'Te': 0.15}" will create a
                        surface where there are 5 Mo atoms for every V, and 0.85 Nb.\n
                        If not specified, the default stoichiometry is set to {'Mo': 1, 'V': 0.3,
                        ' b': 0.15, 'Te': 0.15}''')
    parser.add_argument("-d", "--dimensions",
                        type=lambda d: [int(_) for _ in d.split('x') if len(_) > 0],
                        default=[1, 1, 1],
                        required=False,
                        help='''Specify the number of cells to stitch into the surface (integers),
                        along the x- and y-directions.\n
                        For example: -d 2x2x1 will create a surface containing 4 unit cells, with 2
                        along the x-axis, 2 along the y-axis, and one layer thick.\n
                        If not specified, the default dimensions are a single unit cell producing a
                        1x1x1 surface.
                       ''')
    parser.add_argument("-t", "--template",
                        type=str,
                        default='M1UnitCell.pdb',
                        required=False,
                        help='''Identify the unit cell file to be used to create the surface.\n
                        Note the unit cells are located in the PDB_LIBRARY directory, which defaults
                        to lynx/compounds.\n
                        For example: -t "M1UnitCell.pdb".\n
                        If not specified, the default PDB_LIBRARY/M1UnitCell.pdb is used.
                       ''')
    parser.add_argument("-c", "--crystal_separation",
                        type=float,
                        default=25.0,
                        required=False,
                        help='''Assign a pysical separation (in Angstroems) to the bottom planes of
                        the two crystals corresponding to the top and bottom of the simulation volume
                        within the periodic box.\n
                        Note that this is not the same as the z_box_size, which describes the region
                        available to hydrocarbon molecules in the simulation.\n
                        This value should be larger than the interaction cut-off specified in the
                        forcefield (pair or Coulombic) to prevent the self-interaction of the two
                        surfaces.\n
                        For example: -c 25.0.
                        If not specified, the default value of 2.5 nanometres is used.
                       ''')
    parser.add_argument("-z", "--z_box_size",
                        type=float,
                        default=20.0,
                        required=False,
                        help='''Assign the z-axis size of the simulation (in nm).\n
                        This defines the region available for hydrocarbons to move around in,
                        between the two catalyst plates (region depth = z_box_size
                        - plane_separation - (z_extent * dimension[2])).\n
                        Note that this is not the same as the plane_separation, which describes the
                        physical separation between the bottom layers of the two flipped M1 crystals.\n
                        For example: -z 20.0.\n
                        If not specified, the default value of 20 nanometres is used.
                       ''')
    parser.add_argument("-o", "--organic",
                        type=str,
                        default='ethane.pdb',
                        required=False,
                        help='''Set the hydrocarbon file to use in the simulation.\n
                        Note that the hydrocarbon files should be located in the PDB_LIBRARY directory,
                        which defaults to lynx/compounds.\n
                        For example: -o "ethane.pdb".\n
                        If not specified, the default PDB_LIBRARY/ethane.pdb is used.
                       ''')
    parser.add_argument("-n", "--number_of_organic_mols",
                        type=int,
                        default=200,
                        required=False,
                        help='''Set the number of organic hydrocarbons to be included in the system.\n
                        For example: -n 200.\n
                        If not specified, the default value of 200 hydrocarbons is used.
                       ''')
    parser.add_argument("-f", "--forcefield",
                        type=str,
                        default=None,
                        required=False,
                        help='''Use Foyer to set the forcefield to use when running the simulation.\n
                        Note the forcefields are located in the FF_LIBRARY directory, which defaults to
                        lynx/forcefields.\n
                        For example: -f FF_opls_uff.\n
                        If not specified, the compound will not be saved with forcefield information.
                       ''')
    args = parser.parse_args()
    create_morphology(args)
