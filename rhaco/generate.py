import mbuild as mb
import numpy as np
import argparse
import copy
import os
import re
import zlib
import base64
from rhaco.definitions import PDB_LIBRARY, FF_LIBRARY, FOYER_FF_FORMATS, EXTERNAL_FF_FORMATS, ATOM_MASSES
import xml.etree.cElementTree as ET
from collections import OrderedDict


# Conversion factors
G_TO_AMU = 6.0222E23
CM_TO_NM = 1.0000E07


# Set the defaults for all the required arguments
defaults_dict = {'stoichiometry': {'Mo': 1, 'V': 0.15, 'Nb': 0.13, 'Te': 0.12},
                 'dimensions': [1, 1, 1],
                 'template': 'templateM1.pdb',
                 'reactant_composition': {'C2H6': 1},
                 'reactant_rigid': False,
                 'crystal_separation': 25.0,
                 'z_reactor_size': 20.0,
                 'reactant_num_mol': None,
                 'reactant_density': None,
                 'forcefield': None,
                 'integrate_crystal': False}


def split_argument_into_dictionary(argument):
    """
    Uses the string parsing on the argument to
    split the string into a dictionary.
    Requires:
        argument - string in the from of a dictionary
    Reutrns:
        combine_list - dictionary
    """
    # Create empty dictionary
    dictionary = {}
    # Remove curly brackets
    argument = argument[1:-1]
    # Remove whitespace
    argument = "".join([char for char in argument if char != " "])
    # Identify and split the keys based off of the ',' separator
    argument = argument.split(",")
    # We now have key:val pairs for each element, add these to dict
    for key_val in argument:
        key_val_split = key_val.split(":")
        key = key_val_split[0]
        val = key_val_split[1]
        if ((key[0] == '"') and (key[-1] == '"')) or ((key[0] == "'") and (key[-1] == "'")):
            key = key[1:-1]
        dictionary[key] = float(val)
    return dictionary


class crystal_unit_cell(mb.Compound):
    # This class will contain the unit cell for manipulation and replication
    def __init__(self, template, stoichiometry_dict):
        # Call the mb.Compound initialisation
        super().__init__()
        # Load the unit cell
        mb.load(os.path.join(PDB_LIBRARY, template), compound=self)
        # Replacable atoms in the matrix are assigned as type `X'
        # Note: In both Py2 and Py3, subsequent calls to keys() and values()
        # with no intervening modifications will directly correspond
        # \cite{PyDocumentation}
        atom_types, atom_probs, _ = calculate_probabilities(stoichiometry_dict)
        for particle in self.particles():
            if particle.name == 'X':
                # `Randomly' select an atom type based on the biases given in
                # stoichiometry_dict
                particle.name = np.random.choice(atom_types, p=atom_probs)
        # # Check all the 'X' atom_types got updated
        # assert('X' not in [particle.name for particle in self.particles()])


class crystal_surface(mb.Compound):
    # This class will describe the surface and consist of several
    # crystal_unit_cell instances in a specified dimension
    # Default stoichiometry found in: Nanostructured Catalysts: Selective
    # Oxidations (Hess and Schl\"ogl, 2011, RSC)
    def __init__(self, surface_dimensions, template, stoichiometry_dict,
                 crystal_bonds, x_extent, y_extent, z_extent):
        # Call the mb.Compound initialisation
        super().__init__()
        # OUTER LOOP: Create multiple layers based on the input dimensions
        for z_repeat in range(surface_dimensions[2]):
            # MIDDLE LOOP: Multiply up each x_row to create as many y repeats
            # as specified
            # complete_cell_matrix is required to keep track of new bonds
            # across diagonal elements
            complete_cell_matrix = []
            previous_row = None
            for y_repeat in range(surface_dimensions[1]):
                current_row = []
                # INNER LOOP: First, create as many x repeats as specified
                # Note: Each cell has 159 atoms in it
                previous_cell = None
                for x_repeat in range(surface_dimensions[0]):
                    print("\rAdding " + repr([x_repeat, y_repeat, z_repeat])
                          + " to system...", end=" ")
                    current_cell = crystal_unit_cell(template,
                                                     stoichiometry_dict)
                    current_row.append(current_cell)
                    current_cell.translate([x_repeat * x_extent,
                                            y_repeat * y_extent,
                                            z_repeat * z_extent])
                    self.add(current_cell)
                    if crystal_bonds and (previous_cell is not None):
                        self.add_x_connecting_bonds(previous_cell,
                                                    current_cell)
                    previous_cell = current_cell
                if crystal_bonds and (previous_row is not None):
                    for cell_ID, current_cell_y in enumerate(current_row):
                        previous_cell_y = previous_row[cell_ID]
                        self.add_y_connecting_bonds(previous_cell_y,
                                                    current_cell_y)
                complete_cell_matrix.append(current_row)
                previous_row = current_row
            # Now that the cell_matrix is complete for this layer, there might
            # be a few more bonds to add in.
            # Go across all rows first
            for y_coord in range(surface_dimensions[1]):
                # Now each column
                for x_coord in range(surface_dimensions[0]):
                    # Bonds located across the diagonals (i.e. [0, 0] bonded
                    # to [1, 1]; [0, 1] bonded to [1, 2] etc.)
                    if crystal_bonds and (x_coord + 1 < surface_dimensions[0])\
                       and (y_coord + 1 < surface_dimensions[1]):
                        first_cell = complete_cell_matrix[x_coord][y_coord]
                        second_cell = complete_cell_matrix[
                            x_coord + 1][y_coord + 1]
                        self.add_diagonal_connecting_bonds(first_cell,
                                                           second_cell)
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


class crystal_system(mb.Compound):
    # This class will describe the surface and consist of several
    # crystal_unit_cell instances in a specified dimension
    # Default stoichiometry found in: Nanostructured Catalysts: Selective
    # Oxidations (Hess and Schl\"ogl, 2011, RSC)
    def __init__(self, bottom_crystal, top_crystal, crystal_separation):
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
        # Now shift both crystals in the z-direction away from each other by
        # the (crystal_separation/2.0)
        # Note that crystal_separation is given in Angstroems but currently
        # in nm
        bottom_crystal.translate([0, 0, crystal_separation / 20.0])
        top_crystal.translate([0, 0, -crystal_separation / 20.0])
        # Add both crystal planes to the system
        self.add(bottom_crystal)
        self.add(top_crystal)


class mbuild_template(mb.Compound):
    # This class will contain the mb compound for the reactant
    def __init__(self, template, rigid):
        # Call the mb.Compound initialisation
        super().__init__()
        # Load the unit cell
        mb.load(
            os.path.join(PDB_LIBRARY, ''.join(template.split('.pdb')) + '.pdb'),
            compound=self,
            rigid=rigid,
        )
        atom_types = [atom_type.split('[')[0] for atom_type in
                      self.labels.keys() if '[' in atom_type]
        self.mass = self.get_mass(atom_types)

    def get_mass(self, atom_types):
        mass = 0.0
        for atom_type in set(atom_types):
            try:
                type_mass = ATOM_MASSES[atom_type]
            except KeyError:
                print("***** WARNING *****")
                print("Type", atom_type, "not found in definitions.py."
                      " Assuming a mass of 1.0 AMU (the density specified by"
                      " -rd can no longer be trusted!)")
                type_mass = 1.0
            mass += type_mass * atom_types.count(atom_type)
        return mass


def parse_forcefields(forcefield_string):
    PERMITTED_FF_FORMATS = FOYER_FF_FORMATS + EXTERNAL_FF_FORMATS
    foyer_forcefield_list = []
    external_forcefield_list = []
    if (forcefield_string[0] == "[") and (forcefield_string[-1] == "]"):
        # User has specified a list
        forcefield_string = forcefield_string[1:-1]
    if "," in forcefield_string:
        # Split based on commas
        forcefield_string = forcefield_string.split(",")
    else:
        # Split based on whitespace
        forcefield_string = forcefield_string.split()
    for forcefield in forcefield_string:
        if forcefield.lower() == "none":
            return None
        # Remove any additional whitespace (if split by ", ")
        forcefield = forcefield.strip()
        if ((forcefield[0] == "'") and (forcefield[-1] == "'")) or ((forcefield[0] == '"') and (forcefield[-1] == '"')):
            # Cut out the ' or "
            forcefield = forcefield[1:-1]
        # Check file exists
        forcefield_exists = False
        forcefield_loc = None
        # Loop over permitted file extensions (inc. no extension in case already specified):
        for file_extension in [""] + PERMITTED_FF_FORMATS:
            if len(file_extension) > 0:
                forcefield_w_ext = ".".join([forcefield, file_extension])
            else:
                forcefield_w_ext = str(forcefield)
            # First check the FF library
            forcefield_loc = os.path.join(FF_LIBRARY, forcefield_w_ext)
            forcefield_exists = check_forcefield_exists(forcefield_loc)
            if forcefield_exists:
                break
            # Then check the cwd
            forcefield_exsits = check_forcefield_exists(forcefield_w_ext)
            if forcefield_exists:
                break
        if not forcefield_exists:
            print("***** WARNING *****")
            print("Forcefield", forcefield, "not found with any compatible extension:",
                  repr(PERMITTED_FF_FORMATS), "in either the FF_LIBRARY dir", FF_LIBRARY,
                  "or cwd.")
            print("No forcefield data will be saved at generation.")
            return None
        else:
            FF_file = os.path.split(forcefield_loc)[1]
            # Can't use os.path.splitext here because EAM has two
            # file formats and we need both of them
            FF_format = ".".join(FF_file.split('.')[1:])
            if FF_format in FOYER_FF_FORMATS:
                foyer_forcefield_list.append(forcefield_loc)
            elif FF_format in EXTERNAL_FF_FORMATS:
                external_forcefield_list.append(forcefield_loc)
    if len(foyer_forcefield_list) > 1:
        print("The following forcefields will be implemented with Foyer:", repr(foyer_forcefield_list))
    elif len(foyer_forcefield_list) == 1:
        print("The following forcefield will be implemented with Foyer:", foyer_forcefield_list[0])
    if len(external_forcefield_list) > 1:
        print("The following files will be used as forcefields:", repr(external_forcefield_list))
    elif len(external_forcefield_list) == 1:
        print("The following file will be used as a forcefield:", external_forcefield_list[0])
    return [foyer_forcefield_list, external_forcefield_list]


def check_forcefield_exists(path):
    try:
        with open(path, 'r'):
            pass
        return True
    except FileNotFoundError:
        return False


def parse_reactant_positions(position_string):
    position_coords = []
    position_string = "".join(position_string.split(" "))
    if (position_string[0:2] == "[[") and (position_string[-2:] == "]]"):
        # [[pos1x, pos1y, pos1z], [pos2x, pos2y, pos2z]]
        position_string = position_string[1:-1]
    if "],[" in position_string:
        position_string = position_string.split("],[")
    elif "][" in position_string:
        position_string = position_string.split("][")[1:-1]
    if position_string.count("[") == 1:
        # Only one position specified
        position_coords.append(list(map(float, position_string[1:-1].split(','))))
    else:
        # Multiple positions specified
        position_string[0] = position_string[0][1:]
        position_string[-1] = position_string[-1][:-1]
        for element in position_string:
            position_coords.append(list(map(float, element.split(','))))
    return position_coords


def create_morphology(args):
    output_file = create_output_file_name(args)
    print("Generating first surface (bottom)...")
    surface1 = crystal_surface(args.dimensions, args.template,
                               args.stoichiometry, args.crystal_bonds,
                               args.crystal_x, args.crystal_y, args.crystal_z)
    print("Generating second surface (top)...")
    surface2 = crystal_surface(args.dimensions, args.template,
                               args.stoichiometry, args.crystal_bonds,
                               args.crystal_x, args.crystal_y, args.crystal_z)
    # Now create the system by combining the two surfaces
    system = crystal_system(surface1, surface2, args.crystal_separation)
    # Get the crystal IDs because we're going to need them later so that HOOMD
    # knows not to integrate them.
    crystal_IDs = range(system.n_particles)
    # Now we can populate the box with reactant
    print("Surfaces generated. Generating reactant...")
    reactant_components, reactant_probs, reactant_masses = calculate_probabilities(
        args.reactant_composition, ratio_type='number')
    # Define the regions that the hydrocarbons can go in, so we don't end
    # up with them between layers
    box_top = mb.Box(mins=[-(args.crystal_x * args.dimensions[0]) / 2.0,
                           -(args.crystal_y * args.dimensions[1]) / 2.0,
                           args.crystal_separation / 20.0
                           + (args.crystal_z * args.dimensions[2])],
                     maxs=[(args.crystal_x * args.dimensions[0]) / 2.0,
                           (args.crystal_y * args.dimensions[1]) / 2.0,
                           args.z_reactor_size / 2.0])
    box_bottom = mb.Box(mins=[-(args.crystal_x * args.dimensions[0]) / 2.0,
                              -(args.crystal_y * args.dimensions[1]) / 2.0,
                              -args.z_reactor_size / 2.0],
                        maxs=[(args.crystal_x * args.dimensions[0]) / 2.0,
                              (args.crystal_y * args.dimensions[1]) / 2.0,
                              -args.crystal_separation / 20.0
                              - (args.crystal_z * args.dimensions[2])])
    box_top_vol = np.prod(box_top.maxs - box_top.mins)
    box_bottom_vol = np.prod(box_bottom.maxs - box_bottom.mins)
    reactor_vol = box_top_vol + box_bottom_vol

    # No reactant specified
    if (args.reactant_density is None) and (args.reactant_num_mol is None):
        number_of_reactant_mols = 0
    # Both specified - print error and use the number
    elif (args.reactant_density is not None) and (args.reactant_num_mol is not None):
        print("".join(["Both -rn (", str(args.reactant_num_mol), ") and -rd (",
                       str(args.reactant_density), ") specified. Using -rn ",
                       str(args.reactant_num_mol), "..."]))
        number_of_reactant_mols = args.reactant_num_mol
    # Number specified and not density
    elif (args.reactant_density is None) and (args.reactant_num_mol is not None):
        number_of_reactant_mols = args.reactant_num_mol
    # Density specified but not numbers
    elif (args.reactant_density is not None) and (args.reactant_num_mol is None):
        # Work backwards to come up with how many reactant molecules are needed
        # to get the specified density.
        # Get the average mass for each molecule based on reactant probabilities
        mass_per_n = np.sum([reactant_masses[key] * reactant_probs[index] for index,
                             key in enumerate(reactant_components)])
        # Given the reactor volume and the specified reactant density, calculate
        # the total number of reactant molecules needed.
        # Need to convert from CGS (g/cm^{3} -> AMU/nm^{3})
        reactant_density_conv = args.reactant_density * G_TO_AMU / (CM_TO_NM**3)
        number_of_reactant_mols = int(reactant_density_conv * reactor_vol / mass_per_n)
    if args.reactant_position is not None:
        if len(args.reactant_position) != number_of_reactant_mols:
            print("-rp --reactant_position flag has been specified with length",
                  len(args.reactant_position), "but", str(number_of_reactant_mols),
                  "reactant molecules have been requested!")
            print("Ignoring specified positions and using packmol to randomly place reactant.")
            args.reactant_position = None
        else:
            for _, position in enumerate(args.reactant_position):
                nanoparticle = mbuild_template(
                    reactant_components[0], args.reactant_rigid
                )
                nanoparticle.translate_to(np.array(position))
                system.add(nanoparticle)
    if args.reactant_position is None:
        # Randomly place reactants using packmol
        if number_of_reactant_mols == 1:
            # Only 1 molecule to place, so put it on top of the crystals
            reactant_top = mb.packing.fill_box(
                mbuild_template(reactant_components[0], args.reactant_rigid),
                1, box_top, seed=np.random.randint(0, 2**31 - 1)
            )
            system.add(reactant_top)
        elif number_of_reactant_mols > 1:
            reactant_compounds = []
            n_compounds = []
            for compound_index, reactant_molecule in enumerate(reactant_components):
                reactant_compounds.append(
                    mbuild_template(reactant_molecule, args.reactant_rigid)
                )
                n_compounds.append(int(np.round(
                    reactant_probs[compound_index] * number_of_reactant_mols)))
            # Split the n_compounds to the top and bottom. Note: Top will always have
            # the extra molecule for odd n_compounds
            top_n = list(map(int, np.ceil(np.array(n_compounds) / 2.0)))
            bot_n = list(map(int, np.floor(np.array(n_compounds) / 2.0)))
            reactant_top = mb.packing.fill_box(
                reactant_compounds, top_n, box_top,
                seed=np.random.randint(0, 2**31 - 1)
            )
            reactant_bottom = mb.packing.fill_box(
                reactant_compounds, bot_n, box_bottom,
                seed=np.random.randint(0, 2**31 - 1)
            )
            system.add(reactant_top)
            system.add(reactant_bottom)

    if "M1UnitCell.pdb" in args.template:
        # Check the separation of crystal and reactant that we will use later is
        # correct. Get the set of atom types that produce the crystal (and
        # don't include the base atom type, which we asusme to be oxygen).
        names = [particle.name for particle_ID, particle in
                 enumerate(system.particles()) if (particle_ID in crystal_IDs)
                 and (particle.name != 'O')]
        # Ensure that this is the same as the stoichiometry dictionary keys
        assert(np.array_equal(args.stoichiometry.keys(), set(names)))

    # Generate the morphology box based on the input parameters
    system_box = mb.Box(mins=[-(args.crystal_x * args.dimensions[0]) / 2.0,
                              -(args.crystal_y * args.dimensions[1]) / 2.0,
                              -args.z_reactor_size / 2.0],
                        maxs=[(args.crystal_x * args.dimensions[0]) / 2.0,
                              (args.crystal_y * args.dimensions[1]) / 2.0,
                              args.z_reactor_size / 2.0])
    print("Morphology generated.")
    # Note this logic means a user cannot specify their own FF with the same
    # name as one in our libary!
    if (args.forcefield is None) or ((len(args.forcefield[0]) == 0) and (len(args.forcefield[1]) == 0)):
        print("Saving morphology...")
        system.save(output_file, overwrite=True, box=system_box)
        # Fix the images because mbuild doesn't set them correctly
        morphology = fix_images(output_file)
    else:
        print("Applying forcefield...")
        system.save(output_file, overwrite=True, box=system_box,
                    forcefield_files=args.forcefield[0])
        # Fix the images because mbuild doesn't set them correctly
        morphology = fix_images(output_file)
        # Create a new key here that we can use to tell the simulate.py what
        # forcefield input files it will need to care about (if they aren't
        # already covered by Foyer).
        if len(args.forcefield[1]) > 0:
            morphology["external_forcefields_attrib"] = {"num": str(len(args.forcefield[1]))}
            morphology["external_forcefields_text"] = args.forcefield[1]
    # Identify the crystal atoms in the system by renaming their type to
    # X_<PREVIOUS ATOM TYPE> so we know not to integrate them in HOOMD
    if args.integrate_crystal is False:
        morphology = rename_crystal_types(morphology, crystal_IDs)
    write_morphology_xml(morphology, output_file)
    print("Output generated. Exitting...")


def check_bonds(morphology, bond_dict, box_dims):
    for bond in morphology['bond_text']:
        posn1 = np.array(list(map(float,
                                  morphology['position_text'][int(bond[1])])))
        + (np.array(list(map(float, morphology['image_text'][int(bond[1])])))
           * box_dims)
        posn2 = np.array(list(map(float,
                                  morphology['position_text'][int(bond[2])])))
        + (np.array(list(map(float, morphology['image_text'][int(bond[2])])))
           * box_dims)
        delta_position = posn1 - posn2
        out_of_range = [np.abs(delta_position[axis]) > box_dims[axis] / 2.0
                        for axis in range(3)]
        if any(out_of_range):
            print("Periodic bond found:", bond, "because delta_position =",
                  delta_position, ">=", box_dims, "/ 2.0")
            morphology = move_bonded_atoms(bond[1], morphology, bond_dict,
                                           box_dims)
    return morphology


def zero_out_images(morphology):
    morphology['image_text'] = [['0', '0', '0']]\
        * len(morphology['position_text'])
    morphology['image_attrib'] = {'num': morphology['position_attrib']['num']}
    return morphology


def get_bond_dict(morphology):
    bond_dict = {atom_id: [] for atom_id, atom_type in
                 enumerate(morphology['type_text'])}
    for bond in morphology['bond_text']:
        bond_dict[int(bond[1])].append(int(bond[2]))
        bond_dict[int(bond[2])].append(int(bond[1]))
    return bond_dict


def move_bonded_atoms(central_atom, morphology, bond_dict, box_dims):
    for bonded_atom in bond_dict[central_atom]:
        posn1 = np.array(list(map(float,
                                  morphology['position_text'][central_atom])))
        posn2 = np.array(list(map(float,
                                  morphology['position_text'][bonded_atom])))
        delta_position = posn1 - posn2
        moved = False
        for axis, value in enumerate(delta_position):
            if value > box_dims[axis] / 2.0:
                morphology['position_text'][bonded_atom][axis] = str(
                    posn2[axis] + box_dims[axis])
                moved = True
            if value < -box_dims[axis] / 2.0:
                morphology['position_text'][bonded_atom][axis] = str(
                    posn2[axis] - box_dims[axis])
                moved = True
        if moved:
            morphology = move_bonded_atoms(bonded_atom, morphology, bond_dict,
                                           box_dims)
    return morphology


def load_morphology_xml(xml_file_name):
    morphology_dictionary = OrderedDict()
    with open(xml_file_name, 'r') as xml_file:
        xml_tree = ET.parse(xml_file)
    root = xml_tree.getroot()
    morphology_dictionary['root_tag'] = root.tag
    morphology_dictionary['root_attrib'] = root.attrib
    morphology_dictionary['root_text'] = root.text
    for config in root:
        morphology_dictionary['config_tag'] = config.tag
        morphology_dictionary['config_attrib'] = config.attrib
        morphology_dictionary['config_text'] = config.text
        for child in config:
            if len(child.attrib) > 0:
                morphology_dictionary[child.tag + '_attrib'] = {
                    key.lower(): val for key, val in child.attrib.items()}
            else:
                morphology_dictionary[child.tag + '_attrib'] = {}
            if child.text is not None:
                morphology_dictionary[child.tag + '_text'] = [
                    x.split() for x in child.text.split('\n') if len(x) > 0]
            else:
                morphology_dictionary[child.tag + '_text'] = []
    return morphology_dictionary


def check_wrapped_positions(input_dictionary):
    box_dims = [float(input_dictionary['box_attrib'][axis]) for axis in
                ['lx', 'ly', 'lz']]
    atom_positions = np.array([np.array(list(map(float, _))) for _ in
                               input_dictionary['position_text']])
    atom_images = np.array([np.array(list(map(int, _))) for _ in
                            input_dictionary['image_text']])
    xhi = box_dims[0] / 2.0
    xlo = -box_dims[0] / 2.0
    yhi = box_dims[1] / 2.0
    ylo = -box_dims[1] / 2.0
    zhi = box_dims[2] / 2.0
    zlo = -box_dims[2] / 2.0
    for atom_ID in range(len(atom_positions)):
        while atom_positions[atom_ID][0] > xhi:
            atom_positions[atom_ID][0] -= box_dims[0]
            atom_images[atom_ID][0] += 1
        while atom_positions[atom_ID][0] < xlo:
            atom_positions[atom_ID][0] += box_dims[0]
            atom_images[atom_ID][0] -= 1
        while atom_positions[atom_ID][1] > yhi:
            atom_positions[atom_ID][1] -= box_dims[1]
            atom_images[atom_ID][1] += 1
        while atom_positions[atom_ID][1] < ylo:
            atom_positions[atom_ID][1] += box_dims[1]
            atom_images[atom_ID][1] -= 1
        while atom_positions[atom_ID][2] > zhi:
            atom_positions[atom_ID][2] -= box_dims[2]
            atom_images[atom_ID][2] += 1
        while atom_positions[atom_ID][2] < zlo:
            atom_positions[atom_ID][2] += box_dims[2]
            atom_images[atom_ID][2] -= 1
    input_dictionary['position_text'] = list([list(map(str, _)) for _ in
                                              atom_positions])
    input_dictionary['image_text'] = list([list(map(str, _)) for _ in
                                           atom_images])
    return input_dictionary


def write_morphology_xml(morphology_dictionary, output_file_name):
    # morphology_dictionary is a bunch of keys with the tagnames given for
    # both attributes and text: tag + '_attrib', tag + '_text'
    # The only boilerplate bits are the 'root_tag', 'root_attrib', and
    # 'root_text', which is (obviously) the outmost layer of the xml.
    # Immediately inside are the 'config_tag', 'config_attrib', and
    # 'config_text'. Everything else is a child of config.
    morphology_dictionary = check_wrapped_positions(morphology_dictionary)
    # Build the xml tree.
    root = ET.Element(morphology_dictionary['root_tag'],
                      **morphology_dictionary['root_attrib'])
    root.text = morphology_dictionary['root_text']
    config = ET.Element(morphology_dictionary['config_tag'],
                        **morphology_dictionary['config_attrib'])
    config.text = morphology_dictionary['config_text']
    # Find the remaining elements to make (set is easier here, but a
    # disordered structure, so instead we use lists to keep the order
    # consistent with reading in).
    all_child_tags = ['_'.join(key.split('_')[:-1]) for key in
                      morphology_dictionary.keys()
                      if '_'.join(key.split('_')[:-1]) not in
                      ['root', 'config']]
    child_tags = []
    for tag in all_child_tags:
        # The list comprehension makes two blank entries for some reason and
        # I can't work out why. This will just skip those two entries, as well
        # as make the set.
        if (tag not in child_tags) and (len(tag) > 0):
            child_tags.append(tag)
    for child_tag in child_tags:
        child = ET.Element(child_tag,
                           **morphology_dictionary[child_tag + '_attrib'])
        if child_tag != "external_forcefields":
            data_to_write = '\n'.join(['\t'.join(el) for el in
                                       morphology_dictionary[
                                           child_tag + '_text']])
        else:
            data_to_write = '\n'.join([el for el in morphology_dictionary[child_tag + "_text"]])
        if len(data_to_write) > 0:
            child.text = '\n' + data_to_write + '\n'
        child.tail = '\n'
        config.append(child)
    root.insert(0, config)
    tree = ET.ElementTree(root)
    tree.write(output_file_name, xml_declaration=True, encoding='UTF-8')
    print("XML file written to", str(output_file_name) + "!")


def rename_crystal_types(input_dictionary, AAIDs):
    for atom_index in AAIDs:
        previous_type = input_dictionary['type_text'][atom_index][0]
        input_dictionary['type_text'][atom_index] = ['X_' + previous_type]
    return input_dictionary


def fix_images(file_name):
    print("Fixing the images to ensure everything is wrapped within"
          " the box...")
    morphology = load_morphology_xml(file_name)
    morphology = zero_out_images(morphology)
    bond_dict = get_bond_dict(morphology)
    box_dims = [float(morphology['box_attrib'][axis]) for axis in
                ['lx', 'ly', 'lz']]
    morphology = check_bonds(morphology, bond_dict, box_dims)
    return morphology


def calculate_probabilities(input_dictionary, ratio_type='stoic'):
    '''
    This function takes an input dictionary corresponding to the relative
    ratios of some parameter, then returns normalized probabilities for each
    option that can be used to make choices with appropriate bias.
    '''
    choices = list(input_dictionary.keys())
    number_ratios = np.array(list(input_dictionary.values()))
    mass_dict = None
    if ratio_type == 'number':
        mass_dict = get_masses(input_dictionary.keys())
    probabilities = list(number_ratios / np.sum(number_ratios))
    return choices, probabilities, mass_dict


def get_masses(reactant_names):
    number_ratio = {}
    mass_dict = {}
    # First split out the key names into atoms
    for reactant_name in reactant_names:
        # Consult the mass lookup table
        total_mass = mbuild_template(reactant_name, False).mass
        mass_dict[reactant_name] = total_mass
    # Return dictionary of number ratios
    return mass_dict


def create_output_file_name(args, file_type='hoomdxml'):
    if args.signac is True:
        return 'output.hoomdxml'
    else:
        output_file = "out"
        for (arg_name, arg_val) in sorted(args._get_kwargs()):
            try:
                if isinstance(arg_val, dict):
                    if sorted(arg_val) == sorted(defaults_dict[arg_name]):
                        continue
                else:
                    if (arg_val == defaults_dict[arg_name]):
                        continue
            except KeyError:
                continue
            output_file += "-"
            if arg_name == "stoichiometry":
                output_file += "S_"
                for key, val in arg_val.items():
                    output_file += str(key) + ":" + str(val) + "_"
                output_file = output_file[:-1]
            elif arg_name == "reactant_composition":
                output_file += "RC_"
                for key, val in arg_val.items():
                    output_file += str(key) + ":" + str(val) + "_"
                output_file = output_file[:-1]
            elif arg_name == "dimensions":
                output_file += "D_" + "x".join(list(map(str, arg_val)))
            elif arg_name == "template":
                output_file += "T_" + args.template.split("/")[-1].split(
                    ".")[0]
            elif arg_name == "forcefield":
                if len(args.forcefield[0]) > 0:
                    output_file += "F1"
                    for FF in args.forcefield[0]:
                        output_file += "_" + os.path.split(FF)[1]
                if len(args.forcefield[1]) > 0:
                    if len(args.forcefield[0]) > 0:
                        output_file += "-"
                    output_file += "F2"
                    for FF in args.forcefield[1]:
                        output_file += "_" + os.path.split(FF)[1]
            elif arg_val is False:
                output_file += arg_name[0].upper() + "_Off"
            elif arg_val is True:
                output_file += arg_name[0].upper() + "_On"
            else:
                output_file += arg_name[0].upper() + "_" + str(arg_val)
        return output_file + '.' + file_type


def main():
    parser = argparse.ArgumentParser(prog='rhaco-create-morph',
                                     formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--stoichiometry",
                        type=split_argument_into_dictionary,
                        default={'Mo': 1, 'V': 0.15, 'Nb': 0.13, 'Te': 0.12},
                        required=False,
                        help='''Specify a stoichiometry for the surface.\n
                        Atoms marked as type 'X' in the template will be
                        replaced with atoms of the type given by the keys in
                        the input dictionary, with probabilities determined
                        from the ratios given by the values of the input
                        dictionary.\n
                        For example: -s "{'Mo': 1, 'V': 0.15, 'Nb': 0.13,
                        'Te': 0.12}" will create a surface where there are 5
                        Mo atoms for every V, and 0.85 Nb.\n
                        The above default value is taken from
                        10.1524/zkri.219.3.152.29091\n''')
    parser.add_argument("-d", "--dimensions",
                        type=lambda d: [int(dim) for dim in
                                        d.split('x') if len(dim) > 0],
                        default=[1, 1, 1],
                        required=False,
                        help='''Specify the number of cells to stitch into the
                        surface (integers), along the x- and y-directions.\n
                        For example: -d 2x2x1 will create a surface containing
                        4 unit cells, with 2 along the x-axis, 2 along the
                        y-axis, and one layer thick.\n''')
    parser.add_argument("-t", "--template",
                        type=str,
                        default='M1UnitCell.pdb',
                        required=False,
                        help='''Identify the unit cell file to be used to
                        create the surface.\n
                        Note the unit cells are located in the PDB_LIBRARY
                        directory, which defaults to rhaco/compounds.\n
                        For example: -t "M1UnitCell.pdb".\n
                        If not specified, the default
                        PDB_LIBRARY/M1UnitCell.pdb is used.''')
    parser.add_argument("-c", "--crystal_separation",
                        type=float,
                        default=25.0,
                        required=False,
                        help='''Assign a physical separation (in Angstroems) to
                        the bottom planes of the two crystals corresponding to
                        the top and bottom of the simulation volume within the
                        periodic box.\n
                        Note that this is not the same as the z_reactor_size,
                        which describes the region available to hydrocarbon
                        molecules in the simulation.\n
                        This value should be larger than the interaction
                        cut-off specified in the forcefield (pair or Coulombic)
                        to prevent the self-interaction of the two
                        surfaces.\n
                        For example: -c 25.0.''')
    parser.add_argument("-b", "--crystal_bonds",
                        action='store_true',
                        required=False,
                        help='''Enable the creation of bonds between unit cells
                        of the crystal. This is useful for visualisation, and
                        required when bonds are to be simulated with MD further
                        down the line.\n
                        Note that these bonds are current hard-coded in
                        rhaco/generate.py for the M1 crystal structure. The user
                        can modify these manually if so desired (to place
                        custom bonds), or the user can contact the repository
                        maintainers to express their interest in having an
                        input file that describes the locations of the bonds.
                        ''')
    parser.add_argument("-z", "--z_reactor_size",
                        type=float,
                        default=20.0,
                        required=False,
                        help='''Assign the z-axis size of the simulation
                        (in nm).\n
                        This defines the region available for hydrocarbons to
                        move around in, between the two catalyst plates
                        (region depth = z_reactor_size - plane_separation -
                        (z_extent * dimension[2])).\n
                        Note that this is not the same as the plane_separation,
                        which describes the physical separation between the
                        bottom layers of the two flipped M1 crystals.\n
                        For example: -z 20.0.\n''')
    parser.add_argument("-rc", "--reactant_composition",
                        type=split_argument_into_dictionary,
                        default={'C2H6': 1},
                        required=False,
                        help='''Specify a reactant composition.\n
                        The input proportions are assumed to be by moles, and
                        will be normalized and used to create the reactant,
                        which is incident on the catalyst as it flows into the
                        reactor.\n
                        For example: -rc "{'C2H6': 3, 'O2': 2, 'He': 5}" will
                        set a relative reactant proportion to 0.6:0.4:1
                        respectively by number of molecules.\n
                        Note that keys in the dictionary must be the same as
                        the pdb files located in the PDB_LIBRARY of Rhaco, and
                        the corresponding values are interpreted as the
                        proportion by moles.\n''')
    parser.add_argument("-rr", "--reactant_rigid",
                        type=bool,
                        default=False,
                        required=False,
                        help='''If True, then each reactant molecule will be
                        treated as its own rigid body.''')
    parser.add_argument("-rn", "--reactant_num_mol",
                        type=int,
                        default=None,
                        help='''Set the number of reactant component
                        molecules to be included in the system.\n
                        For example: -rn 1000.\n
                        If unspecified, then no reactant will be
                        used.\n
                        Note that -rn overrides -rd if both are
                        specified.''')
    parser.add_argument("-rd", "--reactant_density",
                        type=float,
                        default=None,
                        help='''Set the density of the reactor
                        fluence in g/cm^{3}.\n
                        For example: -rd 0.05.\n
                        If unspecified, then no reactant will be
                        used.
                        Note that -rn overrides -rd if both are
                        specified.''')
    parser.add_argument("-rp", "--reactant_position",
                        type=parse_reactant_positions,
                        default=None,
                        help='''Set the position of the reactant
                        molecules.\n
                        This only makes sense for a small number of reactant
                        particles (i.e. nanoparticle initializations).
                        For example: -rp [[-50, 0, 50], [50, 0, 50]].\n
                        If unspecified, then reactants will be packed randomly
                        using packmol.
                        ''')
    parser.add_argument("--gecko",
                        action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument("-f", "--forcefield",
                        type=parse_forcefields,
                        default=None,
                        required=False,
                        help='''Use Foyer to set the forcefield to use when
                        running the simulation.\n
                        Note the forcefields are located in the FF_LIBRARY
                        directory, which defaults to rhaco/forcefields.\n
                        For example: -f FF_opls_uff.\n
                        If 'None' specified, the compound will not be saved
                        with forcefield information.''')
    parser.add_argument("-i", "--integrate_crystal",
                        action='store_true',
                        help='''Use HOOMD to integrate the crystal atoms during
                        the molecular dynamics simulation.\n
                        By default, HOOMD will not update the positions and
                        velocities of the crystal atoms during the MD phase
                        (recommended - adding in bond, angle, and dihedral
                        constraints to the crystal surface will dramatically
                        reduce computational efficiency).\n
                        This parameter can be set to allow HOOMD to update the
                        crystal atom positions, providing their behaviour is
                        defined in the forcefield.\n
                        For example: -i.\n
                        If this flag is not passed, the crystal will remain
                        static, but still influence the motion of nearby
                        simulation elements.''')
    parser.add_argument("-sig", "--signac",
                        action='store_true',
                        help='''Change the naming nomenclature to be more
                        specific if signac-flow is not being used (via
                        rhaco-flow).\n
                        By default, output files will be named based on the
                        input parameters used to generate them.\n
                        However, this is not useful when using the signac
                        infrastructure provided by rhaco-flow.\n
                        In this case, signac does the heavy lifting of
                        determining which parameters were used to generate each
                        file, and so only an output.hoomdxml is created.
                        For example: -sig.\n''')
    parser.add_argument("-xx", "--crystal_x",
                        type=float,
                        default=2.148490,
                        required=False,
                        help='''Assign the periodicity of the crystal unit cell
                        in the x direction.\n
                        This defines the 'a' axis spacing between the unit
                        cells in the crystal in nm.\n
                        If unspecified this uses the default for M1 which is
                        2.148490 nm. (Taken from Desanto2006
                        (10.1007/s11244-006-0068-8))\n
                        For example: -xx 2.148490.\n''')
    parser.add_argument("-xy", "--crystal_y",
                        type=float,
                        default=2.664721,
                        required=False,
                        help='''Assign the periodicity of the crystal unit cell
                        in the y direction.\n
                        This defines the 'b' axis spacing between the unit
                        cells in the crystal in nm.\n
                        If unspecified this uses the default for M1 which is
                        2.664721 nm. (Taken from Desanto2006
                        (10.1007/s11244-006-0068-8))\n
                        For example: -xy 2.664721.\n''')
    parser.add_argument("-xz", "--crystal_z",
                        type=float,
                        default=0.400321,
                        required=False,
                        help='''Assign the periodicity of the crystal unit cell
                        in the z direction.\n
                        This defines the 'c' axis spacing between the unit
                        cells in the crystal in nm.\n
                        If unspecified this uses the default for M1 which is
                        0.400321 nm. (Taken from Desanto2006
                        (10.1007/s11244-006-0068-8))\n
                        For example: -xz 0.400321.\n''')
    args = parser.parse_args()
    if args.gecko:
        print(zlib.decompress(base64.decodebytes(
            b"""
            eJyFU0FuxCAMvPOKUS5OpGDujvoTJPOQPL62IWy6zbZWFIHtGYYBgP9C75H+bNxzziT
            CGcgiCkXLHmzIO9knFiMgx0HEKexjONCmkQYb306RAUsz4kcqo1FvAYJLrbdjMubYhq
            X37LaKicYjnUsa/b07IqZ9iVf+hHOv67ltV1fjSajDHlz9/X+IzDHcOWflkI2KSLfwz
            LJs8NT3tywd5Lp93m0GVu75epWjaMThROOrGDE2HIYODT5UVo5D4G2IfgXrytboghyN
            xq7XToXZD/nAkXxLtp7hKJvOxaV6lPgeorJb53r6bWlEbS5I6bCCyPF2LCgsRfcnPq0
            sVfcUpD9xX5RLqSbrN2qWUrMZn29lItf3qH/UEqoIVaRucAzMH56XwNz398D3F9Cy9C
            fxfPvcB2rfSW6mXw==
            """)).decode('utf-8'))
        exit()
    print("The run arguments for this job are:")
    print(print(args))
    create_morphology(args)
