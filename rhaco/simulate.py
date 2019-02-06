import argparse
import copy
import hoomd
import hoomd.md
import hoomd.metal
import hoomd.deprecated
import json
import re
import quaternion
import numpy as np
import xml.etree.cElementTree as ET
from rhaco.definitions import ATOM_MASSES


AVOGADRO = 6.022140857E23
BOLTZMANN = 1.38064852E-23
KCAL_TO_J = 4.184E3
AMU_TO_KG = 1.6605E-27
ANG_TO_M = 1E-10


def set_coeffs(
    file_name,
    system,
    distance_scaling,
    energy_scaling,
    nl_type,
    r_cut,
    groups_list,
    generate_arguments,
):
    """
    Read in the molecular dynamics coefficients exported by Foyer
    """
    coeffs_dict = get_coeffs(file_name, generate_arguments)
    if nl_type == "tree":
        nl = hoomd.md.nlist.tree()
    else:
        nl = hoomd.md.nlist.cell()
    # Ignore any surface-surface interactions
    nl.reset_exclusions(exclusions=["body", "bond", "angle", "dihedral"])
    log_quantities = ["temperature", "pressure", "volume"]
    for quantity in ["potential_energy", "kinetic_energy"]:
        for group in groups_list:
            log_quantities.append("_".join([quantity, group.name])),
    if len(coeffs_dict["pair_coeffs"]) != 0:
        print("Loading LJ pair coeffs...")
        if r_cut is None:
            max_sigma = np.max(
                list(map(float, np.array(coeffs_dict["pair_coeffs"])[:, 2]))
            )
            r_cut = 2.5 * max_sigma
            print("Setting r_cut to 2.5 * max_sigma =", r_cut)
        lj = hoomd.md.pair.lj(r_cut=r_cut, nlist=nl)
        lj.set_params(mode="xplor")
        if "pair_lj_energy" not in log_quantities:
            log_quantities.append("pair_lj_energy")
        for type1 in coeffs_dict["pair_coeffs"]:
            for type2 in coeffs_dict["pair_coeffs"]:
                if (
                    not generate_arguments["integrate_crystal"]
                    and type1[0][:2] == "X_"
                    and type2[0][:2] == "X_"
                ):
                    lj_r_cut = 0.0
                else:
                    lj_r_cut = r_cut
                lj.pair_coeff.set(
                    type1[0],
                    type2[0],
                    epsilon=np.sqrt(type1[1] * type2[1]) / (energy_scaling),
                    sigma=np.sqrt(type1[2] * type2[2]) / (distance_scaling),
                    r_cut=lj_r_cut,
                )
    if len(coeffs_dict["bond_coeffs"]) != 0:
        print("Loading harmonic bond coeffs...")
        harmonic_bond = hoomd.md.bond.harmonic()
        if "bond_harmonic_energy" not in log_quantities:
            log_quantities.append("bond_harmonic_energy")
        for bond in coeffs_dict["bond_coeffs"]:
            harmonic_bond.bond_coeff.set(
                bond[0],
                k=bond[1] / (energy_scaling / distance_scaling ** 2),
                r0=bond[2] / distance_scaling,
            )

    if len(coeffs_dict["angle_coeffs"]) != 0:
        print("Loading harmonic angle coeffs...")
        harmonic_angle = hoomd.md.angle.harmonic()
        if "angle_harmonic_energy" not in log_quantities:
            log_quantities.append("angle_harmonic_energy")
        for angle in coeffs_dict["angle_coeffs"]:
            harmonic_angle.angle_coeff.set(
                angle[0], k=angle[1] / energy_scaling, t0=angle[2]
            )

    if len(coeffs_dict["dihedral_coeffs"]) != 0:
        print("Loading opls dihedral coeffs...")
        harmonic_dihedral = hoomd.md.dihedral.opls()
        if "dihedral_opls_energy" not in log_quantities:
            log_quantities.append("dihedral_opls_energy")
        for dihedral in coeffs_dict["dihedral_coeffs"]:
            harmonic_dihedral.dihedral_coeff.set(
                dihedral[0],
                k1=dihedral[1] / energy_scaling,
                k2=dihedral[2] / energy_scaling,
                k3=dihedral[3] / energy_scaling,
                k4=dihedral[4] / energy_scaling,
            )

    if len(coeffs_dict["external_forcefields"]) != 0:
        for forcefield_loc in coeffs_dict["external_forcefields"]:
            print("Loading external forcefield:", "".join([forcefield_loc, "..."]))
            if forcefield_loc[-7:] == ".eam.fs":
                hoomd.metal.pair.eam(file=forcefield_loc, type="FS", nlist=nl)
                if "pair_eam_energy" not in log_quantities:
                    log_quantities.append("pair_eam_energy")
            elif forcefield_loc[-10:] == ".eam.alloy":
                hoomd.metal.pair.eam(file=forcefield_loc, type="Alloy", nlist=nl)
                if "pair_eam_energy" not in log_quantities:
                    log_quantities.append("pair_eam_energy")
            else:
                print("----==== UNABLE TO PARSE EXTERNAL FORCEFIELD ====----")
                print(
                    forcefield_loc,
                    "is an unspecified file type and will be ignored."
                    " Please code in how to treat this file in rhaco/simulate.py",
                )
    for atomID, atom in enumerate(system.particles):
        atom.mass = coeffs_dict["mass"][atomID]
    # TODO: Support for charges
    # pppmnl = hoomd.md.nlist.cell()
    # pppm = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist = pppmnl)
    # pppm.set_params(Nx=64,Ny=64,Nz=64,order=6,rcut=2.70)
    return system, log_quantities


def get_coeffs(file_name, generate_arguments):
    coeff_dictionary = {
        "pair_coeffs": [],
        "bond_coeffs": [],
        "angle_coeffs": [],
        "dihedral_coeffs": [],
    }
    with open(file_name, "r") as xml_file:
        xml_data = ET.parse(xml_file)
    root = xml_data.getroot()
    for config in root:
        for child in config:
            # First get the masses which are different
            if child.tag == "mass":
                coeff_dictionary["mass"] = [
                    float(mass) for mass in child.text.split("\n") if len(mass) > 0
                ]
            # Now the other coefficients
            elif child.tag in coeff_dictionary.keys():
                if child.text is None:
                    continue
                for line in child.text.split("\n"):
                    if len(line) == 0:
                        continue
                    coeff = line.split()
                    coeff_dictionary[child.tag].append(
                        [coeff[0]] + list(map(float, coeff[1:]))
                    )
                # coeff_dictionary[child.tag] = child.text.split('\n')
    coeff_dictionary["external_forcefields"] = generate_arguments["forcefield"][1]
    return coeff_dictionary


def create_generate_arguments(file_name):
    print("------------============ WARNING ============------------")
    print("Your input file contains no generate arguments.")
    print("This likely means that it was created using Rhaco 1.2 or earlier.")
    print("The following assumptions will be made. If they do not correspond to your"
          " system, then please regenerate your morphologies with rhaco-create-morph"
          " using Rhaco 1.3 or later:")
    print("1) integrate_crystal is False")
    print("2) reactant_rigid is False")
    print("Additionally, forcefields will be read using the old syntax.")
    print("This feature will no longer be supported after Rhaco 1.5 is released.")
    generate_arguments = {
        "integrate_crystal": False,
        "reactant_rigid": False,
    }
    with open(file_name, "r") as xml_file:
        xml_data = ET.parse(xml_file)
    root = xml_data.getroot()
    for config in root:
        generate_arguments["forcefield"] = [[], [
            FF_data
            for FF_data in config.find("external_forcefields").text.split("\n")
            if len(FF_data) > 0
        ]]
        break
    return generate_arguments


def get_generate_arguments(file_name):
    with open(file_name, "r") as xml_file:
        xml_data = ET.parse(xml_file)
    root = xml_data.getroot()
    for config in root:
        try:
            all_args = [
                arg[:-1]
                for arg in config.find("generate_arguments").text.split("\n")
                if len(arg) > 0
            ]
            break
        except AttributeError:
            # Generate arguments does not exist (version < 1.3) 
            return None
    argument_dict = {}
    for argument in all_args:
        key = argument.split(": ")[0]
        val = ": ".join(argument.split(": ")[1:])
        # Many different value types, all as strings currently
        if val == "True":
            argument_dict[key] = True
        elif val == "False":
            argument_dict[key] = False
        elif val == "None":
            argument_dict[key] = None
        elif ("[" in val) or ("{" in val):
            argument_dict[key] = json.loads(val.replace("'", '"'))
        else:
            try:
                argument_dict[key] = int(val)
            except ValueError:
                try:
                    argument_dict[key] = float(val)
                except ValueError:
                    argument_dict[key] = str(val)
    return argument_dict


def rename_crystal_types(snapshot, generate_arguments):
    """
    This function splits the system into atoms to integrate over and atoms to
    not integrate over, based on the specified atom types ('X_<ATOM TYPE>')
    means the atom is part of the catalyst and should not be integrated over.

    Additionally, it removes the X_ from all the atom types so the forcefield
    can be interpreted correctly.
    """
    # Attempt 1: Just rename the types in snapshot.particles.types so that X_A
    # is renamed to just A.
    # Attempt 1.5: This can't be done on an index-by-index basis - the particle
    # types need to be assigned all in one go, so create a new list then update
    # the types data.
    catalyst_type_IDs = []
    new_types = []
    mapping = {}
    for type_index, atom_type in enumerate(snapshot.particles.types):
        if atom_type[:2] == "X_":
            new_atom_type = atom_type[2:]
            catalyst_type_IDs.append(type_index)
        else:
            new_atom_type = atom_type
        if new_atom_type not in new_types:
            new_types.append(new_atom_type)
        mapping[type_index] = new_types.index(new_atom_type)
    catalyst_atom_IDs = []
    for atom_index, type_ID in enumerate(snapshot.particles.typeid):
        if type_ID in catalyst_type_IDs:
            catalyst_atom_IDs.append(atom_index)
    catalyst = hoomd.group.tag_list(name="catalyst", tags=catalyst_atom_IDs)
    gas = hoomd.group.difference(name="gas", a=hoomd.group.all(), b=catalyst)
    # If we're not using rigid bodies in the reactant, then we can assign the surface
    # atoms to all have the same rigid_ID and get a speedup.
    if generate_arguments["reactant_rigid"] is False:
        snapshot.particles.body[catalyst_atom_IDs] = 0
    # If we are specifying an external forcefield (i.e. EAM), then we will need to
    # remove the X_ from the surface atom types otherwise it will break.
    if len(generate_arguments["forcefield"][1]) > 0:
        print("Renaming crystal atoms to remove the X_ for EAM...")
        snapshot.particles.types = new_types
    print("The catalyst group is", catalyst)
    print("The gas group is", gas)
    return snapshot, catalyst, gas


def create_rigid_bodies(file_name, snapshot, generate_arguments):
    # Firstly, return the snapshot if no bodies exist. Return under the following
    # conditions:
    # If integrate_crystal && 0 rigid IDs in system
    # If not integrate_crystal && 1 rigid ID in system
    rigid_IDs = set(snapshot.particles.body)
    # Discard -1 (== 4294967295 in uint32)
    rigid_IDs.discard(4294967295)
    if not generate_arguments["integrate_crystal"]:
        # If we assigned rigid body ID = 0 to the crystal and are not integrating over
        # it, we can remove it (so we don't have to set central particle/rotation etc.,
        # which causes an issue with self-interaction of the surface).
        rigid_IDs.discard(0)
    if len(rigid_IDs) == 0:
        return None, snapshot
    # Mbuild just updates the rigid body number, but rhaco-create-morph creates the
    # central particle and lists it first in the rigid body.
    # For HOOMD to honour the rigid specification, we just need
    # to:
    # 1) Obtain the relative positions of the constituent members of the body from the
    # input xml
    rigid_relative_positions = get_rigid_relative_positions(file_name)
    # 2) Assign a rotation quaternion
    # First, work out what the AAIDs are for the rigid bodies
    rigid_IDs = list(rigid_IDs)
    rigid_body_AAIDs = {}
    for AAID, body in enumerate(snapshot.particles.body):
        # Skip it if it's not a rigid body that we care about
        if body not in rigid_IDs:
            continue
        if body not in rigid_body_AAIDs:
            # This is the first element of a new rigid body
            rigid_body_AAIDs[body] = []
        rigid_body_AAIDs[body].append(AAID)
    # NOTE: This untested for multiple species of rigid bodies in the reactant
    rigid_type_assignments = {}
    all_rigid_body_types = []
    all_rigid_body_positions = []
    rolling_rigid_ID = 0
    for rigid_ID, AAIDs in rigid_body_AAIDs.items():
        # We can use the types to increment the rolling_rigid_ID number
        rigid_body_types = get_rigid_body_types(snapshot.particles, AAIDs)
        # Append our rigid body list with the important rigid body details
        if rigid_body_types not in all_rigid_body_types:
            # NEW RIGID BODY SPECIES
            all_rigid_body_types.append(rigid_body_types)
            # Check the species has the same number of rigid bodies as this one
            assert len(AAIDs) == len(rigid_relative_positions[rolling_rigid_ID])
            # Obtain the relative positions (skip the first element as it's the central
            # particle)
            rigid_body_positions = rigid_relative_positions[rolling_rigid_ID][1:]
            all_rigid_body_positions.append(rigid_body_positions)
            # Calculate the moment of inertia tensor and diagonalize it to get the three
            # components HOOMD needs
            moment_of_inertia = get_moment_of_inertia_tensor(
                rigid_body_positions, rigid_body_types
            )
            # Create the new rigid body type
            rigid_body_type_ID = "_R{}".format(rolling_rigid_ID)
            # Prepare for the next rigid body type (if any)
            rolling_rigid_ID += 1
        rigid_type_assignments[AAIDs[0]] = rigid_body_type_ID
        rotation_matrix = get_rotation_matrix(
            rigid_body_positions, snapshot.particles, AAIDs
        )
        body_rotation = quaternion.from_rotation_matrix(rotation_matrix)
        # Convert the quaternion to the data type that hoomd expects and then
        # assign the rotation to the central rigid_body particle
        snapshot.particles.orientation[AAIDs[0]] = np.array(
            body_rotation.components, dtype=np.float32
        )
        # And assign the moment of inertia we need
        snapshot.particles.moment_inertia[AAIDs[0]] = moment_of_inertia

    # Now assign the rigid atom types so constraint parameters can be set
    for central_particle_ID, atom_type in rigid_type_assignments.items():
        type_id = snapshot.particles.types.index(atom_type)
        snapshot.particles.typeid[central_particle_ID] = type_id
    # Now set the contraints
    rigid = hoomd.md.constrain.rigid()
    for rigid_type_index, rigid_body_types in enumerate(all_rigid_body_types):
        rigid.set_param(
            "_R{}".format(rigid_type_index),
            types=rigid_body_types,
            positions=all_rigid_body_positions[rigid_type_index],
        )

    # **ERROR**: constrain.rigid(): Central particles must have a body tag identical
    # to their contiguous tag.
    tag_is = 0
    change_to = 0
    for AAID, body_tag in enumerate(snapshot.particles.body):
        if AAID == 0:
            previous_body_tag = body_tag
            continue
        if body_tag != previous_body_tag:
            tag_is = body_tag
            change_to = AAID
            previous_body_tag = body_tag
        if body_tag == tag_is:
            snapshot.particles.body[AAID] = change_to

    for rigid_ID, AAIDs in rigid_body_AAIDs.items():
        central_ID = AAIDs[0]
        type_ID = snapshot.particles.typeid[central_ID]
        body_IDs = [snapshot.particles.body[AAID] for AAID in AAIDs[1:]]
        rigid_type_name = rigid_type_assignments[central_ID]
    return rigid, snapshot


def get_rigid_relative_positions(file_name):
    # The rigid relative positions are stored in the input xml:
    with open(file_name, "r") as xml_file:
        xml_tree = ET.parse(xml_file)
    root = xml_tree.getroot()
    for config in root:
        for child in config:
            if "rigid_relative_positions" in child.tag:
                rigid_relative_positions_unsplit = np.array(
                    [
                        np.fromiter(map(float, coords.split("\t")), dtype=np.float32)
                        for coords in child.text.split("\n")
                        if len(coords) > 0
                    ]
                )
    # New rigid bodies in rigid_relative_positions begin with a [0, 0, 0], so split the
    # array up around this element
    split_indices = []
    for index, row in enumerate(rigid_relative_positions_unsplit):
        if np.isclose(row, [0, 0, 0]).all():
            split_indices.append(index)
    split_array = np.split(rigid_relative_positions_unsplit, split_indices)
    # Since the rigid_relative_positions array begins with [0, 0, 0], we can skip the
    # first element of split_array
    rigid_relative_positions = []
    for array in split_array[1:]:
        rigid_relative_positions.append(array)
    return rigid_relative_positions


def get_moment_of_inertia_tensor(positions, types):
    # Calculate the 9 components of the moment_of_inertia_tensor
    # I11 = sum((yi^2 + zi^2) * mi)
    # I22 = sum((xi^2 + zi^2) * mi)
    # I33 = sum((xi^2 + yi^2) * mi)
    # I12 = I21 = -sum(xi * yi * mi)
    # I23 = I32 = -sum(yi * zi * mi)
    # I13 = I31 = -sum(xi * zi * mi)
    I11 = np.sum(
        (positions[:, 1] ** 2 + positions[:, 2] ** 2)
        * np.array([ATOM_MASSES[atom_type] for atom_type in types])
    )
    I22 = np.sum(
        (positions[:, 0] ** 2 + positions[:, 2] ** 2)
        * np.array([ATOM_MASSES[atom_type] for atom_type in types])
    )
    I33 = np.sum(
        (positions[:, 0] ** 2 + positions[:, 1] ** 2)
        * np.array([ATOM_MASSES[atom_type] for atom_type in types])
    )
    I12 = I21 = -np.sum(
        positions[:, 0]
        * positions[:, 1]
        * np.array([ATOM_MASSES[atom_type] for atom_type in types])
    )
    I23 = I32 = -np.sum(
        positions[:, 1]
        * positions[:, 2]
        * np.array([ATOM_MASSES[atom_type] for atom_type in types])
    )
    I13 = I31 = -np.sum(
        positions[:, 0]
        * positions[:, 2]
        * np.array([ATOM_MASSES[atom_type] for atom_type in types])
    )
    MI_tensor = np.array([[I11, I12, I13], [I21, I22, I23], [I31, I32, I33]])
    # Diagonalize the tensor (eigenvalues down lead diagonal)
    eigen = np.linalg.eig(MI_tensor)
    return eigen[0]


def get_rotation_matrix(unrotated_body, particles, AAIDs):
    # Translate rigid body to origin
    CoM = np.array(particles.position[AAIDs[0]])
    rotated_body = np.array([particles.position[AAID] for AAID in AAIDs[1:]]) - CoM
    # First, calculate the basis vectors of the unrotated body
    unrotated_A = unrotated_body[0]
    unrotated_B = unrotated_body[1]
    unrotated_C = unrotated_body[2]
    unrotated_AB = unrotated_B - unrotated_A
    unrotated_AC = unrotated_C - unrotated_A
    unrotated_BC = unrotated_C - unrotated_B
    # Basis vectors are whatever we want, so let's use unrotated_AB = x, normal = z,
    # and np.cross(x, z) = y
    x_axis = unrotated_AB / np.linalg.norm(unrotated_AB)
    normal = np.cross(unrotated_AB, unrotated_AC)
    z_axis = normal / np.linalg.norm(normal)
    y_axis = np.cross(x_axis, z_axis)
    # Normalising shouldn't be necessary but let's be safe.
    y_axis /= np.linalg.norm(y_axis)

    # # Now we have the basis set for our unrotated body. What is the basis set that
    # # describes the rotated body using those same definitions?
    rotated_A = rotated_body[0]
    rotated_B = rotated_body[1]
    rotated_C = rotated_body[2]
    rotated_AB = rotated_B - rotated_A
    rotated_AC = rotated_C - rotated_A
    rotated_BC = rotated_C - rotated_B
    # Same basis vectors as before
    x_prime_axis = rotated_AB / np.linalg.norm(rotated_AB)
    normal_prime = np.cross(rotated_AB, rotated_AC)
    z_prime_axis = normal_prime / np.linalg.norm(normal_prime)
    y_prime_axis = np.cross(x_prime_axis, z_prime_axis)
    y_prime_axis /= np.linalg.norm(y_prime_axis)

    # These vectors form the bases for our two coordinate systems
    original_basis = np.array([x_axis, y_axis, z_axis])
    new_basis = np.array([x_prime_axis, y_prime_axis, z_prime_axis])

    # Given that A' = R A, we now have a coupled system of 9 simultaneous equations we
    # can solve to obtain the rotation matrix
    rotation_matrix = np.linalg.solve(original_basis, new_basis).T
    return rotation_matrix


def get_rigid_body_types(particles, AAIDs):
    types_in_body = []
    for AAID in AAIDs[1:]:
        # Skip the first one because it's type "_R"
        types_in_body.append(particles.types[particles.typeid[AAID]])
    return types_in_body


def get_integrators(args, catalyst, gas, reduced_temperature):
    integrator_list = []
    rigid_gas = hoomd.group.intersection(
        name="rigid_gas", a=hoomd.group.rigid_center(), b=gas
    )
    flexible_gas = hoomd.group.difference(
        name="flexible_gas", a=gas, b=hoomd.group.rigid()
    )
    integrator_list = []
    # Integrate over rigid bodies in the gas
    if len(rigid_gas) > 0:
        integrator_list.append(
            hoomd.md.integrate.nvt(
                group=rigid_gas, tau=args.tau, kT=reduced_temperature
            )
        )
    # Integrate over flexible bodies in the gas
    if len(flexible_gas) > 0:
        integrator_list.append(
            hoomd.md.integrate.nvt(
                group=flexible_gas, tau=args.tau, kT=reduced_temperature
            )
        )
    # Also get the groups:
    groups_list = [group for group in [rigid_gas, flexible_gas] if len(group) > 0]
    return integrator_list, groups_list


def initialize_velocities(snapshot, temperature, gas):
    v = np.random.random((len(gas), 3))
    v -= 0.5
    meanv = np.mean(v, 0)
    meanv2 = np.mean(v ** 2, 0)
    fs = np.sqrt(temperature / meanv2)
    # Shift the velocities such that the average is zero
    v = v - meanv
    # Scale the velocities to match the required temperature
    v *= fs
    # Assign the velocities for this MD phase
    indices_to_modify = [atom.tag for atom in gas]
    for i, index in enumerate(indices_to_modify):
        snapshot.particles.velocity[index] = v[i]
    return snapshot


def main():
    parser = argparse.ArgumentParser(
        prog="rhaco-run-hoomd", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-T",
        "--temperature",
        type=float,
        default=633,
        required=False,
        help="""The desired temperature of the simulation in
                        kelvin (this will be rescaled to produce a reduced
                        temperature that conforms with Foyer's default
                        units (kcal/mol and angstroems).\n""",
    )
    parser.add_argument(
        "-r",
        "--run_time",
        type=float,
        default=1E7,
        required=False,
        help="""The number of timesteps to run the MD
                        simulation for.\n""",
    )
    parser.add_argument(
        "-s",
        "--timestep",
        type=float,
        default=1E-3,
        required=False,
        help="""The integration timestep to use when running
                        the NVT MD simulation.\n""",
    )
    parser.add_argument(
        "-t",
        "--tau",
        type=float,
        default=1E-2,
        required=False,
        help="""The thermostat coupling to use when running
                        the NVT MD simulation.\n""",
    )
    parser.add_argument(
        "-e",
        "--energy_scale_unit",
        type=float,
        default=1.0,
        required=False,
        help="""The energy scaling unit rhaco should use to
                        set the correct temperature, and LJ epsilons in kcal/mol. Default
                        is Foyer's default unit (1.0 kcal/mol). A useful alternative is
                        23.060541945329334 kcal/mol, which corresponds to 1 eV, which is
                        the energy scaling unit for EAM.
                        Be careful with this, it WILL frack everything up.""",
    )
    parser.add_argument(
        "-d",
        "--distance_scale_unit",
        type=float,
        default=1.0,
        required=False,
        help="""The distance scaling unit rhaco should use to
                        set the correct LJ sigmas in angstroems. Default is Foyer's
                        default unit (1.0 angstroem, same as EAM).
                        Be careful with this, it WILL frack everything up.""",
    )
    parser.add_argument(
        "-nl",
        "--nl_type",
        type=str,
        default="tree",
        required=False,
        help='''The neighbour list type to use. Default is tree,
                        other option is "cell"''',
    )
    parser.add_argument(
        "-rc",
        "--r_cut",
        type=float,
        default=None,
        required=False,
        help="""The r_cut value to use in the LJ interactions given
                        in distance_scale_unit. Default = 2.5 * max_lj_sigma""",
    )
    args, file_list = parser.parse_known_args()

    # Foyer gives parameters in terms of kcal/mol for energies and angstroems
    # for distances. Convert these to reduced units for HOOMD using the
    # following conversion, and print a string to inform the user it has been
    # done.
    reduced_temperature = (
        args.temperature * BOLTZMANN * AVOGADRO / (KCAL_TO_J * args.energy_scale_unit)
    )
    timestep_SI = args.timestep * np.sqrt(
        AMU_TO_KG
        * (ANG_TO_M * args.distance_scale_unit) ** 2
        * AVOGADRO
        / (KCAL_TO_J * args.energy_scale_unit)
    )
    print(
        "Using the units of <DISTANCE> =",
        args.distance_scale_unit,
        "Angstroem," " <ENERGY> =",
        args.energy_scale_unit,
        "kcal/mol, and <MASS> = 1 amu," " the input temperature of",
        args.temperature,
        "K corresponds to" " {:.2E}".format(reduced_temperature),
        "in dimensionless HOOMD kT units, and the input timestep",
        args.timestep,
        "corresponds to {:.2E} s.".format(timestep_SI),
    )

    for file_name in file_list:
        hoomd.context.initialize("")
        # hoomd.context.initialize("--notice-level=99", memory_traceback=True)

        system = hoomd.deprecated.init.read_xml(filename=file_name)
        generate_arguments = get_generate_arguments(file_name)
        if generate_arguments is None:
            generate_arguments = create_generate_arguments(file_name)

        # Sort out the rigid bodies (if they exist)
        snapshot = system.take_snapshot()
        rigid, updated_snapshot = create_rigid_bodies(
            file_name, snapshot, generate_arguments
        )
        system.restore_snapshot(updated_snapshot)
        if rigid is not None:
            rigid.validate_bodies()

        # Get the integration groups by ignoring anything that has the X_
        # prefix to the atom type, and rename the types for the forcefield
        snapshot = system.take_snapshot()
        updated_snapshot, catalyst, gas = rename_crystal_types(
            snapshot, generate_arguments
        )
        system.restore_snapshot(updated_snapshot)

        # Create the integrators
        hoomd.md.integrate.mode_standard(dt=args.timestep)
        integrator_list, groups_list = get_integrators(
            args, catalyst, gas, reduced_temperature
        )

        # Apply the forcefield coefficients
        system, log_quantities = set_coeffs(
            file_name,
            system,
            args.distance_scale_unit,
            args.energy_scale_unit,
            args.nl_type,
            args.r_cut,
            groups_list,
            generate_arguments,
        )
        try:
            # Use HOOMD 2.3's randomize_velocities
            for integrator in integrator_list:
                integrator.randomize_velocities(seed=42)
        except AttributeError:
            # Using a previous version of HOOMD - use the old initialization
            # function instead
            snapshot = system.take_snapshot()
            updated_snapshot = initialize_velocities(snapshot, reduced_temperature, gas)
            system.restore_snapshot(updated_snapshot)

        hoomd.dump.gsd(
            filename=".".join(file_name.split(".")[:-1]) + "_traj.gsd",
            period=max([int(args.run_time / 500), 1]),
            group=hoomd.group.all(),
            overwrite=True,
        )
        hoomd.analyze.log(
            filename=".".join(file_name.split(".")[:-1]) + ".log",
            quantities=log_quantities,
            period=max([int(args.run_time / 10000), 1]),
            header_prefix="#",
            overwrite=True,
        )
        ## Now incrementally ramp the charges
        # for chargePhase in range(chargeIncrements + 1):
        #    print("Incrementing charge phase", chargePhase, "of",
        #          chargeIncrements + 1)
        #    for atom in system.particles:
        #        oldCharge = copy.deepcopy(atom.charge)
        #        atom.charge = charges[atom.type] * (chargePhase /
        #                                            float(chargeIncrements))
        #    hoomd.run(chargeTimesteps)

        # hoomd.deprecated.dump.xml(
        #     group=hoomd.group.all(), filename="pre_run.xml", all=True
        # )

        # Get the initial box size dynamically
        hoomd.run_upto(args.run_time)
        hoomd.dump.gsd(
            filename=".".join(file_name.split(".")[:-1]) + "_final.gsd",
            period=None,
            group=hoomd.group.all(),
            overwrite=True,
        )
