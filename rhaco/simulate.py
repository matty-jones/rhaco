import argparse
import copy
import hoomd
import hoomd.md
import hoomd.metal
import hoomd.deprecated
import numpy as np
import quaternion
import xml.etree.cElementTree as ET


AVOGADRO = 6.022140857E23
BOLTZMANN = 1.38064852E-23
KCAL_TO_J = 4.184E3
AMU_TO_KG = 1.6605E-27
ANG_TO_M = 1E-10


def set_coeffs(
    file_name, system, distance_scaling, energy_scaling,
    nl_type, r_cut, groups_list,
):
    '''
    Read in the molecular dynamics coefficients exported by Foyer
    '''
    coeffs_dict = get_coeffs(file_name)
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
    if len(coeffs_dict['pair_coeffs']) != 0:
        print("Loading LJ pair coeffs...")
        if r_cut is None:
            max_sigma = np.max(list(map(
                float, np.array(coeffs_dict['pair_coeffs'])[:,2]
            )))
            r_cut = 2.5 * max_sigma
            print("Setting r_cut to 2.5 * max_sigma =", r_cut)
        lj = hoomd.md.pair.lj(r_cut=r_cut, nlist=nl)
        lj.set_params(mode="xplor")
        if "pair_lj_energy" not in log_quantities:
            log_quantities.append('pair_lj_energy')
        for type1 in coeffs_dict['pair_coeffs']:
            for type2 in coeffs_dict['pair_coeffs']:
                lj.pair_coeff.set(type1[0], type2[0],
                                  epsilon=np.sqrt(type1[1] * type2[1]) / (energy_scaling),
                                  sigma=np.sqrt(type1[2] * type2[2]) / (distance_scaling),
                                 )
    if len(coeffs_dict['bond_coeffs']) != 0:
        print("Loading harmonic bond coeffs...")
        harmonic_bond = hoomd.md.bond.harmonic()
        if "bond_harmonic_energy" not in log_quantities:
            log_quantities.append('bond_harmonic_energy')
        for bond in coeffs_dict['bond_coeffs']:
            harmonic_bond.bond_coeff.set(bond[0],
                                         k=bond[1] / (energy_scaling / distance_scaling**2),
                                         r0=bond[2] / distance_scaling,
                                        )

    if len(coeffs_dict['angle_coeffs']) != 0:
        print("Loading harmonic angle coeffs...")
        harmonic_angle = hoomd.md.angle.harmonic()
        if "angle_harmonic_energy" not in log_quantities:
            log_quantities.append('angle_harmonic_energy')
        for angle in coeffs_dict['angle_coeffs']:
            harmonic_angle.angle_coeff.set(angle[0], k=angle[1] / energy_scaling, t0=angle[2])

    if len(coeffs_dict['dihedral_coeffs']) != 0:
        print("Loading opls dihedral coeffs...")
        harmonic_dihedral = hoomd.md.dihedral.opls()
        if "dihedral_opls_energy" not in log_quantities:
            log_quantities.append('dihedral_opls_energy')
        for dihedral in coeffs_dict['dihedral_coeffs']:
            harmonic_dihedral.dihedral_coeff.set(dihedral[0],
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
                    log_quantities.append('pair_eam_energy')
            elif forcefield_loc[-10:] == ".eam.alloy":
                hoomd.metal.pair.eam(file=forcefield_loc, type="Alloy", nlist=nl)
                if "pair_eam_energy" not in log_quantities:
                    log_quantities.append('pair_eam_energy')
            else:
                print("----==== UNABLE TO PARSE EXTERNAL FORCEFIELD ====----")
                print(forcefield_loc, "is an unspecified file type and will be ignored."
                      " Please code in how to treat this file in rhaco/simulate.py")
    for atomID, atom in enumerate(system.particles):
        atom.mass = coeffs_dict['mass'][atomID]
    # TODO: Support for charges
    #pppmnl = hoomd.md.nlist.cell()
    #pppm = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist = pppmnl)
    #pppm.set_params(Nx=64,Ny=64,Nz=64,order=6,rcut=2.70)
    return system, log_quantities


def get_coeffs(file_name):
    coeff_dictionary = {'pair_coeffs': [], 'bond_coeffs': [],
                        'angle_coeffs': [], 'dihedral_coeffs': [],
                        'external_forcefields': []}
    with open(file_name, 'r') as xml_file:
        xml_data = ET.parse(xml_file)
    root = xml_data.getroot()
    for config in root:
        for child in config:
            # First get the masses which are different
            if child.tag == 'mass':
                coeff_dictionary['mass'] = [float(mass) for mass in
                                            child.text.split('\n') if
                                            len(mass) > 0]
            # Secondly, get the external forcefields, which are also different
            elif child.tag == 'external_forcefields':
                for line in child.text.split('\n'):
                    if len(line) > 0:
                        coeff_dictionary[child.tag].append(line)
            # Now the other coefficients
            elif child.tag in coeff_dictionary.keys():
                if child.text is None:
                    continue
                for line in child.text.split('\n'):
                    if len(line) == 0:
                        continue
                    coeff = line.split()
                    coeff_dictionary[child.tag].append(
                        [coeff[0]] + list(map(float, coeff[1:])))
                # coeff_dictionary[child.tag] = child.text.split('\n')
    return coeff_dictionary


def rename_crystal_types(snapshot):
    '''
    This function splits the system into atoms to integrate over and atoms to
    not integrate over, based on the specified atom types ('X_<ATOM TYPE>')
    means the atom is part of the catalyst and should not be integrated over.

    Additionally, it removes the X_ from all the atom types so the forcefield
    can be interpreted correctly.
    '''
    # Attempt 1: Just rename the types in snapshot.particles.types so that X_A
    # is renamed to just A.
    # Attempt 1.5: This can't be done on an index-by-index basis - the particle
    # types need to be assigned all in one go, so create a new list then update
    # the types data.
    catalyst_type_IDs = []
    new_types = []
    mapping = {}
    for type_index, atom_type in enumerate(snapshot.particles.types):
        if atom_type[:2] == 'X_':
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
    catalyst = hoomd.group.tag_list(name='catalyst', tags=catalyst_atom_IDs)
    gas = hoomd.group.difference(name='gas', a=hoomd.group.all(), b=catalyst)
    # Now use the mapping to remove any duplicate types (needed if the same atom
    # type is present in both the crystal and the reactant)
    snapshot.particles.types = new_types
    for AAID, old_type in enumerate(snapshot.particles.typeid):
        snapshot.particles.typeid[AAID] = mapping[old_type]
    # Finally, add the surface atoms to the same rigid body
    # First get the number of rigid bodies already in the system
    if len(set(snapshot.particles.body)) > 1:
        # Rigid bodies already defined, so 0 is already taken. Increment all rigid
        # bodies by one and then we can set the surface to be zero.
        # Skip rigid body 4294967295 (== -1 for uint32), as these are flexible

        # NOTE: Snapshots don't allow array assignment, gonna do it elementwise instead
        # mask = (snapshot.particles.body != 4294967295).astype(np.uint32)
        # incremented_body = snapshot.particles.body + mask
        # snapshot.particles.body = incremented_body
        for index, rigid_ID in np.ndenumerate(snapshot.particles.body):
            if rigid_ID == 4294967295:
                continue
            else:
                snapshot.particles.body[index] = rigid_ID + 1
    snapshot.particles.body[catalyst_atom_IDs] = 0
    print("The catalyst group is", catalyst)
    print("The gas group is", gas)
    return snapshot, catalyst, gas


def create_rigid_bodies(snapshot):
    # Firstly, return the snapshot if no bodies exist
    # Should be one rigid body for the crystal, and then one (-1) for the flexible
    # reactants == 2 total, or just one (-1) if we are integrating the surface.
    rigid_IDs = set(snapshot.particles.body)
    # Discard -1 (== 4294967295 in uint32)
    rigid_IDs.discard(4294967295)
    if len(set(snapshot.particles.body)) <= 1:
        return snapshot
    # Mbuild just updates the rigid body number, but rhaco-create-morph creates the 
    # central particle and lists it first in the rigid body.
    # For HOOMD to honour the rigid specification, we just need
    # to:
    # 1) Describe the relative positions of the constituent members of the body
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
    # NOTE: This will fail if we have multiple species of rigid bodies
    rigid_type_assignments = {}
    all_rigid_body_types = []
    all_rigid_body_positions = []
    rolling_rigid_ID = 0
    # 2) Assign a rotation quaternion
    for rigid_ID, AAIDs in rigid_body_AAIDs.items():
        rigid_body_types = get_rigid_body_types(
            snapshot.particles, AAIDs
        )
        rigid_body_positions = get_rigid_body_positions(
            snapshot.particles, AAIDs
        )
        # Append our rigid body list with the important rigid body details
        if rigid_body_types not in all_rigid_body_types:
            all_rigid_body_types.append(rigid_body_types)
            all_rigid_body_positions.append(rigid_body_positions)
            rigid_body_type_ID = "R{}".format(rolling_rigid_ID)
            rolling_rigid_ID += 1
        rigid_type_assignments[AAIDs[0]] = rigid_body_type_ID
        euler_a, euler_b, euler_c = get_euler_angles(snapshot.particles, AAIDs)
        euler_a = euler_b = euler_c = 0.0
        body_rotation = quaternion.from_euler_angles(
            euler_a, beta=euler_b, gamma=euler_c
        )
        print(euler_a, euler_b, euler_c)
        print(body_rotation)

        # # Convert the quaternion to the format that hoomd expects
        # # Assign rotation to the central rigid_body particle
        # snapshot.particles.orientation[AAIDs[0]] = np.array(
        #     body_rotation.components, dtype=np.float32
        # )

        # Update the central atom type to include the rigid body species number
        #snapshot.particles.types = snapshot.particles.types + [rigid_body_type_ID]
        #print(snapshot.particles.types)
        #for AAID in AAIDs:
        #    snapshot.particles.typeid[AAID] = len(snapshot.particles.types)
        #print(snapshot.particles.types)
        #exit()

    # Now assign the rigid atom types so constraint parameters can be set
    for central_particle_ID, atom_type in rigid_type_assignments.items():
        type_id = snapshot.particles.types.index(atom_type)
        snapshot.particles.typeid[central_particle_ID] = type_id
    # Now set the contraints
    rigid = hoomd.md.constrain.rigid()
    for rigid_type_index, rigid_body_types in enumerate(all_rigid_body_types):
        rigid.set_param(
            "R{}".format(rigid_type_index),
            types=rigid_body_types,
            positions=all_rigid_body_positions[rigid_type_index],
        )

    # # **ERROR**: constrain.rigid(): Central particles must have a body tag identical
    # # to their contiguous tag.
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
        print(rigid_type_name,
              snapshot.particles.types[type_ID],
              snapshot.particles.body[central_ID],
              snapshot.particles.orientation[central_ID],
              len(body_IDs),
              len(all_rigid_body_positions[int(rigid_type_name.split("R")[-1])]),
             )
    return rigid, snapshot


def get_euler_angles(particles, AAIDs):
    # Translate rigid body to origin
    CoM = np.array(particles.position[AAIDs[0]])
    particle_positions = np.array(
        [particles.position[AAID] for AAID in AAIDs[1:]]
    ) - CoM
    # Use three positions for the rest of the calculation: particle_positions[0] = A,
    # particle_positions[1] = B, particle_positions[2] = C
    pos_A = particle_positions[0]
    pos_B = particle_positions[1]
    pos_C = particle_positions[2]
    vec_AB = pos_B - pos_A
    vec_AC = pos_C - pos_A
    vec_BC = pos_C - pos_B
    # And the unit normal
    normal = np.cross(vec_AB, vec_AC)
    # Normalize to get the "z axis" of the body
    z_axis = normal / np.linalg.norm(normal)
    # Define the "x axis" of the body as the vector starting at C and passing through
    # the midpoint of AB
    x_axis = ((pos_B + pos_A) / 2.0) - pos_C
    x_axis /= np.linalg.norm(x_axis)
    # The y axis is the cross of the x and z axes
    y_axis = np.cross(z_axis, x_axis)
    y_axis /= np.linalg.norm(y_axis)
    # Calculating euler angles using the x-convention:
    alpha = np.arctan2(-z_axis[1], -z_axis[2])
    beta = np.arcsin(z_axis[0])
    gamma = np.arctan2(-y_axis[0], x_axis[0])
    return alpha, beta, gamma


def get_rigid_body_types(particles, AAIDs):
    types_in_body = []
    for AAID in AAIDs[1:]:
        # Skip the first one because it's type "R"
        types_in_body.append(particles.types[particles.typeid[AAID]])
    return types_in_body


def get_rigid_body_positions(particles, AAIDs):
    # Perform coordinate transformation s.t.:
    # 1) Rigid body is at origin
    CoM = np.array(particles.position[AAIDs[0]])
    particle_positions = np.array(
        [particles.position[AAID] for AAID in AAIDs[1:]]
    ) - CoM
    # 2) particle[0], particle[1], particle[2] form plane describing rigid body and
    # has normal [0, 0, 1]
    pos_A = particle_positions[0]
    pos_B = particle_positions[1]
    pos_C = particle_positions[2]
    vec_AB = pos_B - pos_A
    vec_AC = pos_C - pos_A
    vec_BC = pos_C - pos_B
    normal = np.cross(vec_AB, vec_AC)
    normal /= np.linalg.norm(normal)
    # Calculate the rotation matrix required to map the normal onto the baseline
    # (0, 0, 1)
    rotation_matrix = get_rotation_matrix(normal, np.array([0, 0, 1]))
    # Now the rotation matrix is set up, rotate all atoms (positions at origin)
    rotated_positions = np.matmul(rotation_matrix, particle_positions.T).T
    return np.array(rotated_positions, dtype=np.float32)


def get_rotation_matrix(vector1, vector2):
    """
    This function returns the rotation matrix around the origin that maps vector1 to
    vector 2
    """
    cross_product = np.cross(vector1, vector2)
    sin_angle = np.sqrt(
        (
            (cross_product[0] ** 2)
            + ((cross_product[1]) ** 2)
            + ((cross_product[2]) ** 2)
        )
    )
    cos_angle = np.dot(vector1, vector2)
    skew_matrix = np.array(
        [
            [0, -cross_product[2], cross_product[1]],
            [cross_product[2], 0, -cross_product[0]],
            [-cross_product[1], cross_product[0], 0],
        ]
    )
    skew_matrix_squared = skew_matrix @ skew_matrix
    rot_matrix = (
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        + skew_matrix
        + skew_matrix_squared * ((1 - cos_angle) / (sin_angle ** 2))
    )
    return rot_matrix


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
                group=rigid_gas, tau=args.tau, kT=reduced_temperature,
            )
        )
    # Integrate over flexible bodies in the gas
    if len(flexible_gas) > 0:
        integrator_list.append(
            hoomd.md.integrate.nvt(
                group=flexible_gas, tau=args.tau, kT=reduced_temperature,
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
    v = (v - meanv)
    # Scale the velocities to match the required temperature
    v *= fs
    # Assign the velocities for this MD phase
    indices_to_modify = [atom.tag for atom in gas]
    for i, index in enumerate(indices_to_modify):
        snapshot.particles.velocity[index] = v[i]
    return snapshot


def main():
    parser = argparse.ArgumentParser(prog='rhaco-run-hoomd',
                                    formatter_class=argparse.
                                    ArgumentDefaultsHelpFormatter)
    parser.add_argument('-T', '--temperature',
                        type=float,
                        default=633,
                        required=False,
                        help='''The desired temperature of the simulation in
                        kelvin (this will be rescaled to produce a reduced
                        temperature that conforms with Foyer's default
                        units (kcal/mol and angstroems).\n''')
    parser.add_argument('-r', '--run_time',
                        type=float,
                        default=1E7,
                        required=False,
                        help='''The number of timesteps to run the MD
                        simulation for.\n''')
    parser.add_argument('-s', '--timestep',
                        type=float,
                        default=1E-3,
                        required=False,
                        help='''The integration timestep to use when running
                        the NVT MD simulation.\n''')
    parser.add_argument('-t', '--tau',
                        type=float,
                        default=1E-2,
                        required=False,
                        help='''The thermostat coupling to use when running
                        the NVT MD simulation.\n''')
    parser.add_argument('-e', '--energy_scale_unit',
                        type=float,
                        default=1.0,
                        required=False,
                        help='''The energy scaling unit rhaco should use to
                        set the correct temperature, and LJ epsilons in kcal/mol. Default
                        is Foyer's default unit (1.0 kcal/mol). A useful alternative is
                        23.060541945329334 kcal/mol, which corresponds to 1 eV, which is
                        the energy scaling unit for EAM.
                        Be careful with this, it WILL frack everything up.''')
    parser.add_argument('-d', '--distance_scale_unit',
                        type=float,
                        default=1.0,
                        required=False,
                        help='''The distance scaling unit rhaco should use to
                        set the correct LJ sigmas in angstroems. Default is Foyer's
                        default unit (1.0 angstroem, same as EAM).
                        Be careful with this, it WILL frack everything up.''')
    parser.add_argument('-nl', '--nl_type',
                        type=str,
                        default="tree",
                        required=False,
                        help='''The neighbour list type to use. Default is tree,
                        other option is "cell"''')
    parser.add_argument('-rc', '--r_cut',
                        type=float,
                        default=None,
                        required=False,
                        help='''The r_cut value to use in the LJ interactions given
                        in distance_scale_unit. Default = 2.5 * max_lj_sigma''')
    args, file_list = parser.parse_known_args()

    # Foyer gives parameters in terms of kcal/mol for energies and angstroems
    # for distances. Convert these to reduced units for HOOMD using the
    # following conversion, and print a string to inform the user it has been
    # done.
    reduced_temperature = args.temperature * BOLTZMANN * AVOGADRO / (KCAL_TO_J * args.energy_scale_unit)
    timestep_SI = args.timestep * np.sqrt(AMU_TO_KG * (ANG_TO_M * args.distance_scale_unit)**2 * AVOGADRO
                                          / (KCAL_TO_J * args.energy_scale_unit))
    print("Using the units of <DISTANCE> =", args.distance_scale_unit, "Angstroem,"
          " <ENERGY> =", args.energy_scale_unit, "kcal/mol, and <MASS> = 1 amu,"
          " the input temperature of", args.temperature, "K corresponds to"
          " {:.2E}".format(reduced_temperature),
          "in dimensionless HOOMD kT units, and the input timestep",
          args.timestep, "corresponds to {:.2E} s.".format(timestep_SI))

    for file_name in file_list:
        hoomd.context.initialize("")
        # hoomd.context.initialize("--notice-level=99", memory_traceback=True)

        # Get the integration groups by ignoring anything that has the X_
        # prefix to the atom type, and rename the types for the forcefield
        system = hoomd.deprecated.init.read_xml(filename=file_name)
        snapshot = system.take_snapshot()
        updated_snapshot, catalyst, gas = rename_crystal_types(snapshot)
        system.restore_snapshot(updated_snapshot)
        hoomd.deprecated.dump.xml(group=hoomd.group.all(), filename="pre_rigid_bodies.xml", all=True)

        # Sort out any rigid bodies (if they exist)
        snapshot = system.take_snapshot()
        rigid, updated_snapshot = create_rigid_bodies(snapshot)
        system.restore_snapshot(updated_snapshot)
        rigid.validate_bodies()
        hoomd.deprecated.dump.xml(group=hoomd.group.all(), filename="post_rigid_bodies.xml", all=True)

        # Create the integrators
        hoomd.md.integrate.mode_standard(dt=args.timestep);
        integrator_list, groups_list = get_integrators(
            args, catalyst, gas, reduced_temperature
        )

        # Apply the forcefield coefficients
        system, log_quantities = set_coeffs(
            file_name, system, args.distance_scale_unit, args.energy_scale_unit,
            args.nl_type, args.r_cut, groups_list
        )
        try:
            # Use HOOMD 2.3's randomize_velocities
            for integrator in integrator_list:
                integrator.randomize_velocities(seed=42)
        except AttributeError:
            # Using a previous version of HOOMD - use the old initialization
            # function instead
            snapshot = system.take_snapshot()
            updated_snapshot = initialize_velocities(
                snapshot, reduced_temperature, gas
            )
            system.restore_snapshot(updated_snapshot)

        hoomd.deprecated.dump.xml(group=hoomd.group.all(), filename="post_initialization.xml", all=True)
        exit()
        hoomd.dump.gsd(filename=".".join(file_name.split(".")[:-1])
                       + "_traj.gsd", period=max([int(args.run_time/500), 1]),
                       group=hoomd.group.all(), overwrite=True)
        hoomd.analyze.log(filename='.'.join(file_name.split('.')[:-1])
                          + ".log", quantities=log_quantities,
                          period=max([int(args.run_time/10000), 1]),
                          header_prefix='#', overwrite=True)
        ## Now incrementally ramp the charges
        #for chargePhase in range(chargeIncrements + 1):
        #    print("Incrementing charge phase", chargePhase, "of",
        #          chargeIncrements + 1)
        #    for atom in system.particles:
        #        oldCharge = copy.deepcopy(atom.charge)
        #        atom.charge = charges[atom.type] * (chargePhase /
        #                                            float(chargeIncrements))
        #    hoomd.run(chargeTimesteps)

        # Get the initial box size dynamically
        hoomd.run_upto(args.run_time)
        hoomd.dump.gsd(filename=".".join(file_name.split(".")[:-1])
                       + "_final.gsd", period=None,
                       group=hoomd.group.all(), overwrite=True)
