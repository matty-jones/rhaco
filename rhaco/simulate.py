import argparse
import numpy as np
import hoomd
import hoomd.md
import hoomd.metal
import hoomd.deprecated
import xml.etree.cElementTree as ET
import numpy as np


AVOGADRO = 6.022140857E23
BOLTZMANN = 1.38064852E-23
KCAL_TO_J = 4.184E3
AMU_TO_KG = 1.6605E-27
ANG_TO_M = 1E-10


def parse_interactions(omit_string):
    omit_list = []
    if (omit_string[0] == "[") and (omit_string[-1] == "]"):
        omit_string = omit_string[1:-1]
    if "," in omit_string:
        omit_string = "".join(omit_string.split(" "))
        omit_string = omit_string.split(",")
    else:
        omit_string = omit_string.split()
    for interaction in omit_string:
        if (interaction[0] == "'") and (interaction[-1] == "'"):
            omit_list.append(interaction[1:-1])
        elif (interaction[0] == '"') and (interaction[-1] == '"'):
            omit_list.append(interaction[1:-1])
        else:
            omit_list.append(interaction)
    return omit_list


def set_coeffs(
    file_name, system, omit_lj, distance_scaling, energy_scaling,
    nl_type, r_cut
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
    nl.reset_exclusions(exclusions=["body"])
    log_quantities = ['temperature', 'temperature_gas', 'pressure', 'volume',
                      'potential_energy', 'potential_energy_gas', 'kinetic_energy',
                      'kinetic_energy_gas']
    if len(coeffs_dict['pair_coeffs']) != 0:
        print("Loading LJ pair coeffs...")
        if r_cut is None:
            max_sigma = np.max(list(map(
                float, np.array(coeffs_dict['pair_coeffs'])[:,2]
            )))
            r_cut = 2.5 * max_sigma
            print("Setting r_cut to 2.5 * max_lj_sigma =", r_cut)
        lj = hoomd.md.pair.lj(r_cut=r_cut, nlist=nl)
        lj.set_params(mode="xplor")
        if "pair_lj_energy" not in log_quantities:
            log_quantities.append('pair_lj_energy')
        for type1 in coeffs_dict['pair_coeffs']:
            for type2 in coeffs_dict['pair_coeffs']:
                if "-".join([type1[0], type2[0]]) in omit_lj:
                    print("Omitting", "-".join([type1[0], type2[0]]), "0", "0")
                    lj.pair_coeff.set(type1[0], type2[0],
                                      epsilon=0.0,
                                      sigma=0.0)
                else:
                    lj.pair_coeff.set(type1[0], type2[0],
                                      epsilon=np.sqrt(type1[1] * type2[1]) / energy_scaling,
                                      sigma=np.sqrt(type1[2] * type2[2]) / distance_scaling,
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


def rename_types(snapshot):
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
        if new_atom_type in new_types:
            mapping[type_index] = new_types.index(new_atom_type)
        else:
            mapping[type_index] = type_index
            new_types.append(new_atom_type)
    catalyst_atom_IDs = []
    for atom_index, type_ID in enumerate(snapshot.particles.typeid):
        if type_ID in catalyst_type_IDs:
            catalyst_atom_IDs.append(atom_index)
    # Also, add the surface atom to the same rigid body
    snapshot.particles.body[catalyst_atom_IDs] = 0
    catalyst = hoomd.group.tag_list(name='catalyst', tags=catalyst_atom_IDs)
    gas = hoomd.group.difference(name='gas', a=hoomd.group.all(), b=catalyst)
    # Now use the mapping to remove any duplicate types (needed if the same atom
    # type is present in both the crystal and the reactant)
    snapshot.particles.types = new_types
    for AAID, old_type in enumerate(snapshot.particles.typeid):
        snapshot.particles.typeid[AAID] = mapping[old_type]
    print("The catalyst group is", catalyst)
    print("The gas group is", gas)
    return snapshot, catalyst, gas


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
    parser.add_argument('-o', '--omit_lj',
                        type=parse_interactions,
                        default=[],
                        required=False,
                        help='''A list of lj interactions to omit from the
                        rhaco-generated input hoomdxml (useful when using EAM).\n
                        If unspecified, all interactions are considered''')
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
        # Get the integration groups by ignoring anything that has the X_
        # prefix to the atom type, and rename the types for the forcefield
        system = hoomd.deprecated.init.read_xml(filename=file_name)
        snapshot = system.take_snapshot()
        renamed_snapshot, catalyst, gas = rename_types(snapshot)
        # Then, restore the snapshot
        system.restore_snapshot(renamed_snapshot)
        system, log_quantities = set_coeffs(file_name, system, args.omit_lj, args.distance_scale_unit, args.energy_scale_unit,
                                           args.nl_type, args.r_cut)

        hoomd.md.integrate.mode_standard(dt=args.timestep);
        integrator = hoomd.md.integrate.nvt(group=gas, tau=args.tau,
                                            kT=reduced_temperature)
        try:
            # Use HOOMD 2.3's randomize_velocities
            integrator.randomize_velocities(seed=42)
        except AttributeError:
            # Using a previous version of HOOMD - use the old initialization
            # function instead
            snapshot = system.take_snapshot()
            initialized_snapshot = initialize_velocities(snapshot, reduced_temperature, gas)
            system.restore_snapshot(initialized_snapshot)

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
