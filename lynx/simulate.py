import argparse
import numpy as np
import hoomd
import hoomd.md
import hoomd.deprecated
import xml.etree.cElementTree as ET
import numpy as np

AVOGADRO = 6.022140857E23
BOLTZMANN = 1.38064852E-23
KCAL_TO_J = 4.184E3
AMU_TO_KG = 1.6605E-27
ANG_TO_M = 1E-10

def set_coeffs(file_name, system):
    '''
    Read in the molecular dynamics coefficients exported by Foyer
    '''
    coeffs_dict = get_coeffs(file_name)
    # Perylene: Sigma = 3.8 Ang, Epsilon = 0.1217 kCal/mol
    ljnl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=ljnl)
    lj.set_params(mode="xplor")
    for type1 in coeffs_dict['pair_coeffs']:
        for type2 in coeffs_dict['pair_coeffs']:
            lj.pair_coeff.set(type1[0], type2[0],
                              epsilon=np.sqrt(type1[1] * type2[1]),
                              sigma=np.sqrt(type1[2] * type2[2]))

    harmonic_bond = hoomd.md.bond.harmonic()
    for bond in coeffs_dict['bond_coeffs']:
        harmonic_bond.bond_coeff.set(bond[0], k=bond[1], r0=bond[2])

    harmonic_angle = hoomd.md.angle.harmonic()
    for angle in coeffs_dict['angle_coeffs']:
        harmonic_angle.angle_coeff.set(angle[0], k=angle[1], t0=angle[2])

    harmonic_dihedral = hoomd.md.dihedral.opls()
    for dihedral in coeffs_dict['dihedral_coeffs']:
        harmonic_dihedral.dihedral_coeff.set(dihedral[0], k1=dihedral[1],
                                             k2=dihedral[2], k3=dihedral[3],
                                             k4=dihedral[4])

    for atomID, atom in enumerate(system.particles):
        atom.mass = coeffs_dict['mass'][atomID]

    # TODO: Support for charges
    #pppmnl = hoomd.md.nlist.cell()
    #pppm = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist = pppmnl)
    #pppm.set_params(Nx=64,Ny=64,Nz=64,order=6,rcut=2.70)

    return system



def get_coeffs(file_name):
    coeff_dictionary = {'pair_coeffs': [], 'bond_coeffs': [],
                        'angle_coeffs': [], 'dihedral_coeffs': []}
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
            # Now the other coefficients
            elif child.tag in coeff_dictionary.keys():
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
    for type_index, atom_type in enumerate(snapshot.particles.types):
        if atom_type[:2] == 'X_':
            new_types.append(atom_type[2:])
            catalyst_type_IDs.append(type_index)
        else:
            new_types.append(atom_type)
    snapshot.particles.types = new_types
    catalyst_atom_IDs = []
    for atom_index, type_ID in enumerate(snapshot.particles.typeid):
        if type_ID in catalyst_type_IDs:
            catalyst_atom_IDs.append(atom_index)
    catalyst = hoomd.group.tag_list(name='catalyst', tags=catalyst_atom_IDs)
    gas = hoomd.group.difference(name='gas', a=hoomd.group.all(), b=catalyst)
    print("The catalyst group is", catalyst)
    print("The gas group is", gas)
    return snapshot, catalyst, gas


def initialize_velocities(snapshot, temperature):
    v = np.random.random((len(snapshot.particles.velocity), 3))
    v -= 0.5
    meanv = np.mean(v, 0)
    meanv2 = np.mean(v ** 2, 0)
    fs = np.sqrt(temperature / meanv2)
    # Shift the velocities such that the average is zero
    v = (v - meanv)
    # Scale the velocities to match the required temperature
    v *= fs
    # Assign the velocities for this MD phase
    snapshot.particles.velocity[:] = v[:]
    return snapshot


def main():
    parser = argparse.ArgumentParser(prog='lynx-run-hoomd',
                                    formatter_class=argparse.
                                    ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--temperature',
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
    args, file_list = parser.parse_known_args()

    # Foyer gives parameters in terms of kcal/mol for energies and angstroems
    # for distances. Convert these to reduced units for HOOMD using the
    # following conversion, and print a string to inform the user it has been
    # done.
    reduced_temperature = args.temperature * BOLTZMANN * AVOGADRO / KCAL_TO_J
    timestep_SI = args.timestep * np.sqrt(AMU_TO_KG * ANG_TO_M**2 * AVOGADRO
                                          / KCAL_TO_J)
    print("Using the Foyer default units of <DISTANCE> = 1 Angstroem,"
          " <ENERGY> = 1 kcal/mol, and <MASS> = 1 amu, the input temperature"
          " of", args.temperature, "K corresponds to"
          " {:.2E}".format(reduced_temperature),
          "in dimensionless HOOMD kT units, and the input timestep",
          args.timestep, "corresponds to {:.2E} s.".format(timestep_SI))

    for file_name in file_list:
        hoomd.context.initialize("")
        system = hoomd.deprecated.init.read_xml(filename=file_name)
        snapshot = system.take_snapshot()
        # Get the integration groups by ignoring anything that has the X_
        # prefix to the atom type, and rename the types for the forcefield
        renamed_snapshot, catalyst, gas = rename_types(snapshot)
        # Assign the required velocities based on the requested temperature
        initialized_snapshot = initialize_velocities(renamed_snapshot,
                                                     reduced_temperature)
        # Finally, restore the snapshot
        system.restore_snapshot(initialized_snapshot)
        system = set_coeffs(file_name, system)

        hoomd.md.integrate.mode_standard(dt=args.timestep);
        integrator = hoomd.md.integrate.nvt(group=gas, tau=0.1,
                                            kT=reduced_temperature)

        hoomd.dump.gsd(filename=".".join(file_name.split(".")[:-1])
                       + "_traj.gsd", period=max([int(args.run_time/500), 1]),
                       group=hoomd.group.all(), overwrite=True)
        log_quantities = ['temperature', 'pressure', 'volume',
                          'potential_energy', 'kinetic_energy',
                          'pair_lj_energy', 'bond_harmonic_energy',
                          'angle_harmonic_energy', 'dihedral_opls_energy']
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
