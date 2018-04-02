import argparse
import numpy as np
import hoomd
import hoomd.md
import hoomd.deprecated
import xml.etree.cElementTree as ET


def set_coeffs(file_name):
    '''
    Read in the molecular dynamics coefficients exported by Foyer
    '''
    get_coeffs(file_name)
    # Perylene: Sigma = 3.8 Ang, Epsilon = 0.1217 kCal/mol
    ljnl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=ljnl)
    lj.pair_coeff.set('CP', 'CP', epsilon=1.0, sigma=1.0)

    #pppmnl = hoomd.md.nlist.cell()
    #pppm = hoomd.md.charge.pppm(group=hoomd.group.charged(), nlist = pppmnl)
    #pppm.set_params(Nx=64,Ny=64,Nz=64,order=6,rcut=2.70)

    harmonic_bond = hoomd.md.bond.harmonic()
    harmonic_bond.bond_coeff.set('CP-CP', k=30000.0, r0=0.4)

    #harmonic_angle = hoomd.md.angle.harmonic()
    #harmonic_angle.angle_coeff.set('CP-CP-CP', k=380.0, t0=2.09)

    #harmonic_dihedral = hoomd.md.dihedral.opls()
    #harmonic_dihedral.dihedral_coeff.set('CP-CP-CP-CP', k1=0.0, k2=50.0, k3=0.0, k4=0.0)


def get_coeffs(file_name):
    coeff_dictionary = {}
    constraint_types = ['pair', 'bond', 'angle', 'dihedral']
    with open(file_name, 'r') as xml_file:
        xml_data = ET.parse(xml_file)
    for constraint_type in constraint_types:
        data = xml_data.find(constraint_type + '_coeffs').text
        print(data)
        exit()

    #atomDictionary[key] = [list(map(int, x.split(' '))) for x in morphologyConfig.find(key).text.split('\n')[1:-1]]




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
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--temperature',
                        type=float,
                        default=290,
                        required=False,
                        help='''The desired temperature of the simulation in kelvin.\n
                        If unspecified, the temperature is set to the default value of
                        290 K.
                        ''')
    parser.add_argument('-r', '--run_time',
                        type=float,
                        default=1e7,
                        required=False,
                        help='''The number of timesteps to run the MD simulation for.\n
                        If unspecified, the default runtime of 1e7 will be used.
                        ''')
    args, file_list = parser.parse_known_args()

    for file_name in file_list:
        hoomd.context.initialize("")
        system = hoomd.deprecated.init.read_xml(filename=file_name)
        # Set the correct atom Masses
        for atom in system.particles:
            atom.mass = masses[atom.type]

        snapshot = system.take_snapshot()
        # Assign the required velocities based on the requested temperature
        initialized_snapshot = initialize_velocities(snapshot, args.temperature)
        # Finally, restore the snapshot
        system.restore_snapshot(initialized_snapshot)

        set_coeffs(file_name)

        hoomd.md.integrate.mode_standard(dt=0.001);
        carbons = hoomd.group.type('C')
        hydrogens = hoomd.group.type('H')
        hydrocarbons = hoomd.group.union(carbons, hydrogens)
        integrator = hoomd.md.integrate.nvt(group=hydrocarbons, tau=0.1, kT=args.temperature)
        #hoomd.md.nlist.reset_exclusions(exclusions = ['bond', 'angle', 'dihedral', 'body'])

        hoomd.dump.gsd(filename=".".join(file_name.split(".")[:-1]) + "_traj.gsd", period=int(args.run_time/500),
                       overwrite=True)
        logQuantities = ['temperature', 'pressure', 'volume', 'potential_energy', 'kinetic_energy',
                         'pair_lj_energy', 'pair_ewald_energy', 'pppm_energy', 'bond_harmonic_energy',
                         'angle_harmonic_energy', 'dihedral_opls_energy']
        hoomd.analyze.log(filename=fileName.split('.')[0] + ".log", quantities=logQuantities,
                            period=int(args.run_time/10000), header_prefix='#', overwrite=True)
        ## Now incrementally ramp the charges
        #for chargePhase in range(chargeIncrements + 1):
        #    print("Incrementing charge phase", chargePhase, "of", chargeIncrements + 1)
        #    for atom in system.particles:
        #        oldCharge = copy.deepcopy(atom.charge)
        #        atom.charge = charges[atom.type] * (chargePhase / float(chargeIncrements))
        #    hoomd.run(chargeTimesteps)

        # Get the initial box size dynamically
        hoomd.run_upto(args.run_time)
