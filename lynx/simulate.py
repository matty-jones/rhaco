import argparse
import numpy as np
import hoomd
import hoomd.md
import hoomd.deprecated


def set_coeffs():
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
    parser.add_argument("-t", "--temperature", type=float, required=True, help="The temperature of the system")
    args, file_list = parser.parse_known_args()

    #charges = {'CP': 3.880, 'CN': -5.820}
    masses = {'Mo': 1.000, 'Te': 1.000, 'Nb': 1.000, 'V': 1.000, 'C': 1.000, 'H': 1.000}
    run_time = 1e7

    for file_name in file_list:
        hoomd.context.initialize("")
        system = hoomd.deprecated.init.read_xml(filename=file_name)
        # Set the correct atom Masses
        for atom in system.particles:
            atom.mass = masses[atom.type]

        snapshot = system.take_snapshot()
        # Assign the required velocities based on the requested temperature
        initialised_snapshot = initializeVelocities(snapshot, args.temperature)
        # Finally, restore the snapshot
        system.restore_snapshot(initialised_snapshot)

        set_coeffs()

        hoomd.md.integrate.mode_standard(dt=0.001);
        carbons = hoomd.group.type('C')
        hydrogens = hoomd.group.type('H')
        hydrocarbons = hoomd.group.union(carbons, hydrogens)
        integrator = hoomd.md.integrate.nvt(group=hydrocarbons, tau=0.1, kT=args.temperature)
        #hoomd.md.nlist.reset_exclusions(exclusions = ['bond', 'angle', 'dihedral', 'body'])

        hoomd.dump.gsd(filename=".".join(file_name.split(".")[:-1]) + "_traj.gsd", period=int(run_time/500), overwrite=True)
        logQuantities = ['temperature', 'pressure', 'volume', 'potential_energy', 'kinetic_energy',
                         'pair_lj_energy', 'pair_ewald_energy', 'pppm_energy', 'bond_harmonic_energy',
                         'angle_harmonic_energy', 'dihedral_opls_energy']
        hoomd.analyze.log(filename=fileName.split('.')[0] + ".log", quantities=logQuantities,
                            period=int(run_time/10000), header_prefix='#', overwrite=True)
        ## Now incrementally ramp the charges
        #for chargePhase in range(chargeIncrements + 1):
        #    print("Incrementing charge phase", chargePhase, "of", chargeIncrements + 1)
        #    for atom in system.particles:
        #        oldCharge = copy.deepcopy(atom.charge)
        #        atom.charge = charges[atom.type] * (chargePhase / float(chargeIncrements))
        #    hoomd.run(chargeTimesteps)

        # Get the initial box size dynamically
        hoomd.run_upto(run_time)
