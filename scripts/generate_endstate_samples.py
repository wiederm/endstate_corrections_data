# general imports
from openmm import LangevinIntegrator
from endstate_correction.constant import (
    temperature,
)
import numpy as np
from openmm.app import (
    CharmmParameterSet,
    CharmmPsfFile,
    CharmmCrdFile,
    Simulation,
    DCDReporter
)
import endstate_correction
import os
import sys, os
from endstate_correction.constant import zinc_systems, blacklist
from openmmml import MLPotential

########################################################
########################################################
# ------------------- set up system -------------------#
package_path = endstate_correction.__path__[0]
system_idx = int(sys.argv[1])
system_name = zinc_systems[system_idx][0]

if system_name in blacklist:
    print("System in blacklist. Aborting.")
    sys.exit()

env = "vacuum"

print(system_name)
print(env)

# define directory containing parameters
parameter_base = f"{package_path}/data/hipen_data"
# load the charmm specific files (psf, rtf, prm and str files)
psf_file = f"{parameter_base}/{system_name}/{system_name}.psf"
crd_file = f"{parameter_base}/{system_name}/{system_name}.crd"
psf = CharmmPsfFile(psf_file)
crd = CharmmCrdFile(crd_file)
params = CharmmParameterSet(
    f"{parameter_base}/top_all36_cgenff.rtf",
    f"{parameter_base}/par_all36_cgenff.prm",
    f"{parameter_base}/{system_name}/{system_name}.str",
)

# define region that should be treated with the qml
chains = list(psf.topology.chains())
ml_atoms = [atom.index for atom in chains[0].atoms()]
print(f"{ml_atoms=}")
mm_system = psf.createSystem(params=params)
# define system
potential = MLPotential("ani2x")
ml_system = potential.createMixedSystem(
    psf.topology, mm_system, ml_atoms, interpolate=True
)
sim = Simulation(psf.topology, ml_system, LangevinIntegrator(300, 1, 0.001))
##############################################################
# ------------------ Start equilibrium sampling ---------------
# define equilibirum sampling control parameters
run_id = 1
n_samples = 5_000
n_steps_per_sample = 1_000
# path where samples should be stored (will be created if it doesn't exist)
base = f"../data/equilibrium_samples/{system_name}"
os.makedirs(base, exist_ok=True)

nr_lambda_states = 2  # samples equilibrium distribution at both endstates
lambs = np.linspace(0, 1, nr_lambda_states)
# perform sampling for each lambda state

for lamb in lambs:
    print(f"{lamb=}")
    # define where to store samples
    trajectory_file = f"{base}/{system_name}_samples_{n_samples}_steps_{n_steps_per_sample}_lamb_{lamb:.4f}_{env}.dcd"
    print(f"Trajectory saved to: {trajectory_file}")
    # set lambda
    sim.context.setParameter("lambda_interpolate", lamb)
    # set coordinates
    sim.context.setPositions(crd.positions)
    # try to set velocities using openMM, fall back to manual velocity seeding if it fails
    sim.context.setVelocitiesToTemperature(temperature)
    # define DCDReporter
    sim.reporters.append(
        DCDReporter(
            trajectory_file,
            n_steps_per_sample,
        )
    )
    # perform sampling
    sim.step(n_samples * n_steps_per_sample)
    sim.reporters.clear()