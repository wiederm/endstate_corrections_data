# ------------------------------------------------------------ #
# This script performs endstate corrections on the HiPen dataset
# ------------------------------------------------------------ #

# general imports
from openmm import LangevinIntegrator
from openmm.app import (
    CharmmParameterSet,
    CharmmPsfFile,
    CharmmCrdFile,
    PDBFile,
    Simulation,
)
from openmmml import MLPotential
import endstate_correction
from endstate_correction.constant import zinc_systems, blacklist
from endstate_correction.analysis import plot_endstate_correction_results
from endstate_correction.protocol import perform_endstate_correction, Protocol
import mdtraj
import pickle, sys, os

########################################################
########################################################
# ------------------- set up system -------------------#
system_idx = int(sys.argv[1])
system_name = zinc_systems[system_idx][0]
env = "vacuum"

if system_name in blacklist:
    print("System in blacklist. Aborting.")
    sys.exit()

print(f"Perfom endstate correction for {system_name} in {env}")

# define directory containing parameters
parameter_base = f"../hipen_systems/"
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
# write pdb file
with open("temp.pdb", "w") as outfile:
    PDBFile.writeFile(psf.topology, crd.positions, outfile)

# define region that should be treated with the qml
chains = list(psf.topology.chains())
ml_atoms = [atom.index for atom in chains[0].atoms()]
print(f"ML atoms: {ml_atoms=}")
mm_system = psf.createSystem(params=params)
# define system
potential = MLPotential("ani2x")
ml_system = potential.createMixedSystem(
    psf.topology, mm_system, ml_atoms, interpolate=True
)
sim = Simulation(psf.topology, ml_system, LangevinIntegrator(300, 1, 0.001))
########################################################
########################################################
# ------------------- load samples --------------------#
n_samples = 5_000
n_steps_per_sample = 1_000

# define directory containing MM and QML sampling data
traj_base = f"../data/equilibrium_samples/{system_name}"

# load MM samples
mm_samples = []
base = f"{traj_base}/{system_name}_samples_{n_samples}_steps_{n_steps_per_sample}_lamb_0.0000_{env}"
mm_samples = mdtraj.load_dcd(
    f"{base}.dcd",
    top=psf_file,
)[
    int((1_000 / 100) * 20) :
    ]
print(f"Initializing switch from {len(mm_samples)} MM samples")

# load QML samples
qml_samples = []
base = f"{traj_base}/{system_name}_samples_{n_samples}_steps_{n_steps_per_sample}_lamb_1.0000_{env}"
qml_samples = mdtraj.load_dcd(
    f"{base}.dcd",
    top=psf_file,
)[
    int((1_000 / 100) * 20) :
    ]
print(f"Initializing switch from {len(qml_samples)} QML samples")

########################################################
########################################################
# ----------------- perform correction ----------------#

# define the output directory
output_base = f"../data/switching_results/{system_name}/"
os.makedirs(output_base, exist_ok=True)

####################################################
# ---------------- FEP protocol --------------------
####################################################
fep_protocol = Protocol(
    method="FEP",
    sim=sim,
    reference_samples=mm_samples,
    target_samples=qml_samples,
    nr_of_switches=1_000
)

####################################################
# ----------------- NEQ protocol -------------------
####################################################
neq_protocol = Protocol(
    method="NEQ",
    sim=sim,
    reference_samples=mm_samples,
    target_samples=qml_samples,
    nr_of_switches=100,
    neq_switching_length=1_000
    save_endstates=True,
    save_trajs=True,
)

# perform correction
r_fep = perform_endstate_correction(fep_protocol)
r_neq = perform_endstate_correction(neq_protocol)

# save fep and neq results in a pickle file
pickle.dump((r_fep, r_neq), open(f"{output_base}/results.pickle", "wb"))

# plot results
plot_endstate_correction_results(system_name, r_fep, f"{output_base}/results_fep.png")
plot_endstate_correction_results(system_name, r_neq, f"{output_base}/results_neq.png")