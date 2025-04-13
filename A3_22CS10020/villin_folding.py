# File: villin_folding.py
import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.fragment import ConstantLengthFragSet
from pyrosetta.rosetta.protocols.moves import *

# Try importing ClassicFragmentMover
try:
    from pyrosetta.rosetta.protocols.simple_moves import ClassicFragmentMover
    print("Imported ClassicFragmentMover from simple_moves")
except ImportError:
    try:
        from pyrosetta.rosetta.protocols.moves import ClassicFragmentMover
        print("Imported ClassicFragmentMover from moves")
    except ImportError:
        try:
            from pyrosetta.rosetta.protocols.fragment import FragmentMover as ClassicFragmentMover
            print("Using FragmentMover as ClassicFragmentMover fallback")
        except ImportError as e:
            raise ImportError("Failed to import ClassicFragmentMover: " + str(e))

# Import score function
try:
    from pyrosetta import create_score_function
    print("Imported create_score_function from pyrosetta")
except ImportError:
    try:
        from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
        def create_score_function(name):
            return ScoreFunctionFactory.create_score_function(name)
        print("Using ScoreFunctionFactory fallback")
    except ImportError as e:
        raise ImportError("Failed to import score function utilities: " + str(e))

# Initialize PyRosetta
pyrosetta.init(extra_options="-in::file::fullatom -mute all")
print("PyRosetta initialized.")

# Villin headpiece sequence
sequence = "MLSDEDFKAFGMTRSAFANLPLWKQQNLKKEKLLF"

# Step 1: Generate Starting Pose
pose = pose_from_sequence(sequence, "fa_standard")
print("Created fullatom pose from sequence.")

# Step 2: Linearize the Pose
for i in range(1, pose.total_residue() + 1):
    pose.set_phi(i, -150.0)
    pose.set_psi(i, 150.0)
    pose.set_omega(i, 180.0)
print("Linearized pose with extended conformation.")

# Step 3: Convert to Centroid Mode
to_centroid = SwitchResidueTypeSetMover("centroid")
to_centroid.apply(pose)
print("Converted pose to centroid mode.")

# Step 4: Setup MoveMap
movemap = MoveMap()
movemap.set_bb(True)  # Allow all backbone movements
print("Set up MoveMap with full backbone flexibility.")

# Step 5: Load Fragment Libraries and Create Movers
try:
    fragset9 = ConstantLengthFragSet(9, "aat000_09.frag")
    fragset3 = ConstantLengthFragSet(3, "aat000_03.frag")
    print("Loaded 9-mer and 3-mer fragment libraries.")
except FileNotFoundError as e:
    print(f"Error: {e}. Ensure aat000_09.frag and aat000_03.frag are in the directory.")
    raise

mover_9mer = ClassicFragmentMover(fragset9, movemap)
mover_3mer = ClassicFragmentMover(fragset3, movemap)
repeat_9mer = RepeatMover(mover_9mer, 3)
repeat_3mer = RepeatMover(mover_3mer, 3)

seq_mover = SequenceMover()
seq_mover.add_mover(repeat_9mer)
seq_mover.add_mover(repeat_3mer)
print("Set up fragment insertion movers.")

# Step 6: Monte Carlo Sampling
scorefxn = create_score_function("score3")
mc = MonteCarlo(pose, scorefxn, 3.0)
print("Initialized Monte Carlo with score3 and kT=3.0.")

trial_mover = TrialMover(seq_mover, mc)
cycle_mover = RepeatMover(trial_mover, 300)
cycle_mover.apply(pose)
print("Completed 300 cycles of Monte Carlo sampling.")

# Step 7: Recover and Convert Best Decoy
mc.recover_low(pose)
to_fullatom = SwitchResidueTypeSetMover("fa_standard")
to_fullatom.apply(pose)
print("Recovered lowest-energy decoy and converted to fullatom mode.")

# Save output
pose.dump_pdb("villin_predicted.pdb")
print("Saved predicted structure to villin_predicted.pdb.")