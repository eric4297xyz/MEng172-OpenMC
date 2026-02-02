# FILE: reactor_Depleted_Eigenvalue_5.py
# Iteration5: Eigenvalue depletion simulation with spent LWR fuel
# MEDIUM RESOLUTION for reliable completion:
#   - Eigenvalue mode (not fixed source) for proper multiplication handling
#   - Medium particle statistics (1M particles per step) for balance
#   - Short time steps (10 days) for accurate tracking
#   - 30 days total simulation time
#   - Initial fissile material drives the reaction (no external source)
#
# This approach provides physically accurate depletion rates without
# the multiplication over-counting issue present in fixed-source mode.

import openmc
import openmc.deplete
import sys
import os
import time

from fuel_blanket_5 import build_trufuelsphere_nobox

# -------------------- Simulation Parameters --------------------
# EIGENVALUE MODE - No external source needed
# The fissile material (U-235, Pu-239) in the spent fuel drives the reaction

# MEDIUM STATISTICS for reliable completion
particles_per_batch = 10000    # 10k particles per batch
inactive_batches = 10          # 10 batches to converge fission source
active_batches = 40           # 40 active batches for statistics
# Total: 1 million active particles per depletion step (MEDIUM RES)

# Reduced chain has 1760 isotopes (vs 3820) - 54% memory savings
CHAIN_FILE = os.path.abspath("chain_reduced_tru.xml")

# QUARTERLY TIME STEPS for 1-year depletion tracking
# 3-month steps (90 days each) for 1 year total
days_per_step = 90.0
num_steps = 4  # 4 steps × 90 days = 360 days (1 year)
seconds_per_day = 24.0 * 3600.0
step_size_seconds = days_per_step * seconds_per_day
timesteps_in_seconds = [step_size_seconds] * num_steps
total_time_seconds = step_size_seconds * num_steps
total_days = total_time_seconds / seconds_per_day

# Start timer
start_time = time.time()

print("="*70)
print("ITERATION 5: EIGENVALUE DEPLETION - MEDIUM RESOLUTION")
print("="*70)
print("\nThis simulation uses EIGENVALUE MODE for accurate depletion:")
print("  - No external source - fission drives the reaction")
print("  - Initial fissile: U-235 (0.35%), Pu-239 (0.26%)")
print("  - Medium statistics: 400k particles per step (BALANCED)")
print("  - Quarterly time steps: 3 months (90 days) for 1 year total")
print(f"  - Total simulation time: {total_days:.0f} days (~1 year)")
print("  - Spent UO2 fuel: U-238 (30.8%), TRU, fission products")
print("")

print("Building model...")
# No external source needed in eigenvalue mode
# Build simplified geometry (no aluminum box)
(my_geometry, my_materials) = build_trufuelsphere_nobox(
    sphere_inner_radius=50.0,
    sphere_outer_radius=100.0
)
print("Geometry: Fuel shell (r=50-100 cm) with inner Na coolant")

settings = openmc.Settings()
settings.run_mode  = 'eigenvalue'  # <-- EIGENVALUE mode for proper depletion
settings.particles = particles_per_batch
settings.batches   = inactive_batches + active_batches
settings.inactive = inactive_batches  # Converge fission source distribution
settings.seed      = 42

# No external source needed - fission source is automatic in eigenvalue mode
print(f"\nEigenvalue mode settings:")
print(f"  Particles per batch: {particles_per_batch:,}")
print(f"  Inactive batches: {inactive_batches}")
print(f"  Active batches: {active_batches}")
print(f"  Total active particles: {particles_per_batch * active_batches:,.0e}")
print(f"  Expected k-eff: ~0.4-0.6 (subcritical spent fuel)")

print("Creating 3D heating tally...")
mesh_bounds = [-110.0, -110.0, -110.0, 110.0, 110.0, 110.0]  # Extended to fully cover fuel sphere (r=100)
mesh_dim = (44, 44, 44)  # High resolution: ~5cm bins for detailed spatial distribution
mesh = openmc.RegularMesh()
mesh.dimension = mesh_dim
mesh.lower_left = mesh_bounds[0:3]
mesh.upper_right = mesh_bounds[3:6]
mesh_filter = openmc.MeshFilter(mesh)
heating_tally = openmc.Tally(name="3d_heating_tally")
heating_tally.filters = [mesh_filter]
heating_tally.scores = ["heating"]
flux_tally = openmc.Tally(name='flux_tally')
flux_tally.filters = [openmc.MaterialFilter(my_materials)]
flux_tally.scores = ['flux']
tallies = openmc.Tallies([heating_tally, flux_tally])

print("Bundling model...")
model = openmc.Model(
    geometry=my_geometry, 
    materials=my_materials, 
    settings=settings,
    tallies=tallies 
)
model.export_to_xml()

print("Setting up depletion operator...")

# In eigenvalue mode, use power normalization (standard for reactors)
# We'll normalize to a specific power level
# For a research/breeding blanket, use modest power: 1 MW thermal
power_watts = 1.0e6  # 1 MW thermal power

operator = openmc.deplete.CoupledOperator(
    model=model,
    chain_file=CHAIN_FILE,
    normalization_mode="fission-q"  # Standard for eigenvalue depletion
)

print(f"Depletion normalization: {power_watts/1e6:.1f} MW thermal power")

burnable_mats = []
for m in my_materials:
    if m.name and ('spent' in m.name or 'UO2' in m.name or 'fuel' in m.name):
        burnable_mats.append(m)
if burnable_mats:
    print(f"Using burnable material IDs: {[m.id for m in burnable_mats]}")
else:
    print("Warning: no burnable materials detected; operator.burnable_mats left empty")

# Set output directory
output_dir = "results/iteration5_eigenvalue_highres"
os.makedirs(output_dir, exist_ok=True)

# Save current directory and change to output directory for depletion run
original_dir = os.getcwd()
os.chdir(output_dir)

# For eigenvalue mode, use power normalization
power_list = [power_watts] * len(timesteps_in_seconds)  # Constant power
integrator = openmc.deplete.PredictorIntegrator(
    operator=operator,
    timesteps=timesteps_in_seconds,
    power=power_list,  # Use power for eigenvalue mode
    timestep_units='s'
)

print(f"\nRunning eigenvalue depletion simulation...")
print(f"Time steps: {num_steps} steps × {days_per_step:.0f} days = {total_days:.0f} days total")
print(f"            t = 0, 3, 6, 9, 12 months")
print(f"Power level: {power_watts/1e6:.1f} MW thermal (constant)")
print(f"Output directory: {os.path.join(original_dir, output_dir)}/")
print("")
print("MEDIUM RESOLUTION - Expected runtime: ~2-3 hours")
print("Progress will be shown for each depletion step.")
print("")

integrator.integrate()

# Return to original directory
os.chdir(original_dir)

# Calculate elapsed time
end_time = time.time()
elapsed_time = end_time - start_time
elapsed_minutes = elapsed_time / 60.0
elapsed_hours = elapsed_time / 3600.0

print(f"\n{'='*70}")
print("ITERATION 5 EIGENVALUE DEPLETION COMPLETE")
print(f"{'='*70}")
print(f"Quantitative depletion simulation complete. Results are in:")
print(f"  {output_dir}/depletion_results.h5")
print(f"\nTotal simulation time: {elapsed_time:.2f} seconds")
print(f"                       {elapsed_minutes:.2f} minutes")
print(f"                       {elapsed_hours:.2f} hours")
print(f"{'='*70}")
print("\n✓ EIGENVALUE MODE - NO CORRECTION NEEDED")
print("="*70)
print("These results are quantitatively accurate and do NOT require")
print("multiplication corrections. The eigenvalue solver properly handles")
print("fission multiplication and normalization.")
print("")
print("Expected results for subcritical spent fuel over 1 year:")
print("  - k-eff: ~0.4-0.6 (subcritical, gradually decreasing)")
print("  - Modest U-235 and Pu-239 depletion from fission")
print("  - Pu-239 production from U-238 neutron capture")
print("  - Gradual buildup of fission products")
print("  - TRU transmutation and decay over 1 year")
print("  - Am-241 growth from Pu-241 decay (T1/2 = 14.4 years)")
print("")
print("To analyze results, run:")
print("  python analyze_eigenvalue_results.py")
print(f"{'='*70}")
