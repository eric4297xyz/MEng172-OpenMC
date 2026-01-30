# FILE: reactor_Depleted_FixedSource_3.py
# Iteration3: Fixed source depletion simulation with TRU-only fuel (no uranium)

import openmc
import openmc.deplete
import sys
import os
import time

from fuel_blanket_3 import build_trufuelsphere_albox

# --- This block "makes it public" ---
current_file_dir = os.path.dirname(os.path.abspath(__file__))
models_dir = os.path.dirname(current_file_dir)
source_dir = os.path.join(models_dir, 'NeutronSource')
if source_dir not in sys.path:
    sys.path.append(source_dir)
# --- End of block ---
from neutronsource import create_cylindrical_source

# -------------------- Simulation Parameters --------------------
cyl_H = 5.0
cyl_R = 1.0
E_min, E_max = 2.3e6, 2.5e6
source_rate = 1e14  # neutrons per second
particles_per_batch = 1000
num_batches = 20
CHAIN_FILE = "chain_endfb80_pwr.xml"

# Depletion time: 182.5 days (same as Iteration2)
days_per_year = 365.0
hours_per_day = 24.0
seconds_per_hour = 3600.0
time_years = 0.5
time_seconds = time_years * days_per_year * hours_per_day * seconds_per_hour
num_steps = 3
total_time_seconds = time_seconds
step_size = total_time_seconds / num_steps
timesteps_in_seconds = [step_size] * num_steps

# Start timer
start_time = time.time()

print("="*70)
print("ITERATION 3: TRU FUEL (NO URANIUM)")
print("="*70)
print("\nThis simulation tests purified transuranic fuel:")
print("  - All uranium isotopes removed (U233, U234, U235, U236, U237, U238)")
print("  - Contains only transuranics: Pu, Np, Am, Cm")
print("  - Plus fission products and oxygen")
print("  - Tests if TRU can be transmuted without breeding more from U238")
print("")

print("Building model...")
my_source = create_cylindrical_source(
    height=cyl_H,
    radius=cyl_R,
    e_min=E_min,
    e_max=E_max
)
# Set source strength to 1e14 neutrons/second
my_source.strength = source_rate
print(f"Source rate set to {source_rate:.2e} neutrons/second")

(my_geometry, my_materials) = build_trufuelsphere_albox(
    box_side=20.0,
    box_width=20.0,
    sphere_inner_radius=50.0,
    sphere_outer_radius=100.0
)

settings = openmc.Settings()
settings.run_mode  = 'fixed source'  # <-- Fixed source mode
settings.particles = particles_per_batch
settings.batches   = num_batches
settings.seed      = 42
settings.source    = my_source
settings.max_lost_particles = 1000000

print("Creating 3D heating tally...")
mesh_bounds = [-110.0, -110.0, -110.0, 110.0, 110.0, 110.0]  # Extended to fully cover fuel sphere (r=100)
mesh_dim = (22, 22, 22)  # Increased resolution to maintain ~10cm bins
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
operator = openmc.deplete.CoupledOperator(
    model=model,
    chain_file=CHAIN_FILE,
    reduce_chain=False,
    normalization_mode="source-rate"  # Use source-rate for fixed source mode
)

burnable_mats = []
for m in my_materials:
    if m.name and 'TRU' in m.name:
        burnable_mats.append(m)
if burnable_mats:
    print(f"Using burnable material IDs: {[m.id for m in burnable_mats]}")
else:
    print("Warning: no burnable materials detected; operator.burnable_mats left empty")

# Set output directory
output_dir = "results/iteration3_results_tru_fuel"
os.makedirs(output_dir, exist_ok=True)

# Save current directory and change to output directory for depletion run
original_dir = os.getcwd()
os.chdir(output_dir)

# For fixed source mode with source-rate normalization, use source_rate parameter
source_rate_list = [source_rate] * len(timesteps_in_seconds)  # neutrons/s per step
integrator = openmc.deplete.PredictorIntegrator(
    operator=operator,
    timesteps=timesteps_in_seconds,
    source_rates=source_rate_list,  # Use source_rates instead of power
    timestep_units='s'
)

print(f"\nRunning TRU fuel depletion for {time_seconds} seconds in fixed source mode...")
print(f"Duration: {time_seconds / (24*3600):.2f} days ({time_years:.2f} years)")
print(f"Source rate: {source_rate:.2e} neutrons/second")
print(f"Number of steps: {num_steps}")
print(f"Output will be saved to: {os.path.join(original_dir, output_dir)}/")
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
print("ITERATION 3 COMPLETE")
print(f"{'='*70}")
print(f"Depletion simulation complete. Results are in:")
print(f"  {output_dir}/depletion_results.h5")
print(f"\nTotal simulation time: {elapsed_time:.2f} seconds")
print(f"                       {elapsed_minutes:.2f} minutes")
print(f"                       {elapsed_hours:.2f} hours")
print(f"{'='*70}")
print("\nExpected outcome:")
print("  - NO uranium breeding (no U238 present)")
print("  - TRU should decrease via fission and transmutation")
print("  - Fission products should increase")
print("  - Use results_depleted_fixedsource_3.py to analyze")
print(f"{'='*70}")
