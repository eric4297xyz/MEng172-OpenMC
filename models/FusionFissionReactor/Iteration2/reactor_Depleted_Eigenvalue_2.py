# FILE: reactor_U_Eigenvalue.py
# Eigenvalue depletion simulation (explicit k-eff tracking)

import openmc
import openmc.deplete
import sys
import os
import time

#from blanket import build_u238_sphere
from fuel_blanket_2 import build_spentfuelsphere_albox
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
E_min, E_max = 14e6, 14.5e6
u238_sphere_radius = 25.5
particles_per_batch = 1000
num_batches = 20
num_inactive = 5
CHAIN_FILE = "models/FusionFissionReactor/Iteration2/chain_endfb80_pwr.xml"
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

print("Building model...")
my_source = create_cylindrical_source(
    height=cyl_H,
    radius=cyl_R,
    e_min=E_min,
    e_max=E_max
)
(my_geometry, my_materials) = build_spentfuelsphere_albox(
    box_side=20.0,
    box_width=20.0,
    sphere_inner_radius=50.0,
    sphere_outer_radius=100.0
)

settings = openmc.Settings()
settings.run_mode  = 'eigenvalue'  # <-- Explicit eigenvalue mode
settings.particles = particles_per_batch
settings.batches   = num_batches
settings.inactive  = num_inactive
settings.seed      = 42
settings.source    = my_source
settings.max_lost_particles = 1000000

print("Creating 3D heating tally...")
mesh_bounds = [-50.0, -50.0, -50.0, 50.0, 50.0, 50.0]
mesh_dim = (20, 20, 20)
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
    reduce_chain=True,
    reduce_chain_level=5,
    normalization_mode="energy-deposition"  # Use energy-deposition for eigenvalue mode
)

burnable_mats = []
for m in my_materials:
    if m.name and 'UO2' in m.name:
        burnable_mats.append(m)
if burnable_mats:
    print(f"Using burnable material IDs: {[m.id for m in burnable_mats]}")
else:
    print("Warning: no burnable materials detected; operator.burnable_mats left empty")

# Set output directory
output_dir = "results/iteration2_results_depleted_fuel"
os.makedirs(output_dir, exist_ok=True)

# Save current directory and change to output directory for depletion run
original_dir = os.getcwd()
os.chdir(output_dir)

power_list = [1e6] * len(timesteps_in_seconds)  # 1 MW per step
integrator = openmc.deplete.PredictorIntegrator(
    operator=operator,
    timesteps=timesteps_in_seconds,
    power=power_list,
    timestep_units='s'
)

print(f"Running depletion for {time_seconds} seconds in eigenvalue mode...")
print(f"Output will be saved to: {os.path.join(original_dir, output_dir)}/")
integrator.integrate()

# Return to original directory
os.chdir(original_dir)

# Calculate elapsed time
end_time = time.time()
elapsed_time = end_time - start_time
elapsed_minutes = elapsed_time / 60.0
elapsed_hours = elapsed_time / 3600.0

print(f"\nDepletion simulation complete. Results are in '{output_dir}/depletion_results.h5'")
print(f"\n{'='*70}")
print(f"Total simulation time: {elapsed_time:.2f} seconds")
print(f"                       {elapsed_minutes:.2f} minutes")
print(f"                       {elapsed_hours:.2f} hours")
print(f"{'='*70}")
