# FILE: main_simulation.py

import openmc
import openmc.deplete
import sys
import os

from blanket import build_u238_sphere
from U_fuel_blanket import build_fuelsphere_albox
# --- This block "makes it public" ---
# Get the path to the current file's folder (e.g., .../FusionFissionReactor)
current_file_dir = os.path.dirname(os.path.abspath(__file__))
# Get the path to the parent folder (e.g., .../models)
models_dir = os.path.dirname(current_file_dir)
# Get the path to the folder you want to import from (e.g., .../models/NeutronSource)
source_dir = os.path.join(models_dir, 'NeutronSource')

# Add that folder to Python's "search path"
if source_dir not in sys.path:
    sys.path.append(source_dir)
# --- End of block ---

# Import the functions from your other files
from neutronsource import create_cylindrical_source

# -------------------- Simulation Parameters --------------------
cyl_H = 5.0
cyl_R = 1.0
E_min, E_max = 2.3e6, 2.5e6

# Geometry parameters
u238_sphere_radius = 25.5

# --- Transport Settings ---
particles_per_batch = 10_000
num_batches = 3

# --- Depletion Parameters ---
SOURCE_STRENGTH_PER_SEC = 1.0e15  # neutrons / sec

days_per_year = 365.0
hours_per_day = 24.0
seconds_per_hour = 3600.0

time_years = 5.0
time_seconds = time_years * days_per_year * hours_per_day * seconds_per_hour

timesteps_in_seconds = [time_seconds]# Define how many steps you want (e.g., 12 steps to see monthly progress)
num_steps = 12

total_time_seconds = time_seconds # 1 year
step_size = total_time_seconds / num_steps

timesteps_in_seconds = [step_size] * num_steps

# Path to the chain file you downloaded
CHAIN_FILE = "models/FusionFissionReactor/Iteration1/chain_endfb80_pwr.xml"# ---------------------------------------------------------------

print("Building model...")
# 1. Create the Neutron Source
my_source = create_cylindrical_source(
    height=cyl_H,
    radius=cyl_R,
    e_min=E_min,
    e_max=E_max
)

# 2. Create the Geometry and Materials
# (This will use your modified model_definition.py)
#(my_geometry2, my_materials2) = build_u238_sphere(
#    sphere_radius=u238_sphere_radius,
#    world_padding=50.0
#)

(my_geometry, my_materials) = build_fuelsphere_albox(
    box_side=20.0,
    box_width=20.0,
    sphere_inner_radius=50.0,
    sphere_outer_radius=100.0
)

# 3. Define Settings for the transport "snapshot"
settings = openmc.Settings()
settings.run_mode  = 'fixed source'
settings.particles = particles_per_batch
settings.batches   = num_batches
settings.source    = my_source
settings.max_lost_particles = 1000000

# 4. Define a 3D Mesh Tally for Heating
print("Creating 3D heating tally...")
# Define a mesh that covers your 50cm U-238 sphere
mesh_bounds = [-50.0, -50.0, -50.0, 50.0, 50.0, 50.0]
mesh_dim = (80, 80, 80)  # 80x80x80 resolution

mesh = openmc.RegularMesh()
mesh.dimension = mesh_dim
mesh.lower_left = mesh_bounds[0:3]
mesh.upper_right = mesh_bounds[3:6]

mesh_filter = openmc.MeshFilter(mesh)

heating_tally = openmc.Tally(name="3d_heating_tally")
heating_tally.filters = [mesh_filter]
heating_tally.scores = ["heating"]

# B. Flux Tally (FIXED)
flux_tally = openmc.Tally(name='flux_tally')

# FIX: Use MaterialFilter instead of CellFilter.
# This will give you the average flux inside your U238 Materials.
flux_tally.filters = [openmc.MaterialFilter(my_materials)]
flux_tally.scores = ['flux']

tallies = openmc.Tallies([heating_tally, flux_tally])

# 5. Create the main OpenMC model
print("Bundling model...")
model = openmc.Model(
    geometry=my_geometry, 
    materials=my_materials, 
    settings=settings,
    tallies=tallies 
)

# Export the bundled model once so material IDs are consistent between
# the model and the materials file (avoid duplicate registrations).
model.export_to_xml()

print("Setting up depletion operator...")
operator = openmc.deplete.CoupledOperator(
    model=model,
    chain_file=CHAIN_FILE,
    reduce_chain=True,       # Simplifies the chain to only what's needed
    reduce_chain_level=5,
    normalization_mode="source-rate"
)

# Explicitly identify burnable materials (avoid relying on automatic detection
# which can fail when materials were created/exported multiple times).
burnable_mats = []
for m in my_materials:
    if m.name and 'UO2' in m.name:
        burnable_mats.append(int(m.id))

if burnable_mats:
    operator.burnable_mats = burnable_mats
    print(f"Using burnable material IDs: {burnable_mats}")
else:
    print("Warning: no burnable materials detected; operator.burnable_mats left empty")

source_rates_list = [SOURCE_STRENGTH_PER_SEC] * len(timesteps_in_seconds)

integrator = openmc.deplete.PredictorIntegrator(
    operator=operator,
    timesteps=timesteps_in_seconds,
    source_rates=source_rates_list,  # Now both are length 12
    timestep_units='s'
)

print(f"Running depletion for {time_seconds} seconds...")
integrator.integrate()

print("Depletion simulation complete. Results are in 'depletion_results.h5'")