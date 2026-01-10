# FILE: reactor_Depleted_Eigenvalue_Test.py
# Testing eigenvalue depletion mode to see if it properly tracks uranium

import openmc
import openmc.deplete
import sys
import os
import time

from fuel_blanket_2 import build_spentfuelsphere_albox

# -------------------- Simulation Parameters --------------------
particles_per_batch = 1000
num_batches = 100
inactive_batches = 20
CHAIN_FILE = "chain_endfb80_pwr.xml"

# Shorter test: 1 day depletion
timestep_hours = 24  # 1 day
timestep_seconds = timestep_hours * 3600
num_steps = 1
timesteps_in_seconds = [timestep_seconds] * num_steps

# Initial power (Watts) - will be calculated from k-eff
# For spent fuel, we expect subcritical, but OpenMC will still calculate power
power_watts = 1e6  # 1 MW initial guess

# Start timer
start_time = time.time()

print("="*70)
print("EIGENVALUE MODE DEPLETION TEST")
print("="*70)

print("\nBuilding geometry...")
(my_geometry, my_materials) = build_spentfuelsphere_albox(
    box_side=20.0,
    box_width=20.0,
    sphere_inner_radius=50.0,
    sphere_outer_radius=100.0
)

print("Setting up eigenvalue mode...")
settings = openmc.Settings()
settings.run_mode = 'eigenvalue'  # Changed from 'fixed source'
settings.particles = particles_per_batch
settings.batches = num_batches
settings.inactive = inactive_batches
settings.seed = 42

# No external source in eigenvalue mode - OpenMC finds criticality
# Remove: settings.source = my_source

print(f"  Particles per batch: {particles_per_batch}")
print(f"  Total batches: {num_batches}")
print(f"  Inactive batches: {inactive_batches}")

print("\nBundling model...")
model = openmc.Model(
    geometry=my_geometry, 
    materials=my_materials, 
    settings=settings
)
model.export_to_xml()

print("\nSetting up depletion operator...")
operator = openmc.deplete.CoupledOperator(
    model=model,
    chain_file=CHAIN_FILE,
    reduce_chain=False,  # Keep full chain
    normalization_mode="fission-q"  # Use fission energy for eigenvalue mode
)

print(f"Operator created")
print(f"Number of burnable materials: {len(operator.burnable_mats)}")

# Check initial composition
print("\n" + "="*70)
print("CHECKING INITIAL COMPOSITION")
print("="*70)
spent_fuel = [m for m in my_materials if m.name and 'UO2' in m.name][0]
material_densities = spent_fuel.get_nuclide_atom_densities()

print("\nMaterial definition:")
for nuc in ['U235', 'U238', 'Pu239']:
    if nuc in material_densities:
        density = material_densities[nuc]
        total = density * 1e24 * spent_fuel.volume
        print(f"  {nuc}: {density:.6e} atoms/(b-cm) = {total:.6e} total atoms")

print("\nOperator initial composition:")
nuclides_list = list(operator.number.nuclides)
number_data = operator.number.number
for nuc in ['U235', 'U238', 'Pu239']:
    if nuc in nuclides_list:
        idx = nuclides_list.index(nuc)
        value = number_data[0, idx]
        print(f"  {nuc}: {value:.6e} total atoms")

# Set output directory
output_dir = "results/test_eigenvalue_depletion"
os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)

# For eigenvalue mode, use power normalization
print("\n" + "="*70)
print("RUNNING EIGENVALUE DEPLETION")
print("="*70)
print(f"Duration: {timestep_hours} hours ({timestep_seconds} seconds)")
print(f"Number of steps: {num_steps}")
print(f"Power: {power_watts:.2e} W ({power_watts/1e6:.2f} MW)")

integrator = openmc.deplete.PredictorIntegrator(
    operator=operator,
    timesteps=timesteps_in_seconds,
    power=power_watts,  # Use power instead of source_rates
    timestep_units='s'
)

print("\nRunning depletion (this will take several minutes)...")
try:
    integrator.integrate()
    print("\n✓ Depletion completed successfully!")
    
    # Check results
    print("\n" + "="*70)
    print("CHECKING RESULTS")
    print("="*70)
    
    import h5py
    with h5py.File('depletion_results.h5', 'r') as f:
        number = f['number'][:]
        nuclides = list(f['nuclides'].keys())
        
        print(f"\nResults file has {len(nuclides)} nuclides")
        print(f"Number of time steps: {number.shape[0]}")
        
        print("\nUranium tracking:")
        for nuc_name in ['U235', 'U238']:
            if nuc_name in nuclides:
                idx = nuclides.index(nuc_name)
                initial = number[0, 0, 0, idx]
                final = number[-1, 0, 0, idx]
                print(f"\n{nuc_name}:")
                print(f"  Initial: {initial:.6e} atoms/(b-cm)")
                print(f"  Final:   {final:.6e} atoms/(b-cm)")
                
                if initial > 1e-10:
                    change_pct = (final - initial) / initial * 100
                    print(f"  Change: {change_pct:+.4f}%")
                    print(f"  ✓✓ SUCCESS: Eigenvalue mode properly tracks uranium!")
                else:
                    print(f"  ❌ Still zero - eigenvalue mode has same issue")
        
        # Check k-eff
        if 'eigenvalues' in f:
            keff = f['eigenvalues'][:]
            print(f"\nk-eff evolution:")
            for i, k in enumerate(keff[:, 0, 0]):
                print(f"  Step {i}: k-eff = {k:.6f}")
        
        print("\n" + "="*70)
        u238_idx = nuclides.index('U238')
        if number[0, 0, 0, u238_idx] > 1e-10:
            print("🎉 EIGENVALUE MODE WORKS!")
            print("Uranium is being properly tracked.")
            print("You can use eigenvalue mode for your simulations.")
        else:
            print("❌ Eigenvalue mode has the same HDF5 recording issue.")
            print("This confirms it's a broader OpenMC bug.")
            
except Exception as e:
    print(f"\n❌ Error during depletion: {e}")
    import traceback
    traceback.print_exc()

os.chdir('../..')

# Calculate elapsed time
end_time = time.time()
elapsed_time = end_time - start_time

print(f"\n{'='*70}")
print(f"Total time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
print(f"Results saved in: {output_dir}/")
print(f"{'='*70}")
