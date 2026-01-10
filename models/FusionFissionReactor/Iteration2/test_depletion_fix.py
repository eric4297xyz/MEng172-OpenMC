"""
Quick test to verify reduce_chain=False fixes the uranium initialization issue.
Runs 1 short depletion step to check if uranium is properly tracked.
"""
import openmc
import openmc.deplete
import sys
import os

sys.path.append('../NeutronSource')
from neutronsource import create_cylindrical_source
from fuel_blanket_2 import build_spentfuelsphere_albox

print("="*70)
print("TESTING DEPLETION FIX: reduce_chain=False")
print("="*70)

# Build model (same as main simulation)
print("\nBuilding model...")
my_source = create_cylindrical_source(height=5.0, radius=1.0, e_min=2.3e6, e_max=2.5e6)
my_source.strength = 1e14

(my_geometry, my_materials) = build_spentfuelsphere_albox(
    box_side=20.0, box_width=20.0, sphere_inner_radius=50.0, sphere_outer_radius=100.0
)

settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.particles = 500  # Reduced for speed
settings.batches = 10     # Reduced for speed
settings.seed = 42
settings.source = my_source
settings.max_lost_particles = 1000000

model = openmc.Model(geometry=my_geometry, materials=my_materials, settings=settings)

# Create operator with reduce_chain=False (the fix)
print("Creating depletion operator with reduce_chain=False...")
operator = openmc.deplete.CoupledOperator(
    model=model,
    chain_file="chain_endfb80_pwr.xml",
    reduce_chain=False,  # THE FIX
    normalization_mode="source-rate"
)

print(f"Operator created successfully")
print(f"Number of burnable materials: {len(operator.burnable_mats)}")
print(f"Number of nuclides in chain: {len(operator.chain.nuclides)}")

# Check initial uranium composition
print("\n" + "="*70)
print("CHECKING INITIAL COMPOSITION")
print("="*70)

spent_fuel = [m for m in my_materials if m.name and 'UO2' in m.name][0]
material_densities = spent_fuel.get_nuclide_atom_densities()

print("\nMaterial definition (should be correct):")
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
        # operator.number stores total atoms, not density
        value = number_data[0, idx]
        print(f"  {nuc}: {value:.6e} total atoms")
        
        if nuc.startswith('U') and value < 1e20:
            print(f"    ❌ ERROR: {nuc} is essentially zero!")
        else:
            print(f"    ✓ OK")

# Run ONE very short depletion step to test
print("\n" + "="*70)
print("RUNNING SHORT TEST DEPLETION")
print("="*70)

output_dir = "results/test_depletion_fix"
os.makedirs(output_dir, exist_ok=True)
os.chdir(output_dir)

# 1 hour depletion for quick test
timestep = 3600  # 1 hour in seconds
integrator = openmc.deplete.PredictorIntegrator(
    operator=operator,
    timesteps=[timestep],
    source_rates=[1e14],
    timestep_units='s'
)

print(f"Running 1 hour depletion test...")
print(f"(This will take a few minutes)")
try:
    integrator.integrate()
    print("\n✓ Depletion completed successfully!")
    
    # Check results
    print("\n" + "="*70)
    print("CHECKING RESULTS")
    print("="*70)
    
    import h5py
    with h5py.File('depletion_results.h5', 'r') as f:
        number = f['number'][:]  # (n_steps, 1, 1, n_nuclides)
        nuclides = list(f['nuclides'].keys())
        
        print(f"\nResults file has {len(nuclides)} nuclides")
        print(f"Number of time steps: {number.shape[0]}")
        
        for nuc_name in ['U235', 'U238', 'Pu239']:
            if nuc_name in nuclides:
                idx = nuclides.index(nuc_name)
                initial = number[0, 0, 0, idx]
                final = number[-1, 0, 0, idx]
                print(f"\n{nuc_name}:")
                print(f"  Initial: {initial:.6e} atoms/(b-cm)")
                print(f"  Final:   {final:.6e} atoms/(b-cm)")
                
                if nuc_name.startswith('U'):
                    if initial < 1e-10:
                        print(f"  ❌ FAILED: Initial is still zero!")
                    else:
                        print(f"  ✓ SUCCESS: Initial composition preserved!")
                        change_pct = (final - initial) / initial * 100
                        print(f"  Change: {change_pct:+.4f}%")
        
        print("\n" + "="*70)
        if number[0, 0, 0, nuclides.index('U238')] > 1e-10:
            print("✓✓✓ FIX SUCCESSFUL! Uranium is being tracked properly!")
            print("="*70)
            print("\nYou can now run the full simulation with confidence.")
        else:
            print("❌❌❌ FIX FAILED! Uranium still zero.")
            print("="*70)
            print("\nThis may be a deeper OpenMC issue.")
            
except Exception as e:
    print(f"\n❌ Error during depletion: {e}")
    import traceback
    traceback.print_exc()

os.chdir('../..')
print(f"\nTest results saved in: {output_dir}/")
