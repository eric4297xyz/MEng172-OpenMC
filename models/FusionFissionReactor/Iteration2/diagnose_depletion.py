"""
Diagnostic script to understand why uranium isn't being tracked in depletion.
"""
import openmc
import openmc.deplete
import sys
import os

# Add neutronsource to path
sys.path.append('../NeutronSource')
from neutronsource import create_cylindrical_source
from fuel_blanket_2 import build_spentfuelsphere_albox

print("="*70)
print("DEPLETION SETUP DIAGNOSTIC")
print("="*70)

# Build geometry
print("\n1. Building geometry...")
geometry, materials = build_spentfuelsphere_albox(
    box_side=20.0,
    box_width=20.0,
    sphere_inner_radius=50.0,
    sphere_outer_radius=100.0
)

# Find spent fuel material
spent_fuel = None
for mat in materials:
    if mat.name and 'UO2' in mat.name:
        spent_fuel = mat
        break

if spent_fuel is None:
    print("ERROR: No spent fuel material found!")
    sys.exit(1)

print(f"\nSpent fuel material:")
print(f"  ID: {spent_fuel.id}")
print(f"  Name: {spent_fuel.name}")
print(f"  Depletable: {spent_fuel.depletable}")
print(f"  Volume: {spent_fuel.volume} cm³")
print(f"  Density: {spent_fuel.density} atoms/(barn-cm)")
print(f"  Number of nuclides: {len(spent_fuel.nuclides)}")

# Get uranium densities
nuclide_densities = spent_fuel.get_nuclide_atom_densities()
print(f"\nInitial Uranium densities:")
for nuc_name in ['U235', 'U238', 'U236', 'U234']:
    if nuc_name in nuclide_densities:
        density = nuclide_densities[nuc_name]
        total_atoms = density * 1e24 * spent_fuel.volume
        print(f"  {nuc_name}: {density:.6e} atoms/(barn-cm) = {total_atoms:.6e} total atoms")

# Create minimal model
print("\n2. Creating minimal model...")
source = create_cylindrical_source(5.0, 1.0, 2.3e6, 2.5e6)
source.strength = 1e14

settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.particles = 100
settings.batches = 2
settings.source = source
settings.seed = 42

model = openmc.Model(geometry=geometry, materials=materials, settings=settings)

# Test depletion operator WITHOUT reduce_chain
print("\n3. Testing CoupledOperator WITHOUT reduce_chain...")
try:
    operator_no_reduce = openmc.deplete.CoupledOperator(
        model=model,
        chain_file="chain_endfb80_pwr.xml",
        reduce_chain=False,
        normalization_mode="source-rate"
    )
    print(f"  SUCCESS: Operator created")
    print(f"  Number of burnable materials: {len(operator_no_reduce.burnable_mats)}")
    
    # Check what nuclides are being tracked
    if hasattr(operator_no_reduce, 'chain'):
        chain = operator_no_reduce.chain
        print(f"  Chain has {len(chain.nuclides)} nuclides")
        for nuc_name in ['U235', 'U238']:
            if nuc_name in chain.nuclide_dict:
                print(f"    ✓ {nuc_name} in chain")
            else:
                print(f"    ✗ {nuc_name} NOT in chain")
except Exception as e:
    print(f"  ERROR: {e}")

# Test depletion operator WITH reduce_chain
print("\n4. Testing CoupledOperator WITH reduce_chain=True...")
try:
    operator_with_reduce = openmc.deplete.CoupledOperator(
        model=model,
        chain_file="chain_endfb80_pwr.xml",
        reduce_chain=True,
        normalization_mode="source-rate"
    )
    print(f"  SUCCESS: Operator created")
    print(f"  Number of burnable materials: {len(operator_with_reduce.burnable_mats)}")
    
    # Check what nuclides are being tracked
    if hasattr(operator_with_reduce, 'chain'):
        chain = operator_with_reduce.chain
        print(f"  Chain has {len(chain.nuclides)} nuclides")
        for nuc_name in ['U235', 'U238', 'Pu239', 'O16']:
            if nuc_name in chain.nuclide_dict:
                nuc = chain[nuc_name]
                print(f"    ✓ {nuc_name} in chain")
            else:
                print(f"    ✗ {nuc_name} NOT in chain")
    
    # Check initial composition from the operator's number array
    print(f"\n  Checking initial number densities in operator...")
    if hasattr(operator_with_reduce, 'number'):
        initial_number = operator_with_reduce.number
        print(f"  Initial number shape: {initial_number.shape}")
        print(f"  Shape interpretation: (materials, 1, nuclides)")
        
        # Find indices for uranium in the chain
        local_nuclides = list(operator_with_reduce.chain)
        nuclide_names = [str(n) for n in local_nuclides]
        
        print(f"\n  Checking uranium isotopes:")
        for nuc_name in ['U235', 'U238', 'Pu239']:
            if nuc_name in nuclide_names:
                idx = nuclide_names.index(nuc_name)
                value = initial_number[0, 0, idx]
                print(f"    {nuc_name} at index {idx}: {value:.6e} atoms/(barn-cm)")
                if value == 0.0 and nuc_name.startswith('U'):
                    print(f"      ⚠️  WARNING: Uranium has ZERO initial density!")
            else:
                print(f"    {nuc_name} not in operator's nuclide list")
                
except Exception as e:
    print(f"  ERROR: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*70)
print("DIAGNOSIS COMPLETE")
print("="*70)
