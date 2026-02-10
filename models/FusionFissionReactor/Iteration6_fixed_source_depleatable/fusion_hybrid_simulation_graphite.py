#!/usr/bin/env python3
"""
Fusion-Fission Hybrid Reactor Simulation with Fixed-Source Depletion - GRAPHITE MODERATOR STUDY

This script implements fixed-source depletion for a D-D fusion-driven
subcritical fission system using OpenMC's proper tally normalization.

MODIFICATION: Uses graphite as moderator instead of D2O (heavy water) or CO2

System Description:
- External D-D fusion neutron source (2.45 MeV, 10^14 n/s)
- Subcritical spent fuel blanket (k_source ~ 1.1)
- Spherical geometry with graphite moderator (nuclear-grade, 1.7 g/cm³)

Physics Approach:
1. Run fixed-source Monte Carlo transport simulation
2. Extract reaction rates from tallies (normalized by source.strength)
3. Apply nuclide-specific depletion using reaction rates
4. Update material compositions for next time step
5. Repeat for multiple depletion steps

Key Features:
- Correct tally normalization for OpenMC 0.15+ (divide by source.strength)
- Nuclide-specific reaction rates for accurate depletion
- Subcritical neutron multiplication (k_source tracking)
- 3D mesh tallies for spatial visualization
- Graphite moderator with thermal scattering

Author: Graphite Moderator Study
Date: February 2026
"""

import openmc
import openmc.deplete
import numpy as np
import sys
import os
import time
import h5py

from fuel_blanket_5_graphite import build_trufuelsphere_nobox

# ==================== Simulation Parameters ====================

# D-D Fusion Source
source_rate = 1.0e14  # neutrons per second (10^14 n/s typical D-D fusion rate)
source_energy = 2.45e6  # eV (2.45 MeV from D-D fusion)

# Particle statistics (REDUCED FOR QUICK TEST)
particles_per_batch = 5000
num_batches = 30
num_inactive = 10

# Depletion time steps (REDUCED FOR QUICK TEST)
days_per_step = 90.0  # 3 months per step
num_steps = 4  # 4 steps = 1 year to see evolution
seconds_per_day = 24.0 * 3600.0
step_duration_seconds = days_per_step * seconds_per_day

# Chain file
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CHAIN_FILE = os.path.join(SCRIPT_DIR, "chain_reduced_tru.xml")

# Output directory (NEW: Graphite study subdirectory)
output_dir = os.path.join(SCRIPT_DIR, "results_manual_depletion", "graphite_moderator_study")
os.makedirs(output_dir, exist_ok=True)

print("="*70)
print("D-D FUSION-FISSION HYBRID REACTOR SIMULATION - GRAPHITE MODERATOR STUDY")
print("="*70)
print("\nFixed-source depletion with correct OpenMC normalization:")
print(f"  - External D-D fusion source: 2.45 MeV @ {source_rate:.2e} n/s")
print(f"  - Subcritical fission blanket (spent UO₂)")
print(f"  - Graphite moderator (nuclear-grade, 1.7 g/cm³)")
print(f"  - Nuclide-specific reaction tracking")
print(f"\nDepletion: {num_steps} steps × {days_per_step:.0f} days = {num_steps*days_per_step:.0f} days")
print("="*70)

# ==================== Build Geometry ====================

print("\nBuilding geometry with graphite moderator...")
(geometry, materials) = build_trufuelsphere_nobox(
    sphere_inner_radius=50.0,
    sphere_outer_radius=100.0
)
print("Geometry: Fuel shell (r=50-100 cm) with inner graphite moderator")

# Find burnable material
burnable_mat = None
for mat in materials:
    if mat.name and ('spent' in mat.name.lower() or 'uo2' in mat.name.lower() or 'fuel' in mat.name.lower()):
        burnable_mat = mat
        break

if burnable_mat is None:
    print("ERROR: Could not find burnable material!")
    sys.exit(1)

# Mark material as depletable (best practice for fixed-source depletion)
burnable_mat.depletable = True

print(f"Burnable material: {burnable_mat.name} (ID={burnable_mat.id})")
print(f"Initial volume: {burnable_mat.volume:.4e} cm³")

# ==================== Create Point Source ====================

print("\nCreating D-D fusion point source...")
source = openmc.IndependentSource()
source.space = openmc.stats.Point((0, 0, 0))  # Point source at origin
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([source_energy], [1.0])
source.strength = source_rate

# ==================== Setup Transport Settings ====================

settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.particles = particles_per_batch
settings.batches = num_batches
settings.inactive = num_inactive
settings.source = source
settings.seed = 42

print(f"\nTransport settings:")
print(f"  Particles per batch: {particles_per_batch:,}")
print(f"  Active batches: {num_batches - num_inactive}")
print(f"  Total active particles: {particles_per_batch * (num_batches - num_inactive):,}")

# ==================== Create Reaction Rate Tallies ====================

print("\nSetting up reaction rate tallies...")

# Material filter for the fuel
mat_filter = openmc.MaterialFilter([burnable_mat])

# Key fissile/fertile isotopes for accurate depletion
key_nuclides = ['U235', 'U238', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242']

# Create tallies for all important reactions
tallies = openmc.Tallies()

# Nuclide-specific fission reactions (for accurate depletion)
fission_nuclide_tally = openmc.Tally(name='fission_nuclide')
fission_nuclide_tally.filters = [mat_filter]
fission_nuclide_tally.nuclides = key_nuclides  # Specify nuclides directly on tally
fission_nuclide_tally.scores = ['fission', 'nu-fission']
tallies.append(fission_nuclide_tally)

# Nuclide-specific capture (n,gamma)
capture_nuclide_tally = openmc.Tally(name='capture_nuclide')
capture_nuclide_tally.filters = [mat_filter]
capture_nuclide_tally.nuclides = key_nuclides  # Specify nuclides directly on tally
capture_nuclide_tally.scores = ['(n,gamma)']
tallies.append(capture_nuclide_tally)

# Aggregate fission reactions (for power calculation)
fission_tally = openmc.Tally(name='fission')
fission_tally.filters = [mat_filter]
fission_tally.scores = ['fission', 'nu-fission']
tallies.append(fission_tally)

# Aggregate capture (n,gamma)
capture_tally = openmc.Tally(name='capture')
capture_tally.filters = [mat_filter]
capture_tally.scores = ['(n,gamma)']
tallies.append(capture_tally)

# (n,2n) and (n,3n) reactions
n2n_tally = openmc.Tally(name='n2n')
n2n_tally.filters = [mat_filter]
n2n_tally.scores = ['(n,2n)', '(n,3n)']
tallies.append(n2n_tally)

# Flux tally (for verification)
flux_tally = openmc.Tally(name='flux')
flux_tally.filters = [mat_filter]
flux_tally.scores = ['flux']
tallies.append(flux_tally)

print("Tallies created: nuclide-specific fission/capture, aggregate fission/capture, (n,2n), (n,3n), flux")

# Add 3D mesh tallies for spatial visualization
print("\nSetting up 3D mesh tallies for spatial visualization...")
mesh = openmc.RegularMesh()
mesh.lower_left = [-110.0, -110.0, -110.0]
mesh.upper_right = [110.0, 110.0, 110.0]
mesh.dimension = [44, 44, 44]

mesh_filter = openmc.MeshFilter(mesh)

# 3D heating tally
heating_tally_3d = openmc.Tally(name='3d_heating_tally')
heating_tally_3d.filters = [mesh_filter]
heating_tally_3d.scores = ['heating']
tallies.append(heating_tally_3d)

# 3D flux tally
flux_tally_3d = openmc.Tally(name='3d_flux_tally')
flux_tally_3d.filters = [mesh_filter]
flux_tally_3d.scores = ['flux']
tallies.append(flux_tally_3d)

# 3D fission tally
fission_tally_3d = openmc.Tally(name='3d_fission_tally')
fission_tally_3d.filters = [mesh_filter]
fission_tally_3d.scores = ['fission']
tallies.append(fission_tally_3d)

print("3D mesh tallies created: heating, flux, fission (44x44x44 mesh)")

# ==================== Load Depletion Chain ====================

print(f"\nLoading depletion chain from {CHAIN_FILE}...")
chain = openmc.deplete.Chain.from_xml(CHAIN_FILE)
print(f"Chain loaded: {len(chain.nuclides)} nuclides")

# ==================== Depletion Loop ====================

print("\n" + "="*70)
print("STARTING FIXED-SOURCE DEPLETION SIMULATION - GRAPHITE MODERATOR")
print("="*70)

# Storage for results
time_points = [0.0]  # Start at t=0
keff_values = []
compositions = []  # Store material compositions at each step

# Save initial composition
initial_composition = {}
for nuclide, density in burnable_mat.get_nuclide_atom_densities().items():
    initial_composition[nuclide] = density * burnable_mat.volume * 1e24

compositions.append(initial_composition.copy())

start_time = time.time()

for step in range(num_steps):
    print(f"\n{'='*70}")
    print(f"STEP {step+1}/{num_steps}: t = {time_points[-1]/(24*3600):.1f} days")
    print(f"{'='*70}")
    
    # Create model for this step
    model = openmc.Model(
        geometry=geometry,
        materials=materials,
        settings=settings,
        tallies=tallies
    )
    
    # Export and run (keep files in step directory)
    step_dir = os.path.join(output_dir, f"step_{step}")
    os.makedirs(step_dir, exist_ok=True)
    
    model.export_to_xml(step_dir)
    
    print(f"\nRunning transport calculation...")
    model.run(cwd=step_dir)
    
    print("Transport complete. Extracting reaction rates...")
    
    # Open statepoint and extract tallies (use absolute path)
    sp_path = os.path.join(step_dir, f'statepoint.{num_batches}.h5')
    with openmc.StatePoint(sp_path) as sp:
        # Get tallies
        fission_t = sp.get_tally(name='fission')
        capture_t = sp.get_tally(name='capture')
        n2n_t = sp.get_tally(name='n2n')
        flux_t = sp.get_tally(name='flux')
        
        # Extract mean values (reactions per source particle)
        fission_rate_per_source = fission_t.get_values(scores=['fission'])[0, 0, 0]
        nu_fission_rate_per_source = fission_t.get_values(scores=['nu-fission'])[0, 0, 0]
        capture_rate_per_source = capture_t.get_values(scores=['(n,gamma)'])[0, 0, 0]
        n2n_rate_per_source = n2n_t.get_values(scores=['(n,2n)'])[0, 0, 0]
        
        flux_per_source = flux_t.mean[0, 0, 0]
        
        # CRITICAL: In OpenMC 0.15+, tallies in fixed-source mode are internally
        # scaled by source.strength! To get per-particle rates, divide by source.strength
        # Reference: https://github.com/openmc-dev/openmc/pull/1986
        
        # Convert from (per-particle × source.strength) to per-particle
        fission_per_neutron = fission_rate_per_source / source.strength
        nu_fission_per_neutron = nu_fission_rate_per_source / source.strength  
        capture_per_neutron = capture_rate_per_source / source.strength
        n2n_per_neutron = n2n_rate_per_source / source.strength
        flux_per_neutron = flux_per_source / source.strength
        
        print(f"\nPer-source-neutron rates:")
        print(f"  Fission:     {fission_per_neutron:.4e} fissions/neutron")
        print(f"  Nu-fission:  {nu_fission_per_neutron:.4e} neutrons/neutron")
        print(f"  Capture:     {capture_per_neutron:.4e} captures/neutron")
        
        # Now scale to absolute rates (reactions/second)
        fission_rate = fission_per_neutron * source_rate  # fissions/s
        nu_fission_rate = nu_fission_per_neutron * source_rate  # neutrons/s from fission
        capture_rate = capture_per_neutron * source_rate  # captures/s
        n2n_rate = n2n_per_neutron * source_rate  # reactions/s
        flux = flux_per_neutron * source_rate  # n·cm/s (total flux in material)
        
        print(f"\nReaction rates [per second]:")
        print(f"  Fission:     {fission_rate:.4e} fissions/s")
        print(f"  Nu-fission:  {nu_fission_rate:.4e} neutrons/s (from fission)")
        print(f"  Capture:     {capture_rate:.4e} captures/s")
        print(f"  (n,2n):      {n2n_rate:.4e} reactions/s")
        print(f"  Total flux:  {flux:.4e} n·cm/s")
        
        # Calculate k_source (multiplication factor)
        k_source = nu_fission_rate / source_rate
        print(f"\n  k_source (multiplication): {k_source:.4f}")
        print(f"  (Each source neutron produces {k_source:.4f} fission neutrons)")
        
        # Calculate fission power
        energy_per_fission = 200e6 * 1.60218e-19  # 200 MeV in Joules
        power_watts = fission_rate * energy_per_fission
        print(f"\n  Fission power: {power_watts:.2e} W = {power_watts/1e3:.2f} kW = {power_watts/1e6:.4f} MW")
        
    # ==================== PHASE 2: Nuclide-Specific Depletion ====================
    
    print(f"\n{'='*70}")
    print("CALCULATING NUCLIDE-SPECIFIC DEPLETION")
    print(f"{'='*70}")
    print("\nUsing nuclide-specific reaction rate tallies for accurate depletion")
    
    # Get current composition
    current_comp = burnable_mat.get_nuclide_atom_densities()
    current_atoms = {}
    total_atoms = 0
    for nuc, density in current_comp.items():
        atoms = density * burnable_mat.volume * 1e24
        current_atoms[nuc] = atoms
        total_atoms += atoms
    
    print(f"\nCurrent total atoms: {total_atoms:.3e}")
    print(f"Time step: {step_duration_seconds:.2e} seconds ({days_per_step:.0f} days)")
    
    # Extract nuclide-specific reaction rates from tallies
    with openmc.StatePoint(sp_path) as sp:
        fission_nuc_t = sp.get_tally(name='fission_nuclide')
        capture_nuc_t = sp.get_tally(name='capture_nuclide')
    
    # Build nuclide-specific reaction rate dictionary
    nuclide_reactions = {}
    for i, nuc in enumerate(key_nuclides):
        if nuc in current_atoms:
            # Extract reaction rates (these are scaled by source.strength in OpenMC 0.15+)
            # When using tally.nuclides, indexing is [filter_bin, nuclide_bin, score_bin]
            fission_tally_value = fission_nuc_t.get_values(scores=['fission'])[0, i, 0]
            capture_tally_value = capture_nuc_t.get_values(scores=['(n,gamma)'])[0, i, 0]
            
            # Convert from (per-particle × source.strength) to per-particle
            fission_per_neutron = fission_tally_value / source.strength
            capture_per_neutron = capture_tally_value / source.strength
            
            # Convert to absolute rates (reactions/s)
            fission_rate_nuc = fission_per_neutron * source_rate
            capture_rate_nuc = capture_per_neutron * source_rate
            
            nuclide_reactions[nuc] = {
                'fission': fission_rate_nuc,
                'capture': capture_rate_nuc
            }
    
    # Apply depletion using nuclide-specific rates
    new_atoms = current_atoms.copy()
    
    print(f"\nNuclide-specific depletion:")
    for nuc in key_nuclides:
        if nuc in new_atoms and nuc in nuclide_reactions:
            # Get reaction rates for this nuclide
            fission_rate_nuc = nuclide_reactions[nuc]['fission']
            capture_rate_nuc = nuclide_reactions[nuc]['capture']
            
            # Calculate total reactions over time step
            fissions = fission_rate_nuc * step_duration_seconds
            captures = capture_rate_nuc * step_duration_seconds
            
            # Update atom count (simplified: only destruction, not production chains)
            new_atoms[nuc] -= (fissions + captures)
            
            # Prevent negative
            if new_atoms[nuc] < 0:
                new_atoms[nuc] = 0
            
            change = new_atoms[nuc] - current_atoms[nuc]
            pct_change = (change / current_atoms[nuc] * 100) if current_atoms[nuc] > 0 else 0
            if abs(pct_change) > 0.001:  # Only show significant changes
                print(f"  {nuc:8s}: {change:+.3e} atoms ({pct_change:+.3f}%) [F: {fissions:.2e}, C: {captures:.2e}]")
    
    # Show summary of all isotope changes
    print(f"\nTop 10 overall isotope changes:")
    changes = {}
    for nuc in current_atoms:
        if nuc in new_atoms:
            change = new_atoms[nuc] - current_atoms[nuc]
            if abs(change) > 0:
                changes[nuc] = change
    
    sorted_changes = sorted(changes.items(), key=lambda x: abs(x[1]), reverse=True)[:10]
    for nuc, change in sorted_changes:
        pct_change = (change / current_atoms[nuc] * 100) if current_atoms[nuc] > 0 else 0
        print(f"  {nuc:8s}: {change:+.3e} atoms ({pct_change:+.3f}%)")
    
    # Update material composition for next step
    print(f"\nUpdating material composition for next transport step...")
    new_total_density = sum(new_atoms.values()) / (burnable_mat.volume * 1e24)
    burnable_mat.set_density('atom/b-cm', new_total_density)
    
    # Rebuild material with new composition
    # First, clear all nuclides (nuclides is a list of NuclideTuple objects)
    for nuc_tuple in list(burnable_mat.nuclides):
        burnable_mat.remove_nuclide(nuc_tuple.name)
    
    # Add updated nuclides
    for nuc, atoms in new_atoms.items():
        density_atom_per_bcm = atoms / (burnable_mat.volume * 1e24)
        if density_atom_per_bcm > 1e-20:  # Only add if non-negligible
            burnable_mat.add_nuclide(nuc, density_atom_per_bcm)
    
    # Store composition
    compositions.append(new_atoms.copy())
    
    print(f"\nAdvancing time by {days_per_step:.0f} days...")
    time_points.append(time_points[-1] + step_duration_seconds)

elapsed_time = time.time() - start_time

print("\n" + "="*70)
print("SIMULATION COMPLETE - GRAPHITE MODERATOR STUDY")
print("="*70)
print(f"Total time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
print(f"\nResults saved to: {output_dir}/")

print(f"\nCompleted {num_steps} depletion steps")
print(f"Final time: {time_points[-1]/(24*3600):.1f} days")
print("="*70)

# Save results summary
with open(os.path.join(output_dir, 'depletion_summary_graphite.txt'), 'w') as f:
    f.write("D-D FUSION-FISSION HYBRID REACTOR - DEPLETION RESULTS (GRAPHITE MODERATOR)\n")
    f.write("="*70 + "\n\n")
    f.write(f"Moderator: Nuclear-grade Graphite at 1.7 g/cm³\n")
    f.write(f"Source: D-D fusion at {source_energy/1e6:.2f} MeV, {source_rate:.2e} n/s\n")
    f.write(f"Time steps: {num_steps}\n")
    f.write(f"Total time: {time_points[-1]/(24*3600):.1f} days\n")
    f.write(f"\nInitial composition stored\n")
    f.write(f"Final composition stored\n")
    f.write(f"\nKey isotope evolution:\n")
    
    for nuc in ['U235', 'U238', 'Pu239', 'Pu240', 'Pu241']:
        if nuc in compositions[0] and nuc in compositions[-1]:
            initial = compositions[0][nuc]
            final = compositions[-1][nuc]
            change = final - initial
            pct = (change / initial * 100) if initial > 0 else 0
            f.write(f"  {nuc}: {initial:.3e} -> {final:.3e} ({pct:+.2f}%)\n")

print(f"\nSummary saved to: {os.path.join(output_dir, 'depletion_summary_graphite.txt')}")
print("\nNOTE: This simulation uses graphite moderator instead of D2O or CO2")
print("Compare results with other moderator studies in ../results_manual_depletion/")
print("="*70)
