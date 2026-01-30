"""
FILE: spatial_isotope_plots.py
Visualize isotope transmutations and fissions in the x,y plane

This script creates 2D plots showing:
1. Spatial distribution of fission rates in x,y plane
2. Spatial distribution of flux and heating
3. Isotope transmutation changes mapped to spatial locations
4. Combined visualization of nuclear reactions in the reactor geometry
"""

import openmc
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

# ======================== Configuration ========================
RESULTS_DIR = "results/iteration3_results_tru_fuel"
DEPLETION_FILE = f"{RESULTS_DIR}/depletion_results.h5"
OUTPUT_DIR = "spatial_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Reactor geometry parameters (from reactor_Depleted_FixedSource_3.py)
INNER_RADIUS = 50.0  # cm
OUTER_RADIUS = 100.0  # cm
MESH_BOUNDS = [-110.0, -110.0, -110.0, 110.0, 110.0, 110.0]  # Updated to cover full fuel sphere
MESH_DIM = (22, 22, 22)  # Updated mesh dimensions
PLOT_BOUNDS = [-110.0, 110.0, -110.0, 110.0]  # Plot range matches mesh

# ======================== Data Extraction ========================

def load_statepoint_data(statepoint_file, z_slice_index=None):
    """
    Load spatial tally data from a statepoint file.
    
    Args:
        statepoint_file: Path to statepoint .h5 file
        z_slice_index: Index of z-slice to extract (default: middle slice)
    
    Returns:
        dict with heating and flux data for the z-slice
    """
    print(f"Loading {statepoint_file}...")
    sp = openmc.StatePoint(statepoint_file)
    
    # Get heating tally
    heating_tally = sp.get_tally(name='3d_heating_tally')
    heating_data = heating_tally.get_slice(scores=['heating']).mean
    heating_3d = heating_data.reshape(MESH_DIM, order='F')
    
    # Get mesh coordinates
    mesh = heating_tally.filters[0].mesh
    
    # Calculate bin edges from mesh bounds
    nx, ny, nz = MESH_DIM
    x_edges = np.linspace(mesh.lower_left[0], mesh.upper_right[0], nx + 1)
    y_edges = np.linspace(mesh.lower_left[1], mesh.upper_right[1], ny + 1)
    z_edges = np.linspace(mesh.lower_left[2], mesh.upper_right[2], nz + 1)
    
    # Default to middle z-slice
    if z_slice_index is None:
        z_slice_index = MESH_DIM[2] // 2
    
    # Extract 2D slice at z_slice_index
    heating_2d = heating_3d[:, :, z_slice_index]
    
    return {
        'x_edges': x_edges,
        'y_edges': y_edges,
        'heating_2d': heating_2d,
        'z_value': z_edges[z_slice_index],
        'z_index': z_slice_index,
        'heating_3d': heating_3d,
        'mesh_x': x_edges,
        'mesh_y': y_edges,
        'mesh_z': z_edges
    }


def load_depletion_data(h5_file):
    """
    Load depletion results showing isotope evolution.
    
    Returns:
        dict with time, nuclides, and atom counts
    """
    print(f"Loading depletion data from {h5_file}...")
    with h5py.File(h5_file, 'r') as f:
        time = f['time'][:]
        nuclides_group = f['nuclides']
        nuclides = list(nuclides_group.keys())
        number = f['number'][:]  # atoms/(barn·cm) per material zone
        
    return {
        'time': time,
        'nuclides': nuclides,
        'number': number
    }


def calculate_transmutation_changes(depletion_data, volume_cm3=3.6652e6):
    """
    Calculate net transmutation changes for key isotopes.
    
    Args:
        depletion_data: Output from load_depletion_data()
        volume_cm3: Material volume
    
    Returns:
        dict mapping nuclide name to (initial_atoms, final_atoms, change)
    """
    nuclides = depletion_data['nuclides']
    number = depletion_data['number']  # shape: (n_steps, n_materials, n_burnable_nuclides, n_nuclides)
    
    # Convert atom density to total atoms
    # number is in atoms/(barn·cm), multiply by volume and 1e24 to get total atoms
    conversion = volume_cm3 * 1e24
    
    initial_atoms = number[0, 0, 0, :] * conversion
    final_atoms = number[-1, 0, 0, :] * conversion
    
    transmutations = {}
    for i, nuc in enumerate(nuclides):
        initial = initial_atoms[i]
        final = final_atoms[i]
        change = final - initial
        transmutations[nuc] = {
            'initial': initial,
            'final': final,
            'change': change,
            'percent_change': (change / initial * 100) if initial > 1e10 else 0.0
        }
    
    return transmutations


# ======================== Plotting Functions ========================

def plot_heating_distribution(data, step_number, save_path):
    """
    Plot 2D heating distribution in x,y plane.
    """
    fig, ax = plt.subplots(figsize=(10, 9))
    
    heating = data['heating_2d']
    # Mask zeros for better visualization
    heating_masked = np.ma.masked_where(heating < 1e-10, heating)
    
    # Plot with log scale
    im = ax.pcolormesh(data['x_edges'], data['y_edges'], heating_masked.T, 
                       cmap='hot', shading='auto',
                       norm=LogNorm(vmin=heating_masked.min(), vmax=heating_masked.max()))
    
    # Add reactor geometry outlines
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(INNER_RADIUS * np.cos(theta), INNER_RADIUS * np.sin(theta), 
            'c--', linewidth=2, label=f'Inner radius ({INNER_RADIUS} cm)')
    ax.plot(OUTER_RADIUS * np.cos(theta), OUTER_RADIUS * np.sin(theta), 
            'c-', linewidth=2, label=f'Outer radius ({OUTER_RADIUS} cm)')
    
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_title(f'Heating Distribution (Step {step_number}) - z = {data["z_value"]:.1f} cm', 
                 fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    ax.set_xlim(PLOT_BOUNDS[0], PLOT_BOUNDS[1])
    ax.set_ylim(PLOT_BOUNDS[2], PLOT_BOUNDS[3])
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(im, ax=ax, label='Heating (eV/source particle)')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


def plot_fission_proxy(data, step_number, save_path):
    """
    Plot fission rate proxy (heating is proportional to fission rate).
    Since heating ≈ 200 MeV per fission, this gives spatial fission distribution.
    """
    fig, ax = plt.subplots(figsize=(10, 9))
    
    # Convert heating to approximate fission rate
    # Heating is in eV/source_particle, ~200 MeV per fission
    MeV_per_fission = 200
    fission_rate_proxy = data['heating_2d'] / (MeV_per_fission * 1e6)  # convert MeV to eV
    
    # Create radial distance mask to show only fuel sphere region (50-100 cm)
    x_centers = (data['x_edges'][:-1] + data['x_edges'][1:]) / 2
    y_centers = (data['y_edges'][:-1] + data['y_edges'][1:]) / 2
    X_centers, Y_centers = np.meshgrid(x_centers, y_centers, indexing='ij')
    R = np.sqrt(X_centers**2 + Y_centers**2)
    
    # Mask: hide regions outside fuel sphere (r < 50 or r > 100)
    fuel_mask = (R < INNER_RADIUS) | (R > OUTER_RADIUS)
    fission_rate_masked = np.ma.masked_where(
        (fission_rate_proxy < 1e-15) | fuel_mask, 
        fission_rate_proxy
    )
    
    # Use logarithmic scale for better visualization of wide dynamic range
    vmin = max(fission_rate_masked.min(), 1e-15)  # Avoid log(0) issues
    vmax = fission_rate_masked.max()
    
    im = ax.pcolormesh(data['x_edges'], data['y_edges'], fission_rate_masked.T,
                       cmap='plasma', shading='auto',
                       norm=LogNorm(vmin=vmin, vmax=vmax))
    
    # Add reactor geometry
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(INNER_RADIUS * np.cos(theta), INNER_RADIUS * np.sin(theta),
            'w--', linewidth=2, alpha=0.7, label=f'Inner radius ({INNER_RADIUS} cm)')
    ax.plot(OUTER_RADIUS * np.cos(theta), OUTER_RADIUS * np.sin(theta),
            'w-', linewidth=2, alpha=0.7, label=f'Outer radius ({OUTER_RADIUS} cm)')
    
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_title(f'Fission Rate Distribution (Step {step_number}) - z = {data["z_value"]:.1f} cm',
                 fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    ax.set_xlim(PLOT_BOUNDS[0], PLOT_BOUNDS[1])
    ax.set_ylim(PLOT_BOUNDS[2], PLOT_BOUNDS[3])
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, color='white', linewidth=0.5)
    
    cbar = plt.colorbar(im, ax=ax, label='Fissions per source particle (log scale)')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


def plot_radial_profile(data, step_number, save_path):
    """
    Plot radial profile of heating/fission rate.
    Shows how fission varies with distance from center.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    heating_3d = data['heating_3d']
    
    # Create 3D meshgrid with proper bin centers
    x_centers = (data['mesh_x'][:-1] + data['mesh_x'][1:]) / 2
    y_centers = (data['mesh_y'][:-1] + data['mesh_y'][1:]) / 2
    z_centers = (data['mesh_z'][:-1] + data['mesh_z'][1:]) / 2
    
    X, Y, Z = np.meshgrid(x_centers, y_centers, z_centers, indexing='ij')
    
    # Calculate radial distance from center
    R = np.sqrt(X**2 + Y**2 + Z**2)
    
    # Flatten arrays
    r_flat = R.flatten()
    heating_flat = heating_3d.flatten()
    
    # Filter out near-zero values
    mask = heating_flat > 1e-12
    r_filtered = r_flat[mask]
    heating_filtered = heating_flat[mask]
    
    # Create radial bins
    r_bins = np.linspace(0, 120, 60)
    r_centers = (r_bins[:-1] + r_bins[1:]) / 2
    
    # Average heating in each radial bin
    heating_binned = []
    for i in range(len(r_bins)-1):
        mask = (r_filtered >= r_bins[i]) & (r_filtered < r_bins[i+1])
        if mask.any():
            heating_binned.append(np.mean(heating_filtered[mask]))
        else:
            heating_binned.append(0)
    
    ax.plot(r_centers, heating_binned, 'b-o', linewidth=2, markersize=4)
    
    # Mark reactor boundaries
    ax.axvline(INNER_RADIUS, color='red', linestyle='--', linewidth=2, 
               label=f'Inner radius ({INNER_RADIUS} cm)')
    ax.axvline(OUTER_RADIUS, color='green', linestyle='--', linewidth=2,
               label=f'Outer radius ({OUTER_RADIUS} cm)')
    
    ax.set_xlabel('Radial Distance from Center (cm)', fontsize=12)
    ax.set_ylabel('Average Heating (eV/source particle)', fontsize=12)
    ax.set_title(f'Radial Heating Profile (Step {step_number})', fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


def plot_transmutation_summary(transmutations, save_path):
    """
    Plot bar chart showing isotope transmutation changes.
    """
    # Filter for significant changes (TRU isotopes and key fission products)
    significant_changes = {}
    tru_elements = ['Pu', 'Np', 'Am', 'Cm', 'U']
    
    for nuc, data in transmutations.items():
        if data['initial'] > 1e20:  # Significant initial inventory
            # Check if it's TRU or has >5% change
            is_tru = any(nuc.startswith(elem) for elem in tru_elements)
            has_large_change = abs(data['percent_change']) > 5.0
            
            if is_tru or has_large_change:
                significant_changes[nuc] = data
    
    # Sort by absolute change
    sorted_nucs = sorted(significant_changes.items(), 
                        key=lambda x: abs(x[1]['change']), reverse=True)
    
    # Take top 20
    top_nucs = sorted_nucs[:20]
    
    nuclides = [n for n, _ in top_nucs]
    changes = [d['change'] for _, d in top_nucs]
    percent_changes = [d['percent_change'] for _, d in top_nucs]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Absolute change
    colors = ['red' if c < 0 else 'green' for c in changes]
    ax1.barh(nuclides, changes, color=colors, alpha=0.7)
    ax1.set_xlabel('Net Change in Atom Count', fontsize=12)
    ax1.set_ylabel('Isotope', fontsize=12)
    ax1.set_title('Top 20 Isotope Transmutations (Absolute)', fontsize=14, fontweight='bold')
    ax1.axvline(0, color='black', linewidth=0.8)
    ax1.grid(True, alpha=0.3, axis='x')
    ax1.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
    
    # Percent change
    colors2 = ['red' if c < 0 else 'green' for c in percent_changes]
    ax2.barh(nuclides, percent_changes, color=colors2, alpha=0.7)
    ax2.set_xlabel('Percent Change (%)', fontsize=12)
    ax2.set_title('Top 20 Isotope Transmutations (Relative)', fontsize=14, fontweight='bold')
    ax2.axvline(0, color='black', linewidth=0.8)
    ax2.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


def plot_combined_spatial_view(data_list, save_path):
    """
    Create a multi-panel view showing evolution over time.
    """
    n_steps = len(data_list)
    fig, axes = plt.subplots(1, n_steps, figsize=(6*n_steps, 5))
    
    if n_steps == 1:
        axes = [axes]
    
    for i, (data, ax) in enumerate(zip(data_list, axes)):
        heating = data['heating_2d']
        heating_masked = np.ma.masked_where(heating < 1e-10, heating)
        
        im = ax.pcolormesh(data['x_edges'], data['y_edges'], heating_masked.T,
                          cmap='hot', shading='auto',
                          norm=LogNorm(vmin=1e-8, vmax=heating_masked.max()))
        
        # Add geometry
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(INNER_RADIUS * np.cos(theta), INNER_RADIUS * np.sin(theta),
                'c--', linewidth=1.5, alpha=0.8)
        ax.plot(OUTER_RADIUS * np.cos(theta), OUTER_RADIUS * np.sin(theta),
                'c-', linewidth=1.5, alpha=0.8)
        
        ax.set_xlabel('X (cm)')
        ax.set_ylabel('Y (cm)')
        ax.set_title(f'Step {i}', fontweight='bold')
        ax.set_aspect('equal')
        ax.set_xlim(PLOT_BOUNDS[0], PLOT_BOUNDS[1])
        ax.set_ylim(PLOT_BOUNDS[2], PLOT_BOUNDS[3])
        ax.grid(True, alpha=0.2)
    
    # Add single colorbar for all subplots
    fig.colorbar(im, ax=axes, label='Heating (eV/source particle)', 
                 orientation='horizontal', pad=0.1, fraction=0.05)
    
    fig.suptitle(f'Spatial Heating Evolution - z = {data_list[0]["z_value"]:.1f} cm',
                 fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


# ======================== Main Execution ========================

def main():
    print("="*70)
    print("SPATIAL ISOTOPE TRANSMUTATION AND FISSION VISUALIZATION")
    print("="*70)
    
    # Find all statepoint files
    statepoint_files = []
    for i in range(20):  # Check up to 20 steps
        sp_file = f"{RESULTS_DIR}/openmc_simulation_n{i}.h5"
        if os.path.exists(sp_file):
            statepoint_files.append((i, sp_file))
    
    if not statepoint_files:
        print(f"ERROR: No statepoint files found in {RESULTS_DIR}")
        return
    
    print(f"\nFound {len(statepoint_files)} statepoint files")
    
    # Load spatial data for each step
    spatial_data = []
    for step, sp_file in statepoint_files:
        try:
            data = load_statepoint_data(sp_file, z_slice_index=10)  # Middle slice
            spatial_data.append((step, data))
        except Exception as e:
            print(f"Warning: Could not load {sp_file}: {e}")
    
    print(f"Successfully loaded {len(spatial_data)} spatial datasets")
    
    # Load depletion data for transmutations
    if os.path.exists(DEPLETION_FILE):
        depletion_data = load_depletion_data(DEPLETION_FILE)
        transmutations = calculate_transmutation_changes(depletion_data)
        print(f"Loaded depletion data with {len(transmutations)} isotopes")
    else:
        print(f"Warning: Depletion file not found: {DEPLETION_FILE}")
        transmutations = None
    
    # Generate plots
    print("\n" + "="*70)
    print("GENERATING PLOTS")
    print("="*70)
    
    # 1. Individual step plots
    for step, data in spatial_data:
        print(f"\nProcessing Step {step}...")
        
        # Heating distribution
        plot_heating_distribution(data, step, 
                                 f"{OUTPUT_DIR}/heating_step{step}.png")
        
        # Fission rate proxy
        plot_fission_proxy(data, step,
                          f"{OUTPUT_DIR}/fission_rate_step{step}.png")
        
        # Radial profile
        plot_radial_profile(data, step,
                           f"{OUTPUT_DIR}/radial_profile_step{step}.png")
    
    # 2. Combined evolution plot
    if len(spatial_data) > 1:
        print("\nCreating combined evolution plot...")
        data_only = [d for _, d in spatial_data]
        plot_combined_spatial_view(data_only,
                                  f"{OUTPUT_DIR}/evolution_combined.png")
    
    # 3. Transmutation summary
    if transmutations:
        print("\nCreating transmutation summary...")
        plot_transmutation_summary(transmutations,
                                  f"{OUTPUT_DIR}/transmutation_summary.png")
    
    print("\n" + "="*70)
    print("VISUALIZATION COMPLETE")
    print("="*70)
    print(f"\nAll plots saved to: {OUTPUT_DIR}/")
    print(f"\nGenerated {len(os.listdir(OUTPUT_DIR))} plots:")
    for filename in sorted(os.listdir(OUTPUT_DIR)):
        print(f"  - {filename}")
    print("\n" + "="*70)


if __name__ == "__main__":
    main()
