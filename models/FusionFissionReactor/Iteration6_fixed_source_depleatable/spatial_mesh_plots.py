"""
FILE: spatial_mesh_plots.py
Create 2D spatial plots from mesh tallies in depletion results

This script creates:
1. Heating distribution in x-y plane
2. Flux distribution in x-y plane  
3. Fission rate distribution in x-y plane
4. Radial profiles of heating, flux, and fission
5. Multi-panel comparison plots
"""

import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
from PIL import Image

# ======================== Configuration ========================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results_manual_depletion")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "spatial_plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Reactor geometry parameters
INNER_RADIUS = 50.0  # cm
OUTER_RADIUS = 100.0  # cm
MESH_DIM = (44, 44, 44)

# Source parameters for unit conversion
SOURCE_RATE = 1.0e14  # neutrons per second (from fusion_hybrid_simulation.py)

# Calculate mesh cell volume
MESH_BOUNDS = [-110.0, -110.0, -110.0, 110.0, 110.0, 110.0]
x_spacing = (MESH_BOUNDS[3] - MESH_BOUNDS[0]) / MESH_DIM[0]
y_spacing = (MESH_BOUNDS[4] - MESH_BOUNDS[1]) / MESH_DIM[1]
z_spacing = (MESH_BOUNDS[5] - MESH_BOUNDS[2]) / MESH_DIM[2]
CELL_VOLUME_CM3 = x_spacing * y_spacing * z_spacing  # cm³

# ======================== Data Loading ========================

def find_step_directories():
    """Find all step directories."""
    steps = []
    for item in os.listdir(RESULTS_DIR):
        if item.startswith('step_') and os.path.isdir(os.path.join(RESULTS_DIR, item)):
            step_num = int(item.split('_')[1])
            steps.append((step_num, os.path.join(RESULTS_DIR, item)))
    return sorted(steps)


def load_mesh_data(step_dir, z_slice_index=None, convert_units=True):
    """
    Load 3D mesh tally data from statepoint.
    
    Args:
        step_dir: Directory containing statepoint file
        z_slice_index: Index for z-slice (default: middle)
        convert_units: If True, convert to physical units (per cm³ per second)
    
    Returns dict with heating, flux, fission as 2D slices
    """
    statepoint_files = [f for f in os.listdir(step_dir) if f.startswith('statepoint') and f.endswith('.h5')]
    
    if not statepoint_files:
        print(f"No statepoint in {step_dir}")
        return None
    
    sp_file = os.path.join(step_dir, statepoint_files[0])
    print(f"Loading {sp_file}...")
    
    try:
        sp = openmc.StatePoint(sp_file)
        
        # Get mesh tallies
        heating_tally = sp.get_tally(name='3d_heating_tally')
        flux_tally = sp.get_tally(name='3d_flux_tally')
        fission_tally = sp.get_tally(name='3d_fission_tally')
        
        # Extract 3D data
        heating_3d = heating_tally.mean.reshape(MESH_DIM, order='F')
        flux_3d = flux_tally.mean.reshape(MESH_DIM, order='F')
        fission_3d = fission_tally.mean.reshape(MESH_DIM, order='F')
        
        # Convert units from per-source-particle to per cm³ per second
        if convert_units:
            # IMPORTANT: source.strength was set to 1e14 in the simulation,
            # so OpenMC tallies already include that weighting factor.
            # We only need to normalize by cell volume to get volumetric rates.
            
            heating_3d = heating_3d / CELL_VOLUME_CM3  # eV/cm³/s
            flux_3d = flux_3d / CELL_VOLUME_CM3  # n/cm²/s  
            fission_3d = fission_3d / CELL_VOLUME_CM3  # fissions/cm³/s
        
        # Get mesh coordinates
        mesh = heating_tally.filters[0].mesh
        x_edges = np.linspace(mesh.lower_left[0], mesh.upper_right[0], MESH_DIM[0] + 1)
        y_edges = np.linspace(mesh.lower_left[1], mesh.upper_right[1], MESH_DIM[1] + 1)
        z_edges = np.linspace(mesh.lower_left[2], mesh.upper_right[2], MESH_DIM[2] + 1)
        
        # Default to middle z-slice
        if z_slice_index is None:
            z_slice_index = MESH_DIM[2] // 2
        
        # Extract 2D slices at z=0 plane
        heating_2d = heating_3d[:, :, z_slice_index]
        flux_2d = flux_3d[:, :, z_slice_index]
        fission_2d = fission_3d[:, :, z_slice_index]
        
        return {
            'x_edges': x_edges,
            'y_edges': y_edges,
            'z_edges': z_edges,
            'z_index': z_slice_index,
            'z_value': z_edges[z_slice_index],
            'heating_2d': heating_2d,
            'flux_2d': flux_2d,
            'fission_2d': fission_2d,
            'heating_3d': heating_3d,
            'flux_3d': flux_3d,
            'fission_3d': fission_3d
        }
        
    except Exception as e:
        print(f"Error loading mesh data: {e}")
        return None


# ======================== Plotting Functions ========================

def plot_2d_distribution(data, quantity_name, step_num, vmin=None, vmax=None):
    """
    Create 2D heatmap of a quantity in x-y plane.
    """
    x_edges = data['x_edges']
    y_edges = data['y_edges']
    z_value = data['z_value']
    values_2d = data[f'{quantity_name}_2d']
    
    # Create mesh grid for plotting
    X, Y = np.meshgrid(x_edges, y_edges)
    
    # Mask zero values to exclude from color scale
    values_plot = values_2d.T.copy()
    values_masked = np.ma.masked_where(values_plot <= 1e-20, values_plot)
    
    fig, ax = plt.subplots(figsize=(10, 9))
    
    # Plot with log scale, masked values will appear white/transparent
    # Use provided vmin/vmax for constant color scale, or compute from data
    if vmin is None:
        vmin = values_masked.compressed().min()
    if vmax is None:
        vmax = values_masked.max()
    
    im = ax.pcolormesh(X, Y, values_masked, 
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       cmap='viridis', shading='auto')
    
    # Add reactor geometry circles
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(INNER_RADIUS * np.cos(theta), INNER_RADIUS * np.sin(theta), 
            'r--', linewidth=2, label=f'Inner (r={INNER_RADIUS} cm)')
    ax.plot(OUTER_RADIUS * np.cos(theta), OUTER_RADIUS * np.sin(theta), 
            'r--', linewidth=2, label=f'Outer (r={OUTER_RADIUS} cm)')
    
    # Mark center (fusion source location)
    ax.plot(0, 0, 'r*', markersize=15, label='Fusion Source')
    
    ax.set_xlabel('X (cm)', fontsize=12)
    ax.set_ylabel('Y (cm)', fontsize=12)
    ax.set_title(f'{quantity_name.capitalize()} Distribution at z={z_value:.1f} cm\\n'
                 f'Step {step_num} (t={step_num*90:.0f} days)', 
                 fontsize=13, fontweight='bold')
    ax.set_aspect('equal')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Colorbar with physical units (per second, NOT per batch)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    if quantity_name == 'heating':
        cbar.set_label('Heating\n(eV per cm³ per second)', fontsize=11)
    elif quantity_name == 'flux':
        cbar.set_label('Neutron Flux\n(neutrons per cm² per second)', fontsize=11)
    elif quantity_name == 'fission':
        cbar.set_label('Fission Rate\n(fissions per cm³ per second)', fontsize=11)
    
    plt.tight_layout()
    output_file = os.path.join(OUTPUT_DIR, f'{quantity_name}_2d_step{step_num}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_radial_profiles(data, step_num):
    """
    Plot radial profiles of heating, flux, and fission.
    """
    x_edges = data['x_edges']
    y_edges = data['y_edges']
    
    # Calculate radial distance for each mesh cell
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    
    X, Y = np.meshgrid(x_centers, y_centers)
    R = np.sqrt(X**2 + Y**2)
    
    # Get z=0 slice data
    z_idx = data['z_index']
    heating_slice = data['heating_3d'][:, :, z_idx].T
    flux_slice = data['flux_3d'][:, :, z_idx].T
    fission_slice = data['fission_3d'][:, :, z_idx].T
    
    # Flatten arrays
    r_flat = R.flatten()
    heating_flat = heating_slice.flatten()
    flux_flat = flux_slice.flatten()
    fission_flat = fission_slice.flatten()
    
    # Create radial bins
    r_bins = np.linspace(0, 110, 50)
    r_centers = (r_bins[:-1] + r_bins[1:]) / 2
    
    # Average in radial bins
    heating_avg = []
    flux_avg = []
    fission_avg = []
    
    for i in range(len(r_bins)-1):
        mask = (r_flat >= r_bins[i]) & (r_flat < r_bins[i+1])
        if mask.sum() > 0:
            heating_avg.append(heating_flat[mask].mean())
            flux_avg.append(flux_flat[mask].mean())
            fission_avg.append(fission_flat[mask].mean())
        else:
            heating_avg.append(0)
            flux_avg.append(0)
            fission_avg.append(0)
    
    # Create plots
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    fig.suptitle(f'Radial Profiles at Step {step_num} (t={step_num*90:.0f} days)', 
                 fontsize=14, fontweight='bold')
    
    # Heating
    ax = axes[0]
    ax.plot(r_centers, heating_avg, 'o-', linewidth=2, markersize=4, color='red')
    ax.axvline(INNER_RADIUS, color='gray', linestyle='--', alpha=0.5, label='Fuel boundaries')
    ax.axvline(OUTER_RADIUS, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('Heating\n(eV per cm³ per second)', fontsize=10)
    ax.set_title('Heating vs Radius', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_yscale('log')
    
    # Flux
    ax = axes[1]
    ax.plot(r_centers, flux_avg, 'o-', linewidth=2, markersize=4, color='blue')
    ax.axvline(INNER_RADIUS, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(OUTER_RADIUS, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('Neutron Flux\n(neutrons per cm² per second)', fontsize=10)
    ax.set_title('Flux vs Radius', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    # Fission
    ax = axes[2]
    ax.plot(r_centers, fission_avg, 'o-', linewidth=2, markersize=4, color='green')
    ax.axvline(INNER_RADIUS, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(OUTER_RADIUS, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Radius (cm)', fontsize=11)
    ax.set_ylabel('Fission Rate\n(fissions per cm³ per second)', fontsize=10)
    ax.set_title('Fission Rate vs Radius', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    plt.tight_layout()
    output_file = os.path.join(OUTPUT_DIR, f'radial_profiles_step{step_num}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_3_panel_comparison(data, step_num, vmins=None, vmaxs=None):
    """
    Create 3-panel plot comparing heating, flux, and fission.
    
    Args:
        data: Dictionary with mesh data
        step_num: Step number
        vmins: Dict of min values for each quantity (for constant color scale)
        vmaxs: Dict of max values for each quantity (for constant color scale)
    """
    x_edges = data['x_edges']
    y_edges = data['y_edges']
    X, Y = np.meshgrid(x_edges, y_edges)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(f'Spatial Distributions at Step {step_num} (t={step_num*90:.0f} days, z=0 plane)',
                 fontsize=14, fontweight='bold')
    
    quantities = [
        ('heating', 'Heating', 'Reds', 'eV per cm³ per s'),
        ('flux', 'Neutron Flux', 'Blues', 'n per cm² per s'),
        ('fission', 'Fission Rate', 'Greens', 'fissions per cm³ per s')
    ]
    
    # Default vmin/vmax dicts if not provided
    if vmins is None:
        vmins = {}
    if vmaxs is None:
        vmaxs = {}
    
    for ax, (qty, title, cmap, unit) in zip(axes, quantities):
        values = data[f'{qty}_2d'].T.copy()
        values_masked = np.ma.masked_where(values <= 1e-20, values)
        
        # Use provided vmin/vmax for constant color scale, or compute from data
        vmin = vmins.get(qty, values_masked.compressed().min())
        vmax = vmaxs.get(qty, values_masked.max())
        
        im = ax.pcolormesh(X, Y, values_masked,
                          norm=LogNorm(vmin=vmin, vmax=vmax),
                          cmap=cmap, shading='auto')
        
        # Add geometry
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(INNER_RADIUS * np.cos(theta), INNER_RADIUS * np.sin(theta), 
                'k--', linewidth=1.5, alpha=0.7)
        ax.plot(OUTER_RADIUS * np.cos(theta), OUTER_RADIUS * np.sin(theta), 
                'k--', linewidth=1.5, alpha=0.7)
        ax.plot(0, 0, 'r*', markersize=12)
        
        ax.set_xlabel('X (cm)', fontsize=11)
        ax.set_ylabel('Y (cm)', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(unit, fontsize=10)
    
    plt.tight_layout()
    output_file = os.path.join(OUTPUT_DIR, f'3panel_comparison_step{step_num}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


# ======================== GIF Creation ========================

def create_gifs(num_steps):
    """
    Create animated GIFs from the time series of spatial plots.
    
    Args:
        num_steps: Number of timesteps
    """
    print("\n" + "="*70)
    print("CREATING ANIMATED GIFS")
    print("="*70)
    
    gif_types = [
        ('heating_2d', 'Heating Distribution'),
        ('flux_2d', 'Flux Distribution'),
        ('fission_2d', 'Fission Rate Distribution'),
        ('radial_profiles', 'Radial Profiles'),
        ('3panel_comparison', 'Three-Panel Comparison')
    ]
    
    for prefix, description in gif_types:
        # Collect all PNG files for this type
        images = []
        for step in range(num_steps):
            img_path = os.path.join(OUTPUT_DIR, f'{prefix}_step{step}.png')
            if os.path.exists(img_path):
                images.append(Image.open(img_path))
        
        if not images:
            print(f"  Skipping {description} - no images found")
            continue
        
        # Create GIF with 1 second per frame (1000 ms)
        gif_path = os.path.join(OUTPUT_DIR, f'{prefix}_evolution.gif')
        images[0].save(
            gif_path,
            save_all=True,
            append_images=images[1:],
            duration=1000,  # 1 second per frame
            loop=0  # Loop forever
        )
        print(f"  Created: {gif_path} ({len(images)} frames)")
    
    print("\n" + "="*70)


# ======================== Main Execution ========================

def main():
    """
    Main plotting function with constant color scales across timesteps.
    """
    print("="*70)
    print("SPATIAL MESH PLOTS FOR FUSION-FISSION HYBRID REACTOR")
    print("="*70)
    print(f"\nResults: {RESULTS_DIR}")
    print(f"Output: {OUTPUT_DIR}\n")
    
    # Find steps
    steps = find_step_directories()
    print(f"Found {len(steps)} steps: {[s[0] for s in steps]}\n")
    
    if not steps:
        print("ERROR: No step directories found!")
        return
    
    # FIRST PASS: Determine global min/max for constant color scales
    print("\nFirst pass: Determining global color scale ranges...")
    global_mins = {'heating': float('inf'), 'flux': float('inf'), 'fission': float('inf')}
    global_maxs = {'heating': -float('inf'), 'flux': -float('inf'), 'fission': -float('inf')}
    
    all_data = []
    for step_num, step_dir in steps:
        data = load_mesh_data(step_dir)
        if data:
            all_data.append((step_num, step_dir, data))
            for qty in ['heating', 'flux', 'fission']:
                values = data[f'{qty}_2d']
                valid_values = values[values > 1e-20]
                if len(valid_values) > 0:
                    global_mins[qty] = min(global_mins[qty], valid_values.min())
                    global_maxs[qty] = max(global_maxs[qty], valid_values.max())
    
    print("\nGlobal ranges (for constant color scales):")
    print(f"  Heating: {global_mins['heating']:.2e} - {global_maxs['heating']:.2e} eV/cm³/s")
    print(f"  Flux:    {global_mins['flux']:.2e} - {global_maxs['flux']:.2e} n/cm²/s")
    print(f"  Fission: {global_mins['fission']:.2e} - {global_maxs['fission']:.2e} fissions/cm³/s")
    
    # SECOND PASS: Create plots with constant color scales
    print("\nSecond pass: Creating plots with constant color scales...")
    for step_num, step_dir, data in all_data:
        print(f"\n{'='*70}")
        print(f"Processing Step {step_num} (t={step_num*90:.0f} days)")
        print(f"{'='*70}")
        
        # Create plots with constant color scales
        plot_2d_distribution(data, 'heating', step_num, 
                           vmin=global_mins['heating'], vmax=global_maxs['heating'])
        plot_2d_distribution(data, 'flux', step_num,
                           vmin=global_mins['flux'], vmax=global_maxs['flux'])
        plot_2d_distribution(data, 'fission', step_num,
                           vmin=global_mins['fission'], vmax=global_maxs['fission'])
        plot_radial_profiles(data, step_num)
        plot_3_panel_comparison(data, step_num, vmins=global_mins, vmaxs=global_maxs)
    
    # Create animated GIFs from the time series
    create_gifs(len(all_data))
    
    print("\n" + "="*70)
    print("SPATIAL PLOTTING COMPLETE")
    print("="*70)
    print(f"\nAll plots saved to: {OUTPUT_DIR}/")
    print("\nGenerated spatial plots:")
    spatial_files = [f for f in sorted(os.listdir(OUTPUT_DIR)) if any(x in f for x in ['_2d_', 'radial', '3panel'])]
    for f in spatial_files:
        print(f"  - {f}")
    
    print("\nGenerated animations:")
    gif_files = [f for f in sorted(os.listdir(OUTPUT_DIR)) if f.endswith('.gif')]
    for f in gif_files:
        print(f"  - {f}")
    print()


if __name__ == "__main__":
    main()
