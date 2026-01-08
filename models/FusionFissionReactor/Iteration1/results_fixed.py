"""
Post-processing script for depletion simulation results.
Reads depletion_results.h5 directly using h5py and generates summary plots/data.
"""
import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path

# Constants
MEV_PER_FISSION = 200.0  # Average energy per fission
MEV_TO_JOULES = 1.60218e-13
MEV_TO_MJ = MEV_TO_JOULES * 1e-6

RESULTS_FILE = "results/depletion_results.h5"
OUTPUT_DIR = "results"

def read_depletion_results(filename):
    """Read depletion results HDF5 file and return key datasets."""
    with h5py.File(filename, 'r') as f:
        times = f['time'][:]  # shape: (n_steps, 2) = [seconds, ?? ]
        eigenvalues = f['eigenvalues'][:]  # shape: (n_steps, 1, 2) = [value, uncertainty]
        source_rates = f['source_rate'][:]  # shape: (n_steps, 1) = [neutrons/s]
        nuclides = list(f['nuclides'].keys())
        materials = list(f['materials'].keys())
        
        # Read nuclide number densities for material '1' if available
        number_densities = None
        if '1' in f['materials']:
            try:
                mat_group = f['materials']['1']
                if 'nuclides' in mat_group:
                    number_densities = mat_group['nuclides'][:]  # shape: (n_steps, n_nuclides)
            except Exception as e:
                print(f"Warning: could not read number densities: {e}")
        
        # Read number array (nuclide counts or densities)
        number_array = f['number'][:]  # shape: (n_steps, 1, 1, n_nuclides)
        
        return {
            'times': times,
            'eigenvalues': eigenvalues,
            'source_rates': source_rates,
            'nuclides': nuclides,
            'materials': materials,
            'number_array': number_array,
            'number_densities': number_densities,
        }

def extract_burnup_metrics(data):
    """Extract burnup and power metrics from depletion data."""
    times_s = data['times'][:, 0]  # Time in seconds
    eigenvalues = data['eigenvalues'][:, 0, 0]  # k-eff values
    source_rates = data['source_rates'][:, 0]  # neutrons/s
    
    # Calculate power from source rate
    # Assumption: each neutron released comes from ~200 MeV fission energy
    # Power [W] = fission_rate [fissions/s] * energy_per_fission [J/fission]
    # Using source_rate as proxy for fission rate
    power_watts = source_rates * MEV_PER_FISSION * MEV_TO_JOULES
    power_mw = power_watts / 1e6  # Convert to MW
    
    # Calculate cumulative energy
    dt_s = np.diff(times_s)
    dt_s = np.append(dt_s, dt_s[-1] if len(dt_s) > 0 else 0)  # pad last element
    
    # Cumulative neutrons and energy
    cumulative_neutrons = np.cumsum(source_rates * dt_s)
    cumulative_energy_j = np.cumsum(power_watts * dt_s)  # Joules
    cumulative_energy_mwh = cumulative_energy_j / (3.6e9)  # Convert to MWh
    
    return {
        'times_s': times_s,
        'times_h': times_s / 3600.0,
        'times_d': times_s / 86400.0,
        'eigenvalues': eigenvalues,
        'source_rates': source_rates,
        'power_watts': power_watts,
        'power_mw': power_mw,
        'cumulative_neutrons': cumulative_neutrons,
        'cumulative_energy_j': cumulative_energy_j,
        'cumulative_energy_mwh': cumulative_energy_mwh,
    }

def extract_nuclide_inventories(data):
    """Extract key nuclide inventories over time."""
    nuclides = data['nuclides']
    number_array = data['number_array']  # shape: (n_steps, 1, 1, n_nuclides)
    n_steps = number_array.shape[0]
    
    # Flatten to (n_steps, n_nuclides)
    inventory = number_array[:, 0, 0, :]
    
    # Extract key actinides/fission products (with actual values, not forced zeros)
    key_nuclides = ['U235', 'U238', 'Pu239', 'Pu240', 'Pu241', 'Cm244', 'Xe135', 'Sm149']
    results = {}
    
    for nuc in key_nuclides:
        if nuc in nuclides:
            idx = nuclides.index(nuc)
            results[nuc] = inventory[:, idx]
    
    # Also include any nuclides with significant initial or final inventory
    # (helps discover what was actually in the material)
    for i, nuc in enumerate(nuclides):
        initial_val = inventory[0, i]
        final_val = inventory[-1, i]
        # Include if initial > threshold OR final > threshold
        if initial_val > 1e20 or final_val > 1e20:
            if nuc not in results:
                results[nuc] = inventory[:, i]
    
    return results

def generate_summary_report(data, metrics, inventories):
    """Generate a text summary of depletion results."""
    times_s = metrics['times_s']
    eigenvalues = metrics['eigenvalues']
    
    report = []
    report.append("=" * 70)
    report.append("DEPLETION SIMULATION SUMMARY")
    report.append("=" * 70)
    report.append(f"\nSimulation Duration: {times_s[-1]:.2e} seconds ({times_s[-1]/86400:.2f} days)")
    report.append(f"Number of Depletion Steps: {len(times_s)}")
    report.append(f"\nInitial k-eff: {eigenvalues[0]:.6f}")
    report.append(f"Final k-eff: {eigenvalues[-1]:.6f}")
    report.append(f"k-eff Change: {eigenvalues[-1] - eigenvalues[0]:+.6f}")
    
    report.append(f"\nInitial Source Rate: {metrics['source_rates'][0]:.4e} n/s")
    report.append(f"Final Source Rate: {metrics['source_rates'][-1]:.4e} n/s")
    
    report.append(f"\nCumulative Neutrons: {metrics['cumulative_neutrons'][-1]:.4e}")
    
    report.append("\n" + "-" * 70)
    report.append("POWER GENERATION")
    report.append("-" * 70)
    report.append(f"Initial Power: {metrics['power_mw'][0]:.4e} MW ({metrics['power_watts'][0]:.4e} W)")
    report.append(f"Final Power: {metrics['power_mw'][-1]:.4e} MW ({metrics['power_watts'][-1]:.4e} W)")
    report.append(f"Total Energy Generated: {metrics['cumulative_energy_mwh'][-1]:.4e} MWh")
    report.append(f"                        {metrics['cumulative_energy_j'][-1]:.4e} Joules")
    avg_power = np.mean(metrics['power_mw'])
    report.append(f"Average Power Output: {avg_power:.4e} MW")
    
    report.append("\n" + "-" * 70)
    report.append("KEY NUCLIDE INVENTORIES")
    report.append("-" * 70)
    
    # Separate nuclides into categories for better readability
    produced = []  # Initial = 0, Final > 0
    depleted = []  # Initial > 0, decreased
    increased = [] # Initial > 0, increased
    stable = []    # Initial > 0, relatively stable
    
    for nuc, inv in inventories.items():
        initial = inv[0]
        final = inv[-1]
        
        if initial < 1e-10:  # Essentially zero initially
            if final > 1e-10:
                produced.append((nuc, initial, final))
        else:  # Non-zero initially
            rel_change = (final - initial) / initial
            if rel_change < -0.01:  # Decreased by >1%
                depleted.append((nuc, initial, final, 100.0 * rel_change))
            elif rel_change > 0.01:  # Increased by >1%
                increased.append((nuc, initial, final, 100.0 * rel_change))
            else:  # Relatively stable (within ±1%)
                stable.append((nuc, initial, final, 100.0 * rel_change))
    
    if produced:
        report.append("\nProduced during simulation (initially zero):")
        for nuc, initial, final in sorted(produced, key=lambda x: x[2], reverse=True)[:15]:
            report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e}")
    
    if depleted:
        report.append("\nDepleted (decreased from initial):")
        for nuc, initial, final, pct in sorted(depleted, key=lambda x: abs(x[3]), reverse=True)[:15]:
            report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} ({pct:+6.1f}%)")
    
    if increased:
        report.append("\nIncreased from initial:")
        for nuc, initial, final, pct in sorted(increased, key=lambda x: x[3], reverse=True)[:15]:
            report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} ({pct:+6.1f}%)")
    
    if stable:
        report.append("\nRelatively stable (±1%):")
        for nuc, initial, final, pct in sorted(stable, key=lambda x: x[1], reverse=True)[:10]:
            report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} ({pct:+6.1f}%)")
    
    report.append("\n" + "=" * 70)
    
    return "\n".join(report)

def plot_metrics(metrics, output_dir='results'):
    """Generate plots of depletion metrics including power generation."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a 3-panel figure
    fig = plt.figure(figsize=(14, 10))
    
    times_h = metrics['times_h']
    times_d = metrics['times_d']
    eigenvalues = metrics['eigenvalues']
    source_rates = metrics['source_rates']
    power_mw = metrics['power_mw']
    cumulative_energy_mwh = metrics['cumulative_energy_mwh']
    
    # Panel 1: Eigenvalue
    ax1 = plt.subplot(3, 2, 1)
    ax1.plot(times_d, eigenvalues, 'b.-', linewidth=2, markersize=6)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('k-eff')
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Eigenvalue vs. Time')
    
    # Panel 2: Source Rate
    ax2 = plt.subplot(3, 2, 2)
    ax2.semilogy(times_d, source_rates, 'r.-', linewidth=2, markersize=6)
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Source Rate (n/s)')
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Neutron Source Rate vs. Time')
    
    # Panel 3: Power Output (linear)
    ax3 = plt.subplot(3, 2, 3)
    ax3.plot(times_d, power_mw, 'g.-', linewidth=2, markersize=6)
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('Power Output (MW)')
    ax3.grid(True, alpha=0.3)
    ax3.set_title('Power Generation vs. Time')
    
    # Panel 4: Power Output (log scale)
    ax4 = plt.subplot(3, 2, 4)
    ax4.semilogy(times_d, power_mw, 'g.-', linewidth=2, markersize=6)
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Power Output (MW, log scale)')
    ax4.grid(True, alpha=0.3, which='both')
    ax4.set_title('Power Output (Log Scale)')
    
    # Panel 5: Cumulative Energy
    ax5 = plt.subplot(3, 2, 5)
    ax5.semilogy(times_d, cumulative_energy_mwh, 'm.-', linewidth=2, markersize=6)
    ax5.set_xlabel('Time (days)')
    ax5.set_ylabel('Cumulative Energy (MWh, log scale)')
    ax5.grid(True, alpha=0.3, which='both')
    ax5.set_title('Cumulative Energy Generated')
    
    # Panel 6: Summary text
    ax6 = plt.subplot(3, 2, 6)
    ax6.axis('off')
    summary_text = (
        f"Simulation Summary:\n"
        f"Duration: {times_d[-1]:.1f} days\n"
        f"Initial Power: {power_mw[0]:.3e} MW\n"
        f"Final Power: {power_mw[-1]:.3e} MW\n"
        f"Average Power: {np.mean(power_mw):.3e} MW\n"
        f"Total Energy: {cumulative_energy_mwh[-1]:.3e} MWh\n"
        f"k-eff: {eigenvalues[0]:.4f} → {eigenvalues[-1]:.4f}"
    )
    ax6.text(0.1, 0.5, summary_text, fontsize=12, family='monospace',
             verticalalignment='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'depletion_metrics.png'), dpi=150)
    print(f"Saved plot: {os.path.join(output_dir, 'depletion_metrics.png')}")
    plt.close()
    
    # Create a separate power-only detailed plot
    fig_power, ax_power = plt.subplots(figsize=(12, 6))
    ax_power.fill_between(times_d, power_mw, alpha=0.3, label='Power Output')
    ax_power.plot(times_d, power_mw, 'g-', linewidth=2.5, marker='o', markersize=5, label='Power (MW)')
    ax_power.set_xlabel('Time (days)', fontsize=12)
    ax_power.set_ylabel('Power Output (MW)', fontsize=12)
    ax_power.set_title('Detailed Power Generation Profile', fontsize=14, fontweight='bold')
    ax_power.grid(True, alpha=0.3)
    ax_power.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'power_generation.png'), dpi=150)
    print(f"Saved plot: {os.path.join(output_dir, 'power_generation.png')}")
    plt.close()

def main():
    print("Reading depletion results...")
    data = read_depletion_results(RESULTS_FILE)
    print(f"  - {len(data['times'])} depletion steps")
    print(f"  - {len(data['nuclides'])} nuclides tracked")
    print(f"  - Materials: {data['materials']}")
    
    print("\nExtracting burnup metrics...")
    metrics = extract_burnup_metrics(data)
    
    print("Extracting nuclide inventories...")
    inventories = extract_nuclide_inventories(data)
    
    print("\nGenerating summary report...")
    report = generate_summary_report(data, metrics, inventories)
    print(report)
    
    # Save report
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    report_file = os.path.join(OUTPUT_DIR, 'depletion_summary.txt')
    with open(report_file, 'w') as f:
        f.write(report)
    print(f"\nSaved report: {report_file}")
    
    # Generate plots
    print("\nGenerating plots...")
    plot_metrics(metrics, OUTPUT_DIR)
    
    # Save metrics as CSV
    csv_file = os.path.join(OUTPUT_DIR, 'depletion_metrics.csv')
    df = pd.DataFrame({
        'time_s': metrics['times_s'],
        'time_h': metrics['times_h'],
        'time_d': metrics['times_d'],
        'k_eff': metrics['eigenvalues'],
        'source_rate_n_per_s': metrics['source_rates'],
        'power_watts': metrics['power_watts'],
        'power_mw': metrics['power_mw'],
        'cum_neutrons': metrics['cumulative_neutrons'],
        'cum_energy_joules': metrics['cumulative_energy_j'],
        'cum_energy_mwh': metrics['cumulative_energy_mwh'],
    })
    df.to_csv(csv_file, index=False)
    print(f"Saved metrics CSV: {csv_file}")
    
    print("\n✓ Post-processing complete.")

if __name__ == '__main__':
    main()
