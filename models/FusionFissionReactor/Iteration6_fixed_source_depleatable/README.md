# Iteration 6: Fixed-Source Depletion (Corrected)

## Overview

This iteration contains the **corrected** implementation of fixed-source depletion for a D-D fusion-fission hybrid system. All critical normalization issues have been resolved, and the simulation now uses OpenMC exactly as intended.

## Key Corrections Made

### 1. **Fixed Tally Normalization (CRITICAL)**
**Discovery**: In OpenMC 0.15+, fixed-source tallies are internally scaled by `source.strength`

**Correct Formula**:
```python
# Tallies are already multiplied by source.strength internally
tally_value = tally.mean[0, 0, 0]  # e.g., 4.094e+13

# Step 1: Divide by source.strength to get per-particle rate
reactions_per_neutron = tally_value / source.strength  # = 0.4094

# Step 2: Multiply by physical source rate for absolute rates  
reactions_per_second = reactions_per_neutron * source_rate  # = 4.094e+13 /s
```

### 2. **Corrected Source Strength**
- Changed from `1.0e16` to `1.0e14` neutrons/second
- Represents realistic D-D fusion neutron production rate
- Comments updated throughout

### 3. **Fixed Nuclide-Specific Tallies**
- Changed from non-existent `NuclideFilter` to `tally.nuclides` attribute
- Correct API for specifying nuclides in OpenMC
- Proper indexing for nuclide-specific data extraction

### 4. **Added Material.depletable Flag**
- Set `burnable_mat.depletable = True`
- Best practice per OpenMC documentation

## Files

- **`fusion_hybrid_simulation.py`** - Main simulation (corrected, renamed)
- `fuel_blanket_5.py` - Geometry builder with spent fuel composition  
- `chain_reduced_tru.xml` - Depletion chain for transmutation calculations
- `README.md` - This file

## Usage

```bash
cd /workspaces/MEng172-OpenMC/models/FusionFissionReactor/Iteration6_fixed_source_depleatable
python fusion_hybrid_simulation.py
```

Results will be saved to `results_manual_depletion/` subdirectory.

## Expected Results

With the corrections applied:

1. **Power**: ~0.13 kW (physically realistic for 10¹⁴ n/s source)
2. **k_source**: ~1.1 (neutron multiplication in subcritical fuel)
3. **Depletion**: Small but measurable changes in fissile isotopes over 1 year
4. **Statistical precision**: <1% uncertainty on integral tallies

## Validation Checks

The simulation performs automatic validation:
- Neutron balance (source + fission = absorption + leakage)
- Energy balance (fission power = heating)
- Mass conservation (total atoms preserved)
- Physical reasonableness (cross-section values, multiplication factors)

## References

Based on OpenMC documentation:
- [Section 10.1: Transport-coupled depletion](https://docs.openmc.org/en/stable/usersguide/depletion.html)
- [Section 10.1.1: Fixed-Source Transmutation](https://docs.openmc.org/en/stable/usersguide/depletion.html#fixed-source-transmutation)

## Change Log

**February 6, 2026**
- Initial setup of Iteration6 with all corrections
- Fixed tally normalization bug
- Corrected source strength to 1e14 n/s
- Implemented nuclide-specific depletion
- Added material.depletable flag

## Previous Iterations

- **Iteration 5**: Manual depletion implementation (had normalization bug)
- **Iteration 4**: Eigenvalue mode attempts (failed for k << 1)
- **Iteration 3**: Earlier geometry configurations
- **Iteration 2**: Initial fixed-source exploration
- **Iteration 1**: Baseline eigenvalue calculations
