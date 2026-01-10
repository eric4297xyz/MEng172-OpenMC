# Depletion Simulation Diagnosis Summary

## Problem Statement
The depletion results show that all uranium (U-235, U-238) was depleted to near-zero values after only 182.5 days, which is physically unreasonable.

## Key Findings

### ✅ What's Working Correctly:
1. **Material Definition**: The spent fuel material in `fuel_blanket_2.py` is correctly defined with:
   - U-238: 2.19×10⁻² atoms/(barn-cm) = 8.04×10²⁸ total atoms
   - U-235: 2.50×10⁻⁴ atoms/(barn-cm) = 9.15×10²⁶ total atoms
   - Material is properly marked as `depletable=True`
   - Volume is correctly set: 3.67×10⁶ cm³

2. **CoupledOperator Setup**: The depletion operator initializes correctly:
   - Finds 1 burnable material
   - Loads depletion chain with U-235 and U-238 included
   - Initial composition in operator shows correct atom counts

3. **materials.xml**: The exported materials file contains the correct initial composition

### ❌ The Problem:
**The HDF5 depletion results file (`depletion_results.h5`) shows ZERO initial uranium densities:**
```python
U238 densities: [0.0, 1.79e-12, 2.90e-17, 5.18e-22]
U235 densities: [0.0, 1.58e-18, 1.47e-39, 4.97e-40]
```

This means the simulation **started with zero uranium** and somehow produced small amounts during the run (probably from other decay chains).

## Likely Causes

###1. **Depletion Chain Issues with `reduce_chain=True`**
   - While U-235/U-238 are in the reduced chain, there may be issues with reaction cross-sections
   - Fixed source mode + reduce_chain might not properly initialize fissionable material inventories

### 2. **Fixed Source Normalization Issue**
   - Using `normalization_mode="source-rate"` may not properly track material consumption
   - The neutron source rate is constant (1×10¹⁴ n/s), but this doesn't ensure proper depletion tracking

### 3. **First Transport Step Problem**
   - OpenMC may have failed to properly track material changes in the first transport calculation
   - The HDF5 file starts at t=0 with zero uranium, suggesting initialization failed

## Recommended Fixes

### Option 1: Try WITHOUT reduce_chain
```python
operator = openmc.deplete.CoupledOperator(
    model=model,
    chain_file=CHAIN_FILE,
    reduce_chain=False,  # Changed from True
    normalization_mode="source-rate"
)
```

### Option 2: Use Eigenvalue Mode Instead
Fixed source depletion is more experimental. Try eigenvalue mode for spent fuel:
```python
settings.run_mode = 'eigenvalue'
settings.inactive = 50
settings.batches = 100
# Remove source, let it find criticality
```

### Option 3: Manually Initialize the Depletion
After creating the operator, manually verify and set initial composition:
```python
operator = openmc.deplete.CoupledOperator(...)

# Check and print initial composition
print("Initial U-238:", operator.number.number[0, operator.number.nuclides.index('U238')])

# If zero, this confirms the bug and you may need to use OpenMC development version
```

### Option 4: Check OpenMC Version
This might be a bug in OpenMC 0.15.2. Consider:
- Updating to latest OpenMC version
- Reporting the issue to OpenMC developers
- Checking OpenMC GitHub issues for similar problems

## Current Status

The **results script is now correctly annotated** to show:
- Initial composition from material definition (correct values)
- Final composition from HDF5 file (invalid due to simulation bug)
- Clear warnings that the uranium depletion results are not physically meaningful

## Next Steps

1. Re-run simulation with `reduce_chain=False`
2. If that fails, try eigenvalue mode
3. If still failing, this is likely an OpenMC bug that needs to be reported

---
**Generated**: January 10, 2026
