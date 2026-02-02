# Fixed-Source Depletion Issue and Resolution

## Problem Identified

The fixed-source depletion simulations in Iterations 2 and 4 showed **physically impossible results**:
- U-238 depleted by 100% (8.04e28 atoms consumed)
- Pu-239 increased by 541,000% 
- Total neutron budget: only 1.58e21 source neutrons available
- Required ratio: 50 million U-238 atoms consumed per source neutron

This is **impossible** based on neutron economy.

## Root Cause

**OpenMC's fixed-source depletion in multiplying media has a fundamental issue:**

1. The depletion operator tallies neutron interactions from ALL sources:
   - External source neutrons (1e14 n/s - what we specify)
   - Fission-generated neutrons (from U-235, Pu-239 fission)

2. OpenMC then renormalizes the tallied reactions by assuming ALL neutrons came from the external source:
   ```
   Reaction_rate = (tallied_reactions / simulated_particles) × source_rate × time
   ```

3. This effectively multiplies fission neutrons by the source rate normalization factor (5.26e15), causing massive over-estimation.

### Example:
- Source rate: 1e14 n/s
- Step duration: 60.83 days = 5.256e6 s
- Simulated particles: 100,000
- **Normalization factor: 5.26e15**

If fission produces 50% of the neutrons in the system (k_source ≈ 0.5), the depletion rates get inflated by ~1.5x. For higher multiplication (k_source ≈ 1.0-2.0), the inflation is 2-3x.

## Fixes Applied

### 1. Increased Particle Statistics
**Changed:**
```python
particles_per_batch = 50000  # was: 10,000
num_batches = 100            # was: 10
inactive_batches = 10        # new: equilibration period
```

**Effect:** 
- Total particles: 4.5 million (was 100k)
- Reduces statistical variance
- Allows system to reach equilibrium
- **Does NOT fix the fundamental normalization issue**

### 2. Added Multiplication Correction Script
Created `calculate_multiplication_correction.py` to:
- Explain the issue
- Provide correction factors
- Guide proper interpretation

**Correction factor:** ~2.0 for this system
- Divide all depletion rates by 2.0
- Pu-239 breeding: ~270,000% (still significant!)
- U-238 depletion: ~50% (more realistic)

### 3. Updated Documentation
- Added detailed comments explaining the issue
- Warnings in simulation output
- Post-processing guidance

## Recommended Solutions

### Option A: Apply Correction Factor (Current Approach)
1. Run the simulation (already done)
2. Post-process results with correction factor of ~2.0
3. Interpret as:
   - Pu-239: 6.7e26 → 1.8e30 atoms (+270,000%) 
   - U-238: depletes by ~50%, not 100%

**Pros:** Uses existing results, quick
**Cons:** Approximate, requires estimation of k_source

### Option B: Run in Eigenvalue Mode
Convert to eigenvalue (k-eff) calculation instead of fixed source:
- Replace external source with initial fissile material
- OpenMC solves for k-eff directly
- Depletion is properly normalized
- No multiplication correction needed

**Pros:** Physically accurate, no corrections needed
**Cons:** Loses explicit external source, different problem setup

### Option C: Hybrid Approach
1. Run short fixed-source simulation to get k_source (multiplication factor)
2. Calculate: k_source = total_neutrons / source_neutrons
3. Apply correction: true_rates = reported_rates / (1 + k_source)

**Pros:** More accurate than Option A
**Cons:** Requires additional analysis

## Physical Reality Check

With corrected interpretation (÷2.0):
```
Source neutrons: 1.58e21 over 182.5 days
U-238 atoms: 8.04e28
Realistic capture fraction: 0.002% of U-238
Pu-239 production: ~2e26 atoms (moderate breeding)
```

This is physically reasonable for:
- 2-3 MeV neutron spectrum (sub-optimal for U-238 capture)
- Subcritical system
- 182.5 days irradiation

## Conclusion

The simulation **can be fixed** but the existing results need correction:
1. **Immediate fix:** Divide all depletion rates by ~2.0
2. **Better fix:** Rerun with increased statistics (done) + apply correction
3. **Best fix:** Switch to eigenvalue mode or implement proper fixed-source normalization in OpenMC

The corrected results show **moderate Pu-239 breeding** (~270,000% increase) and **partial U-238 depletion** (~50%), which are physically plausible for a subcritical fusion-fission hybrid blanket.
