# Solutions for Fusion-Fission Hybrid Depletion Simulation

## The Challenge

We need to simulate a **fusion-fission hybrid** system:
- **External source**: D-T fusion neutrons at 14.1 MeV (10^14 n/s point source)
- **Multiplying blanket**: Spent UO2 fuel (subcritical, k_eff ≈ 0.515)
- **Accurate depletion**: Track isotope evolution over time

## Problem: Current Approaches Don't Work

### Approach 1: Fixed-Source Depletion (Iterations 2-4)
**Setup:**
- `settings.run_mode = 'fixed source'`
- External D-T source: 10^14 n/s
- `normalization_mode = "source-rate"`

**Issue:** OpenMC's normalization bug
- OpenMC tallies ALL neutrons (source + fission)
- Then multiplies by source_rate normalization
- **Result**: Fission neutrons counted 5×10^15 times → impossible depletion rates
- Pu-239: +541,000% (physically impossible)
- U-238: -100% in 180 days (requires 5×10^7 captures per source neutron)

### Approach 2: Eigenvalue Depletion (Iterations 4-5)
**Setup:**
- `settings.run_mode = 'eigenvalue'`
- No external source
- `normalization_mode = "fission-q"`

**Issue:** No fusion driver
- System driven by fission only (U-235, Pu-239)
- **No 14.1 MeV D-T neutrons**
- Not modeling the fusion-fission hybrid
- Results are accurate for a standalone subcritical assembly, but not our problem

## Solution Options

### Option 1: Manual Depletion Calculation (RECOMMENDED)

Run fixed-source transport WITHOUT depletion, extract reaction rates manually, calculate depletion ourselves.

**Method:**
1. Run OpenMC in fixed-source mode (no depletion)
2. Tally reaction rates (fission, capture, (n,2n), etc.)
3. Extract actual neutron counts (not normalized)
4. Calculate depletion rates externally:
   ```python
   # True reaction rate per second
   reaction_rate = (tally_mean × source_rate × (1 + k_source)) / num_particles
   
   # Where k_source accounts for multiplication
   k_source = (fission_neutrons_produced / source_neutrons)
   ```
5. Use `openmc.deplete` API manually to update compositions

**Pros:**
- Full control over normalization
- Can properly account for fusion source + fission multiplication
- Physically accurate

**Cons:**
- Requires custom depletion script
- More complex than built-in operators
- Need to implement predictor-corrector or similar

**Implementation Strategy:**
```python
# 1. Set up fixed-source transport (no depletion)
settings.run_mode = 'fixed source'
settings.source = dt_point_source  # 14.1 MeV, 10^14 n/s

# 2. Add reaction rate tallies
tallies = [
    ('(n,fission)', material_filter),
    ('(n,gamma)', material_filter),
    ('(n,2n)', material_filter),
    ('(n,alpha)', material_filter),
    # ... all relevant reactions
]

# 3. Run transport only
model.run()

# 4. Extract reaction rates
with openmc.StatePoint('statepoint.h5') as sp:
    fission_rate = sp.get_tally(name='fission').mean
    capture_rate = sp.get_tally(name='capture').mean
    
# 5. Calculate actual rates (accounting for multiplication)
actual_fission_rate = fission_rate * source_rate / num_particles

# 6. Use OpenMC depletion chain to calculate isotope changes
from openmc.deplete import Chain, Integrator
# ... manual integration
```

### Option 2: Fixed-Source with Post-Correction

Use existing fixed-source depletion, but apply correction factor based on k_source.

**Method:**
1. Run fixed-source depletion (as in Iteration 4)
2. Calculate effective multiplication: k_eff/(1-k_eff) ≈ 1.06 for k_eff=0.515
3. Divide all depletion rates by (1 + multiplication_factor)

**Correction:**
```python
# For k_eff = 0.515
multiplication = 0.515 / (1 - 0.515) = 1.06

# True depletion rate
true_rate = openmc_rate / (1 + multiplication) = openmc_rate / 2.06
```

**Pros:**
- Uses existing simulation infrastructure
- Relatively simple post-processing

**Cons:**
- Approximate correction (assumes constant k_eff)
- Still large uncertainties
- Not validated for accuracy

### Option 3: Separate Transport + Depletion

Run two simulations: one for neutronics, one for depletion.

**Method:**
1. **Step 1 (Transport)**: Fixed-source with fusion neutrons
   - Get spatial flux distribution
   - Calculate one-group cross sections

2. **Step 2 (Depletion)**: Use flux from Step 1
   - Calculate reaction rates from φ(r) × Σ
   - Apply to depletion equations
   - Update material compositions

3. **Step 3 (Iterate)**: Repeat with new compositions

**Pros:**
- Decouples neutronics from depletion normalization
- Can use different approximations for each

**Cons:**
- Complex workflow
- Requires manual coupling
- Slower (multiple runs per step)

### Option 4: OpenMC IndependentOperator (If Available)

Check if OpenMC has an `IndependentOperator` for fixed-source without the normalization issue.

**Investigation needed:**
```python
import openmc.deplete
# Check if exists:
hasattr(openmc.deplete, 'IndependentOperator')
```

**If it exists:**
```python
operator = openmc.deplete.IndependentOperator(
    model=model,
    chain_file=CHAIN_FILE,
    normalization_mode=???  # TBD
)
```

**Pros:**
- Built-in solution
- Likely properly handles fixed-source

**Cons:**
- May not exist in OpenMC 0.15.3
- Documentation unclear

### Option 5: Power-Driven Fixed-Source

Instead of source-rate normalization, normalize to power from fission.

**Method:**
1. Run short fixed-source simulation
2. Calculate actual fission power from tallies
3. Use that power for depletion normalization

**Implementation:**
```python
# Short run to get fission power
model.run()

# Extract fission rate and calculate power
fission_energy = 200e6 * 1.60218e-19  # J per fission
fission_rate = get_fission_rate_from_tally()  # fissions/s
power_watts = fission_rate * fission_energy

# Then run depletion with power normalization
operator = openmc.deplete.CoupledOperator(
    model=model,
    chain_file=CHAIN_FILE,
    normalization_mode="fission-q"  # Use fission power
)

integrator = openmc.deplete.PredictorIntegrator(
    operator=operator,
    timesteps=[time_step],
    power=[power_watts],  # Use calculated power
    timestep_units='s'
)
```

**Pros:**
- Uses existing infrastructure
- Normalizes to actual fission power

**Cons:**
- Still doesn't properly account for external source contribution
- May still have issues with multiplication counting

## Recommendation: Proceed with Option 1

**Manual depletion calculation** is the most reliable approach for fusion-fission hybrids.

**Next Steps:**
1. Create `fixed_source_transport.py` - Transport only, no depletion
2. Create `manual_depletion.py` - Extract rates and calculate depletion
3. Create `coupled_simulation.py` - Iterate transport → depletion
4. Validate with neutron balance checks

**Estimated Work:**
- Setup: 2-3 hours
- Testing: 1-2 hours
- Validation: 1-2 hours
- **Total: ~1 day of development**

Would you like me to start implementing Option 1?
