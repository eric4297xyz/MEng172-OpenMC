import openmc
import argparse
import glob
import os

# --- 1. Define materials ---
fuel = openmc.Material(name="UO2 fuel")
fuel.add_element('U', 1, enrichment=0.7)  # natural uranium
fuel.add_element('O', 2)
fuel.add_element('U', 1)
fuel.set_density('g/cm3', 10.0)

water = openmc.Material(name="Water")
water.add_element('H', 2)
water.add_element('O', 1)
water.set_density('g/cm3', 1.0)
#water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([fuel, water])

# --- 2. Geometry ---
# --- 2. Geometry (bounded with an outer vacuum sphere) ---
fuel_region = -openmc.Sphere(r=100, boundary_type='transmission')  # inner sphere (finite surface)
fuel_cell   = openmc.Cell(region=fuel_region, fill=fuel)

# Outer vacuum boundary that closes the model
outer_sphere = openmc.Sphere(r=200.0, boundary_type='vacuum')
mod_region   = +openmc.Sphere(r=100) & -outer_sphere
mod_cell     = openmc.Cell(region=mod_region, fill=water)

# --- 3. Settings ---
settings = openmc.Settings()
# Use temperature interpolation (don’t require windowed-multipole data)
settings.temperature = {'method': 'interpolation'}
settings.batches = 20
settings.inactive = 5
settings.particles = 1000

# --- Source definition ---
# Default: isotropic, monoenergetic point source at the origin (2 MeV)

def build_source(energy_mev: float = 2.0, coords=(0.0, 0.0, 0.0), hemisphere: bool = False, reference_uvw=(0.0, 0.0, 1.0)) -> openmc.Source:
	"""Create an OpenMC point source.

	Parameters
	- energy_mev: float energy in MeV for a monoenergetic source (converted to eV)
	- coords: (x, y, z) point coordinates in cm
	- hemisphere: if True, restrict directions to the hemisphere about reference_uvw
	- reference_uvw: (ux,uy,uz) direction which defines the hemisphere forward direction
	"""
	energy_ev = energy_mev * 1.0e6
	space = openmc.stats.Point(coords)
	# Choose angular distribution: isotropic or hemisphere using PolarAzimuthal
	if hemisphere:
		# mu = cos(theta) in [0,1] corresponds to polar angle theta in [0, pi/2]
		mu = openmc.stats.Uniform(0.0, 1.0)
		phi = openmc.stats.Uniform(0.0, 2.0 * openmc.pi)
		angle = openmc.stats.PolarAzimuthal(mu=mu, phi=phi, reference_uvw=reference_uvw)
	else:
		angle = openmc.stats.Isotropic()
	# Discrete energy distribution for a monoenergetic source
	energy = openmc.stats.Discrete([energy_ev], [1.0])
	return openmc.Source(space=space, angle=angle, energy=energy)

# Parse optional command-line args so you can run different point sources easily
parser = argparse.ArgumentParser(description='Run OpenMC test with a configurable point source')
parser.add_argument('--energy-mev', type=float, default=2.45, help='Monoenergetic source energy in MeV (default: 2.45 MeV)')
parser.add_argument('--x', type=float, default=0.0, help='X coordinate of the point source in cm')
parser.add_argument('--y', type=float, default=0.0, help='Y coordinate of the point source in cm')
parser.add_argument('--z', type=float, default=0.0, help='Z coordinate of the point source in cm')
parser.add_argument('--hemisphere', action='store_false', help='Limit angular distribution to a hemisphere (forward relative to reference_uvw)')
parser.add_argument('--refx', type=float, default=0.0, help='Reference direction x-component for hemisphere (default: 0)')
parser.add_argument('--refy', type=float, default=0.0, help='Reference direction y-component for hemisphere (default: 0)')
parser.add_argument('--refz', type=float, default=1.0, help='Reference direction z-component for hemisphere (default: 1)')
parser.add_argument('--source-rate', type=float, default=1e14, help='Physical source emission rate in neutrons/sec (default: 1e12)')
parser.add_argument('--energy-per-fission-mev', type=float, default=200.0, help='Energy released per fission in MeV (default: 200 MeV)')
args, _ = parser.parse_known_args()

settings.source = build_source(args.energy_mev, (args.x, args.y, args.z), hemisphere=args.hemisphere, reference_uvw=(args.refx, args.refy, args.refz))

# Print mapping from simulated histories to physical neutrons/sec
n_sim = settings.particles * (settings.batches - settings.inactive)
S = args.source_rate
print(f'Total simulated source histories: {n_sim}')
print(f'Physical source rate S = {S:.3e} n/s')
print(f'Each simulated particle represents S / N_sim = {S / n_sim:.6e} n/s')

# --- 4. Export and run ---
# build the geometry and export all XML inputs
geometry = openmc.Geometry([fuel_cell, mod_cell])
materials.export_to_xml()
geometry.export_to_xml()
settings.export_to_xml()

# --- Tallies: add a fission tally for the fuel cell so we can compute power ---
tally = openmc.Tally(name='fission_rate')
tally.filters = [openmc.CellFilter(fuel_cell)]
tally.scores = ['fission']
openmc.Tallies([tally]).export_to_xml()

openmc.run()

# After the run, read the latest statepoint and extract the fission tally
sp_files = glob.glob('statepoint.*.h5')
if sp_files:
	sp = max(sp_files, key=os.path.getmtime)
	print(f'Reading statepoint: {sp}')
	s = openmc.StatePoint(sp)
	try:
		t = s.get_tally(name='fission_rate')
		# mean is per source particle simulated
		mean = t.mean.flatten()[0]
		n_sim = settings.particles * (settings.batches - settings.inactive)
		S = args.source_rate
		fissions_per_s = mean * S
		# energy per fission (MeV -> J)
		e_fiss_J = args.energy_per_fission_mev * 1.0e6 * 1.602176634e-19
		power_W = fissions_per_s * e_fiss_J
		print(f'Fission rate (per particle): {mean:.6e}')
		print(f'Fissions/sec (for S={S:.3e} n/s): {fissions_per_s:.6e} 1/s')
		print(f'Estimated power produced: {power_W:.6e} W = {power_W/1e6:.6f} MW')
	except Exception as e:
		print('Could not read fission tally from statepoint:', e)
else:
	print('No statepoint file found; cannot compute power.')

