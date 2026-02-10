import openmc
import math 

def build_trufuelsphere_nobox(sphere_inner_radius, sphere_outer_radius):
    """
    Builds a geometry with a HOLLOW TRU fuel sphere (shell) at the center,
    surrounded by heavy water (D2O) coolant.
    NO aluminum box in this iteration - simplified geometry.

    This is Iteration4: PURIFIED TRANSURANIC FUEL (no uranium) WITHOUT ALUMINUM BOX
    - Removes all uranium isotopes (U233, U234, U235, U236, U237, U238)
    - Contains only transuranics (Pu, Np, Am, Cm) + fission products + oxygen
    - Density renormalized after uranium removal
    - Simplified geometry: no aluminum box, just fuel shell with coolant

    - Inner radius = 50 cm (default)
    - Outer radius = 100 cm (default)
    - Units: centimeters
    - Returns: (geometry, materials)
    """
    
    # 1) Materials
    
    # Spent fuel with uranium
    mtru_fuel = openmc.Material(1, 'spent UO2')
    mtru_fuel.set_density('atom/b-cm', 7.133315757E-02)
    mtru_fuel.temperature = 300
    
    # Oxygen
    mtru_fuel.add_nuclide('O16', 6.7187968E-01)
    
    # URANIUM ISOTOPES
    mtru_fuel.add_nuclide('U238', 3.0769658E-01)
    mtru_fuel.add_nuclide('U235', 3.4979859E-03)
    mtru_fuel.add_nuclide('U236', 2.1091663E-03)
    mtru_fuel.add_nuclide('U234', 6.4111968E-05)
    mtru_fuel.add_nuclide('U237', 1.7067631E-07)
    
    # TRANSURANICS (Plutonium)
    mtru_fuel.add_nuclide('Pu239', 2.5718196E-03)
    mtru_fuel.add_nuclide('Pu240', 1.0215150E-03)
    mtru_fuel.add_nuclide('Pu241', 6.7015516E-04)
    mtru_fuel.add_nuclide('Pu242', 2.6644975E-04)
    mtru_fuel.add_nuclide('Pu238', 1.3101157E-04)
    
    # TRANSURANICS (Neptunium)
    mtru_fuel.add_nuclide('Np237', 2.7884345E-04)
    
    # TRANSURANICS (Americium)
    mtru_fuel.add_nuclide('Am243', 7.8638873E-05)
    mtru_fuel.add_nuclide('Am241', 3.1275801E-05)
    
    # TRANSURANICS (Curium)
    mtru_fuel.add_nuclide('Cm244', 3.2900525E-05)
    mtru_fuel.add_nuclide('Cm242', 7.3636478E-06)
    mtru_fuel.add_nuclide('Cm245', 1.9560548E-06)
    mtru_fuel.add_nuclide('Cm243', 3.0834446E-07)
    mtru_fuel.add_nuclide('Cm246', 1.9481932E-07)
    
    # FISSION PRODUCTS
    mtru_fuel.add_nuclide('Cs137', 1.0510726E-03)
    mtru_fuel.add_nuclide('Cs133', 9.7963234E-04)
    mtru_fuel.add_nuclide('Tc99', 9.3559531E-04)
    mtru_fuel.add_nuclide('Ru101', 9.1992123E-04)
    mtru_fuel.add_nuclide('Zr93', 8.9218967E-04)
    mtru_fuel.add_nuclide('Mo95', 8.5245707E-04)
    mtru_fuel.add_nuclide('Sr90', 6.7977794E-04)
    mtru_fuel.add_nuclide('Nd143', 6.5286858E-04)
    mtru_fuel.add_nuclide('Nd145', 5.3344636E-04)
    mtru_fuel.add_nuclide('Rh103', 4.8652147E-04)
    mtru_fuel.add_nuclide('Cs135', 4.8545594E-04)
    mtru_fuel.add_nuclide('Pd107', 2.5941176E-04)
    mtru_fuel.add_nuclide('Sm150', 2.2788081E-04)
    mtru_fuel.add_nuclide('I129', 1.4371885E-04)
    mtru_fuel.add_nuclide('Pm147', 1.0269899E-04)
    mtru_fuel.add_nuclide('Eu153', 9.2879588E-05)
    mtru_fuel.add_nuclide('Sm152', 8.8634337E-05)
    mtru_fuel.add_nuclide('Ag109', 8.5141959E-05)
    mtru_fuel.add_nuclide('Sm147', 6.6014495E-05)
    mtru_fuel.add_nuclide('Ru103', 1.9360418E-05)
    mtru_fuel.add_nuclide('Sn126', 1.8018478E-05)
    mtru_fuel.add_nuclide('Nb95', 1.4824839E-05)
    mtru_fuel.add_nuclide('Cl36', 1.4019981E-05)
    mtru_fuel.add_nuclide('Ca41', 1.4019976E-05)
    mtru_fuel.add_nuclide('Ni59', 1.4019973E-05)
    mtru_fuel.add_nuclide('Sm151', 1.3614245E-05)
    mtru_fuel.add_nuclide('Se79', 7.0915868E-06)
    mtru_fuel.add_nuclide('Eu155', 4.7729475E-06)
    mtru_fuel.add_nuclide('Pr143', 2.0561139E-06)
    mtru_fuel.add_nuclide('Sm149', 1.9355989E-06)
    mtru_fuel.add_nuclide('Nd147', 5.1292485E-07)
    mtru_fuel.add_nuclide('Xe133', 9.3481328E-08)
    mtru_fuel.add_nuclide('Gd155', 7.6579955E-08)
    mtru_fuel.add_nuclide('I133', 7.5534463E-08)
    mtru_fuel.add_nuclide('Eu152', 4.5260253E-08)
    mtru_fuel.add_nuclide('Mo99', 1.1472054E-09 / 0.68652)

    # Heavy water coolant (D2O)
    heavy_water = openmc.Material(name='D2O')
    heavy_water.add_nuclide('H2', 2.0)  # Deuterium
    heavy_water.add_nuclide('O16', 1.0)
    heavy_water.set_density('g/cm3', 1.11)  # Heavy water density
    heavy_water.add_s_alpha_beta('c_D_in_D2O')  # Thermal scattering for heavy water

    materials = openmc.Materials([mtru_fuel, heavy_water])

    # 2) Geometry - Simplified: just inner coolant region and fuel shell
    inner_surface_sphere = openmc.Sphere(r=sphere_inner_radius)
    outer_surface_sphere = openmc.Sphere(r=sphere_outer_radius, boundary_type='vacuum')

    # Cells
    inner_coolant_cell = openmc.Cell(
        name='inner_coolant',
        region=-inner_surface_sphere
    )
    inner_coolant_cell.fill = heavy_water
    
    tru_fuel_sphere_cell = openmc.Cell(
        name='TRU_fuel_shell'
    )
    tru_fuel_sphere_cell.region = +inner_surface_sphere & -outer_surface_sphere
    tru_fuel_sphere_cell.fill = mtru_fuel
    
    # 3) VOL: compute spherical shell volume and set it on the MATERIAL (for depletion)
    shell_volume = (4.0/3.0) * math.pi * (sphere_outer_radius**3 - sphere_inner_radius**3)
    mtru_fuel.volume = shell_volume  # Critical for depletion
    tru_fuel_sphere_cell.volume = shell_volume  # Optional; nice for bookkeeping

    root_universe = openmc.Universe(cells=[
        inner_coolant_cell,
        tru_fuel_sphere_cell
    ])
    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()
    materials.export_to_xml()

    return geometry, materials
