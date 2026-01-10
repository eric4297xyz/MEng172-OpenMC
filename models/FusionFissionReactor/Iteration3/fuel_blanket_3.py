import openmc
import math 

def build_trufuelsphere_albox(box_side, box_width, sphere_inner_radius, sphere_outer_radius):
    """
    Builds a geometry with a HOLLOW TRU fuel sphere (shell) at the center,
    surrounded by sodium coolant.
    An aluminum box is placed inside the sphere.

    This is Iteration3: PURIFIED TRANSURANIC FUEL (no uranium)
    - Removes all uranium isotopes (U233, U234, U235, U236, U237, U238)
    - Contains only transuranics (Pu, Np, Am, Cm) + fission products + oxygen
    - Density renormalized after uranium removal

    - Inner radius = 50 cm (default)
    - Outer radius = 100 cm (default)
    - Units: centimeters
    - Returns: (geometry, materials)
    """

    if box_side >= sphere_inner_radius:
        print("box_side >= sphere_inner_radius, box will not fit")
        return
    
    # 1) Materials
    
    # TRU Fuel (purified transuranic fuel - uranium removed)
    # Based on spent UO2 composition but with all U isotopes removed
    # Original density: 7.133315757E-02 atom/b-cm
    # Uranium fraction removed: U238 (0.30769658) + U235 (0.0034979859) + U236 (0.0021091663) + U234 (0.000064111968) + U237 (1.7067631E-07) + U233 (1.2208878E-09)
    # Total U removed: ~0.31348 (31.35% of atoms)
    # New density normalized: 7.133315757E-02 * (1 - 0.31348) = 4.896E-02 atom/b-cm
    
    mtru_fuel = openmc.Material(1, 'TRU fuel (no uranium)')
    mtru_fuel.set_density('atom/b-cm', 4.896E-02)  # Renormalized after U removal
    mtru_fuel.temperature = 300
    
    # Oxygen (renormalized)
    mtru_fuel.add_nuclide('O16', 6.7187968E-01 / 0.68652)  # Scale up to compensate for U removal
    
    # TRANSURANICS (Plutonium)
    mtru_fuel.add_nuclide('Pu239', 2.5718196E-03 / 0.68652)
    mtru_fuel.add_nuclide('Pu240', 1.0215150E-03 / 0.68652)
    mtru_fuel.add_nuclide('Pu241', 6.7015516E-04 / 0.68652)
    mtru_fuel.add_nuclide('Pu242', 2.6644975E-04 / 0.68652)
    mtru_fuel.add_nuclide('Pu238', 1.3101157E-04 / 0.68652)
    mtru_fuel.add_nuclide('Pu244', 9.4590026E-09 / 0.68652)
    
    # TRANSURANICS (Neptunium)
    mtru_fuel.add_nuclide('Np237', 2.7884345E-04 / 0.68652)
    mtru_fuel.add_nuclide('Np239', 3.2856837E-09 / 0.68652)
    
    # TRANSURANICS (Americium)
    mtru_fuel.add_nuclide('Am243', 7.8638873E-05 / 0.68652)
    mtru_fuel.add_nuclide('Am241', 3.1275801E-05 / 0.68652)
    
    # TRANSURANICS (Curium)
    mtru_fuel.add_nuclide('Cm244', 3.2900525E-05 / 0.68652)
    mtru_fuel.add_nuclide('Cm242', 7.3636478E-06 / 0.68652)
    mtru_fuel.add_nuclide('Cm245', 1.9560548E-06 / 0.68652)
    mtru_fuel.add_nuclide('Cm243', 3.0834446E-07 / 0.68652)
    mtru_fuel.add_nuclide('Cm246', 1.9481932E-07 / 0.68652)
    
    # FISSION PRODUCTS (scaled to compensate for U removal)
    mtru_fuel.add_nuclide('Cs137', 1.0510726E-03 / 0.68652)
    mtru_fuel.add_nuclide('Cs133', 9.7963234E-04 / 0.68652)
    mtru_fuel.add_nuclide('Tc99', 9.3559531E-04 / 0.68652)
    mtru_fuel.add_nuclide('Ru101', 9.1992123E-04 / 0.68652)
    mtru_fuel.add_nuclide('Zr93', 8.9218967E-04 / 0.68652)
    mtru_fuel.add_nuclide('Mo95', 8.5245707E-04 / 0.68652)
    mtru_fuel.add_nuclide('Sr90', 6.7977794E-04 / 0.68652)
    mtru_fuel.add_nuclide('Nd143', 6.5286858E-04 / 0.68652)
    mtru_fuel.add_nuclide('Nd145', 5.3344636E-04 / 0.68652)
    mtru_fuel.add_nuclide('Rh103', 4.8652147E-04 / 0.68652)
    mtru_fuel.add_nuclide('Cs135', 4.8545594E-04 / 0.68652)
    mtru_fuel.add_nuclide('Pd107', 2.5941176E-04 / 0.68652)
    mtru_fuel.add_nuclide('Sm150', 2.2788081E-04 / 0.68652)
    mtru_fuel.add_nuclide('I129', 1.4371885E-04 / 0.68652)
    mtru_fuel.add_nuclide('Pm147', 1.0269899E-04 / 0.68652)
    mtru_fuel.add_nuclide('Eu153', 9.2879588E-05 / 0.68652)
    mtru_fuel.add_nuclide('Sm152', 8.8634337E-05 / 0.68652)
    mtru_fuel.add_nuclide('Ag109', 8.5141959E-05 / 0.68652)
    mtru_fuel.add_nuclide('Sm147', 6.6014495E-05 / 0.68652)
    mtru_fuel.add_nuclide('Ru103', 1.9360418E-05 / 0.68652)
    mtru_fuel.add_nuclide('Sn126', 1.8018478E-05 / 0.68652)
    mtru_fuel.add_nuclide('Nb95', 1.4824839E-05 / 0.68652)
    mtru_fuel.add_nuclide('Cl36', 1.4019981E-05 / 0.68652)
    mtru_fuel.add_nuclide('Ca41', 1.4019976E-05 / 0.68652)
    mtru_fuel.add_nuclide('Ni59', 1.4019973E-05 / 0.68652)
    mtru_fuel.add_nuclide('Sm151', 1.3614245E-05 / 0.68652)
    mtru_fuel.add_nuclide('Se79', 7.0915868E-06 / 0.68652)
    mtru_fuel.add_nuclide('Eu155', 4.7729475E-06 / 0.68652)
    mtru_fuel.add_nuclide('Pr143', 2.0561139E-06 / 0.68652)
    mtru_fuel.add_nuclide('Sm149', 1.9355989E-06 / 0.68652)
    mtru_fuel.add_nuclide('Nd147', 5.1292485E-07 / 0.68652)
    mtru_fuel.add_nuclide('Xe133', 9.3481328E-08 / 0.68652)
    mtru_fuel.add_nuclide('Gd155', 7.6579955E-08 / 0.68652)
    mtru_fuel.add_nuclide('I133', 7.5534463E-08 / 0.68652)
    mtru_fuel.add_nuclide('Eu152', 4.5260253E-08 / 0.68652)
    mtru_fuel.add_nuclide('Mo99', 1.1472054E-09 / 0.68652)

    # Aluminum metal
    aluminum = openmc.Material()
    aluminum.add_element('Al', 1)
    aluminum.set_density('g/cm3', 2.7)
    
    # Na coolant 
    Na = openmc.Material()
    Na.add_nuclide('Na23', 1)
    Na.set_density('g/cm3', 0.856)

    materials = openmc.Materials([aluminum, mtru_fuel, Na])

    # 2) Geometry
    Al_box_inner_wall = openmc.model.RectangularParallelepiped(
        -box_side/2, box_side/2, 
        -box_side/2, box_side/2, 
        -box_side/2, box_side/2, 
        boundary_type='transmission'
    )
    Al_box_outer_wall = openmc.model.RectangularParallelepiped(
        -(box_side+box_width)/2, (box_side+box_width)/2, 
        -(box_side+box_width)/2, (box_side+box_width)/2, 
        -(box_side+box_width)/2, (box_side+box_width)/2, 
        boundary_type='transmission'
    )
    inner_surface_sphere = openmc.Sphere(r=sphere_inner_radius)
    outer_surface_sphere = openmc.Sphere(r=sphere_outer_radius, boundary_type='vacuum')

    # Cells
    inner_box_vacuum_cell = openmc.Cell(
        name='inner_Al_vacuum',
        region=-Al_box_inner_wall
    )
    
    al_box_cell = openmc.Cell()
    al_box_cell.region = +Al_box_inner_wall & -Al_box_outer_wall
    al_box_cell.fill = aluminum
    
    outer_box_coolant_cell = openmc.Cell()
    outer_box_coolant_cell.region = +Al_box_outer_wall & -inner_surface_sphere
    outer_box_coolant_cell.fill = Na
    
    tru_fuel_sphere_cell = openmc.Cell()
    tru_fuel_sphere_cell.region = +inner_surface_sphere & -outer_surface_sphere
    tru_fuel_sphere_cell.fill = mtru_fuel
    
    # 3) VOL: compute spherical shell volume and set it on the MATERIAL (for depletion)
    shell_volume = (4.0/3.0) * math.pi * (sphere_outer_radius**3 - sphere_inner_radius**3)
    mtru_fuel.volume = shell_volume  # Critical for depletion
    tru_fuel_sphere_cell.volume = shell_volume  # Optional; nice for bookkeeping

    root_universe = openmc.Universe(cells=[
        inner_box_vacuum_cell,
        al_box_cell, 
        outer_box_coolant_cell, 
        tru_fuel_sphere_cell
    ])
    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()
    materials.export_to_xml()

    return geometry, materials
