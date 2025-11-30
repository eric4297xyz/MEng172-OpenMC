import openmc
import math 

def build_fuelsphere_albox(box_side, box_width,  sphere_inner_radius, sphere_outer_radius,):
    
    
    """
    Builds a geometry with a HOLLOW natural uranium fuel (shell) at the center,
    surrounded by reflector.
    An aluminum box is placed inside the spehere and around the zpinch 

    - Inner radius = 0.5 * outer radius
    - Units: centimeters
    - Returns: (geometry, materials)
    """

    if box_side >= sphere_inner_radius:
        print("box_side >= sphere_inner_radius, box will not fit")
        return
    
    # 1) Materials
    
    #spent_fuel
    # --- 1. Define materials ---
    fuel = openmc.Material(name="UO2 fuel")
    # Add nuclides directly to avoid element expansion (requires cross_sections.xml)
    # Add uranium isotopes explicitly. `Material.add_nuclide` does not accept
    # an `enrichment` keyword, so add U235 and U238 with a simple enrichment
    # fraction. The variable below is the U-235 atom fraction (0.007 = 0.7%).
    enrichment_u235 = 0.007
    total_u = 1.0
    u235_frac = enrichment_u235
    u238_frac = max(0.0, 1.0 - u235_frac)
    fuel.add_nuclide('U235', total_u * u235_frac)
    fuel.add_nuclide('U238', total_u * u238_frac)
    fuel.add_nuclide('O16', 2.0)
    fuel.set_density('g/cm3', 10.0)

    #aluminum metal
    aluminum = openmc.Material(name='Al')
    # Add stable Al-27 nuclide explicitly to avoid element expansion
    aluminum.add_nuclide('Al27', 1.0)
    aluminum.set_density('g/cm3', 2.7)
    #aluminum.temperature = T to be figured out 
    
    #Na coolant 
    Na = openmc.Material(name='Na')
    Na.add_nuclide('Na23', 1.0)
    Na.set_density('g/cm3', 0.856)

    materials = openmc.Materials([aluminum, fuel, Na])

    # 2) Geometry
    Al_box_inner_wall = openmc.model.RectangularParallelepiped(-box_side/2, box_side/2, -box_side/2, box_side/2, -box_side/2, box_side/2, boundary_type='transmission')
    Al_box_outer_wall = openmc.model.RectangularParallelepiped(-(box_side+box_width)/2, (box_side+box_width)/2, -(box_side+box_width)/2, (box_side+box_width)/2, -(box_side+box_width)/2, (box_side+box_width)/2, boundary_type='transmission')
    inner_surface_sphere = openmc.Sphere(r=sphere_inner_radius)
    outer_surface_spheree = openmc.Sphere(r=sphere_outer_radius, boundary_type = 'reflective')

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
    
    fuel_sphere_cell = openmc.Cell()
    fuel_sphere_cell.region= +inner_surface_sphere & -outer_surface_spheree
    fuel_sphere_cell.fill = fuel
    
    # 3) VOL: compute spherical shell volume and set it on the MATERIAL (for depletion)
    shell_volume = (4.0/3.0) * math.pi * (sphere_outer_radius**3 - sphere_inner_radius**3)
    fuel.volume = shell_volume                      # <-- critical for depletion
    fuel_sphere_cell.volume = shell_volume                # optional; nice for bookkeeping

   
    # Do not assign explicit material IDs here; let OpenMC assign IDs when
    # the model is bundled/exported. The calling script should detect the
    # correct material IDs for depletion (e.g., by material name).

    root_universe = openmc.Universe(cells=[inner_box_vacuum_cell,al_box_cell, outer_box_coolant_cell, fuel_sphere_cell])
    geometry = openmc.Geometry(root_universe)
    # Do not export XML here; caller should export once to avoid duplicate
    # OpenMC Material/Geometry registrations (which can cause duplicate ID warnings).

    return geometry, materials
