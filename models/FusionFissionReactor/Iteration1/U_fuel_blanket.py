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
    fuel.add_element('U', 1, enrichment=0.7)  # natural uranium
    fuel.add_element('O', 2)
    fuel.add_element('U', 1)
    fuel.set_density('g/cm3', 10.0)

    #aluminum metal
    aluminum = openmc.Material()
    aluminum.add_element('Al',1)
    aluminum.set_density('g/cm3', 2.7)
    #aluminum.temperature = T to be figured out 
    
    #Na coolant 
    Na = openmc.Material()
    Na.add_nuclide('Na23',1)
    Na.set_density('g/cm3', 0.856)

    materials = openmc.Materials([aluminum, mspentfuel, Na])

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

   
    root_universe = openmc.Universe(cells=[inner_box_vacuum_cell,al_box_cell, outer_box_coolant_cell, fuel_sphere_cell])
    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()
    materials.export_to_xml()

    return geometry, materials
