import openmc
import math

def create_cylindrical_source(height, radius, e_min, e_max):
    """
    Creates and returns an OpenMC Source object for a vertical,
    volumetrically uniform cylindrical source.
    
    The source is centered at (0,0,0) and oriented along the z-axis.
    """
    
    # For uniform-in-volume in a cylinder, p(r) ~ r^1
    r_dist = openmc.stats.PowerLaw(0.0, radius, 1.0)
    
    # Angle in the x-y plane
    phi_dist = openmc.stats.Uniform(0.0, 2.0 * math.pi)
    
    # Uniform distribution along the cylinder's height (z-axis)
    z_dist = openmc.stats.Uniform(-height / 2.0, height / 2.0)
    
    # Create the cylindrical spatial distribution
    space = openmc.stats.CylindricalIndependent(r_dist, phi_dist, z_dist, origin=(0.0, 0.0, 0.0))

    # Isotropic angular distribution
    angle = openmc.stats.Isotropic()
    
    # Uniform energy distribution
    energy = openmc.stats.Uniform(e_min, e_max)
    
    # Combine into a single source object
    source = openmc.Source(space=space, angle=angle, energy=energy)
    
    return source
