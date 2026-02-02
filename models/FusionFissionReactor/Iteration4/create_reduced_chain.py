"""
Create a reduced depletion chain for TRU fuel simulation
This drastically reduces memory usage by only tracking essential isotopes
"""

import openmc.deplete

# Load the full chain
full_chain = openmc.deplete.Chain.from_xml("chain_endfb80_pwr.xml")

# Define isotopes to keep - focus on TRU, key fission products, and absorbers
isotopes_to_keep = [
    # Transuranics (TRU) - the fuel
    'Np237', 'Np238', 'Np239',
    'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'Pu243',
    'Am241', 'Am242', 'Am242_m1', 'Am243',
    'Cm242', 'Cm243', 'Cm244', 'Cm245', 'Cm246',
    
    # Key fission products (neutron poisons)
    'Xe135', 'Sm149', 'Sm151', 'Gd155', 'Gd157',
    'Cs135', 'Cs137', 'Sr90', 'I129', 'Tc99',
    
    # Other important isotopes
    'Zr93', 'Mo95', 'Ru101', 'Pd107', 'Ag109',
    'Nd143', 'Nd145', 'Eu155',
    
    # Oxygen (in fuel)
    'O16', 'O17', 'O18'
]

print(f"Full chain has {len(full_chain.nuclides)} nuclides")
print(f"Reducing to {len(isotopes_to_keep)} essential nuclides")

# Reduce the chain
reduced_chain = full_chain.reduce(isotopes_to_keep)

# Save the reduced chain
reduced_chain.export_to_xml("chain_reduced_tru.xml")

print(f"Reduced chain saved to: chain_reduced_tru.xml")
print(f"Memory savings: ~{100*(1-len(isotopes_to_keep)/len(full_chain.nuclides)):.1f}%")
