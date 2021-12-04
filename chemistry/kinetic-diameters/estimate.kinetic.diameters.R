# Kinetic diameters for molecules are hard to find.
# We need kinetic diameters for a lot of molecules,
# however our model does not strive to make publishable predictions,
# so if necessary we can get by with representative estimates.
# Even a crude way to estimate kinetic diameters would be extremely useful.
#
# As mentioned by Mehio (2014), empirically derived kinetic diameters closely 
# match the minimum diameter of a molecule's isoelectronic density surface 
# at which electron density exceeds 0.0015 e/a0^3, where a0 is the Bohr radius
# This "QM diameter" still requires an expensive software package like TURBOMOLE
# to derive, however Mehio's approach allows us to make guesses for cheap.
# If we assume each element has a constant diameter for its contribution to 
# a molecule's isoelectronic density surface, we can use well known 
# angles and distances within a molecule to determine the minimum diameter 
# of the entire molecule.
#
# We first need to make estimates for these atomic diameters
# A few can be easily derived from the table Mehio provides.
# By knowing that kinetic diameter is the minimum of any possible cross section,
# we can assume some of these atomic diameter estimates are equivalent to the 
# kinetic diameter of a molecule.
#
# NOTE: all our measures of distance are in picometers

He_diameter = 255.7 # direct estimate
H_diameter = 287.7 # from H2
O_diameter = 334.0 # from O2
N_diameter = 357.8 # from N2
C_diameter = 346.9 # mean of the estimates from CO and CO2

# We'll combine this from estimates we can derive from 
# [wikipedia](https://en.wikipedia.org/wiki/Kinetic_diameter):

Ne_diameter = 275 # direct estimate
Ar_diameter = 340 # direct estimate
Kr_diameter = 360 # direct estimate
Xe_diameter = 396 # direct estimate
Cl_diameter = 320 # from both Cl2 and HCl
Br_diameter = 350 # from both Br2 and HBr

# Let's plot them:
# atomic number
D = c(H_diameter, He_diameter, C_diameter, N_diameter, O_diameter, Ne_diameter, Cl_diameter, Ar_diameter, Br_diameter, Kr_diameter, Xe_diameter)
G = c(1,          0,           3,          4,          5,          0,           7,           0,           7,           0,           0          ) 
N = c(1,          2,           6,          7,          8,          10,          17,          18,          35,          36,          54         ) 
row = floor((N-3)/8)+2
plot(row,D, col=G+1); lm(D~row)
plot(G,D, col=row); lm(D~G)
plot(N,D, col=G+1); lm(D~N)

# calculated radius: good fit, may require seperate intercepts for group, 
# but preliminary data suggests we could get by with a single slope
# seems like intercept may be a sinusoid function with 
R = c(53,         31,          67,         56,         48,         38,          79,          71,          94,          88,          108        ) 
plot(R,D, col=G+1); lm(D~R)

# empirical radius: terrible fit
R = c(25,         120,         70,         65,         60,         160,         100,         71,          115,         88,          108        ) 
plot(R,D, col=G+1); lm(D~R)

# van der waals radius: good fit within periodic groups, but may require separate slopes/intercepts for each
R = c(120,        140,         170,        155,        152,        154,         175,         188,         185,         202,         216        ) 
plot(R,D, col=G+1); lm(D~R)

# covalent radius: good fit except for Ne
R = c(38,         32,          77,         75,         73,         69,          99,          97,          114,         110,         130        ) 
plot(R,D, col=G+1); lm(D~R)



# I'd like to find the diameter for F, but I only have the diameter for SF6
# I first need to find the diameter for S, which is easiest to find through H2S.
# H2S is a bent nonlinear molecule so I have to find the relation that derives its
# kinetic diameter from its angles and distances:

get_kinetic_diameter_for_bent_molecule = function(
	bond_length, 
	bond_angle,
	inner_atom_diameter, 
	outer_atom_diameter
) {
	max(inner_atom_diameter, outer_atom_diameter, 
		min(
			bond_length * cos(bond_angle/2) + outer_atom_diameter/2 + inner_atom_diameter/2,
			bond_length * sin(bond_angle/2) + outer_atom_diameter
		)
	)
}

# So working backwards from this:

H2S_diameter = 320
H2S_bond_length = 133.6
H2S_bond_angle = 92.1 * pi/180
S_diameter = ((H2S_diameter - H_diameter/2) - (H2S_bond_length * cos(H2S_bond_angle/2))) * 2

# however we find the diameter we derive for S is unlikely
# because S_diameter << H_diameter. In fact, we even find S_diameter << He_diameter. 

# How about we try SO2 instead?

SO2_diameter = 360
SO2_bond_length = 143.1
SO2_bond_angle = 119 * pi/180
S_diameter = ((SO2_diameter - H_diameter/2) - (SO2_bond_length * cos(SO2_bond_angle/2))) * 2

# we still get ridiculous results when we feed our estimates 