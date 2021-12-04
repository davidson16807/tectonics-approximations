from math import pi, sqrt, exp
import itertools
import plotly.express as px
import numpy as np
import multiprocessing

def clamp(lo, hi, x):
    return max(lo, min(hi, x))

def mix(lo, hi, x):
    return (hi-lo)*x+lo

def step(lo, hi, x):
    return (x-lo)/(lo-hi)


def rho(x,z,r):
    return exp(r-sqrt(x**2+z**2))

def quad(f,lo,hi,n):
    return sum(f(mix(lo,hi,i/n))*(hi-lo)/n for i in range(0,n))

def F(x,z,r,n):
    return quad(lambda xF: rho(xF,z,r), x, x+100, n)

def chapman(x,z,r,F):
    return 1/(rho(x,z,r)/F - abs(x)/sqrt(x**2+z**2))

def x0(z,r):
    return sqrt(max(r**2-z**2, 0))

def row(x,z,r,n):
    return (x,z,r, F(x,z,r,n)) 

def data(n,lx,lz,lr):
    inputs = [
        (x0(z,r)+dx,z,r,n)
        for dx in [mix(0,50,i/(lx-1)) for i in range(lx)]
        for z in [mix(10,510,i/(lz-1)) for i in range(lz)]
        for r in [mix(10,510,i/(lr-1)) for i in range(lr)]
    ]
    return np.array(multiprocessing.Pool().starmap(row, inputs))

data = data(10000,50,1000,20)
np.save('F.npy', data)
data = np.load('F.npy')
x = data[:,0]
z = data[:,1]
r = data[:,2]
F = data[:,3]

def np_rho(x,z,r):
    return np.exp(r-np.sqrt(x**2+z**2))

def np_chapman(x,z,r,F):
    return 1/(np_rho(x,z,r)/F - np.abs(x)/np.sqrt(x**2+z**2))

chapman = np_chapman(x,z,r,F)

def c(*nparrays):
    return np.concatenate([*nparrays])

px.scatter_3d(x=x, y=z, z=np.sqrt(chapman), animation_frame=r, range_x=[10,520], range_y=[0,520], range_z=[0,100], labels={'x':'x','y':'z','z':'1/Ch'}).show()
px.scatter_3d(x=x, y=z, z=chapman, animation_frame=r, range_x=[10,520], range_y=[0,520], range_z=[0,10000], labels={'x':'x','y':'z','z':'Ch'}).show()
px.scatter_3d(x=x, y=z, z=chapman, animation_frame=r, range_x=[10,520], range_y=[0,520], range_z=[0,100], labels={'x':'x','y':'z','z':'Ch'}).show()
px.scatter_3d(x=c(x,x)[0::2,], y=c(z,z)[0::2,], z=c(chapman,(1/(2*(z-x))+1)*sqrt(3.1415/2)*np.sqrt(z)+0.65*x)[0::2,], color=c(0*x,0*x+1)[0::2,], animation_frame=c(r,r)[0::2,], range_x=[10,120], range_y=[0,520], range_z=[0,100], labels={'x':'x','y':'z','z':'Ch'}).show()

SMALL = 1e-20
BIG = 1e20
def sign(x):
    return 1.0 if x>0.0 else 0.0 if x==0.0 else -1.0


#"approx_air_column_density_ratio_through_atmosphere" 
#  calculates the distance you would need to travel 
#  along the surface to encounter the same number of particles in the column. 
#It does this by finding an integral using integration by substitution, 
#  then tweaking that integral to prevent division by 0. 
#All distances are recorded in scale heights.
#"a" and "b" are distances along the ray from closest approach.
#  The ray is fired in the positive direction.
#  If there is no intersection with the planet, 
#  a and b are distances from the closest approach to the upper bound.
#"z2" is the closest distance from the ray to the center of the world, squared.
#"r0" is the radius of the world.
def approx_air_column_density_ratio_through_atmosphere(a, b, z2, r0):
    #GUIDE TO VARIABLE NAMES:
    # "x*" distance along the ray from closest approach
    # "z*" distance from the center of the world at closest approach
    # "r*" distance ("radius") from the center of the world
    # "*0" variable at reference point
    # "*2" the square of a variable
    # "ch" a nudge we give to prevent division by zero, analogous to the Chapman function
    SQRT_HALF_PI = sqrt(pi/2.)
    k = 0.6 #"k" is an empirically derived constant
    x0 = sqrt(max(r0*r0 - z2, SMALL))
    #if obstructed by the world, approximate answer by using a ludicrously large number
    if (a < x0 and -x0 < b and z2 < r0*r0): return BIG
    abs_a  = abs(a)
    abs_b  = abs(b)
    z      = sqrt(z2)
    sqrt_z = sqrt(z)
    ra     = sqrt(a*a+z2)
    rb     = sqrt(b*b+z2)
    ch0    = (1. - 1./(2.*r0)) * SQRT_HALF_PI * sqrt_z + k*x0
    cha    = (1. - 1./(2.*ra)) * SQRT_HALF_PI * sqrt_z + k*abs_a
    chb    = (1. - 1./(2.*rb)) * SQRT_HALF_PI * sqrt_z + k*abs_b
    s0     = min(exp(r0- z),1.) / (x0/r0 + 1./ch0)
    sa     = exp(r0-ra) / max(abs_a/ra + 1./cha, 0.01)
    sb     = exp(r0-rb) / max(abs_b/rb + 1./chb, 0.01)
    return max( sign(b)*(s0-sb) - sign(a)*(s0-sa), 0.0 )


def Hcalculation(v0, v1, y2, zv2, l0, VL, r, beta_sum):
    # For an excellent introduction to what we're try to do here, see Alan Zucconi: 
    #   https://www.alanzucconi.com/2017/10/10/atmospheric-scattering-3/
    # We will be using most of the same terminology and variable names.
    # GUIDE TO VARIABLE NAMES:
    #  Uppercase letters indicate vectors.
    #  Lowercase letters indicate scalars.
    #  Going for terseness because I tried longhand names and trust me, you can't read them.
    #  "*v*"    property of the view ray, the ray cast from the viewer to the object being viewed
    #  "*l*"    property of the light ray, the ray cast from the object to the light source
    #  "y*"     distance from the center of the world to the plane shared by view and light ray
    #  "z*"     distance from the center of the world to along the plane shared by the view and light ray 
    #  "r*"     a distance ("radius") from the center of the world
    #  "h*"     the atmospheric scale height, the distance at which air density reduces by a factor of e
    #  "*2"     the square of a variable
    #  "*0"     property at the start of the raymarch
    #  "*1"     property at the end of the raymarch
    #  "*i"     property during an iteration of the raymarch
    #  "d*"     the change in a property across iterations of the raymarch
    #  "beta*"  a scattering coefficient, the number of e-foldings in light intensity per unit distance
    #  "gamma*" a phase factor, the fraction of light that's scattered in a certain direction
    #  "sigma*" a column density ratio, the density of a column of air relative to surface density
    #  "F*"     fraction of source light that reaches the viewer due to scattering for each color channel
    #  "*_ray"  property of rayleigh scattering
    #  "*_mie"  property of mie scattering
    #  "*_abs"  property of absorption
    # setup variable shorthands
    # express all distances in scale heights 
    # express all positions relative to world origin
    # "beta_*" indicates the rest of the fractional loss.
    # it is dependant on wavelength, and the density ratio, which is dependant on height
    # So all together, the fraction of sunlight that scatters to a given angle is: beta(wavelength) * gamma(angle) * density_ratio(height)
    # number of iterations within the raymarch
    STEP_COUNT = 30
    dv = (v1 - v0) / STEP_COUNT
    vi = 0
    li = 0
    F = 0 # total intensity for each color channel, found as the sum of light intensities for each path from the light source to the camera
    for i in range(STEP_COUNT):
        vi = dv*i + v0
        li = VL*(vi-v0) + l0
        zl2 = max(0, vi*vi + zv2 - li*li)
        # "sigma": columnar density encountered along the entire path, relative to surface density, effectively the distance along the surface needed to obtain a similar column density
        sigma = approx_air_column_density_ratio_through_atmosphere(v0, vi, y2+zv2, r ) + approx_air_column_density_ratio_through_atmosphere(li, 3*r, y2+zl2, r )
        F += exp(r-sqrt(vi*vi+y2+zv2) - beta_sum*sigma) * dv
    return F

def Gcalculation(v0, vi, y2, zv2, l0, VL, r, beta_sum):
    li = VL*(vi-v0) + l0
    zl2 = max(0, vi*vi + zv2 - li*li)
    # "sigma": columnar density encountered along the entire path, relative to surface density, effectively the distance along the surface needed to obtain a similar column density
    sigma = approx_air_column_density_ratio_through_atmosphere(v0, vi, y2+zv2, r ) + approx_air_column_density_ratio_through_atmosphere(li, 3*r, y2+zl2, r )
    return exp(r-sqrt(vi*vi+y2+zv2) - beta_sum*sigma)


appreciable_rgb_fraction = 1e-30

inputs  = [
    (v0, v1, y2, zv2, l0, VL, r, beta_sum) 
    for r in  range(400, 600, 50)
    for y2 in [y*y for y in range(0, r, 50)]
    for zv2 in [zv*zv for zv in range(int(sqrt(r*r-y2)),r, 50)]
    for V0length in range(r, r+10, 1)
    for l0 in [V0length*i/(20-1) for i in range(20)]
    for v0 in [V0length*i/(20-1) for i in range(20)]
    for VL in [i/(20-1) for i in range(0,20,1)]
    for v1 in range(-r, v0, 50)
    for beta_sum in [i/1e6 * 8.5e3 for i in range(0,100, 10)]
]
len(inputs)


def Hrow(v0, v1, y2, zv2, l0, VL, r, beta_sum):
    return (v0, v1, y2, zv2, l0, VL, r, beta_sum, Hcalculation(v0, v1, y2, zv2, l0, VL, r, beta_sum)) 

Htable = np.array(multiprocessing.Pool().starmap(Hrow, inputs))
np.save('H.npy', Htable)

Gcolumn = np.array(multiprocessing.Pool().starmap(Gcalculation, inputs))

def table_function(table):
    def _table_function(*args):
        subtable = table
        for i, arg in enumerate(args):
            if arg is not None:
                values  = np.unique(table[:,i])
                distances = np.abs(values-arg)
                value = values[np.argmin(distances)]
                subtable = subtable[subtable[:,i] == value]
        return subtable
    return _table_function

H=table_function(Htable)

def Hview(r, y2, zv2, V0length, LV0hat, VV0hat, VL, v1):
    '''maps parameter space from a parameterization that can be better understood by users'''
    return H(V0length*VV0hat if V0length and VV0hat else None, v1, 
        y2, zv2, V0length*LV0hat if V0length and LV0hat else None, VL, r)


