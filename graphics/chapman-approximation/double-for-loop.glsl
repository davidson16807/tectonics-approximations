
const float PI = 3.14159265358979323846264338327950288419716939937510;
const float PHI = 1.6180339887;
const float BIG = 1e20;
const float SMALL = 1e-20;

struct maybe_int
{
    int value;
    bool exists;
};
struct maybe_float
{
    float value;
    bool exists;
};
struct maybe_vec2
{
    vec2 value;
    bool exists;
};
struct maybe_vec3
{
    vec3 value;
    bool exists;
};

/*
MENSURATION
this file contains functions for finding perimeters, areas, 
surface areas, and volumes of primitive shapes
*/
float get_surface_area_of_sphere(
    in float radius
) {
    return 4.*PI*radius*radius;
}
//#include "precompiled/academics/math/geometry/point_intersection.glsl"

maybe_vec2 get_bounding_distances_along_ray(in maybe_vec2 distances_along_line){
    return 
      maybe_vec2(
        vec2(
          max(min(distances_along_line.value.x, distances_along_line.value.y), 0.0),
          max(distances_along_line.value.x, distances_along_line.value.y)
        ),
        distances_along_line.exists && max(distances_along_line.value.x, distances_along_line.value.y) > 0.
      );
}
maybe_float get_nearest_distance_along_ray(in maybe_vec2 distances_along_line){
    return 
      maybe_float(
        distances_along_line.value.x < 0.0? distances_along_line.value.y :
        distances_along_line.value.y < 0.0? distances_along_line.value.x :
        min(distances_along_line.value.x, distances_along_line.value.y),
        distances_along_line.exists && max(distances_along_line.value.x, distances_along_line.value.y) > 0.
      );
}
maybe_float get_distance_along_line_to_union(
    in maybe_float shape1,
    in maybe_float shape2
) {
    return maybe_float(
        !shape1.exists ? shape2.value : !shape2.exists ? shape1.value : min(shape1.value, shape2.value),
        shape1.exists || shape2.exists
    );
}
maybe_vec2 get_distances_along_line_to_union(
    in maybe_vec2 shape1,
    in maybe_vec2 shape2
) {
    return maybe_vec2(
        vec2(!shape1.exists ? shape2.value.x : !shape2.exists ? shape1.value.x : min(shape1.value.x, shape2.value.x),
             !shape1.exists ? shape2.value.y : !shape2.exists ? shape1.value.y : max(shape1.value.y, shape2.value.y )),
        shape1.exists || shape2.exists
    );
}
maybe_vec2 get_distances_along_line_to_negation(
    in maybe_vec2 positive,
    in maybe_vec2 negative
) {
    // as long as intersection with positive exists, 
    // and negative doesn't completely surround it, there will be an intersection
    bool exists =
        positive.exists && !(negative.value.x < positive.value.x && positive.value.y < negative.value.y);
    // find the first region of intersection
    float entrance = !negative.exists ? positive.value.x : min(negative.value.y, positive.value.x);
    float exit = !negative.exists ? positive.value.y : min(negative.value.x, positive.value.y);
    // if the first region is behind us, find the second region
    if (exit < 0. && 0. < positive.value.y)
    {
        entrance = negative.value.y;
        exit = positive.value.y;
    }
    return maybe_vec2( vec2(entrance, exit), exists );
}

maybe_vec2 get_distances_along_line_to_intersection(
    in maybe_vec2 shape1,
    in maybe_vec2 shape2
) {
    float x = shape1.exists && shape2.exists ? max(shape1.value.x, shape2.value.x) : 0.0;
    float y = shape1.exists && shape2.exists ? min(shape1.value.y, shape2.value.y ) : 0.0;
    return maybe_vec2(vec2(x,y), shape1.exists && shape2.exists && x < y);
}


/*
A0 line reference
A  line direction, normalized
B0 plane reference
N  plane surface normal, normalized
*/
maybe_float get_distance_along_3d_line_to_plane(
    in vec3 A0,
    in vec3 A,
    in vec3 B0,
    in vec3 N
){
    return maybe_float( -dot(A0 - B0, N) / dot(A, N), abs(dot(A, N)) < SMALL);
}
/*
A0 line reference
A  line direction, normalized
B0 sphere origin
R  sphere radius along each coordinate axis
*/
maybe_vec2 get_distances_along_3d_line_to_sphere(
    in vec3 A0,
    in vec3 A,
    in vec3 B0,
    in float r
){
    float t = dot(B0 - A0, A);
    vec3  At = A0 + A*t - B0;
    float y2 = r*r - dot(At,At);
    float dxr = sqrt(max(y2, SMALL));
    return maybe_vec2(
        vec2(t - dxr, t + dxr),
        y2 > 0.
    );
}


/*
A0 line reference
A  line direction, normalized
B0 cylinder reference
B  cylinder direction, normalized
r  cylinder radius
*/
maybe_vec2 get_distances_along_3d_line_to_infinite_cylinder(
    in vec3 A0,
    in vec3 A,
    in vec3 B0,
    in vec3 B,
    in float r
){
    // INTUITION: simplify the problem by using a coordinate system based around the line and the tube center
    // see closest-approach-between-line-and-cylinder-visualized.scad
    // implementation shamelessly copied from Inigo: 
    // https://www.iquilezles.org/www/articles/intersectors/intersectors.htm
    vec3 D = A0 - B0;
    float BA = dot(B, A);
    float BD = dot(B, D);
    float a = 1.0 - BA * BA;
    float b = dot(D, A) - BD * BA;
    float c = dot(D, D) - BD * BD - r * r;
    float h = sqrt(max(b * b - a * c, 0.0));
    return maybe_vec2(
        vec2((-b + h) / a, (-b - h) / a),
        h > 0.0
    );
}
/*
A0 line reference
A  line direction, normalized
B1 cylinder endpoint 1
B2 cylinder endpoing 2
r  cylinder radius
*/
maybe_vec2 get_distances_along_3d_line_to_cylinder(
    in vec3 A0,
    in vec3 A,
    in vec3 B1,
    in vec3 B2,
    in float r
){
    vec3 B = normalize(B2 - B1);
    maybe_float a1 = get_distance_along_3d_line_to_plane(A0, A, B1, B);
    maybe_float a2 = get_distance_along_3d_line_to_plane(A0, A, B2, B);
    float a_in = min(a1.value, a2.value);
    float a_out = max(a1.value, a2.value);
    maybe_vec2 ends = maybe_vec2(vec2(a_in, a_out), a1.exists || a2.exists);
    maybe_vec2 tube = get_distances_along_3d_line_to_infinite_cylinder(A0, A, B1, B, r);
    maybe_vec2 cylinder = get_distances_along_line_to_intersection(tube, ends);
    // TODO: do we need this line?
    float entrance = max(tube.value.y, a_in);
    float exit = min(tube.value.x, a_out);
    return maybe_vec2(
        vec2(entrance, exit),
        tube.exists && entrance < exit
    );
}
/*
A0 line reference
A  line direction, normalized
B1 capsule endpoint 1
B2 capsule endpoing 2
r  capsule radius
*/
maybe_vec2 get_distances_along_3d_line_to_capsule(
    in vec3 A0,
    in vec3 A,
    in vec3 B1,
    in vec3 B2,
    in float r
){
    maybe_vec2 cylinder = get_distances_along_3d_line_to_cylinder(A0, A, B1, B2, r);
    maybe_vec2 sphere1 = get_distances_along_3d_line_to_sphere(A0, A, B1, r);
    maybe_vec2 sphere2 = get_distances_along_3d_line_to_sphere(A0, A, B2, r);
    maybe_vec2 spheres = get_distances_along_line_to_union(sphere1, sphere2);
    maybe_vec2 capsule = get_distances_along_line_to_union(spheres, cylinder);
    return capsule;
}

/*
A0 point position
B0 sphere origin
r  radius
*/
vec3 get_surface_normal_of_point_near_sphere( in vec3 A0, in vec3 B0 )
{
    return normalize( A0-B0 );
}


const float KELVIN = 1.;
const float KILOGRAM = 1.; // kilograms
const float METER = 1.; // meters
const float SECOND = 1.; // seconds
const float NEWTON = KILOGRAM * METER / (SECOND * SECOND);
const float JOULE = NEWTON * METER;
const float WATT = JOULE / SECOND;
const float EARTH_RADIUS = 6.367e6; // meters
const float ASTRONOMICAL_UNIT = 149597870700.;// meters
const float GLOBAL_SOLAR_CONSTANT = 1361.; // watts/meter^2
const float SOLAR_MASS = 2e30; // kilograms
const float SOLAR_RADIUS = 695.7e6; // meters
const float SOLAR_LUMINOSITY = 3.828e26; // watts
const float SOLAR_TEMPERATURE = 5772.; // kelvin


const float SPEED_OF_LIGHT = 299792458. * METER / SECOND;
const float BOLTZMANN_CONSTANT = 1.3806485279e-23 * JOULE / KELVIN;
const float STEPHAN_BOLTZMANN_CONSTANT = 5.670373e-8 * WATT / (METER*METER* KELVIN*KELVIN*KELVIN*KELVIN);
const float PLANCK_CONSTANT = 6.62607004e-34 * JOULE * SECOND;

// see Lawson 2004, "The Blackbody Fraction, Infinite Series and Spreadsheets"
// we only do a single iteration with n=1, because it doesn't have a noticeable effect on output
float solve_fraction_of_light_emitted_by_black_body_below_wavelength(
    in float wavelength,
    in float temperature
){
    const float iterations = 2.;
    const float h = PLANCK_CONSTANT;
    const float k = BOLTZMANN_CONSTANT;
    const float c = SPEED_OF_LIGHT;
    float L = wavelength;
    float T = temperature;
    float C2 = h*c/k;
    float z = C2 / (L*T);
    float z2 = z*z;
    float z3 = z2*z;
    float sum = 0.;
    float n2=0.;
    float n3=0.;
    for (float n=1.; n <= iterations; n++) {
        n2 = n*n;
        n3 = n2*n;
        sum += (z3 + 3.*z2/n + 6.*z/n2 + 6./n3) * exp(-n*z) / n;
    }
    return 15.*sum/(PI*PI*PI*PI);
}
float solve_fraction_of_light_emitted_by_black_body_between_wavelengths(
    in float lo,
    in float hi,
    in float temperature
){
    return solve_fraction_of_light_emitted_by_black_body_below_wavelength(hi, temperature) -
            solve_fraction_of_light_emitted_by_black_body_below_wavelength(lo, temperature);
}
// This calculates the radiation (in watts/m^2) that's emitted 
// by a single object using the Stephan-Boltzmann equation
float get_intensity_of_light_emitted_by_black_body(
    in float temperature
){
    float T = temperature;
    return STEPHAN_BOLTZMANN_CONSTANT * T*T*T*T;
}
vec3 solve_rgb_intensity_of_light_emitted_by_black_body(
    in float temperature
){
    return get_intensity_of_light_emitted_by_black_body(temperature)
         * vec3(
             solve_fraction_of_light_emitted_by_black_body_between_wavelengths(600e-9*METER, 700e-9*METER, temperature),
             solve_fraction_of_light_emitted_by_black_body_between_wavelengths(500e-9*METER, 600e-9*METER, temperature),
             solve_fraction_of_light_emitted_by_black_body_between_wavelengths(400e-9*METER, 500e-9*METER, temperature)
           );
}
// Rayleigh phase function factor [-1, 1]
float get_fraction_of_rayleigh_scattered_light_scattered_by_angle(
    in float cos_scatter_angle
){
    return 3. * (1. + cos_scatter_angle*cos_scatter_angle)
    / //------------------------
                (16. * PI);
}
// Henyey-Greenstein phase function factor [-1, 1]
// represents the average cosine of the scattered directions
// 0 is isotropic scattering
// > 1 is forward scattering, < 1 is backwards
float get_fraction_of_mie_scattered_light_scattered_by_angle(
    in float cos_scatter_angle
){
    const float g = 0.76;
    return (1. - g*g)
    / //---------------------------------------------
        ((4. + PI) * pow(1. + g*g - 2.*g*cos_scatter_angle, 1.5));
}
// Schlick's fast approximation to the Henyey-Greenstein phase function factor
// Pharr and  Humphreys [2004] equivalence to g above
float approx_fraction_of_mie_scattered_light_scattered_by_angle_fast(
    in float cos_scatter_angle
){
    const float g = 0.76;
    const float k = 1.55*g - 0.55 * (g*g*g);
    return (1. - k*k)
    / //-------------------------------------------
        (4. * PI * (1. + k*cos_scatter_angle) * (1. + k*cos_scatter_angle));
}
/*
"get_fraction_of_microfacets_accessible_to_ray" is Schlick's fast approximation for Smith's function
  see Hoffmann 2015 for a gentle introduction to the concept
  see Schlick (1994) for even more details.
*/
float get_fraction_of_microfacets_accessible_to_ray(
    in float cos_view_angle,
    in float root_mean_slope_squared
){
    float m = root_mean_slope_squared;
    float v = cos_view_angle;
    // float k = m/2.0; return 2.0*v/(v+sqrt(m*m+(1.0-m*m)*v*v)); // Schlick-GGX
    float k = m*sqrt(2./PI); return v/(v*(1.0-k)+k); // Schlick-Beckmann
}
/*
"get_fraction_of_microfacets_with_angle" 
  This is also known as the Beckmann Surface Normal Distribution Function.
  This is the probability of finding a microfacet whose surface normal deviates from the average by a certain angle.
  see Hoffmann 2015 for a gentle introduction to the concept.
  see Schlick (1994) for even more details.
*/
float get_fraction_of_microfacets_with_angle(
    in float cos_angle_of_deviation,
    in float root_mean_slope_squared
){
    float m = root_mean_slope_squared;
    float t = cos_angle_of_deviation;
    float m2 = m*m;
    float t2 = t*t;
    float u = t2*(m2-1.0)+1.0; return m2/(PI*u*u);
    //return exp((t*t-1.)/max(m*m*t*t, 0.1))/max(PI*m*m*t*t*t*t, 0.1);
}
/*
"get_fraction_of_light_reflected_from_facet_head_on" finds the fraction of light that's reflected
  by a boundary between materials when striking head on.
  It is also known as the "characteristic reflectance" within the fresnel reflectance equation.
  The refractive indices can be provided as parameters in any order.
*/
float get_fraction_of_light_reflected_from_facet_head_on(
    in float refractivate_index1,
    in float refractivate_index2
){
    float n1 = refractivate_index1;
    float n2 = refractivate_index2;
    float sqrtF0 = ((n1-n2)/(n1+n2));
    float F0 = sqrtF0 * sqrtF0;
    return F0;
}
/*
"get_rgb_fraction_of_light_reflected_from_facet" returns Fresnel reflectance for each color channel.
  Fresnel reflectance is the fraction of light that's immediately reflected upon striking the surface.
  It is the fraction of light that causes specular reflection.
  Here, we use Schlick's fast approximation for Fresnel reflectance.
  see https://en.wikipedia.org/wiki/Schlick%27s_approximation for a summary 
  see Hoffmann 2015 for a gentle introduction to the concept
  see Schlick (1994) for implementation details
*/
vec3 get_rgb_fraction_of_light_reflected_from_facet(
    in float cos_incident_angle,
    in vec3 characteristic_reflectance
){
    vec3 F0 = characteristic_reflectance;
    float _1_u = 1.-cos_incident_angle;
    return F0 + (1.-F0) * _1_u*_1_u*_1_u*_1_u*_1_u;
}
/*
"get_fraction_of_light_reflected_from_material" is a fast approximation to the Cook-Torrance Specular BRDF.
  It is the fraction of light that reflects from a material to the viewer.
  see Hoffmann 2015 for a gentle introduction to the concept
*/
vec3 get_fraction_of_light_reflected_from_material(
    in float NL, in float NH, in float NV, in float HV,
    in float root_mean_slope_squared,
    in vec3 characteristic_reflectance
){
    float m = root_mean_slope_squared;
    vec3 F0 = characteristic_reflectance;
    return 1.0
        * get_fraction_of_microfacets_accessible_to_ray(NL, m)
        * get_fraction_of_microfacets_with_angle(NH, m)
        * get_fraction_of_microfacets_accessible_to_ray(NV, m)
        * get_rgb_fraction_of_light_reflected_from_facet(HV, F0)
        / max(4.*PI*NV*NL, 0.001);
}
/*
"GAMMA" is the constant that's used to map between 
rgb signals sent to a monitor and their actual intensity
*/
const float GAMMA = 2.2;
/* 
This function returns a rgb vector that quickly approximates a spectral "bump".
Adapted from GPU Gems and Alan Zucconi
from https://www.alanzucconi.com/2017/07/15/improving-the-rainbow/
*/
vec3 get_rgb_intensity_of_rgb_signal(
    in vec3 signal
){
    return vec3(
        pow(signal.x, GAMMA),
        pow(signal.y, GAMMA),
        pow(signal.z, GAMMA)
    );
}
/*
This function returns a rgb vector that best represents color at a given wavelength
It is from Alan Zucconi: https://www.alanzucconi.com/2017/07/15/improving-the-rainbow/
I've adapted the function so that coefficients are expressed in meters.
*/
vec3 get_rgb_signal_of_rgb_intensity(
    in vec3 intensity
){
    return vec3(
        pow(intensity.x, 1./GAMMA),
        pow(intensity.y, 1./GAMMA),
        pow(intensity.z, 1./GAMMA)
    );
}
mat4 get_rotation_matrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    return mat4(oc * axis.x * axis.x + c, oc * axis.x * axis.y - axis.z * s, oc * axis.z * axis.x + axis.y * s, 0.0,
                oc * axis.x * axis.y + axis.z * s, oc * axis.y * axis.y + c, oc * axis.y * axis.z - axis.x * s, 0.0,
                oc * axis.z * axis.x - axis.y * s, oc * axis.y * axis.z + axis.x * s, oc * axis.z * axis.z + c, 0.0,
                0.0, 0.0, 0.0, 1.0);
}
mat4 get_translation_matrix(vec3 offset)
{
    return mat4(1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                offset,1);
}

float F(
    in float x,
    in float y2,
    in float r0
){
    const float SQRT_HALF_PI = sqrt(PI/2.);
    const float k = 0.6; // "k" is an empirically derived constant
    float abs_x  = abs(x);
    float r      = sqrt(x*x+y2);
    float ch     = (1. - 1./(2.*r)) * SQRT_HALF_PI * sqrt(sqrt(y2)) + k*abs_x;
    float F      = exp(r0-r) / (abs_x/r + 1./ch);
    return F;
}

// "approx_air_column_density_ratio_through_atmosphere" 
//   calculates the distance you would need to travel 
//   along the surface to encounter the same number of particles in the column. 
// It does this by finding an integral using integration by substitution, 
//   then tweaking that integral to prevent division by 0. 
// All distances are recorded in scale heights.
// "a" and "b" are distances along the ray from closest approach.
//   The ray is fired in the positive direction.
//   If there is no intersection with the planet, 
//   a and b are distances from the closest approach to the upper bound.
// "z2" is the closest distance from the ray to the center of the world, squared.
// "r0" is the radius of the world.
float approx_air_column_density_ratio_through_atmosphere(
    in float a,
    in float b,
    in float z2,
    in float r0
){
    // GUIDE TO VARIABLE NAMES:
    //  "x*" distance along the ray from closest approach
    //  "z*" distance from the center of the world at closest approach
    //  "r*" distance ("radius") from the center of the world
    //  "*0" variable at reference point
    //  "*2" the square of a variable
    //  "ch" a nudge we give to prevent division by zero, analogous to the Chapman function
    const float SQRT_HALF_PI = sqrt(PI/2.);
    const float k = 0.6; // "k" is an empirically derived constant
    float x0 = sqrt(max(r0*r0 - z2, SMALL));
    return max( sign(b)*(F(x0,z2,r0)-F(b,z2,r0)) - sign(a)*(F(x0,z2,r0)-F(a,z2,r0)), 0.0 );
}

#ifndef PROD
#define ASSERT(test, color) if (!(test)) { return color; }
#else
#define ASSERT(test, color)
#endif

// TODO: multiple scattering events
// TODO: support for light sources from within atmosphere
vec3 get_rgb_fraction_of_distant_light_scattered_by_atmosphere(
    vec3 view_origin, vec3 view_direction, float view_start_length, float view_stop_length,
    vec3 world_position, float world_radius,
    vec3 light_direction, float atmosphere_scale_height,
    vec3 beta_ray, vec3 beta_mie, vec3 beta_abs
){
    // For an excellent introduction to what we're try to do here, see Alan Zucconi: 
    //   https://www.alanzucconi.com/2017/10/10/atmospheric-scattering-3/
    // We will be using most of the same terminology and variable names.
    // GUIDE TO VARIABLE NAMES:
    //  Uppercase letters indicate vectors.
    //  Lowercase letters indicate scalars.
    //  Going for terseness because I tried longhand names and trust me, you can't read them.
    //  "*v*"    property of the view ray, the ray cast from the viewer to the object being viewed
    //  "*l*"    property of the light ray, the ray cast from the object to the light source
    //  "y*"     distance from the center of the world to the plane shared by view and light ray
    //  "z*"     distance from the center of the world to along the plane shared by the view and light ray 
    //  "r*"     a distance ("radius") from the center of the world
    //  "h*"     the atmospheric scale height, the distance at which air density reduces by a factor of e
    //  "*2"     the square of a variable
    //  "*0"     property at the start of the raymarch
    //  "*1"     property at the end of the raymarch
    //  "*i"     property during an iteration of the raymarch
    //  "d*"     the change in a property across iterations of the raymarch
    //  "beta*"  a scattering coefficient, the number of e-foldings in light intensity per unit distance
    //  "gamma*" a phase factor, the fraction of light that's scattered in a certain direction
    //  "sigma*" a column density ratio, the density of a column of air relative to surface density
    //  "F*"     fraction of source light that reaches the viewer due to scattering for each color channel
    //  "*_ray"  property of rayleigh scattering
    //  "*_mie"  property of mie scattering
    //  "*_abs"  property of absorption
    // setup variable shorthands
    // express all distances in scale heights 
    // express all positions relative to world origin
    float h = atmosphere_scale_height;
    float r = world_radius / h;
    vec3 V0 = (view_origin + view_direction * view_start_length - world_position) / h;
    vec3 V1 = (view_origin + view_direction * view_stop_length - world_position) / h;
    vec3 V = view_direction; // unit vector pointing to pixel being viewed
    float v0 = dot(V0,V);
    float v1 = dot(V1,V);
    vec3 L = light_direction; // unit vector pointing to light source
    float VL = dot(V,L);
    // "gamma_*" indicates the fraction of scattered sunlight that scatters to a given angle (indicated by its cosine, A.K.A. "VL").
    // It only accounts for a portion of the sunlight that's lost during the scatter, which is irrespective of wavelength or density
    float gamma_ray = get_fraction_of_rayleigh_scattered_light_scattered_by_angle(VL);
    float gamma_mie = get_fraction_of_mie_scattered_light_scattered_by_angle(VL);
    // "beta_*" indicates the rest of the fractional loss.
    // it is dependant on wavelength, and the density ratio, which is dependant on height
    // So all together, the fraction of sunlight that scatters to a given angle is: beta(wavelength) * gamma(angle) * density_ratio(height)
    vec3 beta_sum = h*(beta_ray + beta_mie + beta_abs);
    vec3 beta_gamma = h*(beta_ray * gamma_ray + beta_mie * gamma_mie);
    // number of iterations within the raymarch
    const float STEP_COUNT = 16.;
    float dv = (v1 - v0) / STEP_COUNT;
    float dl = dv*VL;
    float l0 = dot(V0,L);
    float y  = dot(V0,normalize(cross(V,L)));
    float y2 = y*y;
    float zv2 = dot(V0,V0) - y2 - v0*v0;
    float zl2 = 0.0;
    float vi = 0.0;
    float li = 0.0;
    float dl_dv = 0.0;
    float dzl2_dv = 0.0;
    float sigma; // columnar density encountered along the entire path, relative to surface density, effectively the distance along the surface needed to obtain a similar column density
    vec3 H = vec3(0); // total intensity for each color channel, found as the sum of light intensities for each path from the light source to the camera
    float R = 3.0*r; // some arbitrarily large number, technically the distance to the light source
    for (float i = 0.; i < STEP_COUNT; ++i)
    {
        vi = dv*i + v0;
        li = VL*(vi-v0) + l0;
        zl2 = vi*vi + zv2 - li*li;
        
        dl_dv = VL;
        dzl2_dv = 2.0*vi - 2.0*li*dl_dv;
        
        // do not add if light is obstructed by the planet
        if (li<0.0 && y2+zl2 < r*r){ continue; }
         
        float x0v = sqrt(max(r*r - y2 - zv2, SMALL));
        float x0l = sqrt(max(r*r - y2 - zl2, SMALL));
        sigma = max( sign(vi)*(F(x0v,y2+zv2,r)-F(vi,y2+zv2,r)) - sign(v0)*(F(x0v,y2+zv2,r)-F(v0,y2+zv2,r)), 0.0 )
              + max( sign(R) *(F(x0l,y2+zl2,r)-F(R, y2+zl2,r)) - sign(li)*(F(x0l,y2+zl2,r)-F(li,y2+zl2,r)), 0.0 );
        H += exp(r-sqrt(vi*vi+y2+zv2) - beta_sum*sigma) * dv;
    }
    return H * beta_gamma;
}


vec3 get_rgb_fraction_of_light_transmitted_through_atmosphere(
    in vec3 view_origin, in vec3 view_direction, in float view_start_length, in float view_stop_length,
    in vec3 world_position, in float world_radius, in float atmosphere_scale_height,
    in vec3 beta_ray, in vec3 beta_mie, in vec3 beta_abs
){
    float h = atmosphere_scale_height;
    float r = world_radius / h;
    vec3 V0 = (view_origin + view_direction * view_start_length - world_position) / h;
    vec3 V1 = (view_origin + view_direction * view_stop_length - world_position) / h;
    vec3 V = view_direction; // unit vector pointing to pixel being viewed
    float v0 = dot(V0,V);
    float v1 = dot(V1,V);
    float zv2 = dot(V0,V0) - v0*v0;
    vec3 beta_sum = (beta_ray + beta_mie + beta_abs)*h;
    float sigma = approx_air_column_density_ratio_through_atmosphere(v0,v1,zv2,r);
    return exp(-sigma * beta_sum);
}
mat3 rotate_around_x( const in float angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(1, 0, 0, 0, _cos, -_sin, 0, _sin, _cos);
}



void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    const float SCALE_HEIGHT = 7994.0*METER;
    
    // world
    float earth_radius = 6360e3/SCALE_HEIGHT;      // surface of the world
    float atmo_radius = earth_radius+10.0;   // "top" of the atmosphere
    vec3 earth_origin = vec3(0,0,0);  // center of the world
    
    vec2 mouse = iMouse.xy / iResolution.xy;    
    mat4  projection_matrix_inverse = (mat4(0.5,0,0,0,
                                            0,0.3,0,0,
                                            0,0,0,-50,
                                            0,0,-1,50));
                                            
    vec2  screenspace    = fragCoord/iResolution.xy;
    
    
    bool IS_FIRST_PERSON_POV = screenspace.x > 0.5
    ;
    
    vec2  clipspace = vec2(0.5*(mod(4.0 * screenspace.x, 2.0) - 1.0), 2.0 * screenspace.y - 1.0);
    
    mat4  view_matrix_inverse;
    if(IS_FIRST_PERSON_POV){
        float altitude = 100.0/SCALE_HEIGHT;
        vec3 position = normalize(vec3(0.0,0.0, -1.0))*(earth_radius+altitude);
        vec3 up    = normalize(-position);
        vec3 right = vec3(1,0,0);
        vec3 front = normalize(cross(up, right));
        view_matrix_inverse = 
            get_translation_matrix(-position) *
            get_rotation_matrix(up, 1.2*PI+3.0*PI*mouse.x) * 
            get_rotation_matrix(right,0.1*PI-PI/4.0*mouse.y) * 
            inverse(mat4(vec4(right,0), vec4(up,0), vec4(front,0), vec4(0,0,0,1)));
    } else {
        view_matrix_inverse = 
            get_rotation_matrix(vec3(0,1,0), 6.3*mouse.x) * 
            get_rotation_matrix(vec3(1,0,0), PI/2.0+PI*mouse.y) *
            get_translation_matrix(vec3(0,0,6.0*earth_radius)) *
            mat4(1);
    }
    
    vec3  view_direction = normalize(view_matrix_inverse * projection_matrix_inverse * vec4(clipspace, 1, 1)).xyz;
    vec3  view_origin    = view_matrix_inverse[3].xyz;

    vec3  I =   solve_rgb_intensity_of_light_emitted_by_black_body(SOLAR_TEMPERATURE) 
              * get_surface_area_of_sphere(SOLAR_RADIUS) / get_surface_area_of_sphere(1.*ASTRONOMICAL_UNIT); // intensity of incoming light for each color channel
    vec3  E = vec3(0); // total intensity for each color channel, found as the sum of light intensities for each path from the light source to the camera
    

    // view ray
    vec3 V0 = view_origin;
    vec3 V  = view_direction;
    vec3 Vi = V0;            // point along view ray
    
    
    // light ray
    //vec3 L  = normalize(vec3(1,0,0)); // static
    vec3 L  = (get_rotation_matrix(vec3(0,1,0), 0.87*PI+0.2*PI*sin(1.0*iTime)) * // rotation of the earth
               get_rotation_matrix(vec3(0,0,1), -PI*23.4/180.) * // axial tilt
               vec4(1,0,0,1)).xyz;
    // realistic
    //vec3 L  = (get_rotation_matrix(vec3(0,1,0), 0.1*iTime) * vec4(1,0,0,1)).xyz;

    // "beta_*" is the rest of the fractional loss.
    // it is dependant on wavelength, and the density ratio
    // So all together, the fraction of sunlight that scatters to a given angle is: beta(wavelength) * gamma(angle) * density_ratio
    vec3  beta_ray_air   = vec3(5.20e-6, 1.21e-5, 2.96e-5) * SCALE_HEIGHT;
    vec3  beta_mie_air   = vec3(1e-6)*SCALE_HEIGHT;
    vec3  beta_abs_air   = vec3(0);

    maybe_vec2  air_along_view_ray  = 
        get_bounding_distances_along_ray(
            get_distances_along_line_to_negation(
                get_distances_along_3d_line_to_sphere(V0, V, earth_origin, atmo_radius),
                get_distances_along_3d_line_to_sphere(V0, V, earth_origin, earth_radius)
            ) 
        );
    maybe_float world_along_view_ray = 
        get_nearest_distance_along_ray(
            get_distances_along_3d_line_to_sphere(V0, V, earth_origin, earth_radius)
        );
    
    
    if(!air_along_view_ray.exists && !world_along_view_ray.exists) 
    { 
        // nothing to see here folks, move along
        fragColor = vec4(0.0);
        return;
    }
    
    
    if(world_along_view_ray.exists){
        vec3 Vt = V0+V*(world_along_view_ray.value-0.001);
        maybe_vec2 air_along_light_ray = 
            get_bounding_distances_along_ray(
                get_distances_along_line_to_negation(
                    get_distances_along_3d_line_to_sphere(Vt, L, earth_origin, atmo_radius),
                    get_distances_along_3d_line_to_sphere(Vt, L, earth_origin, earth_radius)
                )
            );
        maybe_vec2 world_along_light_ray = 
                get_distances_along_3d_line_to_sphere(Vt, L, earth_origin, earth_radius);
        
        vec3 I_surface = I;
        if (air_along_light_ray.exists)
        {
            I_surface *= get_rgb_fraction_of_light_transmitted_through_atmosphere(Vt,L, air_along_light_ray.value.x, air_along_light_ray.value.y, earth_origin, earth_radius, 1.0, beta_ray_air, beta_mie_air, beta_abs_air);
        }
        
        vec3 N = get_surface_normal_of_point_near_sphere(Vt, earth_origin);
        vec3 H = normalize(L-V);
        float NV = max(dot(N,-V), 0.);
        float NL = max(dot(N, L), 0.);
        float NH = max(dot(N, H), 0.);
        float HV = max(dot(H,-V), 0.);
        float VL = max(dot(L,-V), 0.);
        vec3 E_surface_ambient = vec3(1e-3); // surface ambient light
        vec3 F_surface_diffuse = vec3(0.0,0.01,0.1); // diffuse surface color
        
        E += E_surface_ambient;
        if (NL > 0.) //!world_along_light_ray.exists 
        {
            const vec3 F0  = vec3(0.04); // NOTE: "0.04" is a representative value for plastics and other diffuse reflectors
            const float m = 3.0;
            vec3 F_reflected = get_fraction_of_light_reflected_from_material(NL,NH,NV, HV, m,F0);
            vec3 E_surface_reflected = I_surface * NL * F_reflected;
            vec3 I_surface_refracted = I_surface * NL * (1. - F_reflected);
            vec3 E_surface_refracted = I_surface_refracted * F_surface_diffuse;
            E += E_surface_reflected + E_surface_refracted;
        }
    }
    if(air_along_view_ray.exists && (air_along_view_ray.value.x < world_along_view_ray.value || !world_along_view_ray.exists))
    {
        vec3 Vt0 = V0+V*(air_along_view_ray.value.x);
        vec3 Vt1 = V0+V*(air_along_view_ray.value.y);
        vec3 I_entrance = I;
        
        
        float v0 = air_along_view_ray.value.x;
        float v1 = air_along_view_ray.value.y;
        E *= get_rgb_fraction_of_light_transmitted_through_atmosphere(V0,V, air_along_view_ray.value.x,  air_along_view_ray.value.y,  earth_origin, earth_radius*0.999, 1.0, beta_ray_air, beta_mie_air, beta_abs_air);
        E += I_entrance * get_rgb_fraction_of_distant_light_scattered_by_atmosphere(
            V0, V, v0, v1, 
            earth_origin, earth_radius, L, 1.0, 
            beta_ray_air, beta_mie_air, beta_abs_air
        );
    }
    
    //fragColor = world_along_view_ray.exists? vec4(vec3(1),1) : vec4(vec3(0),1);
    float exposure_intensity = 15.; // Watts/m^2
    vec3  ldr_tone_map = 1.0 - exp(-E/exposure_intensity);

    fragColor = vec4(get_rgb_signal_of_rgb_intensity(ldr_tone_map), 1);
    //if((E.x==0.)){fragColor.x = 1.0;} else {fragColor.x=0.0;}

}