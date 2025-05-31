/* ephemeris.c - Approximate Ephemeris Calculator using Keplerian Elements */
/* Compatible with K&R C style */

#include <stdio.h>
#include <math.h>
#include <string.h> /* For strcpy in main example */

/* Constants */
#define PI 3.14159265358979323846
#define DEG_TO_RAD (PI / 180.0)
#define RAD_TO_DEG (180.0 / PI)

#define JD_J2000 2451545.0        /* Julian Day of epoch J2000.0 */
#define DAYS_PER_CENTURY 36525.0

#define KEPLER_TOLERANCE 1.0e-9  /* Tolerance for solving Kepler's equation (radians) */
#define MAX_KEPLER_ITERATIONS 10

#define TRUE 1
#define FALSE 0
#define MAX_PLANETS 10 /* Maximum number of planets to load from file */

/* Structures */

/* Orbital elements for a planet at a reference epoch (e.g., J2000.0) */
/* All angles are in degrees, rates are per Julian century */
typedef struct {
    char name[32];
    double a0;          /* Semi-major axis (AU) at epoch */
    double e0;          /* Eccentricity at epoch */
    double I0_deg;      /* Inclination to ecliptic (degrees) at epoch */
    double M0_deg;      /* Mean anomaly (degrees) at epoch */
    double w0_deg;      /* Argument of perihelion (degrees) at epoch */
    double Omega0_deg;  /* Longitude of ascending node (Ω) (degrees) at epoch */

    /* Rates of change (per Julian century of 36525 days) */
    double da_dt;       /* AU/century */
    double de_dt;       /* /century */
    double dI_dt_deg;   /* degrees/century */
    double dM_dt_deg;   /* degrees/century */
    double dw_dt_deg;   /* degrees/century */
    double dOmega_dt_deg;/* degrees/century */
} OrbitalElements;

/* Calculated ephemeris data for a planet */
typedef struct {
    /* Heliocentric ecliptic rectangular coordinates (AU) */
    double x_h_ecl; /* X-coordinate */
    double y_h_ecl; /* Y-coordinate */
    double z_h_ecl; /* Z-coordinate */

    /* Geocentric ecliptic rectangular coordinates (AU) */
    /* Populated only if Earth's heliocentric position is also provided */
    double x_g_ecl;
    double y_g_ecl;
    double z_g_ecl;

    /* Geocentric ecliptic spherical coordinates */
    double lambda_g_deg; /* Ecliptic longitude (degrees) */
    double beta_g_deg;   /* Ecliptic latitude (degrees) */
    double delta_g_au;   /* Distance from Earth (AU) */
} PlanetEphem;


/* Function Prototypes (K&R C style) */
double calculate_julian_day();
double my_fmod(); /* K&R C compatible fmod */
double normalize_angle_deg();
double solve_kepler_equation();
void calculate_planet_helio_ecliptic_coords();
void calculate_geocentric_ecliptic_coords();
int initialize_orbital_elements_from_file(/* char *filename, OrbitalElements elements[], int *num_planets */);

/*
 * my_fmod
 * Calculates x modulo y (the remainder of x/y).
 * The result has the same sign as x.
 * This is a K&R C compatible implementation.
 */
double my_fmod(x, y)
double x;
double y;
{
    double quotient_val;

    if (y == 0.0) {
        /* Standard fmod behavior for y=0 is often a domain error or NaN.
         * For normalize_angle_deg, y is always 360.0, so this won't be hit.
         * Returning x is a simple fallback if an error can't be signaled.
         */
        return x;
    }

    quotient_val = x / y;

    /* K&R C: "When a float or double is converted to an int,
     * the fractional part is truncated." (truncates toward zero)
     */
    return x - (double)((long int)quotient_val) * y;
}

/*
 * calculate_julian_day
 * Computes the Julian Day (JD) for a given Gregorian date and time.
 * Valid for dates after 1582 October 15.
 */
double calculate_julian_day(year, month, day, hour, minute, second)
int year;
int month;
int day;
double hour;
double minute;
double second;
{
    int a, b;
    double jd;
    double day_fraction;

    if (month <= 2) {
        year -= 1;
        month += 12;
    }

    a = (int)(year / 100.0);
    b = 2 - a + (int)(a / 4.0); /* Gregorian calendar correction */

    day_fraction = (hour + minute / 60.0 + second / 3600.0) / 24.0;

    jd = (int)(365.25 * (year + 4716.0)) +
         (int)(30.6001 * (month + 1.0)) +
         day + b - 1524.5 + day_fraction;

    return jd;
}

/*
 * normalize_angle_deg
 * Normalizes an angle in degrees to the range [0, 360).
 */
double normalize_angle_deg(angle_deg)
double angle_deg;
{
    angle_deg = my_fmod(angle_deg, 360.0);
    if (angle_deg < 0.0) {
        angle_deg += 360.0;
    }
    return angle_deg;
}

/*
 * solve_kepler_equation
 * Solves Kepler's equation M = E - e * sin(E) for E (Eccentric Anomaly).
 * M_rad: Mean anomaly in radians.
 * e: Eccentricity.
 * Returns E in radians.
 */
double solve_kepler_equation(M_rad, e)
double M_rad;
double e;
{
    double E_rad;       /* Eccentric anomaly */
    double E_old_rad;
    double dE;
    int iterations;

    /* Initial guess for E */
    E_rad = M_rad + e * sin(M_rad);
    if (e > 0.8) { /* For high eccentricities, a different initial guess might be better */
        E_rad = PI;
    }


    iterations = 0;
    do {
        E_old_rad = E_rad;
        dE = (E_old_rad - e * sin(E_old_rad) - M_rad) / (1.0 - e * cos(E_old_rad));
        E_rad = E_old_rad - dE;
        iterations++;
    } while (fabs(dE) > KEPLER_TOLERANCE && iterations < MAX_KEPLER_ITERATIONS);

    return E_rad;
}

/*
 * calculate_planet_helio_ecliptic_coords
 * Calculates the heliocentric ecliptic coordinates of a planet for a given Julian Day.
 * elements: Orbital elements of the planet at epoch J2000.0.
 * jd: Julian Day for which to calculate the position.
 * ephem: Pointer to PlanetEphem structure to store the results.
 */
void calculate_planet_helio_ecliptic_coords(elements, jd, ephem)
OrbitalElements elements;
double jd;
PlanetEphem *ephem;
{
    double d_days;
    double T_centuries;

    /* Current orbital elements */
    double a_au;
    double e_ecc;
    double I_rad;
    double M_rad;
    double w_rad;       /* Argument of perihelion */
    double Omega_rad;   /* Longitude of ascending node */

    /* Intermediate values */
    double E_rad;       /* Eccentric anomaly */
    double nu_rad;      /* True anomaly */
    double r_au;        /* Heliocentric distance */

    double x_orb, y_orb; /* Coordinates in orbital plane */
    double cos_w, sin_w;
    double cos_Omega, sin_Omega;
    double cos_I, sin_I;

    /* Time parameters */
    d_days = jd - JD_J2000;
    T_centuries = d_days / DAYS_PER_CENTURY;

    /* Calculate orbital elements for the given date */
    a_au = elements.a0 + elements.da_dt * T_centuries;
    e_ecc = elements.e0 + elements.de_dt * T_centuries;
    I_rad = normalize_angle_deg(elements.I0_deg + elements.dI_dt_deg * T_centuries) * DEG_TO_RAD;
    M_rad = normalize_angle_deg(elements.M0_deg + elements.dM_dt_deg * T_centuries) * DEG_TO_RAD;
    w_rad = normalize_angle_deg(elements.w0_deg + elements.dw_dt_deg * T_centuries) * DEG_TO_RAD;
    Omega_rad = normalize_angle_deg(elements.Omega0_deg + elements.dOmega_dt_deg * T_centuries) * DEG_TO_RAD;

    /* Solve Kepler's equation for Eccentric Anomaly E */
    E_rad = solve_kepler_equation(M_rad, e_ecc);

    /* Calculate heliocentric distance r and true anomaly nu */
    r_au = a_au * (1.0 - e_ecc * cos(E_rad));
    nu_rad = atan2(sqrt(1.0 - e_ecc * e_ecc) * sin(E_rad), cos(E_rad) - e_ecc);
    /* Alternative for nu:
     * nu_rad = 2.0 * atan2(sqrt(1.0 + e_ecc) * sin(E_rad / 2.0), sqrt(1.0 - e_ecc) * cos(E_rad / 2.0));
     */


    /* Coordinates in the orbital plane (x_orb pointing to perihelion) */
    x_orb = r_au * cos(nu_rad);
    y_orb = r_au * sin(nu_rad);

    /* Transformation to heliocentric ecliptic coordinates */
    cos_w = cos(w_rad);
    sin_w = sin(w_rad);
    cos_Omega = cos(Omega_rad);
    sin_Omega = sin(Omega_rad);
    cos_I = cos(I_rad);
    sin_I = sin(I_rad);

    ephem->x_h_ecl = x_orb * (cos_w * cos_Omega - sin_w * sin_Omega * cos_I) -
                     y_orb * (sin_w * cos_Omega + cos_w * sin_Omega * cos_I);
    ephem->y_h_ecl = x_orb * (cos_w * sin_Omega + sin_w * cos_Omega * cos_I) +
                     y_orb * (-sin_w * sin_Omega + cos_w * cos_Omega * cos_I);
    ephem->z_h_ecl = x_orb * (sin_w * sin_I) +
                     y_orb * (cos_w * sin_I);
}

/*
 * calculate_geocentric_ecliptic_coords
 * Calculates geocentric ecliptic coordinates from heliocentric positions.
 * planet_helio: Heliocentric ephemeris of the planet.
 * earth_helio: Heliocentric ephemeris of the Earth.
 * planet_geo: Pointer to PlanetEphem structure to store geocentric results for the planet.
 */
void calculate_geocentric_ecliptic_coords(planet_helio, earth_helio, planet_geo)
PlanetEphem planet_helio; /* Heliocentric ephemeris of the planet */
PlanetEphem earth_helio;  /* Heliocentric ephemeris of the Earth */
PlanetEphem *planet_geo;  /* Output: Geocentric ephemeris of the planet */
{
    double xg, yg, zg; /* Geocentric rectangular ecliptic coordinates */
    double lambda_rad, beta_rad;

    /* Copy heliocentric data if it's the same struct */
    if (planet_geo != &planet_helio) {
        planet_geo->x_h_ecl = planet_helio.x_h_ecl;
        planet_geo->y_h_ecl = planet_helio.y_h_ecl;
        planet_geo->z_h_ecl = planet_helio.z_h_ecl;
    }

    xg = planet_helio.x_h_ecl - earth_helio.x_h_ecl;
    yg = planet_helio.y_h_ecl - earth_helio.y_h_ecl;
    zg = planet_helio.z_h_ecl - earth_helio.z_h_ecl;

    planet_geo->x_g_ecl = xg;
    planet_geo->y_g_ecl = yg;
    planet_geo->z_g_ecl = zg;

    /* Rectangular to Spherical (ecliptic longitude, latitude, distance) */
    planet_geo->delta_g_au = sqrt(xg*xg + yg*yg + zg*zg);
    lambda_rad = atan2(yg, xg);
    beta_rad = atan2(zg, sqrt(xg*xg + yg*yg));

    planet_geo->lambda_g_deg = normalize_angle_deg(lambda_rad * RAD_TO_DEG);
    planet_geo->beta_g_deg = beta_rad * RAD_TO_DEG; /* Latitude is -90 to +90 */
}

/*
 * initialize_orbital_elements_from_file
 * Reads orbital elements from a specified file.
 * filename: The path to the data file.
 * elements: Array to store the loaded OrbitalElements.
 * num_planets: Pointer to an integer to store the number of planets read.
 * Returns TRUE on success, FALSE on failure.
 */
int initialize_orbital_elements_from_file(filename, elements, num_planets)
char *filename;
OrbitalElements elements[];
int *num_planets;
{
    FILE *fp;
    char line1[256];
    char line2[256];
    int count = 0;
    char planet_name[32];
    double val_a, val_e, val_I, val_L, val_varpi, val_Omega; /* Epoch values */
    double rate_a, rate_e, rate_I, rate_L, rate_varpi, rate_Omega; /* Rates */

    fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Error opening ephemeris data file");
        return FALSE;
    }

    /* Skip header lines */
    if (fgets(line1, sizeof(line1), fp) == NULL) { fclose(fp); return FALSE; } /* Skip "a e I L Ω varpi" */
    if (fgets(line1, sizeof(line1), fp) == NULL) { fclose(fp); return FALSE; } /* Skip units line */

    while (count < MAX_PLANETS &&
           fgets(line1, sizeof(line1), fp) != NULL &&
           fgets(line2, sizeof(line2), fp) != NULL) {

        if (sscanf(line1, "%s %lf %lf %lf %lf %lf %lf",
                   planet_name, &val_a, &val_e, &val_I, &val_L, &val_varpi, &val_Omega) == 7 &&
            sscanf(line2, "%lf %lf %lf %lf %lf %lf",
                   &rate_a, &rate_e, &rate_I, &rate_L, &rate_varpi, &rate_Omega) == 6) {

            strcpy(elements[count].name, planet_name);
            elements[count].a0 = val_a;
            elements[count].e0 = val_e;
            elements[count].I0_deg = val_I;
            elements[count].Omega0_deg = val_Omega;
            /* M = L - varpi */
            elements[count].M0_deg = normalize_angle_deg(val_L - val_varpi);
            /* w = varpi - Omega */
            elements[count].w0_deg = normalize_angle_deg(val_varpi - val_Omega);

            elements[count].da_dt = rate_a;
            elements[count].de_dt = rate_e;
            elements[count].dI_dt_deg = rate_I;
            elements[count].dOmega_dt_deg = rate_Omega;
            /* dM/dt = dL/dt - dvarpi/dt */
            elements[count].dM_dt_deg = rate_L - rate_varpi;
            /* dw/dt = dvarpi/dt - dOmega/dt */
            elements[count].dw_dt_deg = rate_varpi - rate_Omega;

            count++;
        } else {
            /* Could print a warning for malformed lines if desired */
        }
    }

    fclose(fp);
    *num_planets = count;

    if (count == 0) {
        fprintf(stderr, "No planet data loaded from %s\n", filename);
        return FALSE;
    }

    return TRUE;
}



/* Example Usage */

int main() {
    OrbitalElements earth_elements, mars_elements;
    PlanetEphem earth_ephem, mars_ephem_helio, mars_ephem_geo;
    OrbitalElements all_planets[MAX_PLANETS];
    int num_planets_loaded;
    int i;
    int earth_idx = -1, mars_idx = -1;

    double jd;
    int year, month, day;
    double hour, minute, second;

    /* Initialize orbital elements from file */
    if (!initialize_orbital_elements_from_file("ephemeris_data.txt", all_planets, &num_planets_loaded)) {
        fprintf(stderr, "Failed to load ephemeris data.\n");
        return 1;
    }
    printf("Loaded %d planets from file.\n\n", num_planets_loaded);

    /* Find Earth and Mars in the loaded data */
    /* Note: "EMBary" is Earth-Moon Barycenter, often used for Earth's orbit */
    for (i = 0; i < num_planets_loaded; i++) {
        if (strcmp(all_planets[i].name, "EMBary") == 0) {
            earth_idx = i;
        } else if (strcmp(all_planets[i].name, "Mars") == 0) {
            mars_idx = i;
        }
    }

    if (earth_idx == -1) {
        fprintf(stderr, "Earth (EMBary) data not found in file.\n");
        return 1;
    }
    if (mars_idx == -1) {
        fprintf(stderr, "Mars data not found in file.\n");
        return 1;
    }

    earth_elements = all_planets[earth_idx];
    mars_elements = all_planets[mars_idx];

    /* Date: 2023-10-27 00:00:00 UT */
    year = 2023; month = 10; day = 27;
    hour = 0.0; minute = 0.0; second = 0.0;

    jd = calculate_julian_day(year, month, day, hour, minute, second);
    printf("Julian Day for %d-%02d-%02d %.2f:%.2f:%.2f UT is %.5f\n",
           year, month, day, hour, minute, second, jd);
    printf("--------------------------------------------------\n\n");

    /* Calculate Earth's heliocentric position */
    calculate_planet_helio_ecliptic_coords(earth_elements, jd, &earth_ephem);
    printf("%s Heliocentric Ecliptic Coordinates (J2000.0):\n", earth_elements.name);
    printf("  X: %.6f AU\n", earth_ephem.x_h_ecl);
    printf("  Y: %.6f AU\n", earth_ephem.y_h_ecl);
    printf("  Z: %.6f AU\n\n", earth_ephem.z_h_ecl);

    /* Calculate Mars' heliocentric position */
    calculate_planet_helio_ecliptic_coords(mars_elements, jd, &mars_ephem_helio);
    printf("%s Heliocentric Ecliptic Coordinates (J2000.0):\n", mars_elements.name);
    printf("  X: %.6f AU\n", mars_ephem_helio.x_h_ecl);
    printf("  Y: %.6f AU\n", mars_ephem_helio.y_h_ecl);
    printf("  Z: %.6f AU\n\n", mars_ephem_helio.z_h_ecl);

    /* Calculate Mars' geocentric position */
    calculate_geocentric_ecliptic_coords(mars_ephem_helio, earth_ephem, &mars_ephem_geo);
    printf("%s Geocentric Ecliptic Coordinates (J2000.0):\n", mars_elements.name);
    printf("  Rectangular X: %.6f AU\n", mars_ephem_geo.x_g_ecl);
    printf("  Rectangular Y: %.6f AU\n", mars_ephem_geo.y_g_ecl);
    printf("  Rectangular Z: %.6f AU\n", mars_ephem_geo.z_g_ecl);
    printf("  Ecliptic Longitude (lambda): %.4f deg\n", mars_ephem_geo.lambda_g_deg);
    printf("  Ecliptic Latitude (beta):  %.4f deg\n", mars_ephem_geo.beta_g_deg);
    printf("  Distance (delta):          %.6f AU\n", mars_ephem_geo.delta_g_au);

    return 0;
}
