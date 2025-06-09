/* ephemeris.c - Approximate Ephemeris Calculator using Keplerian Elements */
/* Compatible with K&R C style */

#include <stdio.h>
#include <math.h>
#include <string.h> /* For strcpy in main example */
#include <time.h>   /* For time(), localtime() */
#include "ephemeris.h"

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

    jd = (long int)(365.25 * (year + 4716.0)) +   /* Use long int for large terms */
         (long int)(30.6001 * (month + 1.0)) +  /* Use long int for large terms */
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

/* ... (other definitions and functions) ... */

/*
 * calculate_apparent_az_alt_all_planets
 * Calculates the apparent Azimuth and Altitude of all other celestial bodies
 * as seen by an observer on a specified planet, given their heliocentric ecliptic coordinates.
 * The results are stored in the 'apparent_positions' array.
 *
 * orbital_elements_list: Array of OrbitalElements (for names).
 * all_planet_ephems: Array of PlanetEphem structures, containing current
 *                    HELIOCENTRIC ecliptic coordinates (x_h_ecl, y_h_ecl, z_h_ecl) for all bodies.
 * num_planets: Total number of celestial bodies.
 * observer_planet_idx: Index of the planet where the observer is located.
 * observer_latitude_deg: Observer's latitude on the planet (-90 to +90 degrees).
 * observer_longitude_deg: Observer's longitude on the planet (0 to 360, or -180 to 180 degrees).
 * jd: Julian Day of the observation.
 * apparent_positions: Output array to store calculated Az/Alt and names.
 */
void calculate_apparent_az_alt_all_planets(
    orbital_elements_list,
    all_planet_ephems,
    num_planets,
    observer_planet_idx,
    observer_latitude_deg,
    observer_longitude_deg,
    jd,
    apparent_positions)
    OrbitalElements orbital_elements_list[];
    PlanetEphem all_planet_ephems[]; /* Input: heliocentric ecliptic coords */
    int num_planets;
    int observer_planet_idx;
    double observer_latitude_deg;
    double observer_longitude_deg;
    double jd;
    ApparentSkyPosition apparent_positions[]; /* Output array */
{
    /* K&R C: All variable declarations at the top */
    PlanetEphem observer_helio_coords;
    PlanetEphem target_helio_coords;
    PlanetEphem target_coords_relative_to_observer; /* Stores geocentric ecliptic relative to observer */
    double obs_lat_rad;
    double epsilon_rad; /* Obliquity of the ecliptic */
    double T_cent;      /* Julian centuries since J2000 for GMST calculation */
    double gmst_deg;    /* Greenwich Mean Sidereal Time in degrees */
    double lst_deg;     /* Local Sidereal Time in degrees */ /* K&R: variables must be declared at top */
    double lst_rad;
    int target_idx;

    /* Variables for coordinate transformation for each target */
    double x_rel_ecl, y_rel_ecl, z_rel_ecl; /* Relative ecliptic rectangular */
    double x_eq, y_eq, z_eq;                /* Relative equatorial rectangular */
    double alpha_rad;                       /* Right Ascension */
    double delta_rad;                       /* Declination */
    double H_rad;                           /* Hour Angle */
    double sin_alt, alt_rad, alt_deg;
    double az_y_term, az_x_term, az_rad, az_deg;
    double cos_epsilon, sin_epsilon;
    double dist_eq_proj; /* Projection of distance onto equatorial plane */
    double cos_obs_lat, sin_obs_lat;
    double cos_delta, sin_delta;
    double cos_H, sin_H;

    if (observer_planet_idx < 0 || observer_planet_idx >= num_planets) {
        fprintf(stderr, "Error (calculate_apparent_az_alt): Invalid observer_planet_idx.\n");
        return;
    }

    /* Initialize output array */
    for (target_idx = 0; target_idx < num_planets; ++target_idx) {
        apparent_positions[target_idx].is_valid = FALSE;
        strcpy(apparent_positions[target_idx].name, "");
        apparent_positions[target_idx].azimuth_deg = 0.0;
        apparent_positions[target_idx].altitude_deg = 0.0;
    }


    observer_helio_coords = all_planet_ephems[observer_planet_idx];
    obs_lat_rad = observer_latitude_deg * DEG_TO_RAD;
    epsilon_rad = EPSILON_J2000_DEG * DEG_TO_RAD;

    cos_epsilon = cos(epsilon_rad);
    sin_epsilon = sin(epsilon_rad);
    cos_obs_lat = cos(obs_lat_rad);
    sin_obs_lat = sin(obs_lat_rad);

    /* Calculate Local Sidereal Time (LST) */
    T_cent = (jd - JD_J2000) / DAYS_PER_CENTURY;
    gmst_deg = 280.46061837 + 360.98564736629 * (jd - JD_J2000) +
               0.000387933 * T_cent * T_cent - (T_cent * T_cent * T_cent) / 38710000.0;
    gmst_deg = normalize_angle_deg(gmst_deg);
    lst_deg = normalize_angle_deg(gmst_deg + observer_longitude_deg);
    lst_rad = lst_deg * DEG_TO_RAD;

    for (target_idx = 0; target_idx < num_planets; ++target_idx) {
        if (target_idx == observer_planet_idx) {
            apparent_positions[target_idx].is_valid = FALSE; /* Mark observer as invalid for Az/Alt to self */
            continue; /* Skip observer's own planet */
        }

        target_helio_coords = all_planet_ephems[target_idx];

        /* Get geocentric ecliptic coordinates of target relative to observer */
        /* calculate_geocentric_ecliptic_coords expects first two args by value */
        calculate_geocentric_ecliptic_coords(target_helio_coords, observer_helio_coords,
                                             &target_coords_relative_to_observer);

        x_rel_ecl = target_coords_relative_to_observer.x_g_ecl;
        y_rel_ecl = target_coords_relative_to_observer.y_g_ecl;
        z_rel_ecl = target_coords_relative_to_observer.z_g_ecl;

        /* 1. Transform geocentric ecliptic (relative) to geocentric equatorial (relative) */
        x_eq = x_rel_ecl;
        y_eq = y_rel_ecl * cos_epsilon - z_rel_ecl * sin_epsilon;
        z_eq = y_rel_ecl * sin_epsilon + z_rel_ecl * cos_epsilon;

        /* 2. Calculate Right Ascension (alpha) and Declination (delta) */
        alpha_rad = atan2(y_eq, x_eq);
        dist_eq_proj = sqrt(x_eq * x_eq + y_eq * y_eq);
        delta_rad = atan2(z_eq, dist_eq_proj); /* More robust than asin for delta */

        /* 3. Calculate Hour Angle (H) */
        H_rad = lst_rad - alpha_rad;

        /* Pre-calculate sines and cosines for efficiency and clarity */
        cos_delta = cos(delta_rad);
        sin_delta = sin(delta_rad);
        cos_H = cos(H_rad);
        sin_H = sin(H_rad);

        /* 4. Transform equatorial (RA/Dec/H) to horizontal (Azimuth/Altitude) */
        /* Altitude calculation */
        sin_alt = sin_obs_lat * sin_delta + cos_obs_lat * cos_delta * cos_H;
        alt_rad = asin(sin_alt); /* asin range is -PI/2 to PI/2, correct for altitude */
        alt_deg = alt_rad * RAD_TO_DEG;

        /* Azimuth calculation (Az=0 North, positive Eastward) */
        /* Az = atan2(-cos(Dec)sin(H), sin(Dec)cos(Lat) - cos(Dec)sin(Lat)cos(H)) */
        az_y_term = -cos_delta * sin_H;
        az_x_term = sin_delta * cos_obs_lat - cos_delta * sin_obs_lat * cos_H;
        az_rad = atan2(az_y_term, az_x_term);
        az_deg = normalize_angle_deg(az_rad * RAD_TO_DEG);

        /* Store results in the output array */
        strcpy(apparent_positions[target_idx].name, orbital_elements_list[target_idx].name);
        apparent_positions[target_idx].azimuth_deg = az_deg;
        apparent_positions[target_idx].altitude_deg = alt_deg;
        apparent_positions[target_idx].is_valid = TRUE;
    }
}

/*
 * print_relative_angular_separations
 * Calculates and prints the angular separation between each pair of celestial bodies
 * based on their apparent Azimuth and Altitude from an observer's perspective.
 *
 * apparent_positions: Array of ApparentSkyPosition structures, populated by
 *                     calculate_apparent_az_alt_all_planets.
 * num_planets: Total number of celestial bodies.
 * orbital_elements_list: Array of OrbitalElements (for observer's name).
 * observer_planet_idx: Index of the planet where the observer is located.
 */
void print_relative_angular_separations(
    apparent_positions,
    num_planets,
    orbital_elements_list,
    observer_planet_idx)
    ApparentSkyPosition apparent_positions[];
    int num_planets;
    OrbitalElements orbital_elements_list[];
    int observer_planet_idx;
{
    int i, j;
    double az1_rad, alt1_rad, az2_rad, alt2_rad;
    double cos_angle, angle_rad, angle_deg;

    if (observer_planet_idx < 0 || observer_planet_idx >= num_planets) {
        fprintf(stderr, "Error (print_relative_angular_separations): Invalid observer_planet_idx.\n");
        return;
    }

    printf("\nRelative angular separations as seen from %s:\n", orbital_elements_list[observer_planet_idx].name);
    printf("  Planet 1      Planet 2      Separation (deg)\n");
    printf("  ------------- ------------- ------------------\n");

    for (i = 0; i < num_planets; ++i) {
        if (!apparent_positions[i].is_valid) continue;

        for (j = i + 1; j < num_planets; ++j) {
            if (!apparent_positions[j].is_valid) continue;

            alt1_rad = apparent_positions[i].altitude_deg * DEG_TO_RAD;
            az1_rad  = apparent_positions[i].azimuth_deg  * DEG_TO_RAD;
            alt2_rad = apparent_positions[j].altitude_deg * DEG_TO_RAD;
            az2_rad  = apparent_positions[j].azimuth_deg  * DEG_TO_RAD;

            cos_angle = sin(alt1_rad) * sin(alt2_rad) +
                        cos(alt1_rad) * cos(alt2_rad) * cos(az1_rad - az2_rad);

            /* Clamp cos_angle to [-1, 1] to avoid domain errors with acos due to precision issues */
            if (cos_angle > 1.0) cos_angle = 1.0;
            if (cos_angle < -1.0) cos_angle = -1.0;

            angle_rad = acos(cos_angle);
            angle_deg = angle_rad * RAD_TO_DEG;

            printf("  %-13s %-13s %6.2f\n",
                   apparent_positions[i].name, apparent_positions[j].name, angle_deg);
        }
    }
}

/*
 * get_relative_angular_separations
 * Calculates the angular separation between each pair of celestial bodies
 * based on their apparent Azimuth and Altitude from an observer's perspective.
 * The results are stored in the 'separations_out' array.
 *
 * apparent_positions: Array of ApparentSkyPosition structures.
 * num_planets: Total number of celestial bodies.
 * separations_out: Array to store the calculated angular separations in degrees.
 * max_separations: The maximum number of separations that can be stored in separations_out.
 * Returns the number of separations calculated and stored, or -1 on error.
 */
int get_relative_angular_separations(
    apparent_positions,
    num_planets,
    separations_out,
    max_separations)
    ApparentSkyPosition apparent_positions[];
    int num_planets;
    double separations_out[];
    int max_separations;
{
    /* K&R C: All variable declarations at the top */
    int i, j;
    double az1_rad, alt1_rad, az2_rad, alt2_rad;
    double cos_angle, angle_rad, angle_deg;
    int separations_count;
    int valid_planets_for_pairing[MAX_PLANETS]; /* Store indices of valid planets */
    int num_valid_planets;
    int p1_idx, p2_idx;

    separations_count = 0;
    num_valid_planets = 0;

    /* First, collect indices of all valid planets for pairing */
    for (i = 0; i < num_planets; ++i) {
        if (apparent_positions[i].is_valid) {
            if (num_valid_planets < MAX_PLANETS) {
                valid_planets_for_pairing[num_valid_planets++] = i;
            } else {
                fprintf(stderr, "Error (get_relative_angular_separations): Exceeded internal buffer for valid planets.\n");
                return -1;
            }
        }
    }

    if (num_valid_planets < 2) {
        return 0; /* Not enough valid planets to form a pair */
    }

    for (i = 0; i < num_valid_planets; ++i) {
        p1_idx = valid_planets_for_pairing[i];
        for (j = i + 1; j < num_valid_planets; ++j) {
            p2_idx = valid_planets_for_pairing[j];
            if (separations_count >= max_separations) {
                fprintf(stderr, "Error (get_relative_angular_separations): separations_out array too small.\n");
                return -1; /* Output array full */
            }
            alt1_rad = apparent_positions[p1_idx].altitude_deg * DEG_TO_RAD;
            az1_rad  = apparent_positions[p1_idx].azimuth_deg  * DEG_TO_RAD;
            alt2_rad = apparent_positions[p2_idx].altitude_deg * DEG_TO_RAD;
            az2_rad  = apparent_positions[p2_idx].azimuth_deg  * DEG_TO_RAD;

            cos_angle = sin(alt1_rad) * sin(alt2_rad) +
                        cos(alt1_rad) * cos(alt2_rad) * cos(az1_rad - az2_rad);

            if (cos_angle > 1.0) cos_angle = 1.0;
            if (cos_angle < -1.0) cos_angle = -1.0;

            angle_rad = acos(cos_angle);
            angle_deg = angle_rad * RAD_TO_DEG;

            separations_out[separations_count++] = angle_deg;
        }
    }
    return separations_count;
}

/*
 * print_helio_ecliptic_angular_separations
 * Calculates and prints the angular separation on the ecliptic plane between
 * each pair of planets, as viewed from the Sun (heliocentric perspective).
 *
 * orbital_elements_list: Array of OrbitalElements (for names).
 * all_planet_ephems: Array of PlanetEphem structures, containing current
 *                    HELIOCENTRIC ecliptic coordinates (x_h_ecl, y_h_ecl) for all bodies.
 * num_planets: Total number of celestial bodies.
 */
void print_helio_ecliptic_angular_separations(
    orbital_elements_list,
    all_planet_ephems,
    num_planets)
    OrbitalElements orbital_elements_list[];
    PlanetEphem all_planet_ephems[];
    int num_planets;
{
    /* K&R C: All variable declarations at the top */
    double helio_ecl_lon_deg[MAX_PLANETS]; /* Store heliocentric ecliptic longitude for each planet */
    int i, j;
    double lon_rad;
    double angle_diff_deg;

    if (num_planets <= 1) {
        printf("\nNot enough planets to calculate relative heliocentric ecliptic angles.\n");
        return;
    }

    printf("\nHeliocentric Ecliptic Angular Separations (view from Sun):\n");
    printf("  Planet 1      Planet 2      Separation (deg)\n");
    printf("  ------------- ------------- ------------------\n");

    /* Calculate heliocentric ecliptic longitude for each planet */
    for (i = 0; i < num_planets; ++i) {
        lon_rad = atan2(all_planet_ephems[i].y_h_ecl, all_planet_ephems[i].x_h_ecl);
        helio_ecl_lon_deg[i] = normalize_angle_deg(lon_rad * RAD_TO_DEG);
    }

    /* Calculate and print angular separation between unique pairs */
    for (i = 0; i < num_planets; ++i) {
        for (j = i + 1; j < num_planets; ++j) {
            angle_diff_deg = fabs(helio_ecl_lon_deg[j] - helio_ecl_lon_deg[i]); /* abs difference */
            if (angle_diff_deg > 180.0) { /* Get shortest angle (0-180) */
                angle_diff_deg = 360.0 - angle_diff_deg;
            }
            printf("  %-13s %-13s %6.2f\n",
                   orbital_elements_list[i].name,
                   orbital_elements_list[j].name,
                   angle_diff_deg);
        }
    }
}

/*
 * get_helio_ecliptic_angular_separations
 * Calculates the angular separation on the ecliptic plane between
 * each pair of planets, as viewed from the Sun (heliocentric perspective).
 * The results are stored in the 'separations_out' array.
 *
 * all_planet_ephems: Array of PlanetEphem structures (heliocentric ecliptic coords).
 * num_planets: Total number of celestial bodies.
 * separations_out: Array to store the calculated angular separations in degrees.
 * max_separations: The maximum number of separations that can be stored in separations_out.
 * Returns the number of separations calculated and stored, or -1 on error.
 */
int get_helio_ecliptic_angular_separations(
    all_planet_ephems,
    num_planets,
    separations_out,
    max_separations)
    PlanetEphem all_planet_ephems[];
    int num_planets;
    double separations_out[];
    int max_separations;
{
    /* K&R C: All variable declarations at the top */
    double helio_ecl_lon_deg[MAX_PLANETS];
    int i, j;
    double lon_rad;
    double angle_diff_deg;
    int separations_count;

    separations_count = 0;

    if (num_planets <= 1) {
        return 0; /* Not enough planets to calculate separations */
    }

    /* Calculate heliocentric ecliptic longitude for each planet */
    for (i = 0; i < num_planets; ++i) {
        lon_rad = atan2(all_planet_ephems[i].y_h_ecl, all_planet_ephems[i].x_h_ecl);
        helio_ecl_lon_deg[i] = normalize_angle_deg(lon_rad * RAD_TO_DEG);
    }

    /* Calculate and store angular separation between unique pairs */
    for (i = 0; i < num_planets; ++i) {
        for (j = i + 1; j < num_planets; ++j) {
            if (separations_count >= max_separations) {
                fprintf(stderr, "Error (get_helio_ecliptic_angular_separations): separations_out array too small.\n");
                return -1; /* Output array full */
            }
            angle_diff_deg = fabs(helio_ecl_lon_deg[j] - helio_ecl_lon_deg[i]);
            if (angle_diff_deg > 180.0) {
                angle_diff_deg = 360.0 - angle_diff_deg;
            }
            separations_out[separations_count++] = angle_diff_deg;
        }
    }
    return separations_count;
}

/* Example Usage */

/* Define EPHEMERIS_MAIN_TEST to compile this main function for testing */
#ifdef EPHEMERIS_MAIN_TEST
int main() {
    /* K&R C: All declarations at the top of the block */
    OrbitalElements all_planets[MAX_PLANETS];
    PlanetEphem all_planet_ephems_data[MAX_PLANETS]; /* Store all heliocentric ephems */
    ApparentSkyPosition apparent_sky_pos[MAX_PLANETS]; /* Store Az/Alt for observer */
    int num_planets_loaded;
    int i;

    double jd;
    int year, month, day;
    double hour, minute, second;

    time_t time_now;
    struct tm *local_time_now;

    /* Observer details for Az/Alt calculation (example: Greenwich for Earth) */
    int earth_observer_idx = 2; /* EMBary is 3rd in ephemeris_data.txt (index 2) */
    double observer_lat = 41.77810;  /* Latitude of Naperville, IL (example) */
    double observer_lon = -88.08260; /* Longitude of Naperville, IL (example) */
    
    /* For testing new get_*_separations functions */
    double separations_buffer[MAX_PLANETS * (MAX_PLANETS - 1) / 2]; /* Max possible separations */
    int num_seps_calculated;
    int k; /* Loop counter for printing separations from buffer */
    
    /* Initialize orbital elements from file */
    if (!initialize_orbital_elements_from_file("ephemeris_data.txt", all_planets, &num_planets_loaded)) {
        fprintf(stderr, "Failed to load ephemeris data.\n");
        return 1;
    }
    printf("Loaded %d planets from file.\n\n", num_planets_loaded);

    /* Get current system time */
    time_now = time(NULL);
    local_time_now = localtime(&time_now);

    /* Extract date and time components */
    /* struct tm members: tm_year (years since 1900), tm_mon (0-11), tm_mday (1-31) */
    /* tm_hour (0-23), tm_min (0-59), tm_sec (0-59) */
    year   = local_time_now->tm_year + 1900;
    month  = local_time_now->tm_mon + 1;
    day    = local_time_now->tm_mday;
    hour   = (double)local_time_now->tm_hour;
    minute = (double)local_time_now->tm_min;
    second = (double)local_time_now->tm_sec;

    jd = calculate_julian_day(year, month, day, hour, minute, second);

    printf("Calculating for current date and time: %d-%02d-%02d %02.0f:%02.0f:%02.0f UT\n",
           year, month, day, hour, minute, second);
    printf("Julian Day: %.5f\n", jd);
    printf("-----------------------------------------------------------------------\n");
    printf("Heliocentric Ecliptic Coordinates (J2000.0):\n");

    /* Loop through each loaded planet and calculate its heliocentric position */
    for (i = 0; i < num_planets_loaded; i++) {
        calculate_planet_helio_ecliptic_coords(all_planets[i], jd, &all_planet_ephems_data[i]);
        printf("  %-10s: X: %9.6f AU, Y: %9.6f AU, Z: %9.6f AU\n",
               all_planets[i].name,
               all_planet_ephems_data[i].x_h_ecl,
               all_planet_ephems_data[i].y_h_ecl,
               all_planet_ephems_data[i].z_h_ecl);
    }

    /*
     * Example: Calculate and print apparent Azimuth/Altitude for all planets
     * as seen from Earth (assuming Earth is at index earth_observer_idx).
     */
    if (earth_observer_idx >= 0 && earth_observer_idx < num_planets_loaded) {
        calculate_apparent_az_alt_all_planets(all_planets, all_planet_ephems_data,
                                              num_planets_loaded, earth_observer_idx,
                                              observer_lat, observer_lon, jd,
                                              apparent_sky_pos);

        /* Optional: Print the calculated Az/Alt values */
        printf("\nApparent Azimuth/Altitude from %s (Lat: %.2f, Lon: %.2f) at JD %.5f:\n",
               all_planets[earth_observer_idx].name,
               observer_lat, observer_lon, jd);
        printf("  Target        Azimuth (N=0,E=90)   Altitude\n");
        printf("  ------------- -------------------- --------------------\n");
        for (i = 0; i < num_planets_loaded; i++) {
            if (apparent_sky_pos[i].is_valid) {
                printf("  %-13s Az: %6.2f deg, Alt: %6.2f deg\n",
                       apparent_sky_pos[i].name,
                       apparent_sky_pos[i].azimuth_deg,
                       apparent_sky_pos[i].altitude_deg);
            }
        }

        print_relative_angular_separations(apparent_sky_pos, num_planets_loaded,
                                           all_planets, earth_observer_idx);

        /* Test get_relative_angular_separations */
        num_seps_calculated = get_relative_angular_separations(apparent_sky_pos, num_planets_loaded,
                                                               separations_buffer,
                                                               sizeof(separations_buffer)/sizeof(double));
        if (num_seps_calculated >= 0) {
            printf("\nGot %d relative angular separations (via get_relative_angular_separations):\n", num_seps_calculated);
            for (k = 0; k < num_seps_calculated; k++) {
                printf("  Sep #%d: %.2f deg\n", k + 1, separations_buffer[k]);
            }
        } else {
            printf("\nError getting relative angular separations.\n");
        }
    }

    /* Calculate and print heliocentric ecliptic angular separations */
    print_helio_ecliptic_angular_separations(all_planets, all_planet_ephems_data,
                                             num_planets_loaded);

    /* Test get_helio_ecliptic_angular_separations */
    num_seps_calculated = get_helio_ecliptic_angular_separations(all_planet_ephems_data, num_planets_loaded,
                                                                 separations_buffer,
                                                                 sizeof(separations_buffer)/sizeof(double));
    if (num_seps_calculated >= 0) {
        printf("\nGot %d heliocentric ecliptic angular separations (via get_helio_ecliptic_angular_separations):\n", num_seps_calculated);
        for (k = 0; k < num_seps_calculated; k++) {
            printf("  Sep #%d: %.2f deg\n", k + 1, separations_buffer[k]);
        }
    } else {
        printf("\nError getting heliocentric ecliptic angular separations.\n");
    }
    return 0;
}
#endif /* EPHEMERIS_MAIN_TEST */