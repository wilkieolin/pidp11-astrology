/* ephemeris.h - Public interface for Approximate Ephemeris Calculator */
/* Compatible with K&R C style */

#ifndef EPHEMERIS_H
#define EPHEMERIS_H

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

/* Added constant for J2000.0 obliquity of the ecliptic */
#define EPSILON_J2000_DEG 23.4392911 /* Obliquity of the ecliptic in degrees at J2000.0 */

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

/* Structure to hold apparent Azimuth and Altitude */
typedef struct {
    char name[32];       /* Name of the target body */
    double azimuth_deg;  /* Calculated Azimuth in degrees */
    double altitude_deg; /* Calculated Altitude in degrees */
    int is_valid;        /* TRUE if this entry contains valid Az/Alt data */
} ApparentSkyPosition;


/* Function Prototypes (K&R C style) */
double calculate_julian_day( /* int year, int month, int day, double hour, double minute, double second */ );
double my_fmod( /* double x, double y */ ); /* K&R C compatible fmod */
double normalize_angle_deg( /* double angle_deg */ );
double solve_kepler_equation( /* double M_rad, double e */ );
void calculate_planet_helio_ecliptic_coords( /* OrbitalElements elements, double jd, PlanetEphem *ephem */ );
void calculate_geocentric_ecliptic_coords( /* PlanetEphem planet_helio, PlanetEphem earth_helio, PlanetEphem *planet_geo */ );
int initialize_orbital_elements_from_file( /* char *filename, OrbitalElements elements[], int *num_planets */ );
void calculate_apparent_az_alt_all_planets(
    /* OrbitalElements orbital_elements_list[], */
    /* PlanetEphem all_planet_ephems[], */
    /* int num_planets, */
    /* int observer_planet_idx, */
    /* double observer_latitude_deg, */
    /* double observer_longitude_deg, */
    /* double jd, */
    /* ApparentSkyPosition apparent_positions[] */);
void print_relative_angular_separations(
    /* ApparentSkyPosition apparent_positions[], */
    /* int num_planets, */
    /* OrbitalElements orbital_elements_list[], */
    /* int observer_planet_idx */);
void print_helio_ecliptic_angular_separations(
    /* OrbitalElements orbital_elements_list[], */
    /* PlanetEphem all_planet_ephems[], */
    /* int num_planets */); /* Corrected parameter comment */

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
int get_relative_angular_separations( /* ApparentSkyPosition apparent_positions[], int num_planets, double separations_out[], int max_separations */ );

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
int get_helio_ecliptic_angular_separations( /* PlanetEphem all_planet_ephems[], int num_planets, double separations_out[], int max_separations */ );
#endif /* EPHEMERIS_H */