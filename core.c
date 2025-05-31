#include <stdio.h>
#include <string.h>
#include <math.h> /* For sqrt, pow, atan2, sin, cos. Link with -lm if necessary. */

/* Define M_PI if not already defined by math.h */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* --- Constants --- */
#define G_CONST 6.67430e-11     /* Gravitational constant (m^3 kg^-1 s^-2) */
#define AU 1.495978707e11       /* Astronomical Unit (meters) */
#define DAY_S 86400.0           /* Seconds in a day */
#define NUM_BODIES 10           /* Sun + Mercury to Pluto */

/* Approximate mean radii of celestial bodies (meters) */
/* Made non-const, initialization of global arrays should still work */
double PLANET_RADII_M[NUM_BODIES] = {
    6.957e8,  /* Sun */
    2.4397e6, /* Mercury */
    6.0518e6, /* Venus */
    6.371e6,  /* Earth */
    3.3895e6, /* Mars */
    6.9911e7, /* Jupiter */
    5.8232e7, /* Saturn */
    2.5362e7, /* Uranus */
    2.4622e7, /* Neptune */
    1.1883e6  /* Pluto */
};

/* --- Vector3D Structure and Operations --- */
typedef struct {
    double x, y, z;
} Vector3D;

/* Helper function to initialize Vector3D, as compound literals are not K&R C */
Vector3D V_init(x, y, z)
    double x;
    double y;
    double z;
{
    Vector3D v;
    v.x = x; v.y = y; v.z = z;
    return v;
}

Vector3D V_add(a, b)
    Vector3D a;
    Vector3D b;
{
    Vector3D result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;
    return result;
}

Vector3D V_sub(a, b)
    Vector3D a;
    Vector3D b;
{
    Vector3D result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;
    return result;
}

Vector3D V_scale(v, s)
    Vector3D v;
    double s;
{
    Vector3D result;
    result.x = v.x * s;
    result.y = v.y * s;
    result.z = v.z * s;
    return result;
}

double V_mag_sq(v)
    Vector3D v;
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

double V_dot(a, b)
    Vector3D a;
    Vector3D b;
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/* --- Planet Structure --- */
typedef struct {
    char name[20];
    double mass;           /* kg */
    Vector3D position;     /* m */
    Vector3D velocity;     /* m/s */
    Vector3D acceleration; /* m/s^2 (current net acceleration) */
} Planet;

/* --- File Reading Function --- */
/*
 * Reads planet orbital data from a specified file.
 * File format per planet (values separated by whitespace/newlines):
 *   PlanetName (string, no spaces, max 19 chars)
 *   Mass_kg (double)
 *   PositionX_m PositionY_m PositionZ_m (3 doubles)
 *   VelocityX_mps VelocityY_mps VelocityZ_mps (3 doubles)
 *
 * Returns the number of planets successfully read.
 * Returns -1 if the file cannot be opened.
 * Returns 0 to (max_planets_to_read - 1) if an error occurs during parsing
 * or EOF is reached prematurely.
 */
int read_planet_data_from_file(filename, planets_out, max_planets_to_read)
    char filename[];
    Planet planets_out[];
    int max_planets_to_read;
{
    FILE *fp;
    int i;
    int planets_read_count;
    /* K&R C: all declarations at the top */

    fp = fopen(filename, "r");
    if (fp == NULL) {
        /* perror("Error opening planet data file"); Not available on all K&R systems */
        fprintf(stderr, "Error: Could not open planet data file '%s'\n", filename);
        return -1; /* Indicate file open error */
    }

    planets_read_count = 0;
    for (i = 0; i < max_planets_to_read; ++i) {
        /* fscanf returns the number of items successfully assigned. */
        /* We expect 1 for name, 1 for mass, 3 for position, 3 for velocity. */

        if (fscanf(fp, "%19s", planets_out[i].name) != 1) break;

        if (fscanf(fp, "%lf", &planets_out[i].mass) != 1) break;
        
        if (fscanf(fp, "%lf %lf %lf", &planets_out[i].position.x,
                                      &planets_out[i].position.y,
                                      &planets_out[i].position.z) != 3) break;

        if (fscanf(fp, "%lf %lf %lf", &planets_out[i].velocity.x,
                                      &planets_out[i].velocity.y,
                                      &planets_out[i].velocity.z) != 3) break;
        planets_read_count++;
    }

    fclose(fp);
    return planets_read_count;
}

/* --- File Writing Function --- */
/*
 * Writes planet orbital data to a specified file.
 * Overwrites the file if it exists, or creates it if it does not.
 * File format per planet (values separated by newlines):
 *   PlanetName
 *   Mass_kg
 *   PositionX_m PositionY_m PositionZ_m
 *   VelocityX_mps VelocityY_mps VelocityZ_mps
 *
 * Returns 0 on success.
 * Returns -1 if the file cannot be opened for writing.
 */
int write_planet_data_to_file(filename, planets, num_planets_to_write)
    char filename[];
    Planet planets[];
    int num_planets_to_write;
{
    FILE *fp;
    int i;

    fp = fopen(filename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: Could not open planet data file '%s' for writing.\n", filename);
        return -1; /* Indicate file open error */
    }

    for (i = 0; i < num_planets_to_write; ++i) {
        fprintf(fp, "%s\n", planets[i].name);
        fprintf(fp, "%e\n", planets[i].mass); /* Use %e for scientific notation for mass */
        fprintf(fp, "%.1f %.1f %.1f\n", planets[i].position.x, planets[i].position.y, planets[i].position.z);
        fprintf(fp, "%.2f %.2f %.2f\n", planets[i].velocity.x, planets[i].velocity.y, planets[i].velocity.z);
    }

    fclose(fp);
    return 0; /* Success */
}

/* --- Initialization --- */
/* This function sets up the initial state of the solar system. */
/* It will attempt to load data from data_filename if provided. */
/* Otherwise, it uses hardcoded simplified initial conditions. */
/* Data for specific epochs can be obtained from NASA's JPL HORIZONS system. */
void initialize_solar_system(planets, data_filename)
    Planet planets[];
    char data_filename[]; /* Can be NULL or empty string to use defaults */
{
    /* K&R C: All local variable declarations at the top of the block */
    double r_mercury, r_venus, r_earth, r_mars, r_jupiter, r_saturn, r_uranus, r_neptune, r_pluto;
    int i;
    int planets_loaded_count;

    planets_loaded_count = 0; /* Default to 0, indicating data not loaded from file yet */

    if (data_filename != NULL && data_filename[0] != '\0') {
        planets_loaded_count = read_planet_data_from_file(data_filename, planets, NUM_BODIES);
        if (planets_loaded_count == NUM_BODIES) {
            printf("Successfully loaded %d planets from %s\n", planets_loaded_count, data_filename);
        } else {
            if (planets_loaded_count == -1) { /* File open error already printed by read_planet_data_from_file */
                 fprintf(stderr, "Warning: Using default initial conditions due to file error.\n");
            } else { /* Partial read or format error */
                 fprintf(stderr, "Warning: Failed to load all planets correctly from %s (loaded %d, expected %d). Using default initial conditions.\n",
                    data_filename, planets_loaded_count, NUM_BODIES);
            }
            planets_loaded_count = 0; /* Ensure fallback to defaults */
        }
    } else {
        printf("No data file specified. Using hardcoded default initial conditions for planets.\n");
    }

    if (planets_loaded_count != NUM_BODIES) {
        printf("Initializing planets with hardcoded default values.\n");
        /* Sun (Body 0) */
        strcpy(planets[0].name, "Sun");
        planets[0].mass = 1.989e30;
        planets[0].position = V_init(0.0, 0.0, 0.0);
        planets[0].velocity = V_init(0.0, 0.0, 0.0);

        /* Mercury */
        strcpy(planets[1].name, "Mercury");
        planets[1].mass = 3.301e23;
        r_mercury = 0.387 * AU;
        planets[1].position = V_init(r_mercury, 0.0, 0.0);
        planets[1].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_mercury), 0.0);

        /* Venus */
        strcpy(planets[2].name, "Venus");
        planets[2].mass = 4.867e24;
        r_venus = 0.723 * AU;
        planets[2].position = V_init(r_venus, 0.0, 0.0);
        planets[2].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_venus), 0.0);

        /* Earth */
        strcpy(planets[3].name, "Earth");
        planets[3].mass = 5.972e24;
        r_earth = 1.000 * AU;
        planets[3].position = V_init(r_earth, 0.0, 0.0);
        planets[3].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_earth), 0.0);

        /* Mars */
        strcpy(planets[4].name, "Mars");
        planets[4].mass = 6.417e23;
        r_mars = 1.524 * AU;
        planets[4].position = V_init(r_mars, 0.0, 0.0);
        planets[4].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_mars), 0.0);

        /* Jupiter */
        strcpy(planets[5].name, "Jupiter");
        planets[5].mass = 1.898e27;
        r_jupiter = 5.203 * AU;
        planets[5].position = V_init(r_jupiter, 0.0, 0.0);
        planets[5].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_jupiter), 0.0);

        /* Saturn */
        strcpy(planets[6].name, "Saturn");
        planets[6].mass = 5.683e26;
        r_saturn = 9.537 * AU;
        planets[6].position = V_init(r_saturn, 0.0, 0.0);
        planets[6].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_saturn), 0.0);

        /* Uranus */
        strcpy(planets[7].name, "Uranus");
        planets[7].mass = 8.681e25;
        r_uranus = 19.191 * AU;
        planets[7].position = V_init(r_uranus, 0.0, 0.0);
        planets[7].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_uranus), 0.0);

        /* Neptune */
        strcpy(planets[8].name, "Neptune");
        planets[8].mass = 1.024e26;
        r_neptune = 30.069 * AU;
        planets[8].position = V_init(r_neptune, 0.0, 0.0);
        planets[8].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_neptune), 0.0);

        /* Pluto */
        strcpy(planets[9].name, "Pluto");
        planets[9].mass = 1.303e22;
        r_pluto = 39.482 * AU;
        planets[9].position = V_init(r_pluto, 0.0, 0.0);
        planets[9].velocity = V_init(0.0, sqrt(G_CONST * planets[0].mass / r_pluto), 0.0);
    }

    /* Initialize accelerations to zero for all planets, regardless of loading method */
    for (i = 0; i < NUM_BODIES; ++i) {
        planets[i].acceleration = V_init(0.0, 0.0, 0.0);
    }
}

/* --- Physics Calculations --- */

/* Softening factor to prevent extreme forces at very close distances */
/* and division by zero if planets were to occupy the same point. */
/* Epsilon is a small distance, e.g., 1km = 1000m, so epsilon^2 = 1e6 m^2. */
/* This is a common technique in N-body simulations. */
#define SOFTENING_EPSILON_SQ 1.0e6 /* (1km)^2 */

void calculate_all_accelerations(planets)
    Planet planets[];
{
    /* Local variable declarations */
    int i, j;
    Vector3D r_vec;
    double dist_sq;
    double inv_dist_cubed;
    Vector3D acc_component;

    /* Reset accelerations for all planets */
    for (i = 0; i < NUM_BODIES; ++i) {
        planets[i].acceleration = V_init(0.0, 0.0, 0.0);
    }

    /* Calculate gravitational forces and resulting accelerations */
    /* For each planet i, sum forces from all other planets j */
    for (i = 0; i < NUM_BODIES; ++i) {
        for (j = 0; j < NUM_BODIES; ++j) {
            if (i == j) continue; /* A planet does not exert force on itself */

            r_vec = V_sub(planets[j].position, planets[i].position); /* Vector from i to j */
            dist_sq = V_mag_sq(r_vec);

            /* Acceleration on planet i due to planet j: a_i = G * m_j * (r_j - r_i) / |r_j - r_i|^3 */
            /* We use (dist_sq + SOFTENING_EPSILON_SQ)^(3/2) for the denominator. */
            inv_dist_cubed = 1.0 / pow(dist_sq + SOFTENING_EPSILON_SQ, 1.5);
            
            acc_component = V_scale(r_vec, G_CONST * planets[j].mass * inv_dist_cubed);
            planets[i].acceleration = V_add(planets[i].acceleration, acc_component);
        }
    }
}

/* --- Integration (Updating Kinematics) --- */
void update_kinematics(planets, dt)
    Planet planets[];
    double dt;
{
    int i;
    for (i = 0; i < NUM_BODIES; ++i) {
        /* Using Euler-Cromer (semi-implicit Euler) method: */
        /* v_new = v_old + a_old * dt */
        /* p_new = p_old + v_new * dt */
        planets[i].velocity = V_add(planets[i].velocity, V_scale(planets[i].acceleration, dt));
        planets[i].position = V_add(planets[i].position, V_scale(planets[i].velocity, dt));
    }
}

/* --- Sky Position Calculation --- */
/**
 * Calculates and prints the apparent azimuth and altitude of all other celestial bodies
 * as seen by an observer on a specified planet.
 *
 * Assumes the observer planet's rotation axis is aligned with the global +Z axis,
 * and its prime meridian (0 deg longitude) aligns with the global +X axis.
 *
 * @param planets Array of Planet structures containing their current positions.
 * @param num_planets Total number of celestial bodies in the simulation.
 * @param observer_planet_idx Index of the planet where the observer is located.
 * @param observer_latitude_deg Observer's latitude on the planet (-90 to +90 degrees).
 * @param observer_longitude_deg Observer's longitude on the planet (0 to 360 degrees, or -180 to 180).
 * @param observer_height_m Observer's height above the planet's mean radius (meters).
 * @param planet_radii_m Array of mean radii for each planet in meters.
 */
void calculate_apparent_sky_positions(planets, num_planets,
                                      observer_planet_idx, /* Removed int type from here */
                                      observer_latitude_deg, /* Removed double type from here */
                                      observer_longitude_deg, /* Removed double type from here */
                                      observer_height_m, /* Removed double type from here */
                                      radii_m) /* radii_m was already correctly just the name */
    Planet planets[];
    int num_planets;
    int observer_planet_idx;
    double observer_latitude_deg, observer_longitude_deg, observer_height_m;
    double radii_m[];
{
    Planet* observer_planet;
    double obs_lat_rad;
    double obs_lon_rad;
    Vector3D U_local_coords;
    double total_radius;
    Vector3D observer_offset_from_center;
    Vector3D observer_global_pos;
    Vector3D E_global;
    Vector3D N_global;
    Vector3D U_global;
    int i;
    Planet* target_planet; /* Removed const */
    Vector3D vec_observer_to_target;
    double coord_E, coord_N, coord_U;
    double azimuth_rad, azimuth_deg;
    double horizontal_dist, altitude_rad, altitude_deg;

    if (observer_planet_idx < 0 || observer_planet_idx >= num_planets) {
        fprintf(stderr, "Error: Invalid observer_planet_idx.\n");
        return;
    }

    observer_planet = &planets[observer_planet_idx];
    obs_lat_rad = observer_latitude_deg * M_PI / 180.0;
    obs_lon_rad = observer_longitude_deg * M_PI / 180.0;

    /* Observer's position vector component in planet's local frame (Z=North Pole, X=Prime Meridian) */
    /* U_local_coords is a Unit vector pointing "Up" from planet center to observer */
    U_local_coords.x = cos(obs_lat_rad) * cos(obs_lon_rad);
    U_local_coords.y = cos(obs_lat_rad) * sin(obs_lon_rad);
    U_local_coords.z = sin(obs_lat_rad);

    /* Observer's actual global position */
    /* P_observer = P_planet_center + (R_planet + height) * U_local_coords */
    /* (Since we assumed planet's axes align with global, U_local_coords is also U_global_frame_orientation) */
    total_radius = radii_m[observer_planet_idx] + observer_height_m;
    observer_offset_from_center = V_scale(U_local_coords, total_radius);
    observer_global_pos = V_add(observer_planet->position, observer_offset_from_center);

    /* Define local ENU (East-North-Up) basis vectors at the observer's location */
    /* These are expressed in the global coordinate system due to our simplification. */
    /* E_global: East */
    E_global.x = -sin(obs_lon_rad);
    E_global.y = cos(obs_lon_rad);
    E_global.z = 0.0;

    /* N_global: North */
    N_global.x = -sin(obs_lat_rad) * cos(obs_lon_rad);
    N_global.y = -sin(obs_lat_rad) * sin(obs_lon_rad);
    N_global.z = cos(obs_lat_rad);

    U_global = U_local_coords; /* Up (already calculated as U_local_coords) */

    printf("Apparent sky positions from %s (Lat: %.1f, Lon: %.1f):\n",
           observer_planet->name, observer_latitude_deg, observer_longitude_deg);

    for (i = 0; i < num_planets; ++i) {
        if (i == observer_planet_idx) continue; /* Skip the observer's own planet */

        target_planet = &planets[i];
        vec_observer_to_target = V_sub(target_planet->position, observer_global_pos);

        /* Project this vector onto the local ENU axes */
        coord_E = V_dot(vec_observer_to_target, E_global);
        coord_N = V_dot(vec_observer_to_target, N_global);
        coord_U = V_dot(vec_observer_to_target, U_global);

        azimuth_rad = atan2(coord_E, coord_N);      /* Angle from North, positive towards East */
        azimuth_deg = azimuth_rad * 180.0 / M_PI;
        if (azimuth_deg < 0.0) azimuth_deg += 360.0; /* Normalize to 0-360 */

        horizontal_dist = sqrt(coord_E * coord_E + coord_N * coord_N);
        altitude_rad = atan2(coord_U, horizontal_dist); /* Angle from the EN plane */
        altitude_deg = altitude_rad * 180.0 / M_PI;

        printf("  %-10s: Az: %5.1f deg, Alt: %5.1f deg\n", target_planet->name, azimuth_deg, altitude_deg);
    }
}

/* --- Main Simulation --- */
int main() {
    /* K&R C requires all local variable declarations at the start of the block */
    Planet solar_system[NUM_BODIES];
    double time_step_seconds;
    double total_simulation_time_years;
    double output_interval_days;
    double total_simulation_seconds;
    int num_steps;
    int output_step_interval;
    int step;
    double current_time_days;
    int i;
    char planet_data_filename[] = "planets.dat"; /* Default data file name */
    /* To force use of hardcoded defaults, you could use: */
    /* char* planet_data_filename = NULL; */
    /* Or an empty string: char planet_data_filename[] = ""; */

    /* (1) Intake starting positions (and other properties) */
    initialize_solar_system(solar_system, planet_data_filename);

    /* Simulation parameters */
    time_step_seconds = 1.0 * DAY_S;       /* Time step (e.g., 1 day) */
    total_simulation_time_years = 1.0;     /* Total duration of the simulation */
    output_interval_days = 30.0;           /* How often to print output */

    total_simulation_seconds = total_simulation_time_years * 365.25 * DAY_S;
    num_steps = (int)(total_simulation_seconds / time_step_seconds);
    output_step_interval = (int)((output_interval_days * DAY_S) / time_step_seconds);
    if (output_step_interval == 0) output_step_interval = 1; /* Ensure we output at least once if interval is too small */

    printf("Starting Solar System Simulation\n");
    printf("Time Step: %.2f days (%.0f seconds)\n", time_step_seconds / DAY_S, time_step_seconds);
    printf("Total Simulation Time: %.2f years\n", total_simulation_time_years);
    printf("Outputting positions every %.2f days\n\n", output_interval_days);

    /* Simulation loop */
    for (step = 0; step <= num_steps; ++step) {
        /* (2) Update positions by integrating orbits */
        calculate_all_accelerations(solar_system);
        update_kinematics(solar_system, time_step_seconds);

        /* Output data at specified intervals */
        if (step % output_step_interval == 0) {
            current_time_days = (step * time_step_seconds) / DAY_S;
            printf("--- Time: %.2f days (Year %.3f) ---\n", current_time_days, current_time_days / 365.25);
            for (i = 0; i < NUM_BODIES; ++i) {
                printf("%-10s: Pos (AU): (%.4f, %.4f, %.4f)\tVel (km/s): (%.2f, %.2f, %.2f)\n",
                       solar_system[i].name,
                       solar_system[i].position.x / AU,
                       solar_system[i].position.y / AU,
                       solar_system[i].position.z / AU,
                       solar_system[i].velocity.x / 1000.0,
                       solar_system[i].velocity.y / 1000.0,
                       solar_system[i].velocity.z / 1000.0
                       );
            }

            /* Example: Calculate sky positions from Earth (index 3) at Greenwich (0 lon, ~51.5 lat) */
            if (NUM_BODIES > 3) { /* Ensure Earth is in the simulation */
                calculate_apparent_sky_positions(solar_system, NUM_BODIES, 3, 51.5, 0.0, 0.0, PLANET_RADII_M);
            }
            printf("\n");
        }
    }

    /* (3) Write final state back to file */
    printf("\nWriting final planet states to %s...\n", planet_data_filename);
    if (write_planet_data_to_file(planet_data_filename, solar_system, NUM_BODIES) == 0) {
        printf("Successfully wrote final states.\n");
    } /* Error message handled by write_planet_data_to_file */

    printf("Simulation finished.\n");
    return 0;
}
