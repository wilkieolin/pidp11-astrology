#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES // For M_PI on some compilers (like MSVC)
#include <math.h> // For sqrt, pow. Link with -lm if necessary.

// --- Constants ---
const double G_CONST = 6.67430e-11;     // Gravitational constant (m^3 kg^-1 s^-2)
const double AU = 1.495978707e11;       // Astronomical Unit (meters)
const double DAY_S = 86400.0;           // Seconds in a day
const int NUM_BODIES = 10;              // Sun + Mercury to Pluto

// Approximate mean radii of celestial bodies (meters)
const double PLANET_RADII_M[NUM_BODIES] = {
    6.957e8,  // Sun
    2.4397e6, // Mercury
    6.0518e6, // Venus
    6.371e6,  // Earth
    3.3895e6, // Mars
    6.9911e7, // Jupiter
    5.8232e7, // Saturn
    2.5362e7, // Uranus
    2.4622e7, // Neptune
    1.1883e6  // Pluto
};

// --- Vector3D Structure and Operations ---
typedef struct {
    double x, y, z;
} Vector3D;

Vector3D V_add(Vector3D a, Vector3D b) {
    return (Vector3D){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector3D V_sub(Vector3D a, Vector3D b) {
    return (Vector3D){a.x - b.x, a.y - b.y, a.z - b.z};
}

Vector3D V_scale(Vector3D v, double s) {
    return (Vector3D){v.x * s, v.y * s, v.z * s};
}

double V_mag_sq(Vector3D v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

double V_dot(Vector3D a, Vector3D b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// --- Planet Structure ---
typedef struct {
    char name[20];
    double mass;        // kg
    Vector3D position;  // m
    Vector3D velocity;  // m/s
    Vector3D acceleration; // m/s^2 (current net acceleration)
} Planet;

// --- Initialization ---
// This function sets up the initial state of the solar system.
// IMPORTANT: To meet the requirement of intaking starting positions for a specific date/time,
// you should modify this function to read data from a file or user input.
// Data for specific epochs can be obtained from NASA's JPL HORIZONS system.
// The current implementation uses simplified, approximate initial conditions
// (planets on the x-axis, y-velocities for ~circular orbits) for demonstration.
void initialize_solar_system(Planet planets[]) {
    // Sun (Body 0)
    strcpy(planets[0].name, "Sun");
    planets[0].mass = 1.989e30;
    planets[0].position = (Vector3D){0.0, 0.0, 0.0};
    planets[0].velocity = (Vector3D){0.0, 0.0, 0.0};

    // Mercury
    strcpy(planets[1].name, "Mercury");
    planets[1].mass = 3.301e23;
    double r_mercury = 0.387 * AU;
    planets[1].position = (Vector3D){r_mercury, 0.0, 0.0};
    planets[1].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_mercury), 0.0};

    // Venus
    strcpy(planets[2].name, "Venus");
    planets[2].mass = 4.867e24;
    double r_venus = 0.723 * AU;
    planets[2].position = (Vector3D){r_venus, 0.0, 0.0};
    planets[2].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_venus), 0.0};

    // Earth
    strcpy(planets[3].name, "Earth");
    planets[3].mass = 5.972e24;
    double r_earth = 1.000 * AU;
    planets[3].position = (Vector3D){r_earth, 0.0, 0.0};
    planets[3].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_earth), 0.0};

    // Mars
    strcpy(planets[4].name, "Mars");
    planets[4].mass = 6.417e23;
    double r_mars = 1.524 * AU;
    planets[4].position = (Vector3D){r_mars, 0.0, 0.0};
    planets[4].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_mars), 0.0};

    // Jupiter
    strcpy(planets[5].name, "Jupiter");
    planets[5].mass = 1.898e27;
    double r_jupiter = 5.203 * AU;
    planets[5].position = (Vector3D){r_jupiter, 0.0, 0.0};
    planets[5].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_jupiter), 0.0};

    // Saturn
    strcpy(planets[6].name, "Saturn");
    planets[6].mass = 5.683e26;
    double r_saturn = 9.537 * AU;
    planets[6].position = (Vector3D){r_saturn, 0.0, 0.0};
    planets[6].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_saturn), 0.0};

    // Uranus
    strcpy(planets[7].name, "Uranus");
    planets[7].mass = 8.681e25;
    double r_uranus = 19.191 * AU;
    planets[7].position = (Vector3D){r_uranus, 0.0, 0.0};
    planets[7].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_uranus), 0.0};

    // Neptune
    strcpy(planets[8].name, "Neptune");
    planets[8].mass = 1.024e26;
    double r_neptune = 30.069 * AU;
    planets[8].position = (Vector3D){r_neptune, 0.0, 0.0};
    planets[8].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_neptune), 0.0};

    // Pluto
    strcpy(planets[9].name, "Pluto");
    planets[9].mass = 1.303e22;
    double r_pluto = 39.482 * AU;
    planets[9].position = (Vector3D){r_pluto, 0.0, 0.0};
    planets[9].velocity = (Vector3D){0.0, sqrt(G_CONST * planets[0].mass / r_pluto), 0.0};

    // Initialize accelerations to zero
    for (int i = 0; i < NUM_BODIES; ++i) {
        planets[i].acceleration = (Vector3D){0.0, 0.0, 0.0};
    }
}

// --- Physics Calculations ---
void calculate_all_accelerations(Planet planets[]) {
    // Reset accelerations for all planets
    for (int i = 0; i < NUM_BODIES; ++i) {
        planets[i].acceleration = (Vector3D){0.0, 0.0, 0.0};
    }

    // Calculate gravitational forces and resulting accelerations
    // For each planet i, sum forces from all other planets j
    for (int i = 0; i < NUM_BODIES; ++i) {
        for (int j = 0; j < NUM_BODIES; ++j) {
            if (i == j) continue; // A planet does not exert force on itself

            Vector3D r_vec = V_sub(planets[j].position, planets[i].position); // Vector from i to j
            double dist_sq = V_mag_sq(r_vec);

            // Softening factor to prevent extreme forces at very close distances
            // and division by zero if planets were to occupy the same point.
            // Epsilon is a small distance, e.g., 1km = 1000m, so epsilon^2 = 1e6 m^2.
            // This is a common technique in N-body simulations.
            const double softening_epsilon_sq = 1e6; // (1km)^2

            // Acceleration on planet i due to planet j: a_i = G * m_j * (r_j - r_i) / |r_j - r_i|^3
            // We use (dist_sq + softening_epsilon_sq)^(3/2) for the denominator.
            double inv_dist_cubed = 1.0 / pow(dist_sq + softening_epsilon_sq, 1.5);
            
            Vector3D acc_component = V_scale(r_vec, G_CONST * planets[j].mass * inv_dist_cubed);
            planets[i].acceleration = V_add(planets[i].acceleration, acc_component);
        }
    }
}

// --- Integration (Updating Kinematics) ---
void update_kinematics(Planet planets[], double dt) {
    for (int i = 0; i < NUM_BODIES; ++i) {
        // Using Euler-Cromer (semi-implicit Euler) method:
        // v_new = v_old + a_old * dt
        // p_new = p_old + v_new * dt
        planets[i].velocity = V_add(planets[i].velocity, V_scale(planets[i].acceleration, dt));
        planets[i].position = V_add(planets[i].position, V_scale(planets[i].velocity, dt));
    }
}

// --- Sky Position Calculation ---
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
void calculate_apparent_sky_positions(const Planet planets[], int num_planets,
                                      int observer_planet_idx,
                                      double observer_latitude_deg,
                                      double observer_longitude_deg,
                                      double observer_height_m,
                                      const double radii_m[]) {
    if (observer_planet_idx < 0 || observer_planet_idx >= num_planets) {
        fprintf(stderr, "Error: Invalid observer_planet_idx.\n");
        return;
    }

    const Planet* observer_planet = &planets[observer_planet_idx];
    double obs_lat_rad = observer_latitude_deg * M_PI / 180.0;
    double obs_lon_rad = observer_longitude_deg * M_PI / 180.0;

    // Observer's position vector component in planet's local frame (Z=North Pole, X=Prime Meridian)
    Vector3D U_local_coords; // Unit vector pointing "Up" from planet center to observer
    U_local_coords.x = cos(obs_lat_rad) * cos(obs_lon_rad);
    U_local_coords.y = cos(obs_lat_rad) * sin(obs_lon_rad);
    U_local_coords.z = sin(obs_lat_rad);

    // Observer's actual global position
    // P_observer = P_planet_center + (R_planet + height) * U_local_coords
    // (Since we assumed planet's axes align with global, U_local_coords is also U_global_frame_orientation)
    double total_radius = radii_m[observer_planet_idx] + observer_height_m;
    Vector3D observer_offset_from_center = V_scale(U_local_coords, total_radius);
    Vector3D observer_global_pos = V_add(observer_planet->position, observer_offset_from_center);

    // Define local ENU (East-North-Up) basis vectors at the observer's location
    // These are expressed in the global coordinate system due to our simplification.
    Vector3D E_global; // East
    E_global.x = -sin(obs_lon_rad);
    E_global.y = cos(obs_lon_rad);
    E_global.z = 0.0;

    Vector3D N_global; // North
    N_global.x = -sin(obs_lat_rad) * cos(obs_lon_rad);
    N_global.y = -sin(obs_lat_rad) * sin(obs_lon_rad);
    N_global.z = cos(obs_lat_rad);

    Vector3D U_global = U_local_coords; // Up (already calculated as U_local_coords)

    printf("Apparent sky positions from %s (Lat: %.1f, Lon: %.1f):\n",
           observer_planet->name, observer_latitude_deg, observer_longitude_deg);

    for (int i = 0; i < num_planets; ++i) {
        if (i == observer_planet_idx) continue; // Skip the observer's own planet

        const Planet* target_planet = &planets[i];
        Vector3D vec_observer_to_target = V_sub(target_planet->position, observer_global_pos);

        // Project this vector onto the local ENU axes
        double coord_E = V_dot(vec_observer_to_target, E_global);
        double coord_N = V_dot(vec_observer_to_target, N_global);
        double coord_U = V_dot(vec_observer_to_target, U_global);

        double azimuth_rad = atan2(coord_E, coord_N);      // Angle from North, positive towards East
        double azimuth_deg = azimuth_rad * 180.0 / M_PI;
        if (azimuth_deg < 0.0) azimuth_deg += 360.0; // Normalize to 0-360

        double horizontal_dist = sqrt(coord_E * coord_E + coord_N * coord_N);
        double altitude_rad = atan2(coord_U, horizontal_dist); // Angle from the EN plane
        double altitude_deg = altitude_rad * 180.0 / M_PI;

        printf("  %-10s: Az: %5.1f deg, Alt: %5.1f deg\n", target_planet->name, azimuth_deg, altitude_deg);
    }
}

// --- Main Simulation ---
int main() {
    Planet solar_system[NUM_BODIES];

    // (1) Intake starting positions (and other properties)
    // The function below uses hardcoded values. Modify it to read from a file
    // or prompt the user for data corresponding to a specific date and time.
    initialize_solar_system(solar_system);

    // Simulation parameters
    double time_step_seconds = 1.0 * DAY_S;       // Time step (e.g., 1 day)
    double total_simulation_time_years = 1.0;     // Total duration of the simulation
    double output_interval_days = 30.0;           // How often to print output

    double total_simulation_seconds = total_simulation_time_years * 365.25 * DAY_S;
    int num_steps = (int)(total_simulation_seconds / time_step_seconds);
    int output_step_interval = (int)((output_interval_days * DAY_S) / time_step_seconds);
    if (output_step_interval == 0) output_step_interval = 1; // Ensure we output at least once if interval is too small

    printf("Starting Solar System Simulation\n");
    printf("Time Step: %.2f days (%.0f seconds)\n", time_step_seconds / DAY_S, time_step_seconds);
    printf("Total Simulation Time: %.2f years\n", total_simulation_time_years);
    printf("Outputting positions every %.2f days\n\n", output_interval_days);

    // Simulation loop
    for (int step = 0; step <= num_steps; ++step) {
        // (2) Update positions by integrating orbits
        calculate_all_accelerations(solar_system);
        update_kinematics(solar_system, time_step_seconds);

        // Output data at specified intervals
        if (step % output_step_interval == 0) {
            double current_time_days = (step * time_step_seconds) / DAY_S;
            printf("--- Time: %.2f days (Year %.3f) ---\n", current_time_days, current_time_days / 365.25);
            for (int i = 0; i < NUM_BODIES; ++i) {
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

            // Example: Calculate sky positions from Earth (index 3) at Greenwich (0 lon, ~51.5 lat)
            if (NUM_BODIES > 3) { // Ensure Earth is in the simulation
                calculate_apparent_sky_positions(solar_system, NUM_BODIES, 3, 51.5, 0.0, 0.0, PLANET_RADII_M);
            }
            printf("\n");
        }
    }

    printf("Simulation finished.\n");
    return 0;
}
