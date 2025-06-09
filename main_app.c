/* main_app.c - Main astrology application */
#include <stdio.h>
#include <math.h>           /* For cos, sin, atan2 */
#include <time.h>           /* For time(), localtime() if not using user input */
#include "ephemeris.h"      /* For ephemeris calculations */
#include <string.h>        /* For strlen */
#include "aphorism_utils.h" /* For aphorism generation */

/* Substitute for DBL_MAX from <float.h> if not available */
/* Used for initializing a variable to a very large value before finding a minimum. */
/* 1.0e38 is a common large double literal. */
#define DBL_MAX_SUBSTITUTE 1.0e38

/*
 * generate_astrological_seeds
 * Produces 12 pseudo-random double-precision values based on an input vector.
 * based on an input vector of angles (in degrees).
 * For each of the 12 seeds:
 *  A 52-element slice of the input_vector is used (sliding window).
 *  A point on a unit circle starts at an angle of PI/2 (0,1).
 *  Each of the 52 angles from the slice (converted to radians) is successively
 *  added to the point's current angle on the circle. The angle is normalized
 *  after each addition to the range [0, 2*PI).
 *  After all 52 angles in the slice are processed, the final (x,y) coordinates
 *  of the point are determined from its final accumulated angle.
 *  The astrological seed is atan2(final_y, final_x), normalized to the range [0.0, 1.0].
 *
 * input_vector: Array of double-precision angles in degrees.
 *               Expected to have at least 63 elements for 12 seeds with 52-angle slices.
 * output_seeds: Array to store the 12 generated double-precision seeds (range [0,1]).
 * Returns 1 on success, 0 on failure (e.g., NULL pointers).
 */

/* Local helper to normalize an angle in radians to [0, 2*PI) */
/* Uses my_fmod from ephemeris.c (linked) and PI from ephemeris.h */
static double normalize_angle_rad_local(angle_rad)
    double angle_rad;
{
    double two_pi;
    two_pi = 2.0 * PI; /* PI is from ephemeris.h */

    angle_rad = my_fmod(angle_rad, two_pi); /* my_fmod is from ephemeris.c */
    if (angle_rad < 0.0) {
        angle_rad += two_pi;
    }
    return angle_rad;
}

int generate_astrological_seeds(input_vector, output_seeds)
    double input_vector[]; /* Expected to have at least 63 elements for this logic */
    double output_seeds[]; /* Expected size 12 */
{
    /* K&R C: All variable declarations at the top */
    int seed_idx; /* For iterating 0 to 11 (12 seeds) */
    int angle_in_slice_idx; /* For iterating 0 to 51 (52 angles per slice) */
    double current_point_angle_rad; /* Current angle of the point on the unit circle */
    double input_angle_deg;         /* Angle from input_vector for the current step */
    double shot_angle_relative_rad; /* Angle of "shot" relative to the tangent */
    double tangent_angle_rad;       /* Angle of the tangent at current_point_angle_rad */
    double new_point_angle_rad;     /* Resulting angle on circle after a "shot" */
    double final_x, final_y;
    double raw_seed_angle_rad;

    /* Constants for the algorithm */
    int NUM_SEEDS_TO_PRODUCE = 12;
    int SLICE_LENGTH = 52;
    /* PI and DEG_TO_RAD are defined in ephemeris.h */

    if (input_vector == NULL || output_seeds == NULL) {
        fprintf(stderr, "Error (generate_astrological_seeds): NULL pointer passed.\n");
        return 0; /* Failure */
    }

    for (seed_idx = 0; seed_idx < NUM_SEEDS_TO_PRODUCE; ++seed_idx) {
        current_point_angle_rad = PI / 2.0; /* Initial point (0, 1.0) on unit circle */

        /* Each seed uses a slice of SLICE_LENGTH (52) from input_vector.
         * The slice for seed_idx starts at input_vector[seed_idx].
         */
        for (angle_in_slice_idx = 0; angle_in_slice_idx < SLICE_LENGTH; ++angle_in_slice_idx) {
            /*
             * Get angle from the current slice of input_vector (in degrees).
             * This angle determines the "shot" relative to the tangent.
             */
            input_angle_deg = input_vector[seed_idx + angle_in_slice_idx];
            shot_angle_relative_rad = input_angle_deg * DEG_TO_RAD;

            /* Calculate tangent direction (angle) from current point on circle */
            tangent_angle_rad = current_point_angle_rad + (PI / 2.0);

            /* Calculate the new point's angle on the circle after the "shot" */
            new_point_angle_rad = tangent_angle_rad + shot_angle_relative_rad;
            current_point_angle_rad = normalize_angle_rad_local(new_point_angle_rad);
        }
        /* After processing all 52 angles in the slice, calculate final x, y */
        final_x = cos(current_point_angle_rad);
        final_y = sin(current_point_angle_rad);
        raw_seed_angle_rad = atan2(final_y, final_x); /* Range: -PI to PI */

        /* Normalize to [0.0, 1.0] */
        output_seeds[seed_idx] = (raw_seed_angle_rad + PI) / (2.0 * PI);
    }

    return 1; /* Success */
}

/*
 * select_templates_for_signs
 * Selects an aphorism template for each of the 12 astrological signs/seeds.
 * For each seed, it identifies the template (from the total available pool)
 * whose original index is closest to (seed_value * num_total_templates).
 * Templates are selected without replacement.
 *
 * astrological_seeds: Array of 12 double values (0.0 to 1.0).
 * num_total_templates: The total number of unique aphorism templates available.
 * selected_template_indices_out: Array to store the original indices of the
 *                                selected template for each of the 12 seeds.
 *                                If a template cannot be assigned (e.g., fewer
 *                                than 12 templates available), -1 is stored.
 * Returns the number of seeds for which a template was successfully assigned.
 */
int select_templates_for_signs(astrological_seeds, num_total_templates, selected_template_indices_out)
    double astrological_seeds[]; /* Size 12 */
    int num_total_templates;
    int selected_template_indices_out[]; /* Size 12 */
{
    /* K&R C: All variable declarations at the top */
    int is_template_used[MAX_TEMPLATES]; /* From aphorism_utils.h */
    int i, j;
    double ideal_target_float;
    int best_original_idx;
    double min_abs_diff;
    double current_abs_diff;
    int num_assigned;

    if (num_total_templates > MAX_TEMPLATES) {
        fprintf(stderr, "Error (select_templates_for_signs): num_total_templates (%d) exceeds MAX_TEMPLATES (%d).\n",
                num_total_templates, MAX_TEMPLATES);
        for (i = 0; i < 12; ++i) {
            selected_template_indices_out[i] = -1;
        }
        return 0;
    }
    
    for (j = 0; j < num_total_templates; ++j) {
        is_template_used[j] = APH_FALSE; /* Use APH_FALSE from aphorism_utils.h */
    }
    for (i = 0; i < 12; ++i) { /* Initialize output array */
        selected_template_indices_out[i] = -1;
    }

    num_assigned = 0;
    for (i = 0; i < 12; ++i) { /* For each of the 12 seeds/signs */
        ideal_target_float = astrological_seeds[i] * (double)num_total_templates;
        
        best_original_idx = -1;
        min_abs_diff = DBL_MAX_SUBSTITUTE;

        for (j = 0; j < num_total_templates; ++j) {
            if (!is_template_used[j]) {
                current_abs_diff = fabs((double)j - ideal_target_float);
                if (current_abs_diff < min_abs_diff) {
                    min_abs_diff = current_abs_diff;
                    best_original_idx = j;
                }
            }
        }

        if (best_original_idx != -1) {
            selected_template_indices_out[i] = best_original_idx;
            is_template_used[best_original_idx] = APH_TRUE; /* Use APH_TRUE */
            num_assigned++;
        }
    }
    return num_assigned;
}

/* K&R C style main function */
int main(argc, argv)
    int argc;
    char *argv[];
{
    char aphorism_templates_storage[MAX_TEMPLATES][MAX_TEMPLATE_LEN];
    int num_tpls_read;
    int n_count, v_count, nv_count;
    char filled_aphorism[MAX_TEMPLATE_LEN];    
    static char *test_nouns[] = {"stars", "destiny", NULL}; /* Example for aphorism, K&R: static for init */
    static char *test_verbs[] = {"align", NULL};           /* K&R: static for init */
    static char *test_nv[] = {"cosmic insight", NULL};    /* K&R: static for init */

    /* Example: Using OrbitalElements from ephemeris.h */
    OrbitalElements all_planets[MAX_PLANETS];
    PlanetEphem all_planet_ephems_data[MAX_PLANETS]; /* For heliocentric data */
    ApparentSkyPosition apparent_sky_pos[MAX_PLANETS]; /* For Az/Alt data */
    int num_planets_loaded;
    double jd;
    int year, month, day;
    double hour, minute, second;
#ifndef USE_USER_INPUT_DATETIME
    time_t time_now;
    struct tm *local_time_now;
#endif
    int i; /* Loop counter */

    /* Observer details for apparent sky calculations */
    int observer_planet_idx = 2; /* Default to Earth (EMBary, index 2 in ephemeris_data.txt) */
    double observer_latitude_deg = 41.77810;  /* Example: Naperville, IL latitude */
    double observer_longitude_deg = -88.08260; /* Example: Naperville, IL longitude */

    /* Vector for combined angular separations */
    double combined_angular_separations[MAX_PLANETS * (MAX_PLANETS - 1)]; /* Max possible from both functions */
    int num_relative_seps_actual = 0;
    int num_helio_seps_actual = 0;
    int total_seps_in_vector = 0;
    int current_buffer_offset = 0;
    int max_seps_per_type = MAX_PLANETS * (MAX_PLANETS - 1) / 2; /* Max seps one function type can produce */
    int k; /* Loop counter for printing separations */

    double astrological_seeds[12]; /* For the new function's output */
    int m; /* Loop counter for printing seeds */
    int selected_aphorism_indices[12]; /* For storing indices of selected templates */
    int num_aphorisms_selected_for_signs;

    printf("Welcome to the K&R Astrology Program!\n");

    /* --- Initialize and use ephemeris functions --- */
    printf("\nInitializing Ephemeris Data...\n");
    if (!initialize_orbital_elements_from_file("ephemeris_data.txt", all_planets, &num_planets_loaded)) {
        fprintf(stderr, "Failed to load ephemeris data. Exiting.\n");
        return 1;
    }
    printf("Loaded %d celestial bodies.\n", num_planets_loaded);

#ifdef USE_USER_INPUT_DATETIME
    printf("\nEnter date and time for ephemeris calculation:\n");
    printf("Year (e.g., 2023): ");
    scanf("%d", &year);
    printf("Month (1-12): ");
    scanf("%d", &month);
    printf("Day (1-31): ");
    scanf("%d", &day);
    printf("Hour (0-23, UTC): ");
    scanf("%lf", &hour);
    printf("Minute (0-59): ");
    scanf("%lf", &minute);
    printf("Second (0-59): ");
    scanf("%lf", &second);
    /* Basic validation, can be expanded */
    if (month < 1 || month > 12 || day < 1 || day > 31 ||
        hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0 || second > 59) {
        fprintf(stderr, "Invalid date/time input. Using current system time instead.\n");
        time_now = time(NULL);
        local_time_now = localtime(&time_now);
        year   = local_time_now->tm_year + 1900;
        month  = local_time_now->tm_mon + 1;
        day    = local_time_now->tm_mday;
        hour   = (double)local_time_now->tm_hour;
        minute = (double)local_time_now->tm_min;
        second = (double)local_time_now->tm_sec;
    }
#else
    /* Default: Use current system time */
    time_now = time(NULL);
    local_time_now = localtime(&time_now);
    year   = local_time_now->tm_year + 1900;
    month  = local_time_now->tm_mon + 1;
    day    = local_time_now->tm_mday;
    hour   = (double)local_time_now->tm_hour;
    minute = (double)local_time_now->tm_min;
    second = (double)local_time_now->tm_sec;
#endif
    jd = calculate_julian_day(year, month, day, hour, minute, second);
    printf("\nCalculating for date: %d-%02d-%02d %.0f:%.0f:%.0f UTC\n", year, month, day, hour, minute, second);
    printf("Calculated JD: %f\n", jd);

    /* Calculate heliocentric ecliptic coordinates for all planets */
    printf("Calculating heliocentric positions...\n");
    for (i = 0; i < num_planets_loaded; i++) {
        calculate_planet_helio_ecliptic_coords(all_planets[i], jd, &all_planet_ephems_data[i]);
    }

    /* Validate observer_planet_idx and proceed with observer-dependent calculations */
    if (observer_planet_idx >= num_planets_loaded || observer_planet_idx < 0) {
        fprintf(stderr, "Warning: Default observer index %d is invalid for %d loaded planets.\n", observer_planet_idx, num_planets_loaded);
        if (num_planets_loaded > 0) {
            observer_planet_idx = 0; /* Fallback to the first loaded planet */
            fprintf(stderr, "         Using planet '%s' (index %d) as observer instead.\n", all_planets[observer_planet_idx].name, observer_planet_idx);
        } else {
            fprintf(stderr, "Error: No planets loaded, cannot set an observer or calculate observer-dependent data.\n");
            /* Skip observer-dependent calculations if no planets are loaded */
            goto aphorism_section; /* Jump to aphorism processing */
        }
    }

    /* Calculate apparent Az/Alt for the observer */
    printf("Calculating apparent Az/Alt from %s (Lat:%.2f, Lon:%.2f)...\n",
           all_planets[observer_planet_idx].name, observer_latitude_deg, observer_longitude_deg);
    calculate_apparent_az_alt_all_planets(
        all_planets,
        all_planet_ephems_data,
        num_planets_loaded,
        observer_planet_idx,
        observer_latitude_deg,
        observer_longitude_deg,
        jd,
        apparent_sky_pos);

    /* Get relative angular separations and store in the combined vector */
    printf("Calculating relative angular separations...\n");
    num_relative_seps_actual = get_relative_angular_separations(
        apparent_sky_pos,
        num_planets_loaded,
        combined_angular_separations + current_buffer_offset,
        max_seps_per_type
    );
    if (num_relative_seps_actual >= 0) {
        printf("  Got %d relative angular separations.\n", num_relative_seps_actual);
        current_buffer_offset += num_relative_seps_actual;
    } else {
        fprintf(stderr, "  Error calculating relative angular separations.\n");
        num_relative_seps_actual = 0; /* Ensure count is non-negative for offset logic */
    }

    /* Get heliocentric ecliptic angular separations and append to the combined vector */
    printf("Calculating heliocentric ecliptic angular separations...\n");
    num_helio_seps_actual = get_helio_ecliptic_angular_separations(
        all_planet_ephems_data,
        num_planets_loaded,
        combined_angular_separations + current_buffer_offset,
        max_seps_per_type
    );
    if (num_helio_seps_actual >= 0) {
        printf("  Got %d heliocentric ecliptic angular separations.\n", num_helio_seps_actual);
        current_buffer_offset += num_helio_seps_actual;
    } else {
        fprintf(stderr, "  Error calculating heliocentric ecliptic angular separations.\n");
        num_helio_seps_actual = 0; /* Ensure count is non-negative */
    }
    total_seps_in_vector = current_buffer_offset;

    /* Print the combined vector of angular separations */
    if (total_seps_in_vector > 0) {
        printf("\nCombined Angular Separations Vector (total %d values):\n", total_seps_in_vector);
        for (k = 0; k < total_seps_in_vector; k++) {
            printf("  Value #%d: %.2f deg\n", k + 1, combined_angular_separations[k]);
            if (num_relative_seps_actual > 0 && k == num_relative_seps_actual - 1 && num_helio_seps_actual > 0) {
                 printf("  --- Heliocentric separations start after this value ---\n");
            }
        }
    } else {
        printf("\nNo angular separations were calculated to display.\n");
    }

    /* Generate astrological seeds if enough data is available */
    /* New algorithm needs at least (12 - 1) + 52 = 63 elements */
    if (total_seps_in_vector >= 63) {
        printf("\nGenerating astrological seeds from the first 64 angular separations...\n");
        if (generate_astrological_seeds(combined_angular_separations, astrological_seeds)) {
            printf("Successfully generated 12 astrological seeds:\n");
            for (m = 0; m < 12; ++m) {
                printf("  Seed #%d: %f\n", m + 1, astrological_seeds[m]);
            }
        } else {
            fprintf(stderr, "Failed to generate astrological seeds.\n");
        }
    } else {
        printf("\nNot enough angular separations (%d) to generate astrological seeds (need 63).\n",
               total_seps_in_vector);
    }

aphorism_section: /* Label for goto if ephemeris part is skipped */
    printf("\nReading Aphorism Templates...\n");
    if (read_aphorism_templates("aphorism_templates.txt", aphorism_templates_storage, MAX_TEMPLATES, MAX_TEMPLATE_LEN, &num_tpls_read)) {
        if (num_tpls_read > 0) {
            printf("Successfully read %d aphorism templates.\n", num_tpls_read);

            /* Select aphorism templates for signs if seeds were also generated */
            if (total_seps_in_vector >= 63) { /* Check if seeds were likely generated */
                printf("\nSelecting aphorism templates for 12 signs...\n");
                num_aphorisms_selected_for_signs = select_templates_for_signs(
                    astrological_seeds,
                    num_tpls_read,
                    selected_aphorism_indices
                );
                printf("Selected %d aphorisms for the signs:\n", num_aphorisms_selected_for_signs);
                for (m = 0; m < 12; ++m) {
                    if (selected_aphorism_indices[m] != -1) {
                        printf("  Sign %2d: Template #%2d (\"%.30s%s\")\n",
                               m + 1,
                               selected_aphorism_indices[m],
                               aphorism_templates_storage[selected_aphorism_indices[m]],
                               strlen(aphorism_templates_storage[selected_aphorism_indices[m]]) > 30 ? "..." : "");
                    } else {
                        printf("  Sign %2d: No template assigned.\n", m + 1);
                    }
                }

                /* Now, fill and print the selected aphorisms */
                if (num_aphorisms_selected_for_signs > 0) {
                    printf("\nGenerating aphorisms for the 12 signs:\n");
                    for (m = 0; m < 12; ++m) {
                        if (selected_aphorism_indices[m] != -1) {
                            char* current_template_str = aphorism_templates_storage[selected_aphorism_indices[m]];
                            printf("--- Sign %d (Template #%d) ---\n", m + 1, selected_aphorism_indices[m]);
                            
                            parse_aphorism_template(current_template_str, &n_count, &v_count, &nv_count);
                            printf("  Requires: (N)=%d, (V)=%d, (NV)=%d\n", n_count, v_count, nv_count);

                            if (fill_aphorism_template(current_template_str,
                                                       test_nouns, sizeof(test_nouns)/sizeof(char*),
                                                       test_verbs, sizeof(test_verbs)/sizeof(char*),
                                                       test_nv, sizeof(test_nv)/sizeof(char*),
                                                       filled_aphorism, MAX_TEMPLATE_LEN)) {
                                printf("  Aphorism: %s\n\n", filled_aphorism);
                            } else {
                                printf("  Could not generate aphorism for sign %d.\n\n", m + 1);
                            }
                        } /* No else needed, already printed "No template assigned" above */
                    }
                }
            } else {
                printf("\nSeeds not generated (or not enough angular separations). Cannot select aphorisms for signs.\n");
            }
        } else {
            printf("No aphorism templates found or file empty.\n");
        }
    } else {
        fprintf(stderr, "Failed to read aphorism templates.\n");
    }

    printf("\nProgram finished.\n");
    return 0;
}
