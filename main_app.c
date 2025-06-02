/* main_app.c - Main astrology application */
#include <stdio.h>
#include "ephemeris.h"      /* For ephemeris calculations */
#include "aphorism_utils.h" /* For aphorism generation */

/* K&R C style main function */
int main(argc, argv)
    int argc;
    char *argv[];
{
    char aphorism_templates_storage[MAX_TEMPLATES][MAX_TEMPLATE_LEN];
    int num_tpls_read;
    int n_count, v_count, nv_count;
    char filled_aphorism[MAX_TEMPLATE_LEN];

    /* Example: Using OrbitalElements from ephemeris.h */
    OrbitalElements all_planets[MAX_PLANETS];
    int num_planets_loaded;
    double jd;

    printf("Welcome to the K&R Astrology Program!\n");

    /* --- Initialize and use ephemeris functions --- */
    printf("\nInitializing Ephemeris Data...\n");
    if (!initialize_orbital_elements_from_file("ephemeris_data.txt", all_planets, &num_planets_loaded)) {
        fprintf(stderr, "Failed to load ephemeris data. Exiting.\n");
        return 1;
    }
    printf("Loaded %d celestial bodies.\n", num_planets_loaded);

    jd = calculate_julian_day(2023, 10, 26, 12, 0, 0); /* Example date */
    printf("Calculated JD: %f\n", jd);
    /* ... more calls to ephemeris functions ... */


    /* --- Initialize and use aphorism functions --- */
    printf("\nReading Aphorism Templates...\n");
    if (read_aphorism_templates("aphorism_templates.txt", aphorism_templates_storage, MAX_TEMPLATES, MAX_TEMPLATE_LEN, &num_tpls_read)) {
        if (num_tpls_read > 0) {
            printf("Successfully read %d aphorism templates.\n", num_tpls_read);

            /* Example: Parse the first template */
            parse_aphorism_template(aphorism_templates_storage[0], &n_count, &v_count, &nv_count);
            printf("First template requires: N=%d, V=%d, NV=%d\n", n_count, v_count, nv_count);

            /* Example: Fill the first template (provide some dummy words) */
            char *test_nouns[] = {"stars", "destiny"};
            char *test_verbs[] = {"align"};
            char *test_nv[] = {"cosmic insight"};

            if (fill_aphorism_template(aphorism_templates_storage[0],
                                       test_nouns, 2,
                                       test_verbs, 1,
                                       test_nv, 1,
                                       filled_aphorism, MAX_TEMPLATE_LEN)) {
                printf("Generated Aphorism: %s\n", filled_aphorism);
            } else {
                printf("Could not generate aphorism for the first template.\n");
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
