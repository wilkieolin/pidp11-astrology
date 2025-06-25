/* aphorism_utils.c - Utilities for aphorism template processing */
/* Compatible with K&R C style */

#include <stdio.h>
#include <string.h> /* For strncpy, strlen, strncmp, strcpy, strchr, strtok_r, strdup */
#include <stdlib.h> /* For strtod, free, NULL (strdup also often uses malloc from here) */
#include <time.h>   /* For time(), difftime() */
#include <math.h>   /* For sqrt - assuming this is available in the build environment */
#include "aphorism_utils.h"

#define NUM_ANGLES 49
#define MAX_LINE_LENGTH 1024 /* Reasonably sized buffer for lines */
#define MAX_WORD_LENGTH 256  /* Max expected word length */

/* Substitute for DBL_MAX from <float.h> if not available */
/* Used for initializing a variable to a very large value before finding a minimum. */
/* 1.0e38 is a common large double literal. */
#define DBL_MAX_SUBSTITUTE 1.0e38

/* Define PI if not available (e.g., in strict K&R C without math.h M_PI) */
#define MY_PI 3.14159265358979323846

/* Decay rate for word radius. A value of 8.02e-6 causes the "excess" radius (r-1) to halve roughly every 24 hours. */
#define RADIUS_DECAY_RATE 8.02e-6

/* Amount to increase a word's radius after it has been selected. */
#define RADIUS_INCREASE 1.0

/* Set to 1 to enable debug prints in find_nearest_neighbor, 0 to disable */
#define DEBUG_NEAREST_NEIGHBOR_FLAG 0 /* Set to 1 to enable */

/*
 * read_aphorism_templates
 * Reads aphorism templates from a specified file into a 2D char array.
 */
int read_aphorism_templates(filename, templates_array, max_tpls, max_len, num_tpls_read)
    char *filename;
    char templates_array[][MAX_TEMPLATE_LEN]; /* K&R: Or char (*templates_array)[MAX_TEMPLATE_LEN]; */
    int max_tpls;
    int max_len;
    int *num_tpls_read;
{
    FILE *fp;
    char line_buffer[MAX_TEMPLATE_LEN]; /* Temporary buffer for fgets */
    int count;
    int len; /* For strlen result, K&R strlen returns int */

    count = 0;
    *num_tpls_read = 0; /* Initialize */

    fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Error opening template file");
        return APH_FALSE;
    }

    while (count < max_tpls && fgets(line_buffer, max_len, fp) != NULL) {
        /* Remove trailing newline, if any */
        len = strlen(line_buffer);
        if (len > 0 && line_buffer[len - 1] == '\n') {
            line_buffer[len - 1] = '\0';
            len--; /* Update length */
        }
        /* Also remove trailing carriage return for Windows-style line endings */
        if (len > 0 && line_buffer[len - 1] == '\r') {
            line_buffer[len - 1] = '\0';
        }

        if (strlen(line_buffer) > 0) { /* Only add non-empty lines */
            strncpy(templates_array[count], line_buffer, max_len - 1);
            templates_array[count][max_len - 1] = '\0'; /* Ensure null termination */
            count++;
        }
    }

    fclose(fp);
    *num_tpls_read = count;

    if (count == 0) {
        /* This could be a valid case (empty file), but often indicates an issue. */
        /* fprintf(stderr, "Warning: No templates read from %s\n", filename); */
        /* For this function, returning APH_TRUE for an empty file might be acceptable if 0 templates is valid. */
        /* However, the prompt implies we expect templates, so let's consider 0 an issue for now. */
        /* If an empty set of templates is a valid outcome, this check could be removed or changed. */
    }

    return APH_TRUE; /* Success, even if count is 0, file was processed. */
                     /* Caller should check *num_tpls_read. */
}

/*
 * parse_aphorism_template
 * Parses a template string to count the number of (N), (V), and (NV) placeholders.
 */
void parse_aphorism_template(template_str, num_n, num_v, num_nv)
    char *template_str;
    int *num_n;
    int *num_v;
    int *num_nv;
{
    char *p;
    
    *num_n = 0;
    *num_v = 0;
    *num_nv = 0;

    if (template_str == NULL) {
        return;
    }

    p = template_str;
    while (*p != '\0') {
        if (*p == '(') {
            if (strncmp(p, "(N)", 3) == 0) {
                (*num_n)++;
                p += 3; /* Advance past "(N)" */
            } else if (strncmp(p, "(V)", 3) == 0) {
                (*num_v)++;
                p += 3; /* Advance past "(V)" */
            } else if (strncmp(p, "(NV)", 4) == 0) {
                (*num_nv)++;
                p += 4; /* Advance past "(NV)" */
            } else {
                p++; /* Not a recognized placeholder, advance one char */
            }
        } else {
            p++; /* Not a '(', advance one char */
        }
    }
}

/*
 * fill_aphorism_template
 * Fills a template string by substituting placeholders with provided words.
 */
int fill_aphorism_template(template_str,
                           nouns, num_nouns_avail,
                           verbs, num_verbs_avail,
                           noun_verbs, num_noun_verbs_avail,
                           output_buffer, output_buffer_size)
    char *template_str;
    char *nouns[];
    int num_nouns_avail;
    char *verbs[];
    int num_verbs_avail;
    char *noun_verbs[];
    int num_noun_verbs_avail;
    char *output_buffer;
    int output_buffer_size;
{
    char *p_template; /* Pointer to current position in template_str */
    char *p_out;      /* Pointer to current position in output_buffer */
    int noun_idx;
    int verb_idx;
    int nv_idx;
    int remaining_buffer_size;
    int word_len; /* K&R strlen returns int */

    noun_idx = 0;
    verb_idx = 0;
    nv_idx = 0;

    if (template_str == NULL || output_buffer == NULL || output_buffer_size <= 0) {
        return APH_FALSE;
    }

    p_template = template_str;
    p_out = output_buffer;
    remaining_buffer_size = output_buffer_size;
    *p_out = '\0'; /* Start with an empty string in case template is empty */

    while (*p_template != '\0') {
        if (remaining_buffer_size <= 1) { /* Need space for at least one char + null terminator */
            fprintf(stderr, "Error: Output buffer too small during fill.\n");
            output_buffer[output_buffer_size - 1] = '\0'; /* Ensure null termination */
            return APH_FALSE;
        }

        if (*p_template == '(') {
            if (strncmp(p_template, "(N)", 3) == 0) {
                if (noun_idx < num_nouns_avail && nouns[noun_idx] != NULL) {
                    word_len = strlen(nouns[noun_idx]);
                    if (word_len < remaining_buffer_size) { /* word_len + 1 for null <= remaining_buffer_size */
                        strcpy(p_out, nouns[noun_idx]);
                        p_out += word_len;
                        remaining_buffer_size -= word_len;
                        noun_idx++;
                        p_template += 3;
                    } else {
                        fprintf(stderr, "Error: Output buffer too small for noun '%s'.\n", nouns[noun_idx]);
                        output_buffer[output_buffer_size - 1] = '\0'; return APH_FALSE;
                    }
                } else {
                    fprintf(stderr, "Error: Not enough nouns provided for template.\n");
                    output_buffer[output_buffer_size - 1] = '\0'; return APH_FALSE;
                }
            } else if (strncmp(p_template, "(V)", 3) == 0) {
                if (verb_idx < num_verbs_avail && verbs[verb_idx] != NULL) {
                    word_len = strlen(verbs[verb_idx]);
                    if (word_len < remaining_buffer_size) {
                        strcpy(p_out, verbs[verb_idx]);
                        p_out += word_len;
                        remaining_buffer_size -= word_len;
                        verb_idx++;
                        p_template += 3;
                    } else {
                        fprintf(stderr, "Error: Output buffer too small for verb '%s'.\n", verbs[verb_idx]);
                        output_buffer[output_buffer_size - 1] = '\0'; return APH_FALSE;
                    }
                } else {
                    fprintf(stderr, "Error: Not enough verbs provided for template.\n");
                    output_buffer[output_buffer_size - 1] = '\0'; return APH_FALSE;
                }
            } else if (strncmp(p_template, "(NV)", 4) == 0) {
                if (nv_idx < num_noun_verbs_avail && noun_verbs[nv_idx] != NULL) {
                    word_len = strlen(noun_verbs[nv_idx]);
                    if (word_len < remaining_buffer_size) {
                        strcpy(p_out, noun_verbs[nv_idx]);
                        p_out += word_len;
                        remaining_buffer_size -= word_len;
                        nv_idx++;
                        p_template += 4;
                    } else {
                        fprintf(stderr, "Error: Output buffer too small for noun/verb '%s'.\n", noun_verbs[nv_idx]);
                        output_buffer[output_buffer_size - 1] = '\0'; return APH_FALSE;
                    }
                } else {
                    fprintf(stderr, "Error: Not enough noun/verbs provided for template.\n");
                    output_buffer[output_buffer_size - 1] = '\0'; return APH_FALSE;
                }
            } else { /* Not a recognized placeholder, copy the character */
                *p_out++ = *p_template++;
                remaining_buffer_size--;
            }
        } else { /* Not a '(', copy the character */
            *p_out++ = *p_template++;
            remaining_buffer_size--;
        }
    }

    *p_out = '\0'; /* Null-terminate the output string */
    return APH_TRUE;
}

/**
 * @brief Finds a word in a file, increases its radius, and rewrites the file.
 *
 * This function reads a given file line by line, searching for a specific word.
 * When the word is found, it parses the line to find the radius, adds the
 * specified `radius_increase` to it, and writes the modified line to a
 * temporary file. All other lines are copied verbatim.
 * If the word is found and all operations succeed, the original file is
 * replaced by the temporary file.
 *
 * The file is expected to have lines in the format:
 * word radius angle1 angle2 ...
 *
 * @param word_to_update The word whose radius should be modified.
 * @param radius_increase The amount to add to the current radius.
 * @param filename The path to the space-delimited text file.
 * @return APH_TRUE on success. Returns APH_FALSE if the file cannot be
 *         opened, a temporary file cannot be created, the word is not found,
 *         or a file operation (remove/rename) fails.
 */
int update_word_radius(word_to_update, radius_increase, filename)
    char *word_to_update;
    double radius_increase;
    char *filename;
{
    /* K&R C: All variable declarations at the top of the function block */
    FILE *in_file, *out_file;
    char temp_filename[MAX_WORD_LENGTH]; /* Assuming filename isn't excessively long */
    char line_buffer[MAX_LINE_LENGTH];
    char original_line[MAX_LINE_LENGTH]; /* To preserve original line for writing */
    char current_word[MAX_WORD_LENGTH];
    char *p;
    int token_len;
    int word_found = APH_FALSE;
    int sscanf_ret;
    double old_radius, new_radius;
    char *rest_of_line;

    if (word_to_update == NULL || filename == NULL) {
        return APH_FALSE;
    }

    /* Create a temporary filename */
    strcpy(temp_filename, filename);
    strcat(temp_filename, ".tmp");

    in_file = fopen(filename, "r");
    if (in_file == NULL) {
        perror("update_word_radius: Cannot open input file");
        return APH_FALSE;
    }

    out_file = fopen(temp_filename, "w");
    if (out_file == NULL) {
        perror("update_word_radius: Cannot create temporary file");
        fclose(in_file);
        return APH_FALSE;
    }

    while (fgets(line_buffer, sizeof(line_buffer), in_file) != NULL) {
        strcpy(original_line, line_buffer); /* Keep a pristine copy */
        p = line_buffer;

        /* Parse the first word from the line */
        token_len = 0;
        while (*p != ' ' && *p != '\t' && *p != '\n' && *p != '\r' && *p != '\0' && token_len < MAX_WORD_LENGTH - 1) {
            current_word[token_len++] = *p++;
        }
        current_word[token_len] = '\0';

        /* Check if this is the word to update */
        if (token_len > 0 && strcmp(current_word, word_to_update) == 0) {
            word_found = APH_TRUE;

            while (*p == ' ' || *p == '\t') p++;
            rest_of_line = p;

            sscanf_ret = sscanf(rest_of_line, "%lf", &old_radius);
            if (sscanf_ret != 1) {
                fprintf(stderr, "Warning: Could not parse radius for word '%s'. Writing original line.\n", current_word);
                fputs(original_line, out_file);
                continue;
            }
            new_radius = old_radius + radius_increase;
            while (*rest_of_line != ' ' && *rest_of_line != '\t' && *rest_of_line != '\n' && *rest_of_line != '\r' && *rest_of_line != '\0') {
                rest_of_line++;
            }
            fprintf(out_file, "%s %f%s", current_word, new_radius, rest_of_line);
        } else {
            fputs(original_line, out_file);
        }
    }

    fclose(in_file);
    fclose(out_file);

    if (word_found) {
        if (remove(filename) != 0) { perror("update_word_radius: Error removing original file"); remove(temp_filename); return APH_FALSE; }
        if (rename(temp_filename, filename) != 0) { perror("update_word_radius: Error renaming temporary file"); return APH_FALSE; }
        return APH_TRUE;
    } else {
        remove(temp_filename);
        return APH_FALSE; /* Word not found */
    }
}

/**
 * @brief Applies exponential decay to word radii in a file based on elapsed time.
 *
 * This function reads a timestamp from `access_filename` to determine the
 * time elapsed since the last run. It then reads `word_filename` line by line,
 * calculating a new, smaller radius for each word based on an exponential
 * decay formula. The goal is to smoothly return radii to a baseline of 1.0.
 *
 * The function writes the updated data to a temporary file and then replaces
 * the original `word_filename`. Finally, it updates `access_filename` with the
 * current timestamp.
 *
 * @param word_filename The path to the word data file (e.g., "words/nouns.txt").
 * @param access_filename The path to the file storing the last update timestamp.
 * @return APH_TRUE on success, APH_FALSE on any failure.
 * @note Because this function updates the timestamp file upon completion, calling
 *       it in a loop for multiple word files that share a single timestamp file
 *       will result in only the first file's radii being decayed.
 */
int decay_word_radius(word_filename, access_filename)
    char *word_filename;
    char *access_filename;
{
    /* K&R C: All variable declarations at the top of the function block */
    FILE *in_file, *out_file, *access_file;
    char temp_filename[MAX_WORD_LENGTH];
    char line_buffer[MAX_LINE_LENGTH];
    char original_line[MAX_LINE_LENGTH];
    char current_word[MAX_WORD_LENGTH];
    char *p, *rest_of_line;
    int token_len;
    int sscanf_ret;
    double old_radius, new_radius;
    time_t current_time;
    long last_update_time_long; /* Use long for reading from file */
    time_t last_update_time;
    double time_diff_seconds;

    if (word_filename == NULL || access_filename == NULL) {
        return APH_FALSE;
    }

    /* 1. Get current time and last update time from access file */
    current_time = time(NULL);
    last_update_time = 0; /* Default to 0 */

    access_file = fopen(access_filename, "r");
    if (access_file != NULL) {
        if (fscanf(access_file, "%ld", &last_update_time_long) == 1) {
            last_update_time = (time_t)last_update_time_long;
        } else {
            last_update_time = current_time; /* File exists but is empty/corrupt, so no decay */
        }
        fclose(access_file);
    } else {
        /* File doesn't exist, this is the first run. No decay needed. */
        last_update_time = current_time;
    }

    time_diff_seconds = difftime(current_time, last_update_time);

    /* If no time has passed or time went backwards, no decay is needed. */
    /* We still update the timestamp to the current time and report success. */
    if (time_diff_seconds <= 0) {
        access_file = fopen(access_filename, "w");
        if (access_file == NULL) {
            perror("decay_word_radius: Cannot open access file for writing");
            return APH_FALSE;
        }
        fprintf(access_file, "%ld\n", (long)current_time);
        fclose(access_file);
        return APH_TRUE; /* Success, as no decay was necessary. */
    }

    /* 2. Create a temporary filename for writing updated word data */
    strcpy(temp_filename, word_filename);
    strcat(temp_filename, ".tmp");

    in_file = fopen(word_filename, "r");
    if (in_file == NULL) { perror("decay_word_radius: Cannot open input word file"); return APH_FALSE; }

    out_file = fopen(temp_filename, "w");
    if (out_file == NULL) { perror("decay_word_radius: Cannot create temporary file"); fclose(in_file); return APH_FALSE; }

    /* 3. Process each line of the word file */
    while (fgets(line_buffer, sizeof(line_buffer), in_file) != NULL) {
        strcpy(original_line, line_buffer);
        p = line_buffer;

        token_len = 0;
        while (*p != ' ' && *p != '\t' && *p != '\n' && *p != '\r' && *p != '\0' && token_len < MAX_WORD_LENGTH - 1) { current_word[token_len++] = *p++; }
        current_word[token_len] = '\0';

        if (token_len == 0) { fputs(original_line, out_file); continue; }

        while (*p == ' ' || *p == '\t') p++;
        rest_of_line = p;

        sscanf_ret = sscanf(rest_of_line, "%lf", &old_radius);
        if (sscanf_ret != 1) { fprintf(stderr, "Warning: Could not parse radius for word '%s'. Writing original line.\n", current_word); fputs(original_line, out_file); continue; }

        if (old_radius > 1.0) { new_radius = 1.0 + (old_radius - 1.0) * exp(-RADIUS_DECAY_RATE * time_diff_seconds); } else { new_radius = old_radius; }

        while (*rest_of_line != ' ' && *rest_of_line != '\t' && *rest_of_line != '\n' && *rest_of_line != '\r' && *rest_of_line != '\0') { rest_of_line++; }

        fprintf(out_file, "%s %f%s", current_word, new_radius, rest_of_line);
    }

    fclose(in_file);
    fclose(out_file);

    /* 4. Replace the original file with the temporary file */
    if (remove(word_filename) != 0) { perror("decay_word_radius: Error removing original file"); remove(temp_filename); return APH_FALSE; }
    if (rename(temp_filename, word_filename) != 0) { perror("decay_word_radius: Error renaming temporary file"); return APH_FALSE; }

    /* 5. If all successful, update the access file with the current time */
    access_file = fopen(access_filename, "w");
    if (access_file == NULL) {
        perror("decay_word_radius: Cannot open access file for writing final timestamp");
        return APH_FALSE;
    }
    fprintf(access_file, "%ld\n", (long)current_time);
    fclose(access_file);

    return APH_TRUE;
}

/**
 * @brief Finds the word in a file whose vector is closest to the given input angles.
 *
 * The file is expected to have lines in the format:
 * word radius angle1 angle2 ... angle49
 *
 * Distance is calculated as the cosine distance between the angular vectors,
 * scaled by the word's radius. A larger radius increases the distance, making
 * the word less likely to be chosen. The word with the minimum scaled distance is chosen.
 *
 * @param input_angles An array of NUM_ANGLES doubles representing the target angles.
 * @param filename The path to the space-delimited text file.
 * @return A dynamically allocated string containing the nearest word.
 *         The caller is responsible for freeing this memory.
 *         Returns NULL if an error occurs (e.g., file not found, memory allocation failed)
 *         or if no valid data is found in the file.
 */
char*
find_nearest_neighbor(input_angles, filename)
    double input_angles[NUM_ANGLES]; /* K&R: or double input_angles[] */
    char* filename;
{
    /* K&R C: All variable declarations at the top of the function block */
    FILE *file;
    char line_buffer[MAX_LINE_LENGTH];
    char *nearest_word_str; /* Stores the word with the smallest distance */
    double min_distance; /* Smallest scaled distance found so far */
    char current_word[MAX_WORD_LENGTH];
    double file_angles[NUM_ANGLES];
    double file_radius; /* Radius of the word in spherical coordinates */
    char *p; /* Pointer to walk through the line_buffer */
    char *nl; /* For newline removal */
    int angles_parsed_count;
    int i, j, k; /* Loop counters */
    char token_buffer[MAX_WORD_LENGTH]; /* For individual words/numbers */
    int token_len;
    int sscanf_ret;
    double angle_from_file_rad; /* Temporary variable to store angle read from file in radians */
    double input_angles_deg[NUM_ANGLES]; /* To store input angles in degrees */
    double dot_product;
    double magnitude_input_deg;
    double magnitude_file_deg;
    double similarity;
    double current_cosine_distance; /* Current word's cosine distance to input */
    double current_distance; /* Cosine distance scaled by radius */

    nearest_word_str = NULL; /* Initialize */
    min_distance = DBL_MAX_SUBSTITUTE; /* Initialize with a large value */

    file = fopen(filename, "r");
    if (!file) {
        perror("Error opening word angle data file");
        return NULL;
    }

    /* Convert input_angles to degrees and calculate its magnitude */
    magnitude_input_deg = 0.0;
    for (k = 0; k < NUM_ANGLES; ++k) {
        input_angles_deg[k] = input_angles[k] * (180.0 / MY_PI);
        magnitude_input_deg += input_angles_deg[k] * input_angles_deg[k];
    }
    /* It's important that sqrt is available. If not, this will not compile/link. */
    /* K&R C itself doesn't standardize math.h, but many compilers provide it. */
    /* The use of other ANSI C functions (strtod, strdup) suggests it might be okay. */
    magnitude_input_deg = sqrt(magnitude_input_deg);

    while (fgets(line_buffer, sizeof(line_buffer), file)) {
        p = line_buffer;

        /* Remove newline character if present */
        nl = strchr(line_buffer, '\n');
        if (nl) *nl = '\0';
        nl = strchr(line_buffer, '\r'); /* Handle CRNL too */
        if (nl) *nl = '\0';

        /* 1. Parse the word */
        token_len = 0;
        while (*p != ' ' && *p != '\0' && token_len < MAX_WORD_LENGTH - 1) {
            token_buffer[token_len++] = *p++;
        }
        token_buffer[token_len] = '\0';

        if (token_len == 0) { /* Empty line or line starts with space */
            /* fprintf(stderr, "Warning: Malformed line (missing word): %s\n", line_buffer); */
            continue; /* Skip empty or malformed line */
        }
        strcpy(current_word, token_buffer);

        /* 2. Parse the radius */
        while (*p == ' ') p++; /* Skip leading spaces */
        if (*p == '\0') {
            /* fprintf(stderr, "Warning: Line for word '%s' is missing radius and angles.\n", current_word); */
            continue;
        }
        token_len = 0;
        while (*p != ' ' && *p != '\0' && token_len < MAX_WORD_LENGTH - 1) {
            token_buffer[token_len++] = *p++;
        }
        token_buffer[token_len] = '\0';
        if (token_len == 0) {
            /* fprintf(stderr, "Warning: Line for word '%s' is missing radius and angles.\n", current_word); */
            continue;
        }
        sscanf_ret = sscanf(token_buffer, "%lf", &file_radius);
        if (sscanf_ret != 1) {
            /* fprintf(stderr, "Warning: Could not parse radius for word '%s'.\n", current_word); */
            continue;
        }

        /* 3. Parse the NUM_ANGLES angles */
        angles_parsed_count = 0;
        for (i = 0; i < NUM_ANGLES; ++i) {
            while (*p == ' ') p++; /* Skip leading spaces for the next token */

            if (*p == '\0') { /* End of line, not enough angles */
                break; /* Not enough tokens for all angles */
            }

            token_len = 0;
            while (*p != ' ' && *p != '\0' && token_len < MAX_WORD_LENGTH - 1) {
                token_buffer[token_len++] = *p++;
            }
            token_buffer[token_len] = '\0';

            if (token_len == 0) break; /* No token found (e.g. trailing spaces) */

            sscanf_ret = sscanf(token_buffer, "%lf", &angle_from_file_rad);
            if (sscanf_ret != 1) { /* Check if sscanf successfully assigned one item */
                 angles_parsed_count = -1; /* Mark as error */
                 break;
            }
            file_angles[i] = angle_from_file_rad * (180.0 / MY_PI); /* Convert to degrees */
            angles_parsed_count++;
        }

        if (angles_parsed_count != NUM_ANGLES) {
            /* fprintf(stderr, "Warning: Line for word '%s' did not contain %d valid angles (parsed %d).\n", current_word, NUM_ANGLES, (angles_parsed_count < 0 ? 0 : angles_parsed_count)); */
            continue; /* Skip this line if not all angles were parsed correctly */
        }

        /* 4. Calculate Scaled Cosine Distance */
        dot_product = 0.0;
        magnitude_file_deg = 0.0;

        for (j = 0; j < NUM_ANGLES; ++j) {
            dot_product += input_angles_deg[j] * file_angles[j];
            magnitude_file_deg += file_angles[j] * file_angles[j];
        }
        magnitude_file_deg = sqrt(magnitude_file_deg);

        if (magnitude_input_deg == 0.0 || magnitude_file_deg == 0.0) {
            /* If either vector has zero magnitude, cosine similarity is undefined or 0. */
            /* Treat as dissimilar (similarity = 0, distance = 1). */
            similarity = 0.0;
        } else {
            similarity = dot_product / (magnitude_input_deg * magnitude_file_deg);
        }

        /* Clamp similarity to [-1, 1] to handle potential floating point inaccuracies */
        if (similarity > 1.0) similarity = 1.0;
        if (similarity < -1.0) similarity = -1.0;

        current_cosine_distance = 1.0 - similarity;

        /* Scale distance by radius. A larger radius means the word is "further" away. */
        /* This implements the idea that similarity should decay with radius. */
        current_distance = current_cosine_distance * file_radius;

        /* 5. Update nearest neighbor if this one is closer */
        if (current_distance < min_distance) {
            min_distance = current_distance;
            if (nearest_word_str != NULL) {
                free(nearest_word_str); /* Free previous word if any */
            }
            nearest_word_str = strdup(current_word); /* strdup allocates and copies */
            if (!nearest_word_str) {
                perror("Error allocating memory for nearest_word_str");
                fclose(file);
                return NULL; /* Critical memory error */
            }
            if (DEBUG_NEAREST_NEIGHBOR_FLAG) {
                fprintf(stderr, "[DEBUG_APH_UTIL] New nearest word: '%s' (Distance: %f)\n",
                        nearest_word_str, min_distance);
            }
        }
    }

    fclose(file);

    /* If a nearest word was found, update its radius to make it less likely to be chosen again soon. */
    if (nearest_word_str != NULL) {
        if (!update_word_radius(nearest_word_str, RADIUS_INCREASE, filename)) {
            /* This is a non-fatal warning. The program can continue even if the update fails. */
            fprintf(stderr, "Warning (find_nearest_neighbor): Failed to update radius for word '%s' in file '%s'.\n",
                    nearest_word_str, filename);
        }
    }

    return nearest_word_str; /* Caller must free this */
}

/* Define APHORISM_UTILS_MAIN_TEST to compile this main function for testing */
#ifdef APHORISM_UTILS_MAIN_TEST
int main() {
    char templates[MAX_TEMPLATES][MAX_TEMPLATE_LEN];
    int num_templates_read;
    int i;
    int n_count, v_count, nv_count;
    char filled_aphorism[MAX_TEMPLATE_LEN];

    char *test_nouns[] = {"wisdom", "stars", "future", "code"};
    char *test_verbs[] = {"guides", "illuminates", "reveals", "compiles"};
    char *test_nv[] = {"clarity", "insight", "effort", "debugging"};

    printf("Attempting to read aphorism templates from 'aphorism_templates.txt'...\n");
    if (read_aphorism_templates("aphorism_templates.txt", templates, MAX_TEMPLATES, MAX_TEMPLATE_LEN, &num_templates_read)) {
        printf("Successfully read %d templates.\n\n", num_templates_read);

        for (i = 0; i < num_templates_read; ++i) {
            printf("Template #%d: \"%s\"\n", i + 1, templates[i]);
            parse_aphorism_template(templates[i], &n_count, &v_count, &nv_count);
            printf("  Requires: (N)=%d, (V)=%d, (NV)=%d\n", n_count, v_count, nv_count);

            if (fill_aphorism_template(templates[i],
                                       test_nouns, sizeof(test_nouns)/sizeof(char*),
                                       test_verbs, sizeof(test_verbs)/sizeof(char*),
                                       test_nv, sizeof(test_nv)/sizeof(char*),
                                       filled_aphorism, MAX_TEMPLATE_LEN)) {
                printf("  Filled: \"%s\"\n\n", filled_aphorism);
            } else {
                printf("  Failed to fill template #%d.\n\n", i + 1);
            }
        }
    } else {
        fprintf(stderr, "Failed to read templates or file was empty/not found.\n");
        return 1;
    }

    return 0;
}
#endif /* APHORISM_UTILS_MAIN_TEST */
