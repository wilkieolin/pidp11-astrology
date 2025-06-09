/* aphorism_utils.c - Utilities for aphorism template processing */
/* Compatible with K&R C style */

#include <stdio.h>
#include <string.h> /* For strncpy, strlen, strncmp, strcpy, strchr, strtok_r, strdup */
#include <stdlib.h> /* For strtod, free, NULL (strdup also often uses malloc from here) */
#include "aphorism_utils.h"

#define NUM_ANGLES 49
#define MAX_LINE_LENGTH 1024 /* Reasonably sized buffer for lines */
#define MAX_WORD_LENGTH 256  /* Max expected word length */

/* Substitute for DBL_MAX from <float.h> if not available */
/* Used for initializing a variable to a very large value before finding a minimum. */
/* 1.0e38 is a common large double literal. */
#define DBL_MAX_SUBSTITUTE 1.0e38

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
 * @brief Finds the word in a file whose 49 angular coordinates are closest
 *        to the given input_angles.
 *
 * The file is expected to have lines in the format:
 * word angle1 angle2 ... angle49
 *
 * Distance is calculated as the sum of squared differences between angles.
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
    char *nearest_word_str;
    double min_sq_distance;
    char current_word[MAX_WORD_LENGTH];
    double file_angles[NUM_ANGLES];
    char *p; /* Pointer to walk through the line_buffer */
    char *nl; /* For newline removal */
    int angles_parsed_count;
    int i, j; /* Loop counters */
    double current_sq_distance;
    double diff;
    char token_buffer[MAX_WORD_LENGTH]; /* For individual words/numbers */
    int token_len;
    int sscanf_ret;

    nearest_word_str = NULL; /* Initialize */
    min_sq_distance = DBL_MAX_SUBSTITUTE;

    file = fopen(filename, "r");
    if (!file) {
        perror("Error opening word angle data file");
        return NULL;
    }

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

        /* 2. Parse the NUM_ANGLES angles */
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

            sscanf_ret = sscanf(token_buffer, "%lf", &file_angles[i]);
            if (sscanf_ret != 1) { /* Check if sscanf successfully assigned one item */
                 angles_parsed_count = -1; /* Mark as error */
                 break;
            }
            angles_parsed_count++;
        }

        if (angles_parsed_count != NUM_ANGLES) {
            /* fprintf(stderr, "Warning: Line for word '%s' did not contain %d valid angles (parsed %d).\n", current_word, NUM_ANGLES, (angles_parsed_count < 0 ? 0 : angles_parsed_count)); */
            continue; /* Skip this line if not all angles were parsed correctly */
        }

        /* 3. Calculate squared Euclidean distance */
        current_sq_distance = 0.0;
        for (j = 0; j < NUM_ANGLES; ++j) { /* Use j to avoid conflict with outer i */
            diff = file_angles[j] - input_angles[j];
            current_sq_distance += diff * diff;
        }

        /* 4. Update nearest neighbor if this one is closer */
        if (current_sq_distance < min_sq_distance) {
            min_sq_distance = current_sq_distance;
            if (nearest_word_str != NULL) {
                free(nearest_word_str); /* Free previous word if any */
            }
            nearest_word_str = strdup(current_word); /* strdup allocates and copies */
            if (!nearest_word_str) {
                perror("Error allocating memory for nearest_word_str");
                fclose(file);
                return NULL; /* Critical memory error */
            }
        }
    }

    fclose(file);
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
