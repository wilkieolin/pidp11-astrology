/* aphorism_utils.h - Utilities for aphorism template processing */
/* Compatible with K&R C style */

#ifndef APHORISM_UTILS_H
#define APHORISM_UTILS_H

/* Constants */
#define MAX_TEMPLATES 50     /* Maximum number of templates to read */
#define MAX_TEMPLATE_LEN 256 /* Maximum length of a single template line */

#define APH_TRUE 1
#define APH_FALSE 0

/*
 * read_aphorism_templates
 * Reads aphorism templates from a specified file into a 2D char array.
 * filename: The path to the template file.
 * templates_array: A 2D char array to store the read templates.
 * max_tpls: The maximum number of templates the array can hold (must match first dim of templates_array).
 * max_len: The maximum length of a single template line (must match second dim of templates_array).
 * num_tpls_read: Pointer to an integer to store the number of templates successfully read.
 * Returns APH_TRUE on success, APH_FALSE on file error or if no templates are read.
 */
int read_aphorism_templates( /* char *filename, char templates_array[][MAX_TEMPLATE_LEN], int max_tpls, int max_len, int *num_tpls_read */ );

/*
 * parse_aphorism_template
 * Parses a template string to count the number of (N), (V), and (NV) placeholders.
 * template_str: The null-terminated template string to parse.
 * num_n: Pointer to int to store the count of (N) placeholders.
 * num_v: Pointer to int to store the count of (V) placeholders.
 * num_nv: Pointer to int to store the count of (NV) placeholders.
 */
void parse_aphorism_template( /* char *template_str, int *num_n, int *num_v, int *num_nv */ );

/*
 * fill_aphorism_template
 * Fills a template string by substituting placeholders with provided words.
 * template_str: The null-terminated template string.
 * nouns, verbs, noun_verbs: Arrays of null-terminated strings for placeholders.
 * num_nouns_avail, num_verbs_avail, num_noun_verbs_avail: Counts of available words.
 * output_buffer: Char array to store the resulting filled string.
 * output_buffer_size: The total size of output_buffer.
 * Returns APH_TRUE if successful, APH_FALSE on error (e.g., buffer too small, not enough words).
 */
int fill_aphorism_template( /* char *template_str, char *nouns[], int num_nouns_avail, char *verbs[], int num_verbs_avail, char *noun_verbs[], int num_noun_verbs_avail, char *output_buffer, int output_buffer_size */ );

/*
 * find_nearest_neighbor
 * Finds the word in a file whose 49 angular coordinates are closest
 * to the given input_angles. The constant NUM_ANGLES (49) is defined
 * in aphorism_utils.c.
 * input_angles: An array of doubles representing the target angles.
 * filename: The path to the space-delimited text file.
 * Returns a dynamically allocated string containing the nearest word (caller must free).
 * Returns NULL on error or if no suitable word is found.
 */
char *find_nearest_neighbor( /* double input_angles[], char *filename */ );

#endif /* APHORISM_UTILS_H */
