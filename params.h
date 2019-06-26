#ifndef JW_PARAMS
#define JW_PARAMS

#include <stdio.h>
#include <stdarg.h>

/*!
 * Parameter file error codes.
 */
#define PDATA_SUCCESS 1
#define PDATA_FAILURE 0
#define PDATA_FORMAT -1
#define PDATA_MEM_ERROR -2

/*!
 * Constants that define the data structure.
 */
#define KEYSIZE 80
#define VALSIZE 80
#define MAX_ARRAY_SIZE 2000

/*! 
 * Hash table bucket structure.
 */
typedef struct{
  unsigned int hkey;
  char key[KEYSIZE];
  char value[VALSIZE];
} ht_bucket;

/*!
 * Hash table structure (including size).
 */
typedef struct{
  ht_bucket array[MAX_ARRAY_SIZE];
  size_t array_size;
} pdata;

/*!
 * Enumeration of data types.
 */
typedef enum{ DOUBLE_T, INT_T, STRING_T} pdata_type;

void pdata_initialize( pdata *pd);
int pdata_read_file( pdata *pd, FILE *f);

int pdata_get_var_d( pdata *pd, const char *name, double *val);
int pdata_get_var_i( pdata *pd, const char *name, int *val);
int pdata_get_var_s( pdata *pd, const char *name, char *val);

int pdata_array_length( pdata *pd, const char *name);
int pdata_get_array_d( pdata *pd, const char *name, double *val);
int pdata_get_array_i( pdata *pd, const char *name, int *val);
int pdata_get_array_s( pdata *pd, const char *name, char **val);

int pdata_get_element_d( pdata *pd, const char *name, 
                        unsigned int i, double *val);
int pdata_get_element_i( pdata *pd, const char *name, 
                        unsigned int i, int *val);
int pdata_get_element_s( pdata *pd, const char *name, 
                        unsigned int i, char *val);

int pdata_get_list( pdata *pd, const char *name, int num, ...);

#endif
