#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>

#include "params.h"


/*!
 * The parameter file has list of names and values which we will store
 * in a hash table where the variable names are the keys. 
 */

#define MAXFILL 0.75

/*! 
 * Hash table status codes. 
 */
#define HT_SUCCESS 0
#define HT_FAILURE -1
#define HT_MEM_ERROR 1

/*!
 * Hashing function.
 */
static uint32_t oat_hash( const char *key){
  uint32_t hkey = 0x00000000;
  int i;
  
  for( i=0; key[i] != '\0' && i<KEYSIZE; i++){
    hkey += (uint32_t) key[i];
    hkey += ( hkey << 10 );
    hkey ^= ( hkey >> 6 );
  }
  hkey += ( hkey << 3 );
  hkey ^= ( hkey >> 11 );
  hkey += ( hkey << 15 );

  return hkey;
}


/*!
 * Find in hashtable.
 */
static int hashtab_find( pdata *pd, const char *key, char *value){
  unsigned int i, hkey;

  /* Hash the key */
  hkey = oat_hash( key) % MAX_ARRAY_SIZE;
  i = hkey;

  /* Find first matching key */
  while( strcmp( pd->array[i].key, key) != 0 && 
    pd->array[i].key[0] != '\0' ){
    i=(i+1)%MAX_ARRAY_SIZE;
  }
  if( pd->array[i].key[0] == '\0'){
    return HT_FAILURE;
  }

  strcpy( value, pd->array[i].value);
  return HT_SUCCESS;
}

/*!
 * Add an entry in the hash table.
 */
static void ht_add( pdata *pd, unsigned int index,
        unsigned int hkey, const char *key, const char *value){
  pd->array[index].hkey = hkey;
  strcpy( pd->array[index].key, key);
  strcpy( pd->array[index].value, value);
  pd->array_size++;
}

/*!
 * Shift entries in a hash table down.
 *
 * If the hash table already has an entry at the hashed key `hkey`
 * then we need to shift all the currently filled members down to
 * insert the new value.
 */
static void ht_shift_down( pdata *pd, unsigned int index){
  unsigned int k = 0;

  do{ 
    k++; 
  }while( pd->array[(index+k)%MAX_ARRAY_SIZE].key[0] != '\0');
  while( k>0 ){
    pd->array[(index+k)%MAX_ARRAY_SIZE].hkey = 
      pd->array[(index+k-1)%MAX_ARRAY_SIZE].hkey;
    strcpy( pd->array[(index+k)%MAX_ARRAY_SIZE].key, 
           pd->array[(index+k-1)%MAX_ARRAY_SIZE].key);
    strcpy( pd->array[(index+k)%MAX_ARRAY_SIZE].value, 
           pd->array[(index+k-1)%MAX_ARRAY_SIZE].value);
    k--;
  }
}

/*!
 * Insert into hash table.
 */
static int hashtab_insert( pdata *pd, const char *key, 
                          const char *value){
  unsigned int i, j, hkey;

  if( (double) pd->array_size / (double) MAX_ARRAY_SIZE > MAXFILL ){
    return HT_MEM_ERROR;
  }

  /* Hash the key */
  hkey = oat_hash( key) % MAX_ARRAY_SIZE;
  i = hkey;
  /* Check for direct input */
  if( pd->array[i].key[0] == '\0'){
    ht_add( pd, i, hkey, key, value);
    return HT_SUCCESS;
  }else if( strcmp( pd->array[i].key, key) == 0){
    strcpy( pd->array[i].value, value);
    return HT_SUCCESS;
    
  /* If spot is occupied shift to next spot */
  }else if( pd->array[i].hkey == hkey){
    do{ 
      i=(i+1)%MAX_ARRAY_SIZE; 
    }while( pd->array[i].hkey != hkey &&
      pd->array[i].key[0] != '\0' );
  }
  for( j=0; j<pd->array_size+1; j++){

    /* Direct input if we can */
    if( pd->array[(i+j)%MAX_ARRAY_SIZE].key[0] == '\0'){
      ht_add( pd, (i+j)%MAX_ARRAY_SIZE, hkey, key, value);
      return HT_SUCCESS;
    }
    if( strcmp( pd->array[(i+j)%MAX_ARRAY_SIZE].key, key) == 0){
      strcpy( pd->array[(i+j)%MAX_ARRAY_SIZE].value, value);
      return HT_SUCCESS;
    }
    /* Otherwise shift down and input */
    if( pd->array[(i+j)%MAX_ARRAY_SIZE].hkey != hkey){
      ht_shift_down( pd, (i+j)%MAX_ARRAY_SIZE);
      ht_add( pd, (i+j)%MAX_ARRAY_SIZE, hkey, key, value);
      return HT_SUCCESS;
    }
  }
  return HT_FAILURE;
}


/*!
 * Starting the section of reading the parameter file.
 */
#define VALSTR "var"
#define VALLIST "list"

/*!
 * Initialize the hashtable.
 */
void pdata_initialize( pdata *pd){
  int i;
  pd->array_size = 0;
  for( i=0; i<MAX_ARRAY_SIZE; i++){
    pd->array[i].key[0] = '\0';
  }
}

/*!
 * skip over comment blocks.
 *
 * This function skips over comment blocks in the file. A comment is
 * defined as any line where the first non-whitespace character is
 * `#`.  After the loop the stream pointer points to the first
 * non-whitespace, non-`#` character. 
 */
static int skip_comment( FILE *file){
  char c;
  while( 1 ){
    do{
      c = getc( file);
    }while( isspace( c) && c != EOF );
    switch( c ){
    case '#':
      while( c != '\n' && c != EOF ){
	c = getc( file);
      }
      if( c == EOF ){
	return 1;
      }
      break;
    case EOF:
      return 1;
    default:
      ungetc( c, file);
      return 0;
    }
  }
}

/*!
 * Fill up the pdata data structure from a parameter file.
 *
 * This function loops through the parameter file filling up the pdata
 * data structure. The loop runs by skipping comment blocks. If a
 * comment block (or whitespace) ends the file, it breaks the
 * loop. Inside the loop it reads strings. After any comment blocks
 * the first string should be either the strings `var` or `varlist`
 *
 */
int pdata_read_file( pdata *pd, FILE *file){
  int i, n;
  char str[10], vname[80], vname_iter[87], vval[80];

  /* This like loops over the comment blocks */
  while( skip_comment( file) == 0 ){
    /* we read the next string (limited to 9 characters) */
    n = fscanf( file, "%9s", str);
    /* If we have a read error, we fail with a format error */
    if( strlen( str) == 9 || n != 1 ){
      return PDATA_FORMAT;
    }
    /* If we read a the value keyword read the name and value */
    if( strcmp( str, VALSTR) == 0){
      skip_comment( file);
      n = fscanf( file, "%79s", vname);
      if( strlen( vname) == 79 || n != 1 ){
	      return PDATA_FORMAT;
      }
      skip_comment( file);
      n = fscanf( file, "%79s", vval);
      if( strlen( vval) == 79 || n != 1 ){
	      return PDATA_FORMAT;
      }
      /* Insert into the hash table to lookup by name */
      if( hashtab_insert( pd, vname, vval) == HT_MEM_ERROR ){
	      return PDATA_MEM_ERROR;
      }
      /* If we read the list keyword read the name length */
    }else if( strcmp( str, VALLIST) == 0){
      skip_comment( file);
      n = fscanf( file, "%73s", vname);
      if( strlen( vname) == 73 || n != 1 ){
	      return PDATA_FORMAT;
      }
      i = 0;
      while( 1 ){
        skip_comment( file);
        n = fscanf( file, "%79s", vval);
        if( n == EOF ){
          break;
        }
        if( strlen( vval) == 79 || n != 1 ){
          return PDATA_FORMAT;
        }
        if( strcmp( vval, VALSTR) == 0 ||
           strcmp( vval, VALLIST) == 0
           ){
          fseek( file, -strlen( vval), SEEK_CUR);
          break;
        }
        sprintf( vname_iter, "%s___%03d", vname, i);
        if( hashtab_insert( pd, vname_iter, vval) == HT_MEM_ERROR ){
          return PDATA_MEM_ERROR;
        }
        i++;
      }
    }else{
      return PDATA_FORMAT;
    }
  }
  return PDATA_SUCCESS;
}

/*!
 * Get a single double or integer or string from pdata.
 *
 * These functions are simple table lookups for the key `name`,
 * then copy the string into the appropriate pointer.
 */
int pdata_get_var_d( pdata *pd, const char *name, double *val){
  char str[80];
  if( hashtab_find( pd, name, str) == HT_FAILURE){
    return PDATA_FAILURE;
  };
  *val = atof( str);
  return PDATA_SUCCESS;
}
int pdata_get_var_i( pdata *pd, const char *name, int *val){
  char str[80];
  if( hashtab_find( pd, name, str) == HT_FAILURE){
    return PDATA_FAILURE;
  };
  *val = atoi( str);
  return PDATA_SUCCESS;
}
int pdata_get_var_s( pdata *pd, const char *name, char *val){
  char str[80];
  if( hashtab_find( pd, name, str) == HT_FAILURE){
    return PDATA_FAILURE;
  };
  strcpy( val, str);
  return PDATA_SUCCESS;
}

/*!
 * Get the number of elements stored under the key `name`.
 */
int pdata_array_length( pdata *pd, const char *name){
  int i=0;
  char name_iter[87];
  char str[80];
  sprintf( name_iter, "%s___%03d", name, i);
  while( hashtab_find( pd, name_iter, str) == HT_SUCCESS ){
    i++;
    sprintf( name_iter, "%s___%03d", name, i);
  };
  return i;
}

/*!
 * Get an array of double/integers/strings from pdata.
 *
 * These functions do a series of table lookups and fill up and array
 * pointed to by the appropriate pointer.
 */
int pdata_get_array_d( pdata *pd, const char *name, double *val){
  int i=0;
  char name_iter[87];
  char str[80];
  sprintf( name_iter, "%s___%03d", name, i);
  while( hashtab_find( pd, name_iter, str) == HT_SUCCESS ){
    val[i] = atof( str);
    i++;
    sprintf( name_iter, "%s___%03d", name, i);
  }
  return i;
}
int pdata_get_array_i( pdata *pd, const char *name, int *val){
  int i=0;
  char name_iter[87];
  char str[80];
  sprintf( name_iter, "%s___%03d", name, i);
  while( hashtab_find( pd, name_iter, str) == HT_SUCCESS ){
    val[i] = atoi( str);
    i++;
    sprintf( name_iter, "%s___%03d", name, i);
  }
  return i;
}
int pdata_get_array_s( pdata *pd, const char *name, char **val){
  int i=0;
  char name_iter[80];
  char str[80];
  sprintf( name_iter, "%s___%03d", name, i);
  while( hashtab_find( pd, name_iter, str) == HT_SUCCESS ){
    strcpy( val[i], str);
    i++;
    sprintf( name_iter, "%s___%03d", name, i);
  }
  return i;
}

/*!
 * Get an element from an array of double/integers/strings from pdata.
 *
 * These functions do a series of table lookups and fill up and array
 * pointed to by the appropriate pointer.
 */
int pdata_get_element_d( pdata *pd, const char *name, 
                        unsigned int i, double *val){
  char name_iter[80];
  char str[80];
  sprintf( name_iter, "%s___%03d", name, i);

  if( hashtab_find( pd, name_iter, str) == HT_SUCCESS){
    *val = atof( str);
    return PDATA_SUCCESS;
  }else{
    return PDATA_FAILURE;
  }
}
int pdata_get_element_i( pdata *pd, const char *name, 
                        unsigned int i, int *val){
  char name_iter[80];
  char str[80];
  sprintf( name_iter, "%s___%03d", name, i);

  if( hashtab_find( pd, name_iter, str) == HT_SUCCESS){
    *val = atoi( str);
    return PDATA_SUCCESS;
  }else{
    return PDATA_FAILURE;
  }
}
int pdata_get_element_s( pdata *pd, const char *name, 
                        unsigned int i, char *val){
  char name_iter[80];
  char str[80];
  sprintf( name_iter, "%s___%03d", name, i);

  if( hashtab_find( pd, name_iter, str) == HT_SUCCESS){
    strcpy( val, str);
    return PDATA_SUCCESS;
  }else{
    return PDATA_FAILURE;
  }
}

/*!
 * Read in multiple data of various types.
 *
 * Use the `stdarg.h` macros to make a function that can take a
 * variable number of arguments. The arguments should always come in
 * pairs `data_type` (INT_T, DOUBLE_T, STRING_T), then a pointer to
 * that data type. This can be used to write simple wrappers to fill
 * up custom defined structs.
 */
int pdata_get_list( pdata *pd, const char *name, int num, ...){
  int i, *j;
  double *d;
  char *s;
  char name_iter[80], str[80];
  va_list vl;
 
  /* make sure the last argument is there */
  sprintf( name_iter, "%s___%03d", name, num-1);
  if( hashtab_find( pd, name_iter, str) == HT_FAILURE){
    return PDATA_FAILURE;
  }

  /* Go through the arguments */
  va_start( vl, num);
  for( i=0; i<num; i++){
    sprintf( name_iter, "%s___%03d", name, i);
    if( hashtab_find( pd, name_iter, str) == HT_FAILURE){
      return i;
    };
    switch( va_arg( vl, pdata_type)){
    case DOUBLE_T:
      d = va_arg( vl, double*);
      *d = atof( str);
      break;
    case INT_T:
      j = va_arg( vl, int*);
      *j = atoi( str);
      break;
    case STRING_T:
      s = va_arg( vl, char*);
      strcpy( s, str);
      break;
    default:
      break;
    }
  }
  va_end( vl);

  return PDATA_SUCCESS;
}
