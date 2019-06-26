/**********************************************************************
 * Common header for Polymer program
 * /author Jef Wagner
 * /email wagnerj@union.edu
 * /date 2019-06-15
 * ********************************************************************
 * 
 */
#ifndef JW_POLYMER
#define JW_POLYMER

// Include the fixed size integer types and booleans
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

// Typedef our types from fast versions of given size integer types
typedef uint_fast8_t uint8;
typedef uint_fast16_t uint16;
typedef uint_fast32_t uint32;
typedef uint_fast64_t uint64;

// Define some common mathamtica constants
#define TWOPI 6.2831853071795864769
#define PI 3.1415926535897932385
#define EULER_E 2.7182818284590452354
#define ROOT3 1.7320508075688772935
#define PI_2 1.5707963267948966192
#define ROOT2 1.4142135623730950488
#define PI_3 1.0471975511965977462
#define ROOT3_2 0.86602540378443864676
#define PI_6 0.52359877559829887308
#define ONE_E 0.36787944117144232160
#define PI_180 0.017453292519943295769
// Standard min and max macros (Note: Don't use with ++ or -- notation!)
#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))


/***********************************************************************
 * Define a polymer type
 */
#define DIM 3 // Dimension of the space
#define NB_MAX 3 // Maximum number of bonds for a single monomer

typedef struct {
  double* positions; // position of monomers
  uint8* num_bonds; // number of bonds for each monomer
  uint16* list_bonds; // list of bonds for each monomer
  uint16* bonds; // a list of bonds
} polymer;

// Memory management
polymer* polymer_malloc( uint16 Nm);
void polymer_free(polymer *poly);

int polymer_linear(polymer *poly, uint16 Nm, uint64_t *rng);

// Convience function for finding linked monomers
static inline uint16 linked_monomer( uint16 *bond, uint16 monomer){
	return bond[0]==monomer?bond[1]:bond[0];
}

/***********************************************************************
 * Define structures for the k-d tree
 */

typedef struct kd_node{
	uint16 monomer; // index of monomer in the polymer struct
	uint16 depth; // split axis
	struct kd_node *left; // less than position of monomer along split axis
	struct kd_node *right; // greater than position of monomer along split axis
} kd_node;

// Convenience function for populating a node;
static inline kd_node kd_new_node( uint16 monomer, 
	                               uint16 depth){
	kd_node kdn = {monomer, depth, NULL, NULL};
	return kdn;
}

//
typedef struct{
	double *positions; // position of monomers (same pointer in polymer struct!)
	kd_node *node_array; // block of memory with node array
	uint16 size; // current size of kd_tree
	kd_node *root; // pointer to the root of the kd_tree
	uint16 max_depth;
} kd_tree;

// Memory management
kd_tree* kd_tree_malloc( uint16 max_size);
void kd_tree_free( kd_tree *kdt);
// Populating the kdtree from scratch
void kd_tree_fill( kd_tree *kdt, uint16 size, double *pos);
// Transforming from a tree into a 4xN array for output
void kd_tree_to_array( kd_tree *kdt, uint16 *array);

// kd_node* kd_max( kd_tree *kdt, kd_node *kdn,
// 				 uint16 ax, kd_node **parent);
// kd_node* kd_min( kd_tree *kdt, kd_node *kdn,
// 				 uint16 ax, kd_node **parent);
// kd_node* kd_remove( kd_tree *kdt, kd_node *kdn, kd_node **parent);
// int kd_add( kd_tree *kdt, kd_node *kdn, kd_node *new_node);
// kd_node* kd_find2( kd_tree *kdt, kd_node *node, 
// 				   uint16 monomer, kd_node **parent);
// int kd_check_balance( kd_tree *kdt, kd_node *node);

// kd_node* kd_remove( kd_tree *kdt, kd_node *kdn, kd_node **parent);
// void kd_add( kd_tree *kdt, kd_node *kdn, kd_node *new_node);


// Altering the kdtree by shifting a single position.
void kd_move( kd_tree *kdt, uint16 monomer, double *new_pos);

void kd_ball_search( kd_tree *kdt, kd_node *kdn, 
					 double *pos, double r,
					 uint16 *monomers, uint16 *size);

/********************************************************************
 * Define functions used in the montecarlo moves
 */
void rand1( uint64_t *rng, double r_min, double r_max, 
			double *p0, double *output);
void rand2( uint64_t *rng, double r_min, double r_max,
		    double *p0, double *p1, double *output);

/********************************************************************
 * Define the energy functional
 */
double r_min_calc( const void *params);
double r_max_calc( const void *params);
double r_co_calc( const void *params);

double E_bond( const double *x1, const double *x2, 
			   const void *params);

double E_pair( const double *x1, const double *x2, 
			   const void *params);


/********************************************************************
 * 
 */
typedef struct limits{
	double r_min;
	double r_max;
	double r_co;
} limits;

typedef struct {
  int seed;
  int Nm;
  int num_steps;
  void *enp;
  limits lim; 

  uint64_t rng[2];
  polymer *poly;
  kd_tree *kdtree;
} simulation;

/***********************************************************************
 * Define a input output structures and functions
 */
int read_parameter_file( simulation* sim, 
                         const char* filename, 
                         unsigned int n);

typedef struct file_handle file_handle;

file_handle* open_output_file( simulation *sim, 
							   const char *filename,
							   uint16 num_slices);
int close_output_file( file_handle *ofile);

int write_bond_data( file_handle *ofile, uint16 *bonds);
int write_kdtree_data( file_handle *ofile, uint16 *kdarray);
int write_position_data( file_handle *ofile, double *positions);


/********************************************************************
 * 
 */
double delta_en( simulation *sim, uint16 monomer, double *pos);
double energy( simulation *sim);

void gibbs_step( simulation *sim);

#endif