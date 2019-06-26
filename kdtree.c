#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "polymer.h"
#include "vec3.h"

kd_tree* kd_tree_malloc( uint16 max_size){
	// uint16 inactive_max = num_monomers;
	kd_tree *kdt = (kd_tree *) malloc( sizeof(kd_tree));
	// kdt->inactive = (double *) malloc( DIM*inactive_max*sizeof(double));
	// kdt->inactive_size = 0;
	kdt->node_array = (kd_node *) malloc( max_size*sizeof(kd_node));
	kdt->size = 0;
	kdt->positions = NULL;
	kdt->root = NULL;
	kdt->max_depth = 0;

	return kdt;
}

void kd_tree_free( kd_tree *kdt){
	free(kdt->node_array);
	// free(kdt->inactive);
	free(kdt);
}

///////////////// Median of Medians algorithm


// This function swaps two uint16 values
static inline void swap( uint16 *a, uint16 *b){
	uint16 c = *a;
	*a = *b;
	*b = c;
}

// This function sorts a small list
// and returns the position of the median
uint16 median5( double *pos, uint16 depth,
				uint16 size, uint16 *idx){
	double val0, val1;
	uint16 j, i;
	uint16 ax = depth%DIM;

	for( i=1; i<size; i++){
		j = i;
		val0 = pos[DIM*idx[j-1]+ax];
		val1 = pos[DIM*idx[j]+ax];
		// while( j>0 && val0 > val1){
		// 	swap( &idx[j-1], &idx[j]);
		// 	j--;
		// 	val0 = pos[DIM*idx[j-1]+ax];
		// 	val1 = pos[DIM*idx[j]+ax];
		// }
		while( val0 > val1 ){
			swap( &idx[j-1], &idx[j]);
			j--;
			if( j == 0 ){
				break;
			}
			val0 = pos[DIM*idx[j-1]+ax];
			val1 = pos[DIM*idx[j]+ax];
		}
	}

	return (size-1)/2;
}

// This fuction recursively an approximate median
uint16 median_of_medians( double *pos, uint16 depth,
						  uint16 size, uint16 *idx){
	uint16 i, j, num_groups, final_size;
	uint16 ax = depth%DIM;
	// if less than 5, use median5
	if( size <= 5 ){
		return median5( pos, ax, size, idx);
	}
	// if greater than five, break into groups of 5
	num_groups = (size-1)/5+1;
	for( i=0; i<num_groups-1; i++){
		// find median of each one
		// place those medians at the front
		j = 5*i + median5( pos, ax, 5, &idx[5*i]);
		swap( &idx[i], &idx[j]);
	}
	// The final group might not have 5 elements in it
	i = num_groups-1;
	final_size = size - 5*i;
	j = 5*i + median5( pos, ax, final_size, &idx[5*i]);
	swap( &idx[i], &idx[j]);
	// Take the median of the new group
	return median_of_medians( pos, ax, num_groups, idx);
}

// This function partitions the index list `idx`
// by moving all those less than the pivot to the
// front of the `idx` array, moving all those to
// the right of the pivot to the end, and returning
// the new positions of the pivot
uint16 partition( double *pos, uint16 depth,
			      uint16 size, uint16 *idx, uint16 pivot){
	uint16 i, si=0;
	uint16 ax = depth%DIM;
	double val;
	double pval = pos[DIM*idx[pivot]+ax];

	swap( &idx[pivot], &idx[size-1]);

	for( i=0; i<size-1; i++){
		val = pos[DIM*idx[i]+ax];

		if( val<pval ){
			swap( &idx[si], &idx[i]);
			si++;
		}
	}
	swap( &idx[size-1], &idx[si]);
	return si;
}

uint16 qselect( double *pos, uint16 depth,
				uint16 size, uint16 *idx, 
				uint16 offset, uint16 kth){
	uint16 pivot, new_size, new_offset, new_kth, *new_idx;
	pivot = median_of_medians( pos, depth, size, idx);

	pivot = partition( pos, depth, size, idx, pivot);

	if( kth == pivot ){
		// fprintf( stdout, " piv = %.2f\n", pos[DIM*idx[pivot]+ax]);
		return offset + pivot;
	}
	if( kth < pivot ){
		new_offset = offset;
		new_size = pivot;
		new_idx = idx;
		new_kth = kth;
	}else{
		new_offset = offset+pivot+1;
		new_size = size-pivot-1;
		new_idx = &idx[pivot+1];
		new_kth = kth-pivot-1;
	}
	return qselect( pos, depth, new_size, new_idx, new_offset, new_kth);
}

uint16 qmedian( double *pos, uint16 depth,
				uint16 size, uint16 *idx){
	// uint16 pivot = qselect( pos, ax, size, idx, 0, (size-1)/2);
	// fprintf( stdout, " piv = %.2f\n", pos[DIM*idx[pivot]+ax]);

	return qselect( pos, depth, size, idx, 0, (size-1)/2);
}

void kd_node_fill( kd_tree *kdt, 
				   kd_node *node, uint16 depth, 
				   uint16 size, uint16 *idx){
	uint16 mi, lsize, rsize, *lidx, *ridx;

	mi = qmedian( kdt->positions, depth, size, idx);

	*node = kd_new_node(idx[mi], depth);
	kdt->max_depth = max(kdt->max_depth, depth);

	if( mi != 0){
		lsize = mi;
		lidx = idx;
		node->left = &kdt->node_array[kdt->size];
		kdt->size++;
		kd_node_fill( kdt, node->left, depth+1, lsize, lidx);
	}
	if( mi != size-1){
		rsize = size-mi-1;
		ridx = &idx[mi+1];
		node->right = &kdt->node_array[kdt->size];
		kdt->size++;
		kd_node_fill( kdt, node->right, depth+1, rsize, ridx);
	}
}

void kd_tree_fill( kd_tree *kdt, uint16 size, double *positions){
	kdt->max_depth = 0;
	uint16 *idx = (uint16 *) malloc( size*sizeof(uint16));
	uint16 i;

	for( i=0; i<size; i++){
		idx[i] = i;
	}

	kdt->positions = positions;

	kdt->root = &kdt->node_array[0];
	kdt->size = 1;
	kd_node_fill( kdt, kdt->root, 0, size, idx);

	free( idx);
}

void kd_tree_to_array( kd_tree *kdt, uint16 *array){
	kd_node node, *base;
	uint16 i;

	base = kdt->root;
	for( i=0; i<kdt->size; i++){
		node = kdt->node_array[i];
		array[4*i+0] = node.monomer;
		array[4*i+1] = (uint16) node.depth;
		array[4*i+2] = (uint16) (node.left - base);
		array[4*i+3] = (uint16) (node.right - base);
	}
}

// Find the node, starting with the current node,
// that has the maximum position along an axis
// while also returning the parent of the max node
kd_node* kd_max( kd_tree *kdt, kd_node *kdn,
				 uint16 ax, kd_node **parent){
	double c_max, b_max;
	kd_node *c_node, *b_node;
	kd_node *c_parent, *b_parent;
	// Make sure our axis is in the required range
	ax = ax%DIM;
	// Set the value of the starting node
	// as the maximum
	c_node = kdn;
	c_max = kdt->positions[DIM*c_node->monomer + ax];
	c_parent = *parent;
	// Search the left subtree
	// If the max from the left subtree is larger
	// replace the current value with left value
	if( kdn->left != NULL && ax != kdn->depth%DIM ){
		b_parent = kdn;
		b_node = kd_max( kdt, kdn->left, ax, &b_parent);
		b_max = kdt->positions[DIM*b_node->monomer + ax];
		if( b_max >= c_max){
			c_node = b_node;
			c_parent = b_parent;
			c_max = b_max;
		}
	}
	// Search the right subtree
	if( kdn->right != NULL){
		// If the node axis matches our search axis,
		// the right subtree MUST have a larger value
		if( ax == kdn->depth%DIM ){
			c_parent = kdn;
			c_node = kd_max( kdt, kdn->right, ax, &c_parent);
			c_max = kdt->positions[DIM*c_node->monomer + ax];
		/// Otherwise we compare and replace like above
		}else{
			b_parent = kdn;
			b_node = kd_max( kdt, kdn->right, ax, &b_parent);
			b_max = kdt->positions[DIM*b_node->monomer + ax];
			if( b_max >= c_max){
				c_node = b_node;
				c_parent = b_parent;
				c_max = b_max;
			}
		}
	}

	*parent = c_parent;
	return c_node;
}

kd_node* kd_min( kd_tree *kdt, kd_node *kdn,
				 uint16 ax, kd_node **parent){
	double c_min, b_min;
	kd_node *c_node, *b_node, *c_parent, *b_parent;
	// Make sure our axis is in the required range
	ax = ax%DIM;
	// Set the value of the starting node
	// as the maximum
	c_node = kdn;
	c_min = kdt->positions[DIM*c_node->monomer + ax];
	c_parent = *parent;
	// Search the left subtree
	// If the max from the left subtree is larger
	// replace the current value with left value
	if( kdn->right != NULL && ax != kdn->depth%DIM ){
		b_parent = kdn;
		b_node = kd_min( kdt, kdn->right, ax, &b_parent);
		b_min = kdt->positions[DIM*b_node->monomer + ax];
		if( b_min < c_min){
			c_node = b_node;
			c_parent = b_parent;
			c_min = b_min;
		}
	}
	// Search the right subtree
	if( kdn->left != NULL){
		// If the node axis matches our search axis,
		// the right subtree MUST have a larger value
		if( ax == kdn->depth%DIM ){
			c_parent = kdn;
			c_node = kd_min( kdt, kdn->left, ax, &c_parent);
			c_min = kdt->positions[DIM*c_node->monomer + ax];
		/// Otherwise we compare and replace like above
		}else{
			b_parent = kdn;
			b_node = kd_min( kdt, kdn->left, ax, &b_parent);
			b_min = kdt->positions[DIM*b_node->monomer + ax];
			if( b_min < c_min){
				c_node = b_node;
				c_parent = b_parent;
				c_min = b_min;
			}
		}
	}

	*parent = c_parent;
	return c_node;
} 

kd_node* kd_remove( kd_tree *kdt, kd_node *kdn, kd_node **parent){
	uint16 temp_mon;
	uint16 ax = kdn->depth%DIM;
	kd_node *node;

	if( kdn->left != NULL){
		*parent = kdn;
		node = kd_max(kdt, kdn->left, ax, parent);

		// fprintf( stdout, " %lu <--> %lu  a depth %lu <--> %lu \n",
		//  kdn->monomer, node->monomer, kdn->depth, node ->depth);

		temp_mon = kdn->monomer;
		kdn->monomer = node->monomer;
		node->monomer = temp_mon;
		return kd_remove( kdt, node, parent);
	}else if( kdn->right != NULL){
		*parent = kdn;
		node = kd_min(kdt, kdn->right, ax, parent);

		// fprintf( stdout, " %lu <--> %lu  a depth %lu <--> %lu \n",
		//  kdn->monomer, node->monomer, kdn->depth, node ->depth);

		temp_mon = kdn->monomer;
		kdn->monomer = node->monomer;
		node->monomer = temp_mon;
		return kd_remove( kdt, node, parent);
	}

	int i = 0;
	if( *parent != NULL){
		if( (*parent)->left == kdn){
			(*parent)->left = NULL;
			i += 1;
		}
		if( (*parent)->right == kdn){
			(*parent)->right = NULL;
			i += 1;
		}
	}
	return kdn;
}

/// Adds a node the the kdtree
/// Starts the search at kdn
/// Inserts the node new_node
void kd_add( kd_tree *kdt, kd_node *kdn, kd_node *new_node){
	double *node_pos, *new_pos;
	kd_node *node = kdn;
	uint16 ax;

	node_pos = &kdt->positions[DIM*node->monomer];
	new_pos = &kdt->positions[DIM*new_node->monomer];
	ax = node->depth%DIM;
	if( new_pos[ax] < node_pos[ax]){
		if( node->left == NULL ){
			new_node->depth = node->depth +1;
			node->left = new_node;
			kdt->max_depth = max( kdt->max_depth, new_node->depth);
		}else{
			return kd_add(kdt, node->left, new_node);
		}
	}else{
		if( node->right == NULL ){
			new_node->depth = node->depth +1;
			node->right = new_node;
			kdt->max_depth = max( kdt->max_depth, new_node->depth);
		}else{
			return kd_add(kdt, node->right, new_node);
		}
	}
}

	// while( node->depth <= kdt->max_depth ){
	// 	node_pos = &kdt->positions[DIM*node->monomer];
	// 	ax = node->depth%DIM;
	// 	if( pos[ax] < node_pos[ax] ){
	// 		if( node->left == NULL ){
	// 			new_node->depth = node->depth + 1;
	// 			node->left = new_node;
	// 			kdt->max_depth = max( kdt->max_depth, new_node->depth);
	// 			return 0;
	// 		}else{
	// 			node = node->left;
	// 		}
	// 	}else{
	// 		if( node->right == NULL ){
	// 			new_node->depth = node->depth + 1;
	// 			node->right = new_node;
	// 			kdt->max_depth = max( kdt->max_depth, new_node->depth);
	// 			return 0;
	// 		}else{
	// 			node = node->right;
	// 		}
	// 	}
	// }
	// return -1;
// }

kd_node* kd_find( kd_tree *kdt, kd_node *kdn, uint16 i){
	double *ipos, *npos;
	uint16 ax;

	if( kdn->monomer == i ){
		return kdn;
	}
	ipos = &kdt->positions[DIM*i];
	npos = &kdt->positions[DIM*kdn->monomer];
	ax = kdn->depth%DIM;
	if( ipos[ax] < npos[ax]){
		if( kdn->left == NULL ){
			return NULL;
		}else{
			return kd_find( kdt, kdn->left, i);
		}
	}else{
		if( kdn->right == NULL ){
			return NULL;
		}else{
			return kd_find( kdt, kdn->right, i);
		}
	}
}

kd_node* kd_find2( kd_tree *kdt, kd_node *node, 
				   uint16 monomer, kd_node **parent){
	double *node_pos, *mon_pos;
	uint16 ax;

	if( node->monomer == monomer ){
		return node;
	}

	mon_pos = &kdt->positions[DIM*monomer];
	node_pos = &kdt->positions[DIM*node->monomer];
	ax = node->depth%DIM;
	*parent = node;
	if( mon_pos[ax] < node_pos[ax]){
		if( node->left == NULL ){
			fprintf( stdout, "Lost monomer %lu at depth %lu !\n", monomer, node->depth);
			return NULL;
		}else{
			return kd_find2( kdt, node->left, monomer, parent);
		}
	}else{
		if( node->right == NULL ){
			fprintf( stdout, "Lost monomer %lu at depth %lu !\n", monomer, node->depth);
			return NULL;
		}else{
			return kd_find2( kdt, node->right, monomer, parent);
		}
	}
}

kd_node* kd_comp( kd_tree *kdt, kd_node *node, uint monomer, 
				  double *new_pos, kd_node **parent){
	double *old_pos, *node_pos;
	uint16 ax;

	if( node->monomer == monomer ){
		return NULL;
	}

	old_pos = &kdt->positions[DIM*monomer];
	node_pos = &kdt->positions[DIM*node->monomer];
	ax = node->depth%DIM;
	if( (new_pos[ax] < node_pos[ax] && old_pos[ax] >= node_pos[ax]) ||
		(new_pos[ax] >= node_pos[ax] && old_pos[ax] < node_pos[ax])){
		return kd_find2( kdt, node, monomer, parent);
	}

	*parent = node;
	if( new_pos[ax] < node_pos[ax] && old_pos[ax] < node_pos[ax] ){
		if( node->left == NULL ){
			fprintf( stdout, "Lost Something!\n");
			return NULL;
		}else{
			return kd_comp( kdt, node->left, monomer, new_pos, parent);
		}
	}else{
		if( node->right == NULL ){
			fprintf( stdout, "Lost Something!\n");
			return NULL;
		}else{
			return kd_comp( kdt, node->right, monomer, new_pos, parent);
		}
	}
}


void kd_move( kd_tree *kdt, uint16 monomer, double *new_pos){
	double *pos;
	kd_node *node, *parent;
	
	parent = NULL;
	node = kd_find2( kdt, kdt->root, monomer, &parent);
	// fprintf( stdout, " move monomer %lu at depth %lu \n", node->monomer, node->depth);

	node = kd_remove( kdt, node, &parent);

	pos = &kdt->positions[DIM*node->monomer];
	pos[0] = new_pos[0];
	pos[1] = new_pos[1];
	pos[2] = new_pos[2];

	kd_add( kdt, kdt->root, node);
	// }
}

/// Moves a monomer with index `monomer` to a new positions `new_pos`
/// If the kd_tree is still valid after the move, we simply move the
/// monomer and leave the kdtree alone. If it is not valid, we remove
/// the node and re-add the new node
void kd_move2( kd_tree *kdt, uint16 monomer, double *new_pos){
	double *pos, *node_pos, *b_pos;
	kd_node *node;//, *parent;
	kd_node *b_node, *b_parent;
	uint16 ax, temp_mon;

	pos = &kdt->positions[DIM*monomer];
	// parent = NULL;
	node = kd_find( kdt, kdt->root, monomer);	
	// node = kdt->root;
	// while( node->monomer != monomer ){
	// 	node_pos = &kdt->positions[DIM*node->monomer];
	// 	ax = node->depth%DIM;
	// 	if( pos[ax] < node_pos[ax]){
	// 		parent = node;
	// 		node = node->left;
	// 	}else{
	// 		parent = node;
	// 		node = node->right;
	// 	}
	// }

	kd_node *t1, *t2;

	ax = node->depth%DIM;
	// b_parent = parent;
	if( new_pos[ax] < pos[ax]){
		// Move the monomer
		memcpy( pos, new_pos, DIM*sizeof(double));
		// FIX // FIX // FIX //
		if( node->left != NULL ){
			b_parent = node;
			b_node = kd_max( kdt, node->left, ax, &b_parent);
			b_pos = &kdt->positions[DIM*b_node->monomer];
			if( b_pos[ax] >= new_pos[ax]){
				// Swap the monomers
				fprintf( stdout, " %lu <--> %lu  a depth %lu <--> %lu \n",
				 node->monomer, b_node->monomer, node->depth, b_node ->depth);
				temp_mon = node->monomer;
				node->monomer = b_node->monomer;
				b_node->monomer = temp_mon;
				// remove the node we just swapped out
				t1 = b_node;
				b_node = kd_remove( kdt, b_node, &b_parent);
				if( b_node->monomer != monomer ){
					fprintf( stdout, "removed isn't correct!\n");
					temp_mon += 1;
				}
				t2 = kd_find( kdt, kdt->root, b_node->monomer);
				if( t2 != NULL ){
					fprintf( stdout, "didn't remove node!\n");
					temp_mon += 1;
				}
				// add the node back in
				kd_add(kdt, kdt->root, b_node);
				t2 = kd_find( kdt, kdt->root, b_node->monomer);
				if( t2 == NULL ){
					fprintf( stdout, "didn't add node correctly!\n");
					temp_mon += 1;
				}
			}
		}
	}else{
		// Move the monomer
		memcpy( pos, new_pos, DIM*sizeof(double));
		// FIX // FIX // FIX //
		if( node->right != NULL ){
			b_parent = node;
			b_node = kd_min( kdt, node->right, ax, &b_parent);
			b_pos = &kdt->positions[DIM*b_node->monomer];
			if( b_pos[ax] < new_pos[ax]){
				// Swap the monomers
				fprintf( stdout, " %lu <--> %lu  a depth %lu <--> %lu \n",
				 node->monomer, b_node->monomer, node->depth, b_node ->depth);
				temp_mon = node->monomer;
				node->monomer = b_node->monomer;
				b_node->monomer = temp_mon;
				// remove the node we just swapped out
				t1 = b_node;
				b_node = kd_remove( kdt, b_node, &b_parent);
				if( b_node->monomer != monomer ){
					fprintf( stdout, "removed isn't correct!\n");
					temp_mon += 1;
				}
				t2 = kd_find( kdt, kdt->root, b_node->monomer);
				if( t2 != NULL ){
					fprintf( stdout, "didn't remove node!\n");
					temp_mon += 1;
				}
				// add the node back in
				kd_add(kdt, kdt->root, b_node);
				t2 = kd_find( kdt, kdt->root, b_node->monomer);
				if( t2 == NULL ){
					fprintf( stdout, "didn't add node correctly!\n");
					temp_mon += 1;
				}
			}
		}
	}
}

int kd_check_left( kd_tree *kdt, kd_node *node, double *pos, uint16 ax){
	int t0, tl, tr;
	t0 = (kdt->positions[DIM*node->monomer+ax] < pos[ax]);
	if( node->left == NULL ){
		tl = 1;
	}else{
		tl = kd_check_left( kdt, node->left, pos, ax);
	}
	if( node->right == NULL ){
		tr = 1;
	}else{
		tr = kd_check_left( kdt, node->right, pos, ax);
	}
	return t0 && tl && tr;
}

int kd_check_right( kd_tree *kdt, kd_node *node, double *pos, uint16 ax){
	int t0, tl, tr;
	t0 = (kdt->positions[DIM*node->monomer+ax] >= pos[ax]);
	if( node->left == NULL ){
		tl = 1;
	}else{
		tl = kd_check_right( kdt, node->left, pos, ax);
	}
	if( node->right == NULL ){
		tr = 1;
	}else{
		tr = kd_check_right( kdt, node->right, pos, ax);
	}
	return t0 && tl && tr;
}

int kd_check_balance( kd_tree *kdt, kd_node *node){
	int tl, tr, bl, br;
	double *pos = &kdt->positions[DIM*node->monomer];
	uint16 ax = node->depth%DIM;
	if( node->left == NULL ){
		tl = 1;
		bl = 1; 
	}else{
		tl = kd_check_left( kdt, node->left, pos, ax);
		bl = kd_check_balance( kdt, node->left);
	}
	if( node->right == NULL ){
		tr = 1;
		br = 1;
	}else{
		tr = kd_check_right( kdt, node->right, pos, ax);
		br = kd_check_balance( kdt, node->right);
	}
	return tl && tr && bl && br;
}

// This searches for all monomers within a distance r from
// the positions pos. The list of monomers are written to 
// an array `neibors` of `size` elements.
void kd_ball_search( kd_tree *kdt, kd_node *kdn, 
					 double *pos, double r,
					 uint16 *neighbors, uint16 *size){
	kd_node *node = kdn;
	double *node_pos = &kdt->positions[DIM*node->monomer];
	uint16 ax = node->depth%DIM;
	if( dist2(pos,node_pos) < r*r ){
		neighbors[*size] = node->monomer;
		*size += 1;
	}
	if( pos[ax]-r < node_pos[ax] && node->left != NULL){
		kd_ball_search(kdt, node->left, pos, r, neighbors, size);
	}
	if( pos[ax]+r >= node_pos[ax] && node->right != NULL){
		kd_ball_search(kdt, node->right, pos, r, neighbors, size);
	}
}

// /***********************************
//  * Implement a stack
//  */
// typedef struct{
// 	uint16 *a;
// 	uint16 size;
// 	uint16 max;
// } stack;

// stack *stack_malloc( uint16 max){
// 	stack *s = (stack *) malloc( sizeof(stack));
// 	s->a = (uint16 *) malloc( max*sizeof(stack));
// 	s->size = 0;
// 	s->max = max;
// 	return s;
// }

// void stack_free( stack *s){
// 	free(s->a);
// 	free(s);
// }

// void stack_grow( stack *s){
// 	s = (stack *) realloc( s, 2*s->max);
// 	s->max *= 2;
// }

// void st_push( stack *s, uint16 i){
// 	if( s->size == s->max ){
// 		stack_grow(s);
// 	}
// 	s->a[s->size] = i;
// 	s->size++;
// }

// uint16 st_pop( stack *s){
// 	s->size--;
// 	return s->a[s->size];
// }