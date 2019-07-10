#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "polymer.h"
#include "energy.h"
#include "params.h"
#include <hdf5.h>

#define ND 4
#define NI 3
int read_parameter_file( simulation* sim, 
  const char* filename, unsigned int n){
  // Variable Declarations
  FILE *f;
  int nd[ND], sd[ND], id[ND];
  double temp_d;
  double *var_d[ND];
  int ni[NI], si[NI], ii[NI];
  int temp_i;
  int *var_i[NI];
  int k, tmp_n, n_tot;
  // Initialize params data structure with all zeros */
  pdata pd = { {{0, "\0", "\0"}}, 0};
  // List out the names of the parameters
  const char *var_d_names[] = {
    "k",
    "a",
    "ep",
    "sig"
  };
  const char *var_i_names[] = {
    "seed",
    "Nm",
    "num_steps"
  };
  //List out the poniters
  en_params *enp = (en_params *) sim->enp;
  var_d[0] = &(enp->k);
  var_d[1] = &(enp->a);
  var_d[2] = &(enp->ep);
  var_d[3] = &(enp->sig);
  var_i[0] = &(sim->seed);
  var_i[1] = &(sim->Nm);
  var_i[2] = &(sim->num_steps);

  //Default values
  enp->k = 30.0;
  enp->a = 1.5;
  enp->ep = 1.0;
  enp->sig = 1.0;
  sim->seed = (unsigned int) time(NULL);
  sim->Nm = 256;
  sim->num_steps = 1000;

  // Open the file
  f = fopen( filename, "r");
  if( f == NULL ){
    fprintf( stderr, "Could not open the parameter file %s.\n", filename);
    return 1;
  }
  // Read in the parameter file
  if( pdata_read_file( &pd, f) == PDATA_FORMAT ){
    fprintf( stderr, "Could not parse the parameter file %s.\n", filename);
    return 2;
  }
  fclose(f);

  // Find the number of elements for each parameter
  for( k=0; k<ND; k++){
    nd[k] = 0;
    sd[k] = pdata_get_var_d( &pd, var_d_names[k], &temp_d);
    if( sd[k] == PDATA_SUCCESS ){
      nd[k] = 1;
    }else{
      nd[k] = pdata_array_length( &pd, var_d_names[k]);
    }
  }
  for( k=0; k<NI; k++){
    ni[k] = 0;
    si[k] = pdata_get_var_i( &pd, var_i_names[k], &temp_i);
    if( si[k] == PDATA_SUCCESS){
      ni[k] = 1;
    }else{
      ni[k] = pdata_array_length( &pd, var_i_names[k]);
    }
  }
  // From that get the total possible combinations
  n_tot = 1;
  for( k=0; k<ND; k++){
  	// fprintf( stdout, "nd[%d] = %d \n", k, nd[k]);
    if( nd[k] != 0 ){
      n_tot *= nd[k];
    }
  }
  for( k=0; k<NI; k++){
  	// fprintf( stdout, "ni[%d] = %d \n", k, ni[k]);
    if( ni[k] != 0 ){
      n_tot *= ni[k];
    }
  }
  // Check that n is less than n_tot
  if( n >= n_tot){
    fprintf( stderr, "Simulation number of out bounds\n");
    fprintf( stderr, "simulation_number: %d\n", n);
    fprintf( stderr, "total number of simulations: %d\n", n_tot);
    fprintf( stderr, "simulation number must be at least one less than\n");
    fprintf( stderr, " the total number of simulations.\n");
    return 3;
  }
  // Find the index for each parameter
  tmp_n = n;
  for( k=0; k<ND; k++){
    if( nd[k] != 0 ){
      id[k] = tmp_n%nd[k];
      tmp_n /= nd[k];
    }else{
      id[k] = -1;
    }
  }
  for( k=0; k<NI; k++){
    if( ni[k] != 0 ){
      ii[k] = tmp_n%ni[k];
      tmp_n /= ni[k];
    }else{
      ii[k] = -1;
    }
  }
  // Read each parameter
  for( k=0; k<ND; k++){
    if( nd[k] != 0){
      if( sd[k] == PDATA_SUCCESS){
        pdata_get_var_d( &pd, var_d_names[k], var_d[k]);
      }else{
        pdata_get_element_d( &pd, var_d_names[k], id[k], var_d[k]);
      }
    }
  }
  for( k=0; k<NI; k++){
    if( ni[k] != 0){
      if( si[k] == PDATA_SUCCESS){
        pdata_get_var_i( &pd, var_i_names[k], var_i[k]);
      }else{
        pdata_get_element_i( &pd, var_i_names[k], ii[k], var_i[k]);
      }
    }
  }

  return 0;
}

struct file_handle{
	char filename[100];
	hid_t outfile;
	hsize_t current_size;
};

file_handle* open_output_file( simulation *sim, const char *filename,
							   uint16 num_slices){
	hid_t outfile;
	hid_t space, data;
	hsize_t dims[3];

	// Create the file
	outfile = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// Create the position dataset in the file
	dims[0] = num_slices;
	dims[1] = sim->Nm,
	dims[2] = 3;

	space = H5Screate_simple( 3, dims, NULL);
    data = H5Dcreate( outfile, "positions", H5T_IEEE_F64LE, 
                      space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return NULL;
	}
	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return NULL;		
	}

	dims[0] = num_slices;
	dims[1] = sim->Nm-1;
	dims[2] = 2;
	
	space = H5Screate_simple( 3, dims, NULL);
    data = H5Dcreate( outfile, "bonds", H5T_STD_U16LE, space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return NULL;
	}
	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return NULL;	
	}

	// Create the kdtree dataset in the file
	dims[0] = num_slices;
	dims[1] = sim->Nm;
	dims[2] = 4;

	space = H5Screate_simple( 3, dims, NULL);
    data = H5Dcreate( outfile, "kdtree", H5T_STD_U16LE, space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return NULL;
	}
	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return NULL;	
	}

	// Create the energy dataset
	dims[0] = num_slices;

	space = H5Screate_simple( 1, dims, NULL);
	data = H5Dcreate( outfile, "energy", H5T_IEEE_F64LE, space,
					  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return NULL;
	}
	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return NULL;	
	}

	file_handle *ofile = (file_handle *) malloc(sizeof(file_handle));
	strcpy( ofile->filename, filename);
	ofile->outfile = outfile;
	ofile->current_size = 0;
	return ofile;
}

int close_output_file( file_handle *ofile){
	if( H5Fclose( ofile->outfile)<0 ){
		fprintf(stderr, "Error closing HDF5 file/n");
		return -1;
	}
	free(ofile);
	return 0;
}

int write_att( file_handle *ofile, char *name, enum datatype type, void *val){
	hid_t file, data, space, attr;

	double *f_val;
	uint16 *u_val;

	file = ofile->outfile;

	data = H5Dopen2(file, "/positions", H5P_DEFAULT);
	if( data < 0 ){
		fprintf(stderr, "Error opening HDF5 dataset/n");
		return -1;
	}

	space = H5Screate(H5S_SCALAR);

	if( type == T_DOUBLE){
		attr = H5Acreate2( data, name, H5T_IEEE_F64LE, space,
					   H5P_DEFAULT, H5P_DEFAULT);
		f_val = (double *) val;
		H5Awrite(attr, H5T_NATIVE_DOUBLE, f_val);
	}else{
		attr = H5Acreate2( data, name, H5T_STD_U16LE, space,
					   H5P_DEFAULT, H5P_DEFAULT);
		u_val = (uint16 *) val;
		H5Awrite(attr, H5T_NATIVE_UINT_FAST16, u_val);
	}

	if( H5Sclose( space ) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Aclose( attr ) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return -1;
	}
	return 0;
}

int write_2D_slab( file_handle *ofile, const char *name, int set_num, 
				   enum datatype type, void *array){
	hid_t file, data, space, memspace;

	hsize_t dims[3], memdims[2];
	hsize_t maxdims[3];

	double *f_array;
	uint16 *u_array;

	file = ofile->outfile;

	data = H5Dopen2(file, name, H5P_DEFAULT);
	if( data < 0 ){
		fprintf(stderr, "Error opening HDF5 dataset/n");
		return -1;
	}

	space = H5Dget_space( data);
	if( space < 0 ){
		fprintf(stderr, "Error opening HDF5 data space/n");
		return -1;
	}

	// rank = H5Sget_simple_extent_ndims( space);
	H5Sget_simple_extent_dims(space, dims, maxdims);

	memdims[0] = dims[1];
	memdims[1] = dims[2];
	memspace = H5Screate_simple( 2, memdims, NULL);	

	hsize_t offset[3], count[3];
	hsize_t stride[] = {1,1,1};
	hsize_t block[] = {1,1,1};
	offset[0] = set_num;
	offset[1] = 0;
	offset[2] = 0;
	count[0] = 1;
	count[1] = dims[1];
	count[2] = dims[2];

	H5Sselect_hyperslab( space, H5S_SELECT_SET,
						 offset, stride, count, block);

	if( type == T_DOUBLE){
		f_array = (double *) array;
		H5Dwrite( data, H5T_NATIVE_DOUBLE, memspace, 
				  space, H5P_DEFAULT, f_array);
	}else{
		u_array = (uint16 *) array;
		H5Dwrite( data, H5T_NATIVE_UINT_FAST16, memspace, 
				  space, H5P_DEFAULT, u_array);
	}

	if( H5Sclose( memspace) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return -1;
	}
	return 0;
}

int write_val( file_handle *ofile, char *name, int set_num, 
			   enum datatype type, void *val){
	hid_t file, data, space, memspace;

	hsize_t pos = set_num;

	double *f_val;
	uint16 *u_val;

	file = ofile->outfile;

	data = H5Dopen2(file, name, H5P_DEFAULT);
	if( data < 0 ){
		fprintf(stderr, "Error opening HDF5 dataset/n");
		return -1;
	}

	space = H5Dget_space( data);
	if( space < 0 ){
		fprintf(stderr, "Error opening HDF5 data space/n");
		return -1;
	}

	memspace = H5Screate(H5S_SCALAR);

	H5Sselect_elements( space, H5S_SELECT_SET, 1, &pos);

	if( type == T_DOUBLE){
		f_val = (double *) val;
		H5Dwrite( data, H5T_NATIVE_DOUBLE, memspace, 
				  space, H5P_DEFAULT, f_val);
	}else{
		u_val = (uint16 *) val;
		H5Dwrite( data, H5T_NATIVE_UINT_FAST16, memspace, 
				  space, H5P_DEFAULT, u_val);
	}

	if( H5Sclose( memspace) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return -1;
	}
	return 0;
}

void write_parameters( file_handle *ofile, simulation *sim){
	write_att( ofile, "seed", T_UINT, &sim->seed);
	write_att( ofile, "num_monomers", T_UINT, &sim->Nm);
	write_att( ofile, "num_steps", T_UINT, &sim->num_steps);

	en_params *enp = (en_params *) sim->enp;
	write_att( ofile, "FENE_k", T_DOUBLE, &enp->k);
	write_att( ofile, "FENE_a", T_DOUBLE, &enp->a);
	write_att( ofile, "TLJ_ep", T_DOUBLE, &enp->ep);
	write_att( ofile, "TLJ_sig", T_DOUBLE, &enp->sig);

	write_att( ofile, "r_max", T_DOUBLE, &sim->lim.r_max);
	write_att( ofile, "r_min", T_DOUBLE, &sim->lim.r_min);
	write_att( ofile, "r_cutoff", T_DOUBLE, &sim->lim.r_co);
}

void write_polymer_state( file_handle *ofile, simulation *sim, double en){
	uint16 *kdarray = (uint16 *) malloc( 4*sim->Nm*sizeof(uint16));
	kd_tree_to_array(sim->kdtree, kdarray);

	write_2D_slab( ofile, "/positions", ofile->current_size, T_DOUBLE, sim->poly->positions);
	write_2D_slab( ofile, "/bonds", ofile->current_size, T_UINT, sim->poly->bonds);
	write_2D_slab( ofile, "/kdtree", ofile->current_size, T_UINT, kdarray);
	write_val( ofile, "/energy", ofile->current_size, T_DOUBLE, &en);

	ofile->current_size++;
	free(kdarray);
}

// int write_kdtree_data( file_handle *ofile, uint16 *kdarray){
// 	hid_t file, data, space;

// 	file = ofile->outfile;
// 	data = H5Dopen2(file, "/kdtree", H5P_DEFAULT);
// 	if( data < 0 ){
// 		fprintf(stderr, "Error opening HDF5 dataset/n");
// 		return -1;
// 	}
// 	space = H5Dget_space( data);
// 	if( space < 0 ){
// 		fprintf(stderr, "Error opening HDF5 data space/n");
// 		return -1;
// 	}

// 	H5Dwrite( data, H5T_NATIVE_UINT_FAST16, H5S_ALL, space, H5P_DEFAULT, kdarray);

// 	if( H5Sclose( space) < 0 ){
// 		fprintf(stderr, "Error closing HDF5 data space./n");
// 		return -1;
// 	}
// 	if( H5Dclose( data) < 0 ){
// 		fprintf(stderr, "Error closing HDF5 dataset/n");
// 		return -1;
// 	}
// 	return 0;
// }


// int write_bond_data( file_handle *ofile, uint16 *bonds){
// 	hid_t file, data, space;

// 	// int rank;
// 	// hsize_t dims[10];
// 	// hsize_t maxdims[10];

// 	file = ofile->outfile;
// 	data = H5Dopen2(file, "/bonds", H5P_DEFAULT);
// 	if( data < 0 ){
// 		fprintf(stderr, "Error opening HDF5 dataset/n");
// 		return -1;
// 	}
// 	space = H5Dget_space( data);
// 	if( space < 0 ){
// 		fprintf(stderr, "Error opening HDF5 data space/n");
// 		return -1;
// 	}

// 	// rank = H5Sget_simple_extent_ndims( space);
// 	// fprintf( stdout, "The bonds dataset has rank %d.\n", rank);
// 	// H5Sget_simple_extent_dims(space, dims, maxdims);
// 	// fprintf( stdout, "    size %llu x %llu \n", dims[0], dims[1]);
// 	// fprintf( stdout, " maxsize %llu x %llu \n\n", maxdims[0], maxdims[1]);

// 	H5Dwrite( data, H5T_NATIVE_UINT_FAST16, H5S_ALL, space, H5P_DEFAULT, bonds);

// 	if( H5Sclose( space) < 0 ){
// 		fprintf(stderr, "Error closing HDF5 data space./n");
// 		return -1;
// 	}
// 	if( H5Dclose( data) < 0 ){
// 		fprintf(stderr, "Error closing HDF5 dataset/n");
// 		return -1;
// 	}
// 	return 0;
// }

//enum datatype{ T_UINT, T_DOUBLE};

// int write_int_slab( file_handle *ofile, char *name, 
// 					int set_num, uint16 *data){
// 	hid_t file, data, space, memspace;

// 	hsize_t dims[3], memdims[2];
// 	hsize_t maxdims[3];

// 	file = ofile->outfile;

// 	data = H5Dopen2(file, name, H5P_DEFAULT);
// 	if( data < 0 ){
// 		fprintf(stderr, "Error opening HDF5 dataset/n");
// 		return -1;
// 	}

// 	space = H5Dget_space( data);
// 	if( space < 0 ){
// 		fprintf(stderr, "Error opening HDF5 data space/n");
// 		return -1;
// 	}

// 	// rank = H5Sget_simple_extent_ndims( space);
// 	H5Sget_simple_extent_dims(space, dims, maxdims);

// 	memdims[0] = dims[1];
// 	memdims[1] = dims[2];
// 	memspace = H5Screate_simple( 2, memdims, NULL);	

// 	hsize_t offset[3], count[3];
// 	hsize_t stride[] = {1,1,1};
// 	hsize_t block[] = {1,1,1};
// 	offset[0] = set_num;
// 	offset[1] = 0;
// 	offset[2] = 0;
// 	count[0] = 1;
// 	count[1] = dims[1];
// 	count[2] = dims[2];

// 	H5Sselect_hyperslab( space, H5S_SELECT_SET,
// 						 offset, stride, count, block);

// 	H5Dwrite( data, H5T_NATIVE_UINT_FAST16, memspace, space, H5P_DEFAULT, data);

// 	if( H5Sclose( space) < 0 ){
// 		fprintf(stderr, "Error closing HDF5 data space./n");
// 		return -1;
// 	}

// 	if( H5Dclose( data) < 0 ){
// 		fprintf(stderr, "Error closing HDF5 dataset/n");
// 		return -1;
// 	}
// 	return 0;
// }

int write_position_data( file_handle *ofile, double *positions){
	hid_t file, data, space, memspace;

//	int rank;
	hsize_t dims[3], memdims[2];
	hsize_t maxdims[3];

	file = ofile->outfile;
	data = H5Dopen2(file, "/positions", H5P_DEFAULT);
	if( data < 0 ){
		fprintf(stderr, "Error opening HDF5 dataset/n");
		return -1;
	}

	space = H5Dget_space( data);
	if( space < 0 ){
		fprintf(stderr, "Error opening HDF5 data space/n");
		return -1;
	}

	// rank = H5Sget_simple_extent_ndims( space);
	H5Sget_simple_extent_dims(space, dims, maxdims);

	memdims[0] = dims[1];
	memdims[1] = dims[2];
	memspace = H5Screate_simple( 2, memdims, NULL);

	hsize_t offset[3], count[3];
	hsize_t stride[] = {1,1,1};
	hsize_t block[] = {1,1,1};
	offset[0] = ofile->current_size;
	offset[1] = 0;
	offset[2] = 0;
	count[0] = 1;
	count[1] = dims[1];
	count[2] = dims[2];

	H5Sselect_hyperslab( space, H5S_SELECT_SET,
						 offset, stride, count, block);

	H5Dwrite( data, H5T_NATIVE_DOUBLE, memspace, space, H5P_DEFAULT, positions);

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return -1;
	}
	return 0;
}

int write_kdtree_data( file_handle *ofile, uint16 *kdarray){
	hid_t file, data, space, memspace;

//	int rank;
	hsize_t dims[3], memdims[2];
	hsize_t maxdims[3];

	file = ofile->outfile;
	data = H5Dopen2(file, "/kdtree", H5P_DEFAULT);
	if( data < 0 ){
		fprintf(stderr, "Error opening HDF5 dataset/n");
		return -1;
	}

	space = H5Dget_space( data);
	if( space < 0 ){
		fprintf(stderr, "Error opening HDF5 data space/n");
		return -1;
	}

	// rank = H5Sget_simple_extent_ndims( space);
	H5Sget_simple_extent_dims(space, dims, maxdims);

	memdims[0] = dims[1];
	memdims[1] = dims[2];
	memspace = H5Screate_simple( 2, memdims, NULL);

	hsize_t offset[3], count[3];
	hsize_t stride[] = {1,1,1};
	hsize_t block[] = {1,1,1};
	offset[0] = ofile->current_size;
//	ofile->current_size++;
	offset[1] = 0;
	offset[2] = 0;
	count[0] = 1;
	count[1] = dims[1];
	count[2] = dims[2];

	H5Sselect_hyperslab( space, H5S_SELECT_SET,
						 offset, stride, count, block);

	H5Dwrite( data, H5T_NATIVE_UINT_FAST16, memspace, space, H5P_DEFAULT, kdarray);

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return -1;
	}
	return 0;
}

int write_bond_data( file_handle *ofile, uint16 *bonds){
	hid_t file, data, space, memspace;

//	int rank;
	hsize_t dims[3], memdims[2];
	hsize_t maxdims[3];

	file = ofile->outfile;
	data = H5Dopen2(file, "/bonds", H5P_DEFAULT);
	if( data < 0 ){
		fprintf(stderr, "Error opening HDF5 dataset/n");
		return -1;
	}

	space = H5Dget_space( data);
	if( space < 0 ){
		fprintf(stderr, "Error opening HDF5 data space/n");
		return -1;
	}

	// rank = H5Sget_simple_extent_ndims( space);
	H5Sget_simple_extent_dims(space, dims, maxdims);

	memdims[0] = dims[1];
	memdims[1] = dims[2];
	memspace = H5Screate_simple( 2, memdims, NULL);

	hsize_t offset[3], count[3];
	hsize_t stride[] = {1,1,1};
	hsize_t block[] = {1,1,1};
	offset[0] = ofile->current_size;
//	ofile->current_size++;
	offset[1] = 0;
	offset[2] = 0;
	count[0] = 1;
	count[1] = dims[1];
	count[2] = dims[2];

	H5Sselect_hyperslab( space, H5S_SELECT_SET,
						 offset, stride, count, block);

	H5Dwrite( data, H5T_NATIVE_UINT_FAST16, memspace, space, H5P_DEFAULT, bonds);

	if( H5Sclose( space) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return -1;
	}

	if( H5Dclose( data) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return -1;
	}
	return 0;
}
