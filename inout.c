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

	// Create the file
	outfile = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// Create the position dataset in the file
	hid_t pspace, pdata;
	int prank = 3;
	hsize_t pdims[3] = {num_slices, sim->Nm, 3};

	pspace = H5Screate_simple( prank, pdims, NULL);
    pdata = H5Dcreate( outfile, "positions", 
                       H5T_IEEE_F64LE, pspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if( H5Sclose( pspace) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return NULL;
	}
	if( H5Dclose( pdata) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return NULL;		
	}

	// Create the bonds dataset in the file
	hid_t bspace, bdata;
	int brank = 3;
	hsize_t bdims[3] = {num_slices, sim->Nm-1, 2};

	bspace = H5Screate_simple( brank, bdims, NULL);
    bdata = H5Dcreate( outfile, "bonds", 
                       H5T_STD_U16LE, bspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if( H5Sclose( bspace) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return NULL;
	}
	if( H5Dclose( bdata) < 0 ){
		fprintf(stderr, "Error closing HDF5 dataset/n");
		return NULL;	
	}

	// Create the kdtree dataset in the file
	hid_t kspace, kdata;
	int krank = 3;
	hsize_t kdims[3] = {num_slices, sim->Nm, 4};

	kspace = H5Screate_simple( krank, kdims, NULL);
    kdata = H5Dcreate( outfile, "kdtree", 
                       H5T_STD_U16LE, kspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	if( H5Sclose( kspace) < 0 ){
		fprintf(stderr, "Error closing HDF5 data space./n");
		return NULL;
	}
	if( H5Dclose( kdata) < 0 ){
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
	ofile->current_size++;
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
