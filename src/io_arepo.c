#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "io_arepo.h"

void dump_data_block(hid_t file, char * block_name, float *data, int n_points, int n_cols){
  hid_t dataspace, status, dataset, plist;
  hsize_t dims[2];


  dims[0] = n_points;
  dims[1] = n_cols;
  dataspace = H5Screate_simple(2, dims, NULL); 

  plist     = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist, 2, dims);
  dataset = H5Dcreate(file, block_name, 
		      H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);
  status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, data);
  
  status = H5Dclose(dataset);
  status = H5Pclose(plist);
  status = H5Sclose(dataspace);
}

float * load_data_block(hid_t file, char * block_name, int *n_points, int *n_cols){
  hid_t dset, filespace, status_n, memspace;
  int rank;
  hsize_t dims[2];
  long long n_items;
  float *data;

  dset = H5Dopen2(file, block_name, H5P_DEFAULT);
  fprintf(stdout, "getting dataset %s\n", block_name);
  filespace = H5Dget_space (dset);

  rank = H5Sget_simple_extent_ndims (filespace);
  fprintf(stdout, "rank is: %d\n", rank);

  status_n = H5Sget_simple_extent_dims (filespace, dims, NULL);
  *n_points = dims[0];
  *n_cols = dims[1];
  if(rank==1){
    dims[1] = 1;
    *n_cols = 1;
  }
  fprintf(stdout, "dimesions are: %d %d\n", (int)(dims[0]), (int)(dims[1]));


  n_items = (long long)dims[0] * (long long)dims[1];
  if(!(data=malloc(sizeof(float) * n_items))){
    fprintf(stderr, "problem with data allocation\n");
  }
  memspace = H5Screate_simple (rank,dims,NULL);  
  status_n = H5Dread (dset, H5T_NATIVE_FLOAT, memspace, filespace,
		    H5P_DEFAULT, data);
  return data;
}
