#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#define DM_COORDINATES "/PartType1/Coordinates"
#define DM_VELOCITIES "/PartType1/Velocities"
#define GAS_GRADIENT "/PartType0/VelocityGradient"
#define USAGE "./snap2web.x snapfilename"

/*
  Program to read AREPO snapshots and compute for each cell:
  - Eigenvalues and Eigenvectors of the velocitiy shear tensor.
  - The vorticity vector.
  - The helicity scalar.

  It takes as an input the fileame to the AREPO snapshot. It is assumed
  that this snapshot includes a file named VelocityGradient. This is usually
  constructed by postprocessing the AREPO snapshot first.
  
  Author: Jaime Forero-Romero (Uniandes)
  Written: 12-June-2014
  Modified: 

*/

float * load_data_block(hid_t file, char * block_name, int *n_points, int *n_cols);
int main(int argc, char **argv){
  hid_t file;
  float *gas_gradient;
  int n_points, n_cols;
  float shear

  if(argc!=2){
    fprintf(stderr, "USAGE:%s\n", USAGE);
    exit(1);
  }

  file = H5Fopen (argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  fprintf(stdout, "opening file %s\n", argv[1]);

  gas_gradient = load_data_block(file, GAS_GRADIENT, &n_points, &n_cols);  
 
  return 0;
}

float * load_data_block(hid_t file, char * block_name, int *n_points, int *n_cols){
  hid_t dset, filespace, status_n, memspace;
  int rank;
  hsize_t dims[2];
  long long n_items;
  float *data;

  /*READ DM coordinates*/
  dset = H5Dopen2(file, block_name, H5P_DEFAULT);
  fprintf(stdout, "getting dataset %s\n", block_name);
  filespace = H5Dget_space (dset);

  rank = H5Sget_simple_extent_ndims (filespace);
  fprintf(stdout, "rank is: %d\n", rank);

  status_n = H5Sget_simple_extent_dims (filespace, dims, NULL);
  fprintf(stdout, "dimesions are: %d %d\n", (int)(dims[0]), (int)(dims[1]));

  *n_points = dims[0];
  *n_cols = dims[1];

  n_items = (long long)dims[0] * (long long)dims[1];
  if(!(data=malloc(sizeof(float) * n_items))){
    fprintf(stderr, "problem with data allocation\n");
  }
  memspace = H5Screate_simple (rank,dims,NULL);  
  status_n = H5Dread (dset, H5T_NATIVE_FLOAT, memspace, filespace,
		    H5P_DEFAULT, data);
  return data;
}

/*
  group = H5Gopen(file, "/Header", H5P_DEFAULT);
  attr = H5Aopen(group, "BoxSize", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_FLOAT, &box);

  fprintf(stdout, "BoxSize %f\n", box);




  dset = H5Dopen2(file, DM_VELOCITIES,H5P_DEFAULT);
  fprintf(stdout, "getting dataset %s\n", DM_VELOCITIES);
  filespace = H5Dget_space (dset);

  rank = H5Sget_simple_extent_ndims (filespace);
  fprintf(stdout, "rank is: %d\n", rank);

  status_n = H5Sget_simple_extent_dims (filespace, dims, NULL);
  fprintf(stdout, "dimesions are: %d %d\n", (int)(dims[0]), (int)(dims[1]));

  n_items = (int)dims[0] * (int)dims[1];
  if(!(data_vel=malloc(sizeof(float) * n_items))){
    fprintf(stderr, "problem with data allocation\n");
  }
  memspace = H5Screate_simple (rank,dims,NULL);  
  status_n = H5Dread (dset, H5T_NATIVE_FLOAT, memspace, filespace,
		    H5P_DEFAULT, data_vel);



  for(i=0;i<10;i++){
    fprintf(stdout, "%f %f %f\n", 
	    data_pos[i*3 + 0 ], data_pos[i*3 + 1], data_pos[i*3 + 2]);
  }

  if(!(binout=fopen(argv[2], "w"))){
    fprintf(stdout, "problem opening file %s\n", argv[3]);
    exit(1);
  }

  fwrite(&n_items, sizeof(int), 1, binout);

  fwrite(&(data_pos[0]), sizeof(float), n_items, binout);
  fwrite(&(data_vel[0]), sizeof(float), n_items, binout);
  
  fclose(binout);
  return 0;
}

*/

