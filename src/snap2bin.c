#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#define DM_COORDINATES "/PartType1/Coordinates"
#define DM_VELOCITIES "/PartType1/Velocities"
int main(int argc, char **argv){
  hid_t           file, dset, dcpl, filespace;    /* Handles */
  hid_t dsetv_id, status_n, memspace;
  H5D_layout_t    layout;
  long long npoints;
  float *data_pos;
  float *data_vel;
  int rank;
  hsize_t dims[2];
  int i;
  int n_items;
  FILE *binout;

  file = H5Fopen (argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  fprintf(stdout, "opening file %s\n", argv[1]);

  dset = H5Dopen (file, "Heade");
  
  exit(1);


  /*READ DM coordinates*/
  dset = H5Dopen (file, DM_COORDINATES);
  fprintf(stdout, "getting dataset %s\n", DM_COORDINATES);
  filespace = H5Dget_space (dset);

  rank = H5Sget_simple_extent_ndims (filespace);
  fprintf(stdout, "rank is: %d\n", rank);

  status_n = H5Sget_simple_extent_dims (filespace, dims, NULL);
  fprintf(stdout, "dimesions are: %d %d\n", (int)(dims[0]), (int)(dims[1]));

  n_items = (int)dims[0] * (int)dims[1];
  if(!(data_pos=malloc(sizeof(float) * n_items))){
    fprintf(stderr, "problem with data allocation\n");
  }
  memspace = H5Screate_simple (rank,dims,NULL);  
  status_n = H5Dread (dset, H5T_NATIVE_FLOAT, memspace, filespace,
		    H5P_DEFAULT, data_pos);

  /*read DM velocities*/
  dset = H5Dopen (file, DM_VELOCITIES);
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


