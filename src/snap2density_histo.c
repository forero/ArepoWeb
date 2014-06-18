#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "io_arepo.h"
#include "array.h"
#define USAGE "./snap2histo.x snapfilename fileout min_dens max_dens n_points"

/*
  Program to read a full arepo snapshot and create a hisogram of the 
  logarithm of the gas density.

  The output is written in ascii format
*/

int main(int argc, char **argv){
  float *gas_density;
  hid_t file, group, attr, status_n;
  int n_points, n_cols, n_files;
  char filein[512];
  char fileout[512];
  int i_file, j;
  double *dens_array_edges;
  double *dens_histo;
  double delta_dens;
  double dummy;
  double min_dens, max_dens;
  int n_points_histo;

  if(argc!=6){
    fprintf(stderr, "USAGE:%s\n", USAGE);
    exit(1);    
  }

  min_dens =atof(argv[3]);
  max_dens =atof(argv[4]);
  n_points_histo =atoi(argv[5]);

  strcpy(fileout, argv[2]);
  fprintf(stderr, "outputs go to %s\n", fileout);

  /*generate the density array*/
  dens_array_edges = gen_linspace(min_dens, max_dens, n_points_histo+1, &delta_dens);

  /*initialize the histogram*/
  dens_histo = gen_linspace(0.0, 0.0, n_points_histo, &dummy);  

  /*Read the first file of the snapshot*/
  sprintf(filein, "%s.0.hdf5", argv[1]);
  file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
  group = H5Gopen(file, "/Header", H5P_DEFAULT);
  attr = H5Aopen(group, "NumFilesPerSnapshot", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_INT, &n_files);
  H5Fclose (file);
  fprintf(stdout, "%d\n", n_files);

  
  /*loop over the files*/
  for(i_file=0;i_file<n_files;i_file++){
    sprintf(filein, "%s.%d.hdf5", argv[1], i_file);
    file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
    gas_density = load_data_block(file, GAS_DENSITY, &n_points, &n_cols);  
    H5Fclose (file);
    
    /*we want to have the density in log scale*/
    for(j=0;j<n_points;j++){
      gas_density[j] = log10(gas_density[j]);
    }
    
    /*update the hisrtogram according to the data in the array*/
    add_to_histo(dens_histo, dens_array_edges, n_points_histo, gas_density, n_points);
    
    free(gas_density);
  }

  /*write the density to disk*/
  dump_histo(fileout, dens_histo, dens_array_edges, n_points_histo);


  return 0;
}




