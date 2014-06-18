#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "io_arepo.h"
#include "array.h"
#define USAGE "./snap2histo.x snapfilename webfile fileout min_length max_length max_lambda_th"



/*
  This program calculates the mass and volume filling fractions for all the
  kinds of environments within a range of scale length for the cells. 
  The scale lenght is estimated as the radius of the sphere with the same volume.

  The output is written in ascii format.
*/

int main(int argc, char **argv){
  hid_t file, group, attr, status_n;
  int n_points, n_cols, n_files;
  int i_file, j;
  float *eigenvalues;
  float *gas_density;
  float *gas_volume;
  double min_len;
  double max_len;
  double max_lambda_th;
  char filein[512];
  char fileout[512];
  char longfileout[512];
  double *lambda_th_edges;
  double *volume_void;
  double *volume_sheet;
  double *volume_filament;
  double *volume_peak;
  double *mass_void;
  double *mass_sheet;
  double *mass_filament;
  double *mass_peak;
  double delta_dens;
  int n_lambda = 10;
  double dummy;

  if(argc!=7){
    fprintf(stderr, "USAGE: %s\n", USAGE);
    exit(1);
  }
  
  min_len = atof(argv[4]);
  max_len = atof(argv[5]);
  max_lambda_th = atof(argv[6]);

  /*generate the array for the lambda threshold values*/
  lambda_th_edges = gen_linspace(0.0, max_lambda_th, n_lambda+1, &delta_dens);

  /*initialize all the arrays to get the volumes and masses*/
  mass_void = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  mass_sheet = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  mass_filament = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  mass_peak = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  volume_void = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  volume_sheet = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  volume_filament = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  volume_peak = gen_linspace(0.0, 0.0, n_lambda, &dummy);

  strcpy(fileout, argv[3]);
  fprintf(stderr, "outputs go to %s\n", fileout);

  /*Read the first file of the snapshot*/
  sprintf(filein, "%s.0.hdf5", argv[1]);
  file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
  group = H5Gopen(file, "/Header", H5P_DEFAULT);
  attr = H5Aopen(group, "NumFilesPerSnapshot", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_INT, &n_files);
  H5Fclose (file);
  fprintf(stdout, "%d\n", n_files);


  /*loop over the snapshots and the eigenvalues*/
  for(i_file=0;i_file<n_files;i_file++){
    /*density*/
    sprintf(filein, "%s.%d.hdf5", argv[1], i_file);
    file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
    gas_density = load_data_block(file, GAS_DENSITY, &n_points, &n_cols);  
    H5Fclose (file);

    /*volume*/
    sprintf(filein, "%s.%d.hdf5", argv[1], i_file);
    file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
    gas_volume = load_data_block(file, GAS_VOLUME, &n_points, &n_cols);  
    H5Fclose (file);

    /*eigenvalues*/
    sprintf(filein, "%s.%d.hdf5", argv[2], i_file);
    file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
    eigenvalues = load_data_block(file, EIGENVALUES, &n_points, &n_cols);  
    H5Fclose (file);

    /*convert the density into the mass*/
    for(j=0;j<n_points;j++){
      gas_density[j] = gas_density[j] * gas_volume[j];
    }


    add_to_environments(volume_void, lambda_th_edges, gas_volume, eigenvalues, n_lambda, n_points, VOID,
			gas_volume, min_len, max_len);    
    add_to_environments(volume_sheet, lambda_th_edges, gas_volume, eigenvalues, n_lambda, n_points, SHEET,
			gas_volume, min_len, max_len);    
    add_to_environments(volume_filament, lambda_th_edges, gas_volume, eigenvalues, n_lambda, n_points, FILAMENT,
			gas_volume, min_len, max_len);    
    add_to_environments(volume_peak, lambda_th_edges, gas_volume, eigenvalues, n_lambda, n_points, PEAK,
			gas_volume, min_len, max_len);
  }

  /*normalize to get the total values*/
  for(j=0;j<n_lambda;j++){
    
    dummy = volume_void[j] + volume_sheet[j] + volume_filament[j] + volume_peak[j];
    volume_void[j] /= dummy;
    volume_sheet[j] /= dummy;
    volume_filament[j] /= dummy;
    volume_peak[j] /= dummy;
    
  }

  /*dump to file*/

  sprintf(longfileout,"%s_max_lambda_%.1e_len_%.1f_%.1f_void_VFF.dat", fileout, max_lambda_th, min_len, max_len);
  dump_histo(longfileout, volume_void, lambda_th_edges, n_lambda);
  sprintf(longfileout,"%s_max_lambda_%.1e_len_%.1f_%.1f_sheet_VFF.dat", fileout, max_lambda_th, min_len, max_len);
  dump_histo(longfileout, volume_sheet, lambda_th_edges, n_lambda);
  sprintf(longfileout,"%s_max_lambda_%.1e_len_%.1f_%.1f_filament_VFF.dat", fileout, max_lambda_th, min_len, max_len);
  dump_histo(longfileout, volume_filament, lambda_th_edges, n_lambda);
  sprintf(longfileout,"%s_max_lambda_%.1e_len_%.1f_%.1f_peak_VFF.dat", fileout, max_lambda_th, min_len, max_len);
  dump_histo(longfileout, volume_peak, lambda_th_edges, n_lambda);

   
  return 0;
}


