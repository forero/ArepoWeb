#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "io_arepo.h"
#include "array.h"
#include "struct.h"
#define USAGE "./snap2lambda.x snapfilename webfile fileout min_length max_length sigma_lambda"



/*
  This program calculates the histogram for the eigenvalues and its trace
  in cells within a range of scale lengths.

  The scale lenght is estimated as the radius of the sphere with the same volume.

  The output is written in ascii format.
*/

int main(int argc, char **argv){
  hid_t file, group, attr, status_n;
  int n_points, n_cols, n_files;
  int i_file;
  float *eigenvalues;
  float *gas_volume;
  double min_len;
  double max_len;
  char filein[512];
  char fileout[512];
  char longfileout[512];
  double *lambda1_th_edges;
  double *lambda2_th_edges;
  double *lambda3_th_edges;
  double *histo_lambda1;
  double *histo_lambda2;
  double *histo_lambda3;
  double *histo_trace;
  float sigma_lambda;
  double median, variance;
  int n_lambda = 100;
  double dummy;

  if(argc!=7){
    fprintf(stderr, "USAGE: %s\n", USAGE);
    exit(1);
  }
  
  min_len = atof(argv[4]);
  max_len = atof(argv[5]);
  sigma_lambda = atof(argv[6]);

  strcpy(fileout, argv[3]);
  fprintf(stderr, "outputs go to %s\n", fileout);

  /*read the first eigenvalue files*/
  sprintf(filein, "%s.%d.hdf5", argv[2], 10);
  file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
  eigenvalues = load_data_block(file, EIGENVALUES, &n_points, &n_cols);  
  H5Fclose (file);  

  /*read the volume*/
  sprintf(filein, "%s.%d.hdf5", argv[1], 10);
  file = H5Fopen (filein, H5F_ACC_RDONLY, H5P_DEFAULT);
  gas_volume = load_data_block(file, GAS_VOLUME, &n_points, &n_cols);  
  H5Fclose (file);

  /*generate the interval for the lambda values from the variance and the mean*/
  //  get_median_variance_cut(eigenvalues, n_points, n_cols, EIGEN1, &median, &variance, gas_volume, min_len, max_len);
  lambda1_th_edges = gen_linspace(-3.0*sigma_lambda, +3.0*sigma_lambda, n_lambda+1, &dummy);

  //  get_median_variance_cut(eigenvalues, n_points, n_cols, EIGEN2, &median, &variance, gas_volume, min_len, max_len);
  lambda2_th_edges = gen_linspace(-3.0*sigma_lambda, +3.0*sigma_lambda, n_lambda+1, &dummy);

  //  get_median_variance_cut(eigenvalues, n_points, n_cols, EIGEN3, &median, &variance, gas_volume, min_len, max_len);
  lambda3_th_edges = gen_linspace(-3.0*sigma_lambda, +3.0*sigma_lambda, n_lambda+1, &dummy);

  free(eigenvalues);
  free(gas_volume);

  /*initialize all the arrays to get the volumes and masses*/
  histo_lambda1 = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  histo_lambda2 = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  histo_lambda3 = gen_linspace(0.0, 0.0, n_lambda, &dummy);
  histo_trace = gen_linspace(0.0, 0.0, n_lambda, &dummy);


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

    add_to_histo_cut(histo_lambda1, lambda1_th_edges, n_lambda, eigenvalues, n_points, 3, EIGEN1,
		     gas_volume, min_len, max_len);    
    add_to_histo_cut(histo_lambda2, lambda2_th_edges, n_lambda, eigenvalues, n_points, 3, EIGEN2,
		     gas_volume, min_len, max_len);    
    add_to_histo_cut(histo_lambda3, lambda3_th_edges, n_lambda, eigenvalues, n_points, 3, EIGEN3,
		     gas_volume, min_len, max_len);    

    free(gas_volume);
    free(eigenvalues);
  }

  /*dump to file*/

  sprintf(longfileout,"%s_len_%.1f_%.1f_lambda1_histo.dat", fileout, min_len, max_len);
  dump_histo(longfileout, histo_lambda1, lambda1_th_edges, n_lambda);
  sprintf(longfileout,"%s_len_%.1f_%.1f_lambda2_histo.dat", fileout, min_len, max_len);
  dump_histo(longfileout, histo_lambda2, lambda2_th_edges, n_lambda);
  sprintf(longfileout,"%s_len_%.1f_%.1f_lambda3_histo.dat", fileout, min_len, max_len);
  dump_histo(longfileout, histo_lambda3, lambda3_th_edges, n_lambda);
   
  return 0;
}


