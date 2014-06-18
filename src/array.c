#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io_arepo.h"
#include "struct.h"


void add_to_environments(double *total_env, double *lambda_th, float *cargo, float *eigenvalues, int n_lambda, int n_points, int env_kind, float *volume, float min_len, float max_len){
  int i, j;
  float my_len;

  for(i=0;i<n_points;i++){
    my_len = pow((3.0*volume[i]/(4.0*PI)), 1.0/3.0);
    if((my_len>min_len) && (my_len<max_len)){
      for(j=0;j<n_lambda;j++){      
	if(env_kind==VOID){
	  if(eigenvalues[i*3 + EIGEN1] < lambda_th[j]){
	    total_env[j] +=cargo[j];
	  }
	}else if(env_kind==SHEET){
	  if((eigenvalues[i*3 + EIGEN2] < lambda_th[j]) && (eigenvalues[i*3 + EIGEN1]>lambda_th[j])){
	    total_env[j] +=cargo[j];
	  }
	}else if(env_kind==FILAMENT){
	  if((eigenvalues[i*3 + EIGEN3] < lambda_th[j]) && (eigenvalues[i*3 + EIGEN2]>lambda_th[j])){
	    total_env[j] +=cargo[j];
	  }
	}else if(env_kind==PEAK){
	  if(eigenvalues[i*3 + EIGEN3]>lambda_th[j]){
	    total_env[j] +=cargo[j];
	  }
	}else{
	  fprintf(stderr, "ERROR: uknown kind of environment\n");
	  exit(1);
	}
      }
    }
  }
}

void dump_histo(char *filename, double*histo, double*edges, int n_histo){
  FILE *out;
  int i;
  double dummy;
  if(!(out=fopen(filename, "w"))){
    fprintf(stderr, "Problem opening file %s\n", filename);
    exit(1);
  }
  
  for(i=0;i<n_histo;i++){
    //    dummy = 0.5*(edges[i] + edges[i+1]);
    dummy = edges[i];
    fprintf(out, "%f %f\n", dummy, histo[i]);
  }  
  fclose(out);
}

double*gen_linspace(double min_val, double max_val, int n_points, double*delta_val){
  double*data;
  double delta;
  int i;
  if(!(data=malloc(sizeof(double )*n_points))){
    fprintf(stderr, "allocation problem\n");
    exit(1);
  }
  
  delta = (max_val - min_val)/(n_points - 1);
  
  for(i=0;i<n_points;i++){
    data[i] = min_val + i*delta;
  }
  *delta_val = delta;
  return data;
}

void  add_to_histo(double *histo, double*edges, int n_histo, float *data, int n_data){
  int i,j; 
  int in;
  for(i=0;i<n_data;i++){
    in = 0;
    j = 0;
    do{
      if((data[i]<edges[j+1])&&(data[i]>edges[j])){
	histo[j] += 1.0;
	in=1;
      }
      j++;
    }while((j<n_histo) && (in==0));
  }

}
