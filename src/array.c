#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io_arepo.h"
#include "struct.h"

void get_min_max(float *array, int n_points, double *min_val, double *max_val){
  double min=1E10, max=-1E10;
  int i;
  for(i=0;i<n_points;i++){
    if(array[i]<min){
      min = array[i];
    }      
    if(array[i]>max){
      max = array[i];
    }
  }
  *min_val = min;
  *max_val = max;
}

void get_min_max_cut(float *array, int n_points, int n_cols, int this_column , double *min_val, double *max_val, float *volume, float min_len, float max_len){
  double min=1E10, max=-1E10;
  int i;
  float my_len;
  
  for(i=0;i<n_points;i++){
    my_len = pow((3.0*volume[i]/(4.0*PI)), 1.0/3.0);      
      if((my_len>min_len) && (my_len<max_len)){      
	if(array[i*n_cols + this_column]<min){
	  min = array[i*n_cols + this_column];
	}      
	if(array[i*n_cols + this_column]>max){
	  max = array[i*n_cols + this_column];
	}
      }
    
  }
  *min_val = min;
  *max_val = max;
}

void get_median_variance_cut(float *array, int n_points, int n_cols, int this_column , double *median_out, double *variance_out, float *volume, float min_len, float max_len){
  double median=0.0;
  double variance=0.0;
  int i;
  float my_len;
  float point;
  int first_point=0;  
  float median_old;

  i=0;
  do{
    my_len = pow((3.0*volume[i]/(4.0*PI)), 1.0/3.0);      
    if((my_len>min_len) && (my_len<max_len)){      
      median = array[i*n_cols + this_column];
      first_point=1;
    }    
    i++;
  }while(i<n_points && first_point==0);

  for(i=0;i<n_points;i++){    
    my_len = pow((3.0*volume[i]/(4.0*PI)), 1.0/3.0);      
    if((my_len>min_len) && (my_len<max_len)){      
      point = array[i*n_cols + this_column];
      first_point++;
      median_old = median;
      median = median + (point - median)/first_point;
      variance = variance + (point - median)*(point - median_old);
    }    
  }
  variance = sqrt(variance/(first_point-1));
  fprintf(stdout, "median variance %f %f, total points %d\n", median, variance, first_point);
  *median_out = median;
  *variance_out = variance;
}

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
  
  fprintf(stdout,"min max in linear : %f %f\n", min_val, max_val);

  delta = (max_val - min_val)/(n_points - 1);
  
  for(i=0;i<n_points;i++){
    data[i] = min_val + i*delta;
  }
  *delta_val = delta;
  return data;
}

void  add_to_histo(double *histo, double*edges, int n_histo, float *data, int n_data, int n_cols, int this_column){
  int i,j; 
  int in;
  for(i=0;i<n_data;i++){
    in = 0;
    j = 0;
    do{
      if((data[i*n_cols + this_column]<edges[j+1])&&(data[i*n_cols + this_column]>edges[j])){
	histo[j] += 1.0;
	in=1;
      }
      j++;
    }while((j<n_histo) && (in==0));
  }

}


void  add_to_histo_cut(double *histo, double*edges, int n_histo, float *data, int n_data, int n_cols, int this_column,
		       float *volume, float min_len, float max_len){
  int i,j; 
  int in;
  float my_len;

  for(i=0;i<n_data;i++){
    my_len = pow((3.0*volume[i]/(4.0*PI)), 1.0/3.0);
    if((my_len>min_len) && (my_len<max_len)){
      in = 0;
      j = 0;
      do{
	if((data[i*n_cols + this_column]<edges[j+1])&&(data[i*n_cols + this_column]>edges[j])){
	  histo[j] += 1.0;	  
	  in=1;
	}
	j++;
      }while((j<n_histo) && (in==0));
    }
  }
}
