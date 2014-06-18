#include <stdio.h>
#include <stdlib.h>

void dump_histo(char *filename, double*histo, double*edges, int n_histo){
  FILE *out;
  int i;
  double dummy;
  if(!(out=fopen(filename, "w"))){
    fprintf(stderr, "Problem opening file %s\n", filename);
    exit(1);
  }
  
  for(i=0;i<n_histo;i++){
    dummy = 0.5*(edges[i] + edges[i+1]);
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
