#ifndef __ARRAY_H__
#define __ARRAY_H__

void dump_histo(char *filename, double *histo, double *edges, int n_histo);
double *gen_linspace(double min_val, double max_val, int n_points, double *delta_val);
void  add_to_histo(double *histo, double *edges, int n_histo, float *data, int n_data);
#endif
