#ifndef __ARRAY_H__
#define __ARRAY_H__

void add_to_environments(double *total_env, double *lambda_th, float *cargo, float *eigenvalues, int n_lambda, int n_points, int env_kind, float *volume, float min_len, float max_len);
void dump_histo(char *filename, double *histo, double *edges, int n_histo);
double *gen_linspace(double min_val, double max_val, int n_points, double *delta_val);
void  add_to_histo(double *histo, double *edges, int n_histo, float *data, int n_data);
#endif
