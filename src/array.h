#ifndef __ARRAY_H__
#define __ARRAY_H__

void get_min_max(float *array, int n_points, double *min_val, double *max_val);
void add_to_environments(double *total_env, double *lambda_th, float *cargo, float *eigenvalues, int n_lambda, int n_points, int env_kind, float *volume, float min_len, float max_len);
void dump_histo(char *filename, double *histo, double *edges, int n_histo);
double *gen_linspace(double min_val, double max_val, int n_points, double *delta_val);
void  add_to_histo(double *histo, double*edges, int n_histo, float *data, int n_data, int n_cols, int this_column);
void  add_to_histo_cut(double *histo, double*edges, int n_histo, float *data, int n_data, int n_cols, int this_column,
		       float *volume, float min_len, float max_len);
void get_min_max_cut(float *array, int n_points, int ncols, int this_column, double *min_val, double *max_val, float *volume, float min_len, float max_len);
void get_median_variance_cut(float *array, int n_points, int n_cols, int this_column , double *median, double *variance, float *volume, float min_len, float max_len);
#endif
