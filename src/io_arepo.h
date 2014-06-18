#ifndef __IO_AREPO__
#define __IO_AREPO__

#define DM_COORDINATES "/PartType1/Coordinates"
#define DM_VELOCITIES "/PartType1/Velocities"
#define GAS_GRADIENT "/PartType0/VelocityGradient"
#define GAS_VELOCITIES "/PartType0/Velocities"
#define GAS_DENSITY "/PartType0/Density"
#define EIGENVECTORS "/ShearTensor/Eigenvectors"
#define EIGENVALUES "/ShearTensor/Eigenvalues"
#define VORTICITY "/ShearTensor/Vorticity"
#define HELICITY "/ShearTensor/Helicity"
float * load_data_block(hid_t file, char * block_name, int *n_points, int *n_cols);float *malloc_data(int n_points, int n_cols);
void dump_data_block(hid_t file, char * block_name, float *data, int n_points, int n_cols);

#endif
