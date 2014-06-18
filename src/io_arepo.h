#include <hdf5.h>
#ifndef __IO_AREPO__
#define __IO_AREPO__

#define DM_COORDINATES "/PartType1/Coordinates"
#define DM_VELOCITIES "/PartType1/Velocities"
#define GAS_GRADIENT "/PartType0/VelocityGradient"
#define GAS_VELOCITIES "/PartType0/Velocities"
#define GAS_DENSITY "/PartType0/Density"
#define GAS_VOLUME "/PartType0/Volume"
#define EIGENVECTORS "/ShearTensor/Eigenvectors"
#define EIGENVALUES "/ShearTensor/Eigenvalues"
#define VORTICITY "/ShearTensor/Vorticity"
#define HELICITY "/ShearTensor/Helicity"
#define VOID 0
#define SHEET 1
#define FILAMENT 2
#define PEAK 3
#define EIGEN1 0
#define EIGEN2 1
#define EIGEN3 2

float * load_data_block(hid_t file, char * block_name, int *n_points, int *n_cols);
float *malloc_data(int n_points, int n_cols);
void dump_data_block(hid_t file, char * block_name, float *data, int n_points, int n_cols);

#endif
