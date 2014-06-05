#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv){
  hid_t           file, space, dset, dcpl;    /* Handles */


  file = H5Fopen (argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  fprintf(stdout, "opening file %s\n", argv[1]);

  return 0;
}


