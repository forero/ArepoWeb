#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define DM_COORDINATES "/PartType1/Coordinates"
#define DM_VELOCITIES "/PartType1/Velocities"
#define GAS_GRADIENT "/PartType0/VelocityGradient"
#define GAS_VELOCITIES "/PartType0/Velocities"
#define EIGENVECTORS "/ShearTensor/Eigenvectors"
#define EIGENVALUES "/ShearTensor/Eigenvalues"
#define VORTICITY "/ShearTensor/Vorticity"
#define HELICITY "/ShearTensor/Helicity"
#define USAGE "./snap2web.x snapfilename fileout"
#define XX  0
#define XY  1
#define XZ  2
#define YX  3
#define YY  4
#define YZ  5
#define ZX  6
#define ZY  7
#define ZZ  8
#define PI  3.141592653589793238462643383279502884197169

#define HUBBLE  0.10 /*Hubble constant in AREPO units*/
/*
  Program to read AREPO snapshots and compute for each cell:
  - Eigenvalues and Eigenvectors of the velocitiy shear tensor.
  - The vorticity vector.
  - The helicity scalar.

  It takes as an input the fileame to the AREPO snapshot. It is assumed
  that this snapshot includes a file named VelocityGradient. This is usually
  constructed by postprocessing the AREPO snapshot first.

  The output is written in HDF5 format.
  
  Author: Jaime Forero-Romero (Uniandes)
  Written: 12-June-2014
  Modified: 

*/

int Diagonalise3x3(double *M,double *L,double *vec);
float * load_data_block(hid_t file, char * block_name, int *n_points, int *n_cols);float *malloc_data(int n_points, int n_cols);
void dump_data_block(hid_t file, char * block_name, float *data, int n_points, int n_cols);
int main(int argc, char **argv){
  char fileout[512];
  hid_t file, grp;
  float *gas_gradient;
  float *gas_velocity;
  float *gas_eigenvalues;
  float *gas_eigenvectors;
  float *gas_vorticity;
  float *gas_helicity;
  int n_points, n_cols;
  double shear_tensor[9];
  double eigenvalues[3];
  double eigenvectors[9];
  float vorticity[3];
  float helicity;
  int i,j;

  if(argc!=3){
    fprintf(stderr, "USAGE:%s\n", USAGE);
    exit(1);    
  }

  /*output filename*/
  strcpy(fileout, argv[2]);

  /*load the input data*/
  file = H5Fopen (argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
  fprintf(stdout, "opening file %s\n", argv[1]);
  gas_gradient = load_data_block(file, GAS_GRADIENT, &n_points, &n_cols);  
  gas_velocity = load_data_block(file, GAS_VELOCITIES, &n_points, &n_cols);  
  H5Fclose (file);

  /*save memory for the outputs*/
  gas_eigenvalues = malloc_data(n_points, 3);
  gas_eigenvectors = malloc_data(n_points, 9);
  gas_vorticity = malloc_data(n_points, 3);
  gas_helicity = malloc_data(n_points, 1);

  for(i=0;i<n_points;i++){   
    shear_tensor[XX] = -0.5*(gas_gradient[i*9 + XX] + gas_gradient[i*9 + XX]); 
    shear_tensor[XY] = -0.5*(gas_gradient[i*9 + YY] + gas_gradient[i*9 + YY]); 
    shear_tensor[XZ] = -0.5*(gas_gradient[i*9 + ZZ] + gas_gradient[i*9 + ZZ]); 
    shear_tensor[XY] = -0.5*(gas_gradient[i*9 + XY] + gas_gradient[i*9 + YX]); 
    shear_tensor[XZ] = -0.5*(gas_gradient[i*9 + XZ] + gas_gradient[i*9 + YZ]); 
    shear_tensor[YZ] = -0.5*(gas_gradient[i*9 + YZ] + gas_gradient[i*9 + ZY]); 
    shear_tensor[YX] = shear_tensor[XY];
    shear_tensor[ZX] = shear_tensor[XZ];
    shear_tensor[ZY] = shear_tensor[YZ];

    /*Get the Eigenvalues and eigenvectors from the shear tensor*/
    Diagonalise3x3(shear_tensor, eigenvalues,eigenvectors);

    /*Computes the vorticity*/
    vorticity[0] = gas_gradient[i*9 +ZY] - gas_gradient[i*9 + YZ];
    vorticity[1] = gas_gradient[i*9 +XZ] - gas_gradient[i*9 + ZX];
    vorticity[2] = gas_gradient[i*9 +YX] - gas_gradient[i*9 + XY];

    /*compute the helicity*/
    helicity = 0.0;
    for(j=0;j<3;j++){
      helicity+= vorticity[j] * gas_velocity[i*3 + j];
    }

    /*fill the arrays*/
    gas_helicity[i] = helicity;

    for(j=0;j<3;j++){
      gas_eigenvalues[3*i + j] = eigenvalues[j];
      gas_vorticity[3*i + j] = vorticity[j];
    }
    for(j=0;j<9;j++){
      gas_eigenvectors[9*i + j] = eigenvectors[j];
    }
  } 

  /*Dump the data to disk*/
  file = H5Fcreate(fileout, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  grp = H5Gcreate(file, "/ShearTensor", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dump_data_block(file, EIGENVECTORS, gas_eigenvectors, n_points, 9);
  dump_data_block(file, EIGENVALUES, gas_eigenvalues, n_points, 3);
  dump_data_block(file, VORTICITY, gas_vorticity, n_points, 3);
  dump_data_block(file, HELICITY, gas_helicity, n_points, 1);

  H5Fclose (file);  
  return 0;
}

float *malloc_data(int n_points, int n_cols){
  float *data;
  
  if(!(data=malloc(sizeof(float) * n_points * n_cols))){
    fprintf(stderr, "problem with data allocation\n");
    exit(1);
  }
  return data;
}

void dump_data_block(hid_t file, char * block_name, float *data, int n_points, int n_cols){
  hid_t dataspace, status, dataset, plist;
  hsize_t dims[2];


  dims[0] = n_points;
  dims[1] = n_cols;
  dataspace = H5Screate_simple(2, dims, NULL); 

  plist     = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist, 2, dims);
  dataset = H5Dcreate(file, block_name, 
		      H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, plist, H5P_DEFAULT);
  status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, data);


  status = H5Dclose(dataset);
  status = H5Pclose(plist);
  status = H5Sclose(dataspace);
}

float * load_data_block(hid_t file, char * block_name, int *n_points, int *n_cols){
  hid_t dset, filespace, status_n, memspace;
  int rank;
  hsize_t dims[2];
  long long n_items;
  float *data;

  /*READ DM coordinates*/
  dset = H5Dopen2(file, block_name, H5P_DEFAULT);
  fprintf(stdout, "getting dataset %s\n", block_name);
  filespace = H5Dget_space (dset);

  rank = H5Sget_simple_extent_ndims (filespace);
  fprintf(stdout, "rank is: %d\n", rank);

  status_n = H5Sget_simple_extent_dims (filespace, dims, NULL);
  fprintf(stdout, "dimesions are: %d %d\n", (int)(dims[0]), (int)(dims[1]));

  *n_points = dims[0];
  *n_cols = dims[1];

  n_items = (long long)dims[0] * (long long)dims[1];
  if(!(data=malloc(sizeof(float) * n_items))){
    fprintf(stderr, "problem with data allocation\n");
  }
  memspace = H5Screate_simple (rank,dims,NULL);  
  status_n = H5Dread (dset, H5T_NATIVE_FLOAT, memspace, filespace,
		    H5P_DEFAULT, data);
  return data;
}

/*
  group = H5Gopen(file, "/Header", H5P_DEFAULT);
  attr = H5Aopen(group, "BoxSize", H5P_DEFAULT);
  status_n = H5Aread(attr, H5T_NATIVE_FLOAT, &box);

  fprintf(stdout, "BoxSize %f\n", box);




  dset = H5Dopen2(file, DM_VELOCITIES,H5P_DEFAULT);
  fprintf(stdout, "getting dataset %s\n", DM_VELOCITIES);
  filespace = H5Dget_space (dset);

  rank = H5Sget_simple_extent_ndims (filespace);
  fprintf(stdout, "rank is: %d\n", rank);

  status_n = H5Sget_simple_extent_dims (filespace, dims, NULL);
  fprintf(stdout, "dimesions are: %d %d\n", (int)(dims[0]), (int)(dims[1]));

  n_items = (int)dims[0] * (int)dims[1];
  if(!(data_vel=malloc(sizeof(float) * n_items))){
    fprintf(stderr, "problem with data allocation\n");
  }
  memspace = H5Screate_simple (rank,dims,NULL);  
  status_n = H5Dread (dset, H5T_NATIVE_FLOAT, memspace, filespace,
		    H5P_DEFAULT, data_vel);



  for(i=0;i<10;i++){
    fprintf(stdout, "%f %f %f\n", 
	    data_pos[i*3 + 0 ], data_pos[i*3 + 1], data_pos[i*3 + 2]);
  }

  if(!(binout=fopen(argv[2], "w"))){
    fprintf(stdout, "problem opening file %s\n", argv[3]);
    exit(1);
  }

  fwrite(&n_items, sizeof(int), 1, binout);

  fwrite(&(data_pos[0]), sizeof(float), n_items, binout);
  fwrite(&(data_vel[0]), sizeof(float), n_items, binout);
  
  fclose(binout);
  return 0;
}

*/



//Diagonalise 3x3 symetric matrix M, L are the eigenvalues (L[0]>=L[1]>=L[2]) and vec the eigenvectors
//returns 0 if OK, true if degenerate
//                   (M0 M1 M2)
//M = (M0 .. M5) --> (M1 M3 M4)
//                   (M2 M4 M5)
int Diagonalise3x3(double *M,double *L,double *vec)
{
    double a,b,c,u,v,w;
    double B,C,D;
    double p,q,r;
    double phi;
    double tmp;
    int i;
    
    a=M[0];u=M[1];v=M[2];
    b=M[3];c=M[5];w=M[4];

    //if ((a==b)&&(b==c)&&(c==0)) return -1;

    //Compute characteristic polynome
    B=-(a+b+c);
    C=a*b+a*c+b*c-u*u-v*v-w*w;
    D=u*u*c+v*v*b+w*w*a-a*b*c-2*u*v*w;

    //solve 3rd degree equation
    p=(3.*C-B*B)/9.;
    q=(2.*B*B*B/27.-B*C/3.+D)/2.;
    
    if (q*q+p*p*p>=0) 
    {
      L[0]=L[1]=L[2]=0;
      return -1;
    }
    
   
    r=((q<0)?-1:1)*sqrt(fabs(p));
    phi = acos(q/(r*r*r))/3.;
	
    r*=2;//B/=3;
    L[0]=-r*cos(phi)-B/3.;
    L[1]=r*cos(PI/3-phi)-B/3.;
    L[2]=r*cos(PI/3+phi)-B/3.;

    
    
    if ((L[2])>(L[1])) {tmp=L[1];L[1]=L[2];L[2]=tmp;}
    if ((L[1])>(L[0])) {tmp=L[1];L[1]=L[0];L[0]=tmp;}
    if ((L[2])>(L[1])) {tmp=L[1];L[1]=L[2];L[2]=tmp;}
    
    //Some special case ...

    //matrix is already diagonal
    if ((u==0)&&(v==0)&&(w==0))
    {
	vec[0]=0;vec[1]=0;vec[2]=0;
	//associate each axis (x,y and z) with it s eigenvalues
	for (i=0;i<3;i++)
	{
	    if (fabs(L[i])>1.E-8)
	    {
		vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;
		if (fabs((a-L[i])/L[i])<1.E-4) vec[3*i+0]=1;
		else if (fabs((b-L[i])/L[i])<1.E-4) vec[3*i+1]=1;
		else if (fabs((c-L[i])/L[i])<1.E-4) vec[3*i+2]=1;
	    }
	    else {L[i]=0;vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;}
	}
    }
    
    //We already have 1 eigenvalue (matrix is block diagonal)
    else if (((u==0)&&(v==0))||((v==0)&&(w==0))||((u==0)&&(w==0)))
    {
      double tmpa,tmpb,tmpu;
	int tmpd;
	//Finds which axis are already eigenvector (x,y or z)
	vec[0]=0;vec[1]=0;vec[2]=0;
	for (i=0;i<3;i++)
	{
	    if (fabs(L[i])>1.E-8)
	    {
		vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;
		//Check if we already have this eigenvector
		if (fabs((a-L[i])/L[i])<1.E-4) vec[3*i+0]=1;
		else if (fabs((b-L[i])/L[i])<1.E-4) vec[3*i+1]=1;
		else if (fabs((c-L[i])/L[i])<1.E-4) vec[3*i+2]=1;
		//no ...
		//diagonalise 2x2 block ...
		else
		{
		    if ((u==0)&&(v==0)) {tmpa=b;tmpb=c;tmpu=w;tmpd=1;}
		    if ((v==0)&&(w==0)) {tmpa=a;tmpb=b;tmpu=u;tmpd=0;}
		    if ((u==0)&&(w==0)) {tmpa=a;tmpb=c;tmpu=v;tmpd=2;}
		    if ((tmpa==tmpb)&&(tmpb==0)) 
			B=C=0;
		    else
		    {
			B=(L[i]+tmpu-tmpa)/(L[i]+tmpu-tmpb);
			C=1./sqrt(1+B*B);
		    }
		    if (tmpd==2)
		    {
			vec[3*i]=C;
			vec[3*i+2]=B*C;
		    }
		    else
		    {
			vec[3*i+tmpd]=C;
			vec[3*i+tmpd+1]=B*C;
		    }
		}
	    }
	    else {L[i]=0;vec[3*i+0]=0;vec[3*i+1]=0;vec[3*i+2]=0;}
	}
    }

    //General case ...
    else for (i=0;i<3;i++)
    {

	B=L[i]-b-u*u/(L[i]-a);
	B=(w+u*v/(L[i]-a))/B;
	
	C=u/(L[i]-a)*B+v/(L[i]-a);
	D=1./sqrt(B*B+C*C+1.);

	vec[3*i+2]=D;
	vec[3*i+1]=B*D;
	vec[3*i]=C*D;
    }
    
    return 0;
}

