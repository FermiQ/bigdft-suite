#include <iostream>
#include <stdio.h>
#include "cufft.h"
#include "cuda.h"
#include "cuda_runtime_api.h"
 
#define DOUBLE

#ifdef DOUBLE
#define Complex  cufftDoubleComplex
#define Real double
#define Transform CUFFT_Z2Z
#define TransformExec cufftExecZ2Z
#else
#define Complex  cufftComplex
#define Real float
#define Transform CUFFT_C2C
#define TransformExec cufftExecC2C
#endif

#define TILE_DIM  8

// synchronize blocks
extern "C" void synchronize_() {
 
  cudaThreadSynchronize();
}

// allocate device memory
extern "C" void cudamalloc_(int *size, Real **d_data) {

  cudaMalloc((void**)d_data, sizeof(Real)*(*size));
  if( cudaGetLastError() != cudaSuccess)
      printf("allocate error\n");
}

extern "C" void cudafree_(Real **d_data) {

  cudaFree(*d_data);
}

// set device memory
extern "C" void reset_gpu_data_(int *size, Real* h_data, Real **d_data){

 cudaMemcpy(*d_data, h_data, sizeof(Real)*(*size),
         cudaMemcpyHostToDevice);
 if( cudaGetLastError() != cudaSuccess)
      printf("transfer error\n");

}

// read device memory
extern "C" void get_gpu_data_(int *size, Real *h_data, Real **d_data) {

 cudaMemcpy(h_data, *d_data, sizeof(Real)*(*size),
         cudaMemcpyDeviceToHost);
 if (cudaGetLastError() != cudaSuccess)
        printf("transfer back error\n");
}


// transpose
__global__ void transpose(Complex *idata, Complex *odata,
        int width, int height)
{
  __shared__ Complex tile[TILE_DIM][TILE_DIM+1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*(width);
  int xIndex1 = blockIdx.y * TILE_DIM + threadIdx.x;
  int yIndex1 = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex1 + (yIndex1)*height;

  if (xIndex < width && yIndex < height)
      tile[threadIdx.y][threadIdx.x] = idata[index_in];
    __syncthreads();

  if (xIndex1 < height && yIndex1 < width) {
      odata[index_out] = tile[threadIdx.x][threadIdx.y];
  }
}

// transpose together with spread operation
__global__ void transpose_spread(Complex *idata, Complex *odata, 
	int width, int height, int bign_h)
{
  __shared__ Complex tile[TILE_DIM][TILE_DIM+1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*(width);
  int xIndex1 = blockIdx.y * TILE_DIM + threadIdx.x;
  int yIndex1 = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex1 + (yIndex1)*height;
  int div = index_out / bign_h;
  int mod = index_out % bign_h;
  index_out = div * (bign_h << 1) + mod+bign_h;
  int plus = -bign_h;

  if (xIndex < width && yIndex < height)
      tile[threadIdx.y][threadIdx.x] = idata[index_in];
    __syncthreads();

  if (xIndex1 < height && yIndex1 < width) {
      odata[index_out] = tile[threadIdx.x][threadIdx.y];
    #ifdef DOUBLE
      odata[index_out + plus] = make_double2(0., 0.);
    #else
      odata[index_out + plus] = make_float2(0.f, 0.f);
    #endif
  }
}

// transpose together with inverse spread operation
__global__ void transpose_spread_i(Complex *idata, Complex *odata,
        int width, int height, int bign_h)
{
  __shared__ Complex tile[TILE_DIM][TILE_DIM+1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*(width);
  int xIndex1 = blockIdx.y * TILE_DIM + threadIdx.x;
  int yIndex1 = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex1 + (yIndex1)*height;
  int div = index_in / bign_h;
  int mod = index_in % bign_h;
  index_in = div * (bign_h << 1) + mod;

  if (xIndex < width && yIndex < height)
      tile[threadIdx.y][threadIdx.x] = idata[index_in];
    __syncthreads();

  if (xIndex1 < height && yIndex1 < width)
      odata[index_out] = tile[threadIdx.x][threadIdx.y];
}

// spread operation
__global__ void spread(Real* src, unsigned int spitch, Real* dst, unsigned int dpitch)
{
   unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
   unsigned int tid = threadIdx.x;
 
   Real res = (tid >= spitch) ? src[bid * spitch + tid-spitch] : 0.0;
   if( tid < dpitch) {
	dst[bid * dpitch + tid] = res;
   }
}

// inverse spread operation
__global__ void spread_i(Real* src, unsigned int spitch, Real* dst, unsigned int dpitch)
{
   unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
   unsigned int tid = threadIdx.x;

   Real res = src[bid * dpitch + tid];
   if( tid < dpitch) dst[bid * spitch + tid] = res;
}

// spread operation for 2nd dim
__global__ void spread_y(Complex* src, Complex* dst)
{
   unsigned int tid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
   unsigned int tid1 = (blockIdx.y * gridDim.x * 2 + blockIdx.x) * blockDim.x + threadIdx.x;

   Complex res =  src[tid];
   dst[tid1 + blockDim.x*gridDim.x] = res;
   int plus = -gridDim.x;
#ifdef DOUBLE
   dst[tid1 + blockDim.x*(gridDim.x + plus)] = make_double2(0., 0.);
#else
   dst[tid1 + blockDim.x*(gridDim.x + plus)] = make_float2(0.f, 0.f);
#endif
}

// inverse spread operation for 2nd dim
__global__ void spread_y_i(Complex* src, Complex* dst)
{
   unsigned int tid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
   unsigned int tid1 = (blockIdx.y * gridDim.x * 2 + blockIdx.x) * blockDim.x + threadIdx.x;

   Complex res =  src[tid1];
   dst[tid] = res;
}

// multiply with potential
__global__ void multiply_kernel(int nx, int ny, int nz, Complex *d_data, Real *d_kernel, Real scal) {

 int tj = threadIdx.x;
 int td = blockDim.x;

 int blockData = (nx*ny*nz)/(gridDim.x*gridDim.y);

 int jj = (blockIdx.y*gridDim.x + blockIdx.x)*blockData;

 for (int k=0; k<blockData/td; k++) {
     d_data[jj + tj+ k*td].x *= d_kernel[jj + tj+ k*td]*scal;
     d_data[jj + tj+ k*td].y *= d_kernel[jj + tj+ k*td]*scal;
 }

}

// scale
__global__ void scale_kernel(int nx, int ny, int nz, Complex *d_data, Real mult) {

 int tj = threadIdx.x;
 int td = blockDim.x;

 int blockData = (nx*ny*nz)/(gridDim.x*gridDim.y);

 int jj = (blockIdx.y*gridDim.x + blockIdx.x)*blockData;

 for (int k=0; k<blockData/td; k++) {
     d_data[jj + tj+ k*td].x *= mult;
     d_data[jj + tj+ k*td].y *= mult;
 }

}

__global__ void zero(int nx, int ny, int nz, Real *z) {

        int tj = threadIdx.x;
        int td = blockDim.x;

	int blockData = (nx*ny*nz)/(gridDim.x*gridDim.y);

        int jj = ((blockIdx.y)*gridDim.x + (blockIdx.x))*blockData;

        for (int k=0; k<blockData/td; k++) {
        	z[jj + tj+ k*td] = 0.0;
        }
}

__global__ void copy_0(int nx, int ny, int nz, Real *in, Real *out) {

	int tj = threadIdx.x;
        int td = blockDim.x;

        int jj =  (blockIdx.y*nx*ny/4 + blockIdx.x*nx/2);
        int jj1 =  ((blockIdx.y+nz/2)*nx*ny + (blockIdx.x+ny/2)*nx);


        out[jj1+tj+td] = in[jj+tj];

}

__global__ void copy(int nx,int ny,int nz, Real *in, Real *out) {

        int tj = threadIdx.x;
        //int td = blockDim.x;

        int jj =  (blockIdx.y*nx*ny/4 + blockIdx.x*nx/2);
        int jj1 =  ((blockIdx.y)*nx*ny + (blockIdx.x)*nx);

        out[jj+tj] = in[jj1+tj];
}

/************ 1D transform *************/

extern "C" void cuda_1d_plan_(int *NX_p, int *Nbatch_p,
                 cufftHandle *plan) {

 int NX = *NX_p;
 int Nbatch = *Nbatch_p;

 int n1d[3]= {NX, 1, 1};

 if(cufftPlanMany(plan,  1, n1d,
              NULL, 1, NX,
              NULL, 1, NX, Transform, Nbatch) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 //cufftPlan1d(plan, NX, Transform, Nbatch );

}

extern "C" void cuda_1d_forward_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in 1D forward transform\n");
   }

}

extern "C" void cuda_1d_inverse_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in 1D inverse transform\n");
   }

}

/************ 2D transform *************/

extern "C" void cuda_2d_plan_(int *NX_p, int *NY_p, int *Nbatch_p,
                 cufftHandle *plan) {

 int NX = *NX_p;
 int NY = *NY_p;
 int Nbatch = *Nbatch_p;

 int n1d[3]= {NX, NY, 1};

 if(cufftPlanMany(plan,  1, n1d,
              NULL, 1, NX*NY,
              NULL, 1, NX*NY, Transform, Nbatch) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

}

extern "C" void cuda_2d_forward_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in 2D forward transform\n");
   }

}

extern "C" void cuda_2d_inverse_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in 2D inverse transform\n");
   }

}

/************ 3D transform *************/

extern "C" void cuda_3d_plan_(int *NX_p, int *NY_p, int *NZ_p,
                 cufftHandle *plan) {

 int NX = *NX_p;
 int NY = *NY_p;
 int NZ = *NZ_p;

 int n[3] = { NZ, NY, NX };
 if(cufftPlanMany(plan, 3, n,
              NULL, 1, NX*NY*NZ,
              NULL, 1, NX*NY*NZ, Transform, 1) != CUFFT_SUCCESS)
      printf("Error creating plan\n");
}

extern "C" void cuda_3d_forward_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in 3D forward transform\n");
   }

}

extern "C" void cuda_3d_inverse_(int *NX_p, int *NY_p, int *NZ_p ,cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   int NX = *NX_p;
   int NY = *NY_p;
   int NZ = *NZ_p;

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in 3D inverse transform\n");
   }

   // scale kernel paramters
   int nThreads = NX;
   dim3 nBlocks(NY,NZ,1);

   scale_kernel <<< nBlocks, nThreads >>> (NX,NY,NZ,*d_data2, 1.0/double(NX*NY*NZ));
}

/************ 3D Poisson Solver for periodic boundary *************/

extern "C" void cuda_3d_psolver_per_plan_(int *NX_p, int *NY_p, int *NZ_p,
                 cufftHandle *plan, cufftHandle *plan1) {

 int NX = *NX_p;
 int NY = *NY_p;
 int NZ = *NZ_p;

 int n[3] = { NZ, NY, NX };
 if(cufftPlanMany(plan, 3, n,
              NULL, 1, NX*NY*NZ,
              NULL, 1, NX*NY*NZ, CUFFT_D2Z, 1) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 if(cufftPlanMany(plan1, 3, n,
              NULL, 1, NX*NY*NZ,
              NULL, 1, NX*NY*NZ, CUFFT_Z2D, 1) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

}


extern "C" void cuda_3d_psolver_per_(int *NX_p, int *NY_p, int *NZ_p,cufftHandle *plan,
             cufftHandle *plan1, Complex **d_data, Complex **d_data2, Real **d_kernel, Real *scal_p,
	     int *geo_p) {

 int NX = *NX_p;
 int NY = *NY_p;
 int NZ = *NZ_p;

 int geo = *geo_p;

 Real scal = *scal_p;

 // multiply kernel paramters
 int nThreads = NX/2+1;
 dim3 nBlocks(NY,NZ,1);

 // copy kernel paramters
 int nthreads = NX/2;
 dim3 nblocks(NY/2,NZ/2,1);

 Complex* dst = *d_data;
 Complex* src = *d_data2;

   if (geo==0) {
    src = *d_data;
    dst = *d_data2;
    zero <<< nblocks, nthreads >>> (NX,NY,NZ, (Real*)dst);
    copy_0 <<< nblocks, nthreads  >>> (NX,NY,NZ, (Real*)src, (Real*)dst);
   }

   // Forward FFT

   if( cufftExecD2Z(*plan, (Real*)dst, src)!= CUFFT_SUCCESS){
      printf("error in PSper forward transform\n");
   }

   // multiply with kernel

   multiply_kernel <<< nBlocks, nThreads >>> (NX/2+1,NY,NZ,src,*d_kernel,scal);

   // Inverse FFT

   if( cufftExecZ2D(*plan1, src, (Real*)dst)!= CUFFT_SUCCESS){
      printf("error in PSper inverse transform\n");
   }

   if (geo==0)
     copy <<< nblocks, nthreads >>> (NX,NY,NZ, (Real*)dst, (Real*)src);
}

/************ 3D Poisson Solver for general boundary *************/

extern "C" void cuda_3d_psolver_general_plan_(int *NX_p, int *NY_p, int *NZ_p,
                 cufftHandle *plan1, cufftHandle *plan1_, cufftHandle *plan2,
                 cufftHandle *plan3, cufftHandle *plan3_, int *switch_alg,
		 int *geo1_p, int *geo2_p, int *geo3_p) {

 int NX = *NX_p;
 int NY = *NY_p;
 int NZ = *NZ_p;

 //int geo1 = *geo1_p;
 int geo2 = *geo2_p;
 int geo3 = *geo3_p;

 int n1d[3]= {1, 1, 1};

 int ysize = NY/2 + geo2 * NY/2;
 int zsize = NZ/2 + geo3 * NZ/2;

 n1d[0] = NX;
 if(cufftPlanMany(plan1,  1, n1d,
              NULL, 1, NX,
              NULL, 1, NX, CUFFT_D2Z, ysize*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 if(cufftPlanMany(plan1_,  1, n1d,
              NULL, 1, NX,
              NULL, 1, NX, CUFFT_Z2D, ysize*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 n1d[0] = NY;
 if(cufftPlanMany(plan2,  1, n1d,
              NULL, 1, NY,
              NULL, 1, NY, Transform, (NX/2+1)*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 n1d[0] = NZ;
 if(cufftPlanMany(plan3,  1, n1d,
              NULL, 1, NZ,
              NULL, 1, NZ, Transform, (NX/2+1)*NY) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 int nPrimeSize = 17;
 int primeSize[] = {92,104,116,124,136,148,152,164,172,184,188,204,208,220,228,232,248};

 *switch_alg = 0;
 for (int p=0; p<nPrimeSize; p++)
   if (NZ == primeSize[p]) {
     *switch_alg = 1;
     break;
   }

 n1d[0] = NZ;

 int inembed[1];
 int onembed[1];
 inembed[0] = 1;
 onembed[0] = 1;
 if(cufftPlanMany(plan3_,  1, n1d,
              inembed, NY, 1,
              onembed, NY, 1, Transform, NY) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

}

extern "C" void cuda_3d_psolver_general_(int *NX_p, int *NY_p, int *NZ_p,
          cufftHandle *plan1, cufftHandle *plan1_, cufftHandle *plan2,
          cufftHandle *plan3, cufftHandle *plan3_,
          Complex **d_data, Complex **d_data2, Real **d_kernel, int *switch_alg,
          int *geo1_p, int *geo2_p, int *geo3_p, Real *scal_p) {

 int NX = *NX_p;
 int NY = *NY_p;
 int NZ = *NZ_p;

 Real scal = *scal_p;

 int geo1 = *geo1_p;
 int geo2 = *geo2_p;
 int geo3 = *geo3_p;

 int ysize=NY/2+geo2*NY/2;
 int zsize=NZ/2+geo3*NZ/2;

 // transpose kernel parameters
 dim3 grid((NX+TILE_DIM-1)/TILE_DIM,(ysize*zsize+TILE_DIM-1)/TILE_DIM,1);
 dim3 threads(TILE_DIM,TILE_DIM,1);

 // spread kernel parameters
 dim3 nblocks(zsize,ysize,1);

 // multiply kernel paramters
 int nThreads = NX/2+1;
 dim3 nBlocks(NZ,NY,1);

 Complex* dst = *d_data;
 Complex* src = *d_data2;

 // X transform 

   if (geo1==0) {
     src = *d_data;
     dst = *d_data2;
     spread<<<nblocks, NX>>>((Real*)src, NX/2, (Real*)dst, NX);
   }

   if( cufftExecD2Z(*plan1, (Real*)dst, src)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 1\n");
   }

   if (geo2==0) {
     transpose_spread <<< grid, threads >>>(src, dst,NX/2+1,ysize*zsize,NY/2);
   } else {
     transpose <<< grid, threads >>>(src, dst,NX/2+1,ysize*zsize);
   }

   // Y transform
   if( TransformExec(*plan2, dst, src, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 2\n");
   }

  // Z transform, on entire cube
  if (!(*switch_alg)) {
   grid.x = (NY+TILE_DIM-1)/TILE_DIM;
   grid.y = ((NX/2+1)*zsize+TILE_DIM-1)/TILE_DIM;

   if (geo3==0) {
     transpose_spread <<< grid, threads >>>(src, dst,NY,(NX/2+1)*NZ/2,NZ/2);
   } else {
     transpose <<< grid, threads >>>(src, dst,NY,(NX/2+1)*NZ);
   }

   if( TransformExec(*plan3, dst, src, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 3\n");
   }
  }
  else {
   if (geo3==0) {
      nblocks.x=zsize;
      nblocks.y=NX;
      spread_y<<<nblocks, NY>>>(src, dst);
   }

   for(int k=0; k<NX; ++k){
     if( TransformExec(*plan3_, dst, src, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 3\n");
     }
     src += NY*NZ;
     dst += NY*NZ;
   }

   src -= NX*NY*NZ;
   dst -= NX*NY*NZ;
  }

  // multiply with kernel

  multiply_kernel <<< nBlocks, nThreads >>> (NX/2+1,NY,NZ,src,*d_kernel,scal);

  // inverse transform

  // Z transform, on entire cube 
  if (!(*switch_alg)) {
   if( TransformExec(*plan3, src, dst, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 1\n");
   }

   grid.x = (zsize*(NX/2+1)+TILE_DIM-1)/TILE_DIM;
   grid.y = (NY+TILE_DIM-1)/TILE_DIM;

   if (geo3==0) {
     transpose_spread_i <<< grid, threads >>>(dst,src,NZ/2*(NX/2+1),NY,NZ/2);
   } else {
     transpose <<< grid, threads >>>(dst, src,NZ*(NX/2+1),NY);
   }

  }
  else {

   for(int k=0; k<NX; ++k){
     if( TransformExec(*plan3_, src, dst, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 3\n");
     }
     src += NY*NZ;
     dst += NY*NZ;
   }

   src -= NX*NY*NZ;
   dst -= NX*NY*NZ;

   if (geo3==0)
      spread_y_i<<<nblocks, NY>>>(dst, src);
  }

  // Y transform

   if( TransformExec(*plan2, src, dst, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 2\n");
   }

   grid.x = (ysize*zsize+TILE_DIM-1)/TILE_DIM;
   grid.y = (NX/2+1+TILE_DIM-1)/TILE_DIM;

   if (geo2==0) {
      transpose_spread_i <<< grid, threads >>>(dst, src,ysize*zsize,NX/2+1, NY/2);
   } else
      transpose <<< grid, threads >>>(dst, src,ysize*zsize,NX/2+1);

   // X transform

   if( cufftExecZ2D(*plan1_, src, (Real*)dst)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 3\n");
   }

   nblocks.x=zsize;
   nblocks.y=ysize;
   if (geo1==0) {
      spread_i<<<nblocks, NX/2>>>((Real*)dst, NX/2, (Real*)src, NX);
   }

   //scale_kernel<<< nBlocks, nThreads >>> (NX/2+1,NY,NZ,(Real*)dst,scal); 
}


