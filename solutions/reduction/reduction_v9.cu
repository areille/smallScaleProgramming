// Copyright 2018, Cranfield University
// All rights reserved
// Author: Salvatore Filippone

#include <iostream>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_timer.h>

const int VERSION=9;

// Computes the reduction using the CPU.
double Reduction(int n, const double* x) {
  double result = 0.0f;
  for (int i = 0; i < n; ++i) {
    result += x[i] ;
  }
  return result;
}
// Memory management for device side
const int BLOCKS_PER_MP = 32; // Sufficiently large for memory transaction hiding
int thread_block = 0; // Must be a power of 2 >= 64
int max_blocks=0;   // Blocks in a grid
int red_sz=0;       // Size of reduction buffer
double *o_data=NULL, *d_res_data=NULL, *h_res_data=NULL;
static struct cudaDeviceProp *prop=NULL;
void reduce_alloc_wrk()
{
  int mpCnt;
  if (prop == NULL) {
    if ((prop=(struct cudaDeviceProp *) malloc(sizeof(struct cudaDeviceProp)))==NULL) {
      fprintf(stderr,"CUDA Error gpuInit3: not malloced prop\n");
      return;
    }
    cudaSetDevice(0); // BEWARE: you may have more than one device
    cudaGetDeviceProperties(prop,0); 
  }
  if (thread_block <= 0)
    std::cerr << "thread_block must be a power of 2 between 64 and 1024" << std::endl;

  if (max_blocks == 0) {
    mpCnt      = prop->multiProcessorCount;
    max_blocks = mpCnt*BLOCKS_PER_MP;
    // Enough to do the second-level reduction
    red_sz     = (max_blocks+thread_block-1)/thread_block;
    //std::cerr << mpCnt << ' '<<max_blocks << ' '<<thread_block<< std::endl;
  }
  if (o_data == NULL) cudaMalloc(&o_data,max_blocks*sizeof(double));
  if (d_res_data == NULL) cudaMalloc(&d_res_data,(red_sz)*sizeof(double));
  if (h_res_data == NULL) h_res_data = (double *)malloc((red_sz)*sizeof(double));
}

// Fully unrolled. Assuming thread_block >= 64
__device__ void warpReduce(volatile double *sdata, int tid) {
  sdata[tid] += sdata[tid + 32];
  sdata[tid] += sdata[tid + 16];
  sdata[tid] += sdata[tid +  8];
  sdata[tid] += sdata[tid +  4];
  sdata[tid] += sdata[tid +  2];
  sdata[tid] += sdata[tid +  1];
}

template <unsigned int THD> __global__ void reduce(int n, double *g_idata, double *g_odata) {
  extern __shared__ double sdata[];
  // each thread loads one element from global to shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
  unsigned int gridSize = blockDim.x * gridDim.x;
  sdata[tid] = 0.0;
  while (i<n) {
    sdata[tid] += g_idata[i] ;
    i += gridSize;
  }
  __syncthreads();
  // do reduction in shared mem
  if (THD >= 1024){ if (tid < 512) { sdata[tid] += sdata[tid + 512]; }  __syncthreads();  }
  if (THD >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; }  __syncthreads();  }
  if (THD >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; }  __syncthreads();  }
  if (THD >= 128) { if (tid < 64)  { sdata[tid] += sdata[tid +  64]; }  __syncthreads();  }
  // write result for this block to global mem
  if (tid < 32) warpReduce(sdata,tid);             
  if (tid == 0) g_odata[blockIdx.x] += sdata[0];

}


void do_gpu_reduce(int n, double *g_idata, double *g_odata)
{
  const int shmem_size   = thread_block*sizeof(double);
  int   nblocks = ((n + thread_block - 1) / thread_block);
  if (nblocks > max_blocks) nblocks = max_blocks;

  switch(thread_block) {
  case 1024:
    reduce<1024><<<nblocks,1024,shmem_size,0>>>(n,g_idata,g_odata); break;
  case 512:
    reduce<512><<<nblocks,512,shmem_size,0>>>(n,g_idata,g_odata); break;
  case 256:
    reduce<256><<<nblocks,256,shmem_size,0>>>(n,g_idata,g_odata); break;
  case 128:
    reduce<128><<<nblocks,128,shmem_size,0>>>(n,g_idata,g_odata); break;
  case 64:
    reduce<64><<<nblocks,64,shmem_size,0>>>(n,g_idata,g_odata); break;
  default:
    std::cerr << "thread_block must be a power of 2 between 64 and 1024" << std::endl;
  }
  return;
  
}

double gpu_reduce(int n, double *d_v)
{
  reduce_alloc_wrk();
  cudaMemset((void *) o_data, 0, max_blocks*sizeof(double));
  cudaMemset((void *)d_res_data,0,(red_sz)*sizeof(double));

  do_gpu_reduce(n, d_v, o_data);
  do_gpu_reduce(max_blocks,o_data,d_res_data);
  cudaError_t err = cudaMemcpy(h_res_data, d_res_data,
			       red_sz*sizeof(double), cudaMemcpyDeviceToHost);
  return(Reduction(red_sz,h_res_data));
  
}
  


// Returns a random number from range [0, 1).
double rand_double() {
  return static_cast<double>(rand()) / RAND_MAX;
}

int main(int argc, char** argv) {

  if (argc < 3) { 
    std::cerr << "Usage: " <<argv[0] << " N  Threads_per_block" << std::endl;
    exit(1);
  }

  int N = atoi(argv[1]);
  thread_block = atoi(argv[2]);
  switch(thread_block) {
  case 1024:
  case 512:
  case 256:
  case 128:
  case 64:
    break;
  default: 
    std::cerr << "thread_block must be a power of 2 between 64 and 1024" << std::endl;
    exit(1);
  }

  double bdwdth;
  double *h_x=(double *) malloc(N*sizeof(double));
  double *d_x;
  
  srand(time(0));

  for (int i=0; i<N; i++)
    h_x[i]=rand_double();
  
  cudaError_t err=cudaMalloc((void **)&d_x,(N*sizeof(double)));
  err = cudaMemcpy(d_x, h_x, N*sizeof(double), cudaMemcpyHostToDevice);
  reduce_alloc_wrk();

  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  std::cout << "Testing reduction algorithm " << VERSION << "  on a DOUBLE vector of size:  " << N << std::endl;
  // Calculate the reduction on the host.
  timer->start();
  double cpu_sum = Reduction(N, h_x);

  timer->stop();
  std::cout << "CPU reduction:            " << cpu_sum
      << " " << timer->getTime() << " ms. " << std::endl;


// ------------ GPU reduction
  timer->reset();
  timer->start();
  
  double gpu_sum = gpu_reduce(N,d_x);

  timer->stop();
  bdwdth = ((double)N*sizeof(double))/timer->getTime();
  bdwdth *= 1.e-6;
  std::cout << "GPU reduction:            " << gpu_sum
      << " " << timer->getTime() << " ms. " << std::endl;
  
  std::cout << "Relative difference:            " << abs(gpu_sum-cpu_sum)/gpu_sum << std::endl;
  std::cout << "Measured bandwidth:             " << bdwdth << " GB/s" << std::endl;
  

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;

  checkCudaErrors(cudaDeviceReset());
  return 0;
}
