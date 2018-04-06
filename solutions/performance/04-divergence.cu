// Copyright 2014, Cranfield University
// All rights reserved
// Author: Michał Czapiński (mczapinski@gmail.com)
//
// A program demonstrating the impact of warp divergence on the overall performance.
// It can only be observed when enough blocks are run (typicaly more than 10000).

#include <iostream>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_timer.h>

// Size of each block of threads on the GPU.
const int BLOCK_DIM = 32;

// Number of blocks to run.
const dim3 GRID_DIM(1024, 1024);

__device__ float result[BLOCK_DIM];

/*! GPU kernel executing the given number of divergent threads -
 *  (divergentThreads + 1) execution paths.
 *
 * 	\param divergentThreads number of threads diverging from the main execution path.
 */
__global__ void WarpDivergence(int divergentThreads) {
  int tid = threadIdx.x % warpSize;
  if (tid < divergentThreads) {
    switch (tid) {
      case 0:
        result[tid] = expf(0.1f);
        break;
      case 1:
        result[tid] = expf(0.2f);
        break;
      case 2:
        result[tid] = expf(0.3f);
        break;
      case 3:
        result[tid] = expf(0.4f);
        break;
      case 4:
        result[tid] = expf(0.5f);
        break;
      case 5:
        result[tid] = expf(0.6f);
        break;
      case 6:
        result[tid] = expf(0.7f);
        break;
      case 7:
        result[tid] = expf(0.8f);
        break;
      case 8:
        result[tid] = expf(0.9f);
        break;
      case 9:
        result[tid] = expf(0.11f);
        break;
      case 10:
        result[tid] = expf(0.12f);
        break;
      case 11:
        result[tid] = expf(0.13f);
        break;
      case 12:
        result[tid] = expf(0.14f);
        break;
      case 13:
        result[tid] = expf(0.15f);
        break;
      case 14:
        result[tid] = expf(0.16f);
        break;
      case 15:
        result[tid] = expf(0.17f);
        break;
      case 16:
        result[tid] = expf(0.18f);
        break;
      case 17:
        result[tid] = expf(0.19f);
        break;
      case 18:
        result[tid] = expf(0.21f);
        break;
      case 19:
        result[tid] = expf(0.22f);
        break;
      case 20:
        result[tid] = expf(0.23f);
        break;
      case 21:
        result[tid] = expf(0.24f);
        break;
      case 22:
        result[tid] = expf(0.25f);
        break;
      case 23:
        result[tid] = expf(0.26f);
        break;
      case 24:
        result[tid] = expf(0.27f);
        break;
      case 25:
        result[tid] = expf(0.28f);
        break;
      case 26:
        result[tid] = expf(0.29f);
        break;
      case 27:
        result[tid] = expf(0.31f);
        break;
      case 28:
        result[tid] = expf(0.32f);
        break;
      case 29:
        result[tid] = expf(0.33f);
        break;
      case 30:
        result[tid] = expf(0.34f);
        break;
      case 31:
        result[tid] = expf(0.35f);
        break;
    }
  } else {
    result[tid] = expf(1.0f);
  }
}

int main(int argc, char** argv) {
  checkCudaErrors(cudaSetDevice(gpuGetMaxGflopsDeviceId()));

  StopWatchInterface* timer = NULL;
  sdkCreateTimer(&timer);

  // This is to "wake up" the GPU and avoid spurious overhead on the first time measurement.
  checkCudaErrors(cudaDeviceSynchronize());

  for (int i = 0; i < BLOCK_DIM; ++i) {
    timer->reset();
    timer->start();

    WarpDivergence<<< GRID_DIM, BLOCK_DIM >>>(i);
    checkCudaErrors(cudaDeviceSynchronize());

    timer->stop();
    std::cout << (i + 1) << " execution path(s): " << timer->getTime() << " ms" << std::endl;
  }

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;
  return 0;
}
