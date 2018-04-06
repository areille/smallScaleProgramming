// Copyright 2014, Cranfield University
// All rights reserved
// Author: Michał Czapiński (mczapinski@gmail.com)
//
// A program demonstrating the use of shared memory.
// It is based on a GPU reduction example from a presentation by
// Mark Haris "Optimizing Parallel Reduction in CUDA",
// http://developer.download.nvidia.com/assets/cuda/files/reduction.pdf

#include <iostream>

#include <cuda_runtime.h>  // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers

// Number of elements in the vector to reduce.
const int N = 128 * 1024 * 1024;

// Size of each block of threads on the GPU (must be a power of 2!)
const dim3 BLOCK_DIM = 256;

/*! Sequential CPU implementation of a global reduction (sum).
 *
 *  \param n number of integers in input array.
 * 	\param input array of integers to add up.
 * 	\return sum of all the integers in input array.
 */
int CpuReduction(int n, const int* input) {
  int sum = 0;
  for (int i = 0; i < n; ++i) {
    sum += input[i];
  }
  return sum;
}

/*! GPU implementation of a partial reduction (sum)
 *  Only integers in the same block are summed up.
 *
 *  \param n number of integers in input array.
 *  \param input array of integers to add up.
 * 	\param output array for the partial sums for each block.
 */
__global__ void GpuReduction(int n, const int* input, int* output) {
  int tid = threadIdx.x;
  int bid = blockIdx.y * gridDim.x + blockIdx.x;  // Block index
  int gid = bid * blockDim.x + tid;  // Index in global input array

  // Dynamic shared memory declaration. Size is determined at the kernel invocation.
  extern __shared__ int aux[];

  // Load the data from global to shared memory.
  // Don't read outside the allocated global memory. Adding 0 doesn't change the result.
  aux[tid] = (gid < n ? input[gid] : 0);
  __syncthreads();  // Wait for all the threads to read their data.

  // Add blockDim.x numbers in shared memory and store the result in the first element.
  for (int s = 1; s < blockDim.x; s *= 2) {
    if (tid % (2 * s) == 0) {
      aux[tid] += aux[tid + s];
    }
    __syncthreads();
  }

  // Need to write only one output value per block. Would also work if all the threads
  // write the output value, but it could be slower. Why?
  if (tid == 0) {
    output[bid] = aux[0];
  }
}

int main(int argc, char** argv) {

// ----------------------- Host memory initialisation ----------------------- //

  srand(time(0));
  int* numbers = new int[N];
  for (int i = 0; i < N; ++i) {
    numbers[i] = 1 + rand() % 10;
  }

// ---------------------- Device memory initialisation ---------------------- //

  int* dev_old = 0;
  checkCudaErrors(cudaMalloc((void**) &dev_old, N * sizeof(int)));
  checkCudaErrors(cudaMemcpy(dev_old, numbers, N * sizeof(int), cudaMemcpyHostToDevice));

  int* dev_new = 0;
  checkCudaErrors(cudaMalloc((void**) &dev_new, N * sizeof(int)));

// --------------------- Calculations for CPU implementation ---------------- //

  // Create the CUDA SDK timer.
  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  timer->start();
  int result = CpuReduction(N, numbers);

  timer->stop();
  std::cout << "CPU: " << timer->getTime() << " ms, sum = " << result << std::endl;

// ---------------------- Kernel launch specifications ---------------------- //

  // The GPU kernel performs a partial reduction only and needs to be run
  // multiple times, therefore grid_dim will need to be updated.
  dim3 grid_dim;
  const size_t shmem_size = BLOCK_DIM.x * sizeof(int);

// --------- Calculations for GPU implementation with bank conflicts -------- //

  timer->reset();
  timer->start();

  // current_n holds the number of elements remaining to be reduced.
  int current_n = N;
  while (current_n > 1) {
    int blocks_required = (current_n - 1) / BLOCK_DIM.x + 1;
    grid_dim.x = static_cast<int>(ceil(sqrt(blocks_required)));
    grid_dim.y = ((blocks_required - 1) / grid_dim.x) + 1;

    GpuReduction<<< grid_dim, BLOCK_DIM, shmem_size >>>(current_n, dev_old, dev_new);

    std::swap(dev_new, dev_old);
    current_n = blocks_required;
  }
  checkCudaErrors(cudaDeviceSynchronize());
  timer->stop();

  // Copy the result (single int value) from the device.
  checkCudaErrors(cudaMemcpy(&result, dev_old, sizeof(int), cudaMemcpyDeviceToHost));
  std::cout << "GPU: " << timer->getTime() << " ms, sum = " << result << std::endl;

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;

  checkCudaErrors(cudaFree(dev_old));
  checkCudaErrors(cudaFree(dev_new));

  delete[] numbers;
  return 0;
}
