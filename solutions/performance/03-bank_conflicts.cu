// Copyright 2014, Cranfield University
// All rights reserved
// Author: Michał Czapiński (mczapinski@gmail.com)
//
// A program demonstrating the impact of shared memory bank conflicts on
// the overall performance. It is based on a GPU reduction example from a
// presentation by Mark Haris "Optimizing Parallel Reduction in CUDA",
// http://developer.download.nvidia.com/assets/cuda/files/reduction.pdf

#include <iostream>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_timer.h>

// Numer of elements in the vector to reduce.
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

/*! GPU implementation of a partial reduction (sum) - only integers in the same
 *  block are summed up. Number of shared memory banks used decreases as the
 *  stride length s increases, causing bank conflicts.
 *
 *  \param n number of integers in input array.
 *  \param input array of integers to add up.
 * 	\param output array for the partial sums for each block.
 */
__global__ void GpuReductionConflicts(int n, const int* input, int* output) {
  int tid = threadIdx.x;
  int bid = blockIdx.y * gridDim.x + blockIdx.x;
  int gid = bid * blockDim.x + tid;

  extern __shared__ int aux[];

  // Don't read outside allocated global memory. Adding 0 doesn't change the result.
  aux[tid] = (gid < n ? input[gid] : 0);
  __syncthreads();

  // Add blockDim.x numbers in shared memory and store the result in the first element.
  for (int s = 1; s < blockDim.x; s *= 2) {
    int index = 2 * s * tid;

    // The number of threads performing reduction decreases, as the stride increases.
    if (index < blockDim.x) {
      aux[index] += aux[index + s];
    }
    __syncthreads();
  }

  if (tid == 0) {
    output[bid] = aux[0];
  }
}

/*! GPU implementation of a partial reduction (sum) - only integers in the same
 *  block are summed up. No shared memory bank conflicts occur.
 *
 *  \param n number of integers in input array.
 * 	\param input array of integers to add up.
 * 	\param output array for the partial sums for each block.
 */
__global__ void GpuReductionNoConflicts(int n, const int* input, int* output) {
  int tid = threadIdx.x;
  int bid = blockIdx.y * gridDim.x + blockIdx.x;
  int gid = bid * blockDim.x + tid;

  extern __shared__ int aux[];

  // Don't read outside allocated global memory. Adding 0 doesn't change the result.
  aux[tid] = (gid < n ? input[gid] : 0);
  __syncthreads();

  // >> operator is a bit-wise shift to the right, i.e. '>> 1' is integer division by two.
  for (int s = (blockDim.x >> 1); s > 0; s >>= 1) {
    if (tid < s) {
      // Threads access consecutive addresses in shared memory, i.e. no bank conflict occurs
      aux[tid] += aux[tid + s];
    }
    __syncthreads();
  }

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

  checkCudaErrors(cudaSetDevice(gpuGetMaxGflopsDeviceId()));

  int* dev_old = 0;
  checkCudaErrors(cudaMalloc((void**) &dev_old, N * sizeof(int)));
  checkCudaErrors(cudaMemcpy(dev_old, numbers, N * sizeof(int), cudaMemcpyHostToDevice));

  int* dev_new = 0;
  checkCudaErrors(cudaMalloc((void**) &dev_new, N * sizeof(int)));

// --------------------- Calculations for CPU implementation ---------------- //

  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  timer->start();
  int result = CpuReduction(N, numbers);

  timer->stop();
  std::cout << "CPU:                     " << timer->getTime() << " ms, sum = " << result << std::endl;

// ---------------------- Kernel launch specifications ---------------------- //

  // The GPU kernel performs a partial reduction only and needs to be run
  // multiple times, therefore grid_dim will need to be updated.
  dim3 grid_dim;
  const size_t shmem_size = BLOCK_DIM.x * sizeof(int);

// --------- Calculations for GPU implementation with bank conflicts -------- //

  timer->reset();
  timer->start();

  int current_n = N;
  while (current_n > 1) {
    int blocks_required = (current_n - 1) / BLOCK_DIM.x + 1;
    grid_dim.x = static_cast<int>(ceil(sqrt(blocks_required)));
    grid_dim.y = ((blocks_required - 1) / grid_dim.x) + 1;

    GpuReductionConflicts
        <<< grid_dim, BLOCK_DIM, shmem_size >>>(current_n, dev_old, dev_new);

    std::swap(dev_new, dev_old);
    current_n = blocks_required;
  }
  checkCudaErrors(cudaDeviceSynchronize());
  timer->stop();

  checkCudaErrors(cudaMemcpy(&result, dev_old, sizeof(int), cudaMemcpyDeviceToHost));
  std::cout << "GPU (bank conflicts):    " << timer->getTime() << " ms, sum = " << result << std::endl;

// ------- Calculations for GPU implementation without bank conflicts ------- //

  // Start with the same set of numbers on the GPU.
  checkCudaErrors(cudaMemcpy(dev_old, numbers, N * sizeof(int), cudaMemcpyHostToDevice));

  timer->reset();
  timer->start();

  current_n = N;
  while (current_n > 1) {
    int blocks_required = (current_n - 1) / BLOCK_DIM.x + 1;
    grid_dim.x = static_cast<int>(ceil(sqrt(blocks_required)));
    grid_dim.y = ((blocks_required - 1) / grid_dim.x) + 1;

    GpuReductionNoConflicts
        <<< grid_dim, BLOCK_DIM, shmem_size >>>(current_n, dev_old, dev_new);

    std::swap(dev_new, dev_old);
    current_n = blocks_required;
  }
  checkCudaErrors(cudaDeviceSynchronize());
  timer->stop();

  checkCudaErrors(cudaMemcpy(&result, dev_old, sizeof(int), cudaMemcpyDeviceToHost));
  std::cout << "GPU (no bank conflicts): " << timer->getTime() << " ms, sum = " << result << std::endl;

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;

  checkCudaErrors(cudaFree(dev_old));
  checkCudaErrors(cudaFree(dev_new));
  delete[] numbers;

  // Required for cuda-memcheck to detect leaks correctly.
  checkCudaErrors(cudaDeviceReset());
  return 0;
}
