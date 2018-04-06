// Copyright 2014, Cranfield University
// All rights reserved
// Author: Michał Czapiński (mczapinski@gmail.com)
//
// A program demonstrating the impact of non-coalesced global memory reads on
// the GPU kernel performance. It contains CPU and GPU implementations
// using two data layouts: array of structures, and structure of arrays.

#include <iostream>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_timer.h>

// Size of each block of threads on the GPU.
const dim3 BLOCK_DIM = 192;

// Number of blocks to run.
const int BLOCKS_COUNT = 8 * 1024;

// Number of vectors.
const int N = BLOCKS_COUNT * BLOCK_DIM.x;

// Dimensionality of each vector.
const int DIM = 8;

// TODO Play with DIM value (1-32) - higher values increase the stride in
// array of structures implementation.
// How does it affect the CPU performance? How does it affect the GPU performance?

// Tolerance used to compare floating-point numbers.
const float EPS = 1e-3;

#define POW2(x) (x) * (x)

/*! Calculates 2-norms of all the vectors in an array (CPU implementation
 *  using array of structures).
 *
 *  \param n number of vectors in an array.
 *  \param dim dimensionality of vectors in an array.
 *  \param vecs vectors values in array of structures format (n * dim elements).
 *  \param norms vector for the computed 2-norms.
 */
void CpuNorm2Aos(int n, int dim, const float* vecs, float* norms) {
  for (int i = 0; i < n; ++i) {
    norms[i] = 0.0f;
    for (int d = 0; d < dim; ++d) {
      norms[i] += POW2(vecs[i * dim + d]);
    }
    norms[i] = sqrtf(norms[i]);
  }
}

/*! Method calculates 2-norms of all the vectors in an array (CPU implementation
 *  using structure of arrays).
 *
 *  \param n number of vectors in an array.
 *  \param dim dimensionality of vectors in an array.
 *  \param vecs vectors values in structure of arrays format (n*dim elements).
 *  \param norms vector for the computed 2-norms.
 */
void CpuNorm2Soa(int n, int dim, const float* vecs, float* norms) {
  for (int i = 0; i < n; ++i) {
    norms[i] = 0.0f;
    for (int d = 0; d < dim; ++d) {
      norms[i] += POW2(vecs[d * n + i]);
    }
    norms[i] = sqrtf(norms[i]);
  }
}

/*! GPU kernel calculating 2-norms of all the vectors in an array (using global
 *  memory and array of structures layout).
 *
 *  \param n number of vectors in an array.
 *  \param dim dimensionality of vectors in an array.
 *  \param vecs vectors values in structure of arrays format (n*dim elements).
 *  \param norms vector for the computed 2-norms.
 */
__global__ void GpuNorm2Aos(int n, int dim, const float* vecs, float* norms) {
  // Using 2D grid of blocks enables handling data of arbitrary size.
	int gid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	if (gid >= n) return;

	float norm = 0.0f;
	for (int i = 0; i < dim; ++i) {
	  // The first element is at grid * dim position, then the stride between elements is equal to one.
		float tmp = vecs[gid * dim + i];
		norm += POW2(tmp);
	}
	norms[gid] = sqrtf(norm);
}

/*! GPU kernel calculating 2-norms of all the vectors in an array (using global
 *  memory and structure of arrays layout).
 *
 *  \param n number of vectors in an array.
 *  \param dim dimensionality of vectors in an array.
 *  \param vecs vectors values in structure of arrays format (n*dim elements).
 *  \param norms vector for the computed 2-norms.
 */
__global__ void GpuNorm2Soa(int n, int dim, const float* vecs, float* norms) {
  // Using 2D grid of blocks enables handling data of arbitrary size.
	int gid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	if (gid >= n) return;

	float norm = 0.0f;
	for (int i = 0; i < dim; ++i) {
	  // The first element is at gid position, then the stride between elements is equal to n.
		float tmp = vecs[gid + i * n];
		norm += POW2(tmp);
	}
	norms[gid] = sqrtf(norm);
}

void CheckGpuResult(int n, const float* cpu_norms, const float* gpu_norms) {
  int errors = 0;
  for (int i = 0; i < n; ++i) {
    if (fabs(gpu_norms[i] - cpu_norms[i]) >= EPS) {
      std::cerr << i << "\t" << gpu_norms[i] << " != " << cpu_norms[i] << std::endl;
      ++errors;
    }
  }

  if (errors) {
    std::cerr << "ERROR! " << errors << " wrong values have been encountered." << std::endl;
  } else {
    std::cout << "SUCCESS! The GPU result is the same as the CPU result." << std::endl;
  }
}

int main(int argc, char** argv) {

// ----------------------- Host memory initialisation ----------------------- //

  float* vecs = new float[DIM * N];
  float* cpu_norms = new float[N];
  float* gpu_norms = new float[N];

  srand(time(0));
  for (int i = 0; i < DIM * N; ++i) {
    vecs[i] = static_cast<float>(rand()) / RAND_MAX;
  }

// ---------------------- Device memory initialisation ---------------------- //

  checkCudaErrors(cudaSetDevice(gpuGetMaxGflopsDeviceId()));

  float* dev_vecs = 0;
  checkCudaErrors(cudaMalloc((void**) &dev_vecs, DIM * N * sizeof(float)));
  checkCudaErrors(cudaMemcpy(dev_vecs, vecs, DIM * N * sizeof(float), cudaMemcpyHostToDevice));

  float* dev_norms = 0;
  checkCudaErrors(cudaMalloc((void**) &dev_norms, N * sizeof(float)));

// ---------------------- Kernel launch specifications ---------------------- //

  // Make grid dimensions as close to square as possible.
  dim3 grid_dim;
  grid_dim.x = static_cast<int>(ceil(sqrt(BLOCKS_COUNT)));
  grid_dim.y = ((BLOCKS_COUNT - 1) / grid_dim.x) + 1;

// ------ Calculations for array of structures (global memory version) ------ //

  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  std::cout << "Array of structures" << std::endl;

  timer->start();
  CpuNorm2Aos(N, DIM, vecs, cpu_norms);

  timer->stop();
  std::cout << "CPU: " << timer->getTime() << " ms" << std::endl;

  timer->reset();
  timer->start();
  GpuNorm2Aos<<< grid_dim, BLOCK_DIM >>>(N, DIM, dev_vecs, dev_norms);
  checkCudaErrors(cudaDeviceSynchronize());

  timer->stop();
  std::cout << "GPU: " << timer->getTime() << " ms" << std::endl;

  // Getting results from the GPU.
  checkCudaErrors(cudaMemcpy(gpu_norms, dev_norms, N * sizeof(float), cudaMemcpyDeviceToHost));
  CheckGpuResult(N, cpu_norms, gpu_norms);

// ------ Calculations for structure of arrays (global memory version) ------ //

  std::cout << std::endl << "Structure of arrays" << std::endl;

  timer->reset();
  timer->start();
  CpuNorm2Soa(N, DIM, vecs, cpu_norms);

  timer->stop();
  std::cout << "CPU: " << timer->getTime() << " ms" << std::endl;

  timer->reset();
  timer->start();
  GpuNorm2Soa<<< grid_dim, BLOCK_DIM >>>(N, DIM, dev_vecs, dev_norms);
  checkCudaErrors(cudaDeviceSynchronize());

  timer->stop();
  std::cout << "GPU: " << timer->getTime() << " ms" << std::endl;

  // Getting results from the GPU.
  checkCudaErrors(cudaMemcpy(gpu_norms, dev_norms, N * sizeof(float), cudaMemcpyDeviceToHost));
  CheckGpuResult(N, cpu_norms, gpu_norms);

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;

  checkCudaErrors(cudaFree(dev_vecs));
  checkCudaErrors(cudaFree(dev_norms));

  delete[] vecs;
  delete[] cpu_norms;
  delete[] gpu_norms;

  // Required for cuda-memcheck to detect leaks correctly.
  checkCudaErrors(cudaDeviceReset());
  return 0;
}
