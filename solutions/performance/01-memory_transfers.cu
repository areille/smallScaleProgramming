// Copyright 2014, Cranfield University
// All rights reserved
// Author: Michał Czapiński (mczapinski@gmail.com)
//
// Demonstrates the use of page-locked host memory (to speed up CPU-GPU transfers).
// Also demonstrates handling of 2D arrays: allocating row-aligned arrays (to meet the
// global memory coalescing requirements), and the use of cudaMemcpy2D routine.

#include <iostream>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_timer.h>

// Matrix dimensions.
const int COLS = 4000;
const int ROWS = 4096;

// TODO Play with different COLS values.
// What happens if it naturaly meets the alignment requirements (multiple of 128)?
// What happens if it does not? What happens if COLS becomes small? Why?

int main(int argc, char** argv) {

// ----------------------- Host memory initialisation ----------------------- //

  // Allocate h_A in page-locked memory.
  // It is accessed in the same way as memory allocated with malloc or new[].
  float* h_A = 0;
  checkCudaErrors(cudaHostAlloc((void**) &h_A, ROWS * COLS * sizeof(float), cudaHostAllocDefault));

  // h_B is allocated in pageable memory.
  float* h_B = new float[ROWS * COLS];

  srand(time(0));
  for (int row = 0; row < ROWS; ++row) {
    for (int col = 0; col < COLS; ++col) {
      int idx = row * COLS + col;
      h_A[idx] = static_cast<float>(rand()) / RAND_MAX;
      h_B[idx] = static_cast<float>(rand()) / RAND_MAX;
    }
  }

// ---------------- Device memory allocation (no alignment) ----------------- //

  checkCudaErrors(cudaSetDevice(gpuGetMaxGflopsDeviceId()));

  // Allocate global memory on the GPU.
  float *d_A, *d_B;

  const size_t BYTES = ROWS * COLS * sizeof(float);
  checkCudaErrors(cudaMalloc((void**) &d_A, BYTES));
  checkCudaErrors(cudaMalloc((void**) &d_B, BYTES));

// ----------------- Upload data to the GPU (no alignment) ------------------ //

  std::cout << "Uploading linear data to the GPU (cudaMemcpy)" << std::endl;

  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  timer->start();
  checkCudaErrors(cudaMemcpy(d_A, h_A, BYTES, cudaMemcpyHostToDevice));

  timer->stop();
  std::cout << "Page-locked memory: " << timer->getTime() << " ms ("
      << 1e-6 * BYTES / timer->getTime() << " GB/s)" << std::endl;

  timer->reset();
  timer->start();
  checkCudaErrors(cudaMemcpy(d_B, h_B, BYTES, cudaMemcpyHostToDevice));

  timer->stop();
  std::cout << "Pageable memory:    " << timer->getTime() << " ms ("
      << 1e-6 * BYTES / timer->getTime() << " GB/s)" << std::endl << std::endl;

// --------------- Download data from the GPU (no alignment) ---------------- //

  std::cout << "Downloading linear data from the GPU (cudaMemcpy)" << std::endl;

  timer->reset();
  timer->start();
  checkCudaErrors(cudaMemcpy(h_A, d_A, BYTES, cudaMemcpyDeviceToHost));
  timer->stop();

  std::cout << "Page-locked memory: " << timer->getTime() << " ms ("
        << 1e-6 * BYTES / timer->getTime() << " GB/s)" << std::endl;

  timer->reset();
  timer->start();
  checkCudaErrors(cudaMemcpy(h_B, d_B, BYTES, cudaMemcpyDeviceToHost));
  timer->stop();

  std::cout << "Pageable memory:    " << timer->getTime() << " ms ("
      << 1e-6 * BYTES / timer->getTime() << " GB/s)" << std::endl << std::endl;

// --------------- Device memory allocation (with alignment) ---------------- //

  checkCudaErrors(cudaFree(d_A));
  checkCudaErrors(cudaFree(d_B));

  size_t pitch = 0;
  checkCudaErrors(cudaMallocPitch((void**) &d_A, &pitch, COLS * sizeof(float), ROWS));
  checkCudaErrors(cudaMallocPitch((void**) &d_B, &pitch, COLS * sizeof(float), ROWS));
  pitch /= sizeof(float);  // Pitch in float elements rather than in bytes.
  // Because both arrays have the same width (COLS * sizeof(float) bytes), pitch will be the same.

  std::cout << "Matrix has " << COLS << " columns per row." << std::endl;
  std::cout << "To meet the row alignment requirements " << pitch
      << " elements per row have been allocated." << std::endl << std::endl;

// ---------------- Upload data to the GPU (with alignment) ----------------- //

  std::cout << "Uploading row-aligned data to the GPU (cudaMemcpy2D)" << std::endl;

  timer->reset();
  timer->start();
  checkCudaErrors(cudaMemcpy2D(d_A, pitch * sizeof(float), h_A, COLS * sizeof(float),
                               COLS * sizeof(float), ROWS, cudaMemcpyHostToDevice));

  timer->stop();
  std::cout << "Page-locked memory: " << timer->getTime() << " ms ("
      << 1e-6 * BYTES / timer->getTime() << " GB/s)" << std::endl;

  timer->reset();
  timer->start();
  checkCudaErrors(cudaMemcpy2D(d_B, pitch * sizeof(float), h_B, COLS * sizeof(float),
                               COLS * sizeof(float), ROWS, cudaMemcpyHostToDevice));

  timer->stop();
  std::cout << "Pageable memory:    " << timer->getTime() << " ms ("
      << 1e-6 * BYTES / timer->getTime() << " GB/s)" << std::endl << std::endl;

// -------------- Download data from the GPU (with alignment) --------------- //

  std::cout << "Downloading row-aligned data from the GPU (cudaMemcpy2D)" << std::endl;

  timer->reset();
  timer->start();
  checkCudaErrors(cudaMemcpy2D(h_A, COLS * sizeof(float), d_A, pitch * sizeof(float),
                               COLS * sizeof(float), ROWS, cudaMemcpyDeviceToHost));

  timer->stop();
  std::cout << "Page-locked memory: " << timer->getTime() << " ms ("
      << 1e-6 * BYTES / timer->getTime() << " GB/s)" << std::endl;

  timer->reset();
  timer->start();
  checkCudaErrors(cudaMemcpy2D(h_B, COLS * sizeof(float), d_B, pitch * sizeof(float),
                               COLS * sizeof(float), ROWS, cudaMemcpyDeviceToHost));

  timer->stop();
  std::cout << "Pageable memory:    " << timer->getTime() << " ms ("
      << 1e-6 * BYTES / timer->getTime() << " GB/s)" << std::endl;

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;

  checkCudaErrors(cudaFree(d_A));
  checkCudaErrors(cudaFree(d_B));

  // Page-locked memory has to be deallocated using CUDA API.
  checkCudaErrors(cudaFreeHost(h_A));
  delete[] h_B;

  // Required for cuda-memcheck to detect leaks correctly.
  checkCudaErrors(cudaDeviceReset());
  return 0;
}
