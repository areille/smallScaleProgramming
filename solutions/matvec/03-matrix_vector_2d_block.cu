// 
// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
//

// Computes matrix-vector product. Matrix A is in row-major order
// i.e. A[i, j] is stored in i * ncols + j element of the vector.
//

#include <iostream>

#include <cuda_runtime.h>  // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h>  // For CUDA SDK timers


// What is a good initial guess for XBD and YBD (both
// greater than 1) ?
// After you get the code to work, experiment with different sizes
// to find the best possible performance
// Note: For meaningful time measurements you need sufficiently large matrices.
#define XBD 64
#define YBD 16
const dim3 BLOCK_DIM(XBD,YBD);

// Simple CPU implementation of matrix addition.
void CpuMatrixVector(int rows, int cols, const float* A, const float* x, float* y) {
  for (int row = 0; row < rows; ++row) {
    float t=0.0;
    for (int col = 0; col < cols; ++col) {
      int idx = row * cols + col;
      t += A[idx] * x[col];
    }
    y[row] = t;
  }
}

// GPU implementation of matrix_vector product
// using a block of threads for each block of rows. 
__global__ void gpuMatrixVector(int rows, int cols, const float* A,
				const float* x, float* y) {
  __shared__ float ax[YBD][XBD];
  int tr     = threadIdx.y;
  int tc     = threadIdx.x;
  int row    = blockIdx.x*blockDim.y + tr;
  ax[tr][tc] = 0.0;
  if (row < rows) {
    // Starting address of indexing within matrix A
    int idxm = row*cols + tc;
    int ic   = tc; 
    float t  = 0.0;
    for ( ; ic<cols; ic += XBD) {
      t += A[idxm]*x[ic];
      idxm += XBD;
    }
    if (ic < cols) {
      t += A[idxm]*x[ic];
    }
    ax[tr][tc] = t;
  }
  __syncthreads();
  for (int s=XBD/2; s >32; s >>=1){
    if (tc<s)
      ax[tr][tc] += ax[tr][tc+s]; 
    __syncthreads();
  }
  for (int s=min(32,XBD/2); s >0; s >>=1){
    if (tc<s)
      ax[tr][tc] += ax[tr][tc+s]; 
  }

  if ((tc == 0)&&(row<rows))
    y[row] = ax[tr][tc];
  
}

int main(int argc, char** argv) {

  if (argc < 3) {
    fprintf(stderr,"Usage: %s  rows cols\n",argv[0]);
  }
  int nrows=atoi(argv[1]);
  int ncols=atoi(argv[2]);
  
  
  // ----------------------- Host memory initialisation ----------------------- //

  float* h_A = new float[nrows * ncols];
  float* h_x = new float[ncols];
  float* h_y = new float[nrows];
  float* h_y_d = new float[nrows];

  srand(123456);
  for (int row = 0; row < nrows; ++row) {
    for (int col = 0; col < ncols; ++col) {
      int idx = row * ncols + col;
      h_A[idx] = 100.0f * static_cast<float>(rand()) / RAND_MAX;
    }
    h_y[row] = 0.0;
  }
  for (int col = 0; col < ncols; ++col) {
    h_x[col] = 100.0f * static_cast<float>(rand()) / RAND_MAX;
  }

  std::cout << "Matrix-vector product: 2D thread block version " << std::endl;
  std::cout << "Test case: " << nrows  << " x " << ncols << std::endl;
// ---------------------- Device memory initialisation ---------------------- //

  float *d_A, *d_x, *d_y;

  checkCudaErrors(cudaMalloc((void**) &d_A, nrows * ncols * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**) &d_x, ncols * sizeof(float)));
  checkCudaErrors(cudaMalloc((void**) &d_y, nrows * sizeof(float)));

  // Copy matrices from the host (CPU) to the device (GPU).
  checkCudaErrors(cudaMemcpy(d_A, h_A, nrows * ncols * sizeof(float), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_x, h_x,  ncols * sizeof(float), cudaMemcpyHostToDevice));

  // ------------------------ Calculations on the CPU ------------------------- //
  float flopcnt=2.e-6*nrows*ncols;
  
  // Create the CUDA SDK timer.
  StopWatchInterface* timer = 0;
  sdkCreateTimer(&timer);

  timer->start();
  CpuMatrixVector(nrows, ncols, h_A, h_x, h_y);

  timer->stop();
  float cpuflops=flopcnt/ timer->getTime();
  std::cout << "  CPU time: " << timer->getTime() << " ms." << " GFLOPS " << cpuflops << std::endl;
  
// ------------------------ Calculations on the GPU ------------------------- //

  // Calculate the dimension of the grid of blocks. A 1D grid suffices.
  const dim3 GRID_DIM((nrows - 1+ BLOCK_DIM.y)/ BLOCK_DIM.y  ,1);

  timer->reset();
  timer->start();
  gpuMatrixVector<<<GRID_DIM, BLOCK_DIM >>>(nrows, ncols, d_A, d_x, d_y);
  checkCudaErrors(cudaDeviceSynchronize());

  timer->stop();
  float gpuflops=flopcnt/ timer->getTime();
  std::cout << "  GPU time: " << timer->getTime() << " ms." << " GFLOPS " << gpuflops<<std::endl;

  // TODO Download the resulting vector d_y from the device and store it in h_y_d.
  checkCudaErrors(cudaMemcpy(h_y_d, d_y, nrows*sizeof(float),cudaMemcpyDeviceToHost));


  // Now let's check if the results are the same.
  float reldiff = 0.0f;
  float diff = 0.0f;
  
  for (int row = 0; row < nrows; ++row) {
    float maxabs = std::max(std::abs(h_y[row]),std::abs(h_y_d[row]));
    if (maxabs == 0.0) maxabs=1.0;
    reldiff = std::max(reldiff, std::abs(h_y[row] - h_y_d[row])/maxabs);
    diff = std::max(diff, std::abs(h_y[row] - h_y_d[row]));
  }
  std::cout << "Max diff = " << diff << "  Max rel diff = " << reldiff << std::endl;
  // Rel diff should be as close as possible to unit roundoff; float
  // corresponds to IEEE single precision, so unit roundoff is
  // 1.19e-07
  // 

// ------------------------------- Cleaning up ------------------------------ //

  delete timer;

  checkCudaErrors(cudaFree(d_A));
  checkCudaErrors(cudaFree(d_x));
  checkCudaErrors(cudaFree(d_y));

  delete[] h_A;
  delete[] h_x;
  delete[] h_y;
  delete[] h_y_d;
  return 0;
}
