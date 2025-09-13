/*
 * GPUTools.hpp
 *
 *  Created on: 21.03.2022
 *      Author: Fabian Brandt-Tumescheit, fabratu
 */

#ifndef NETWORKIT_AUXILIARY_GPU_TOOLS_HPP_
#define NETWORKIT_AUXILIARY_GPU_TOOLS_HPP_

#include <iostream>
#include <string>
#include <vector>

#ifdef __CUDACC__
#include <cuda.h>
#endif // __CUDACC__

namespace Aux {

/**
 * Tools to deal with GPU-computing
 */
namespace GPUTools {

/**
 * @brief Tries to initialize CUDA-capable GPU-device. Can also be used to to runtime testing in
 * order to create branches for CPU/GPU-computing.
 * @param dev Device id, which should be used for computation. Default: 0
 */
bool initAndTestCUDA(int dev = 0);

/**
 * @brief Finishes current GPU-based computation and transfer resulting data to host.
 *
 */
void synchronizeCUDA();

/**
 * @brief Resets a current CUDA-context and device-session.
 *
 */
void resetCUDA();

/**
 * @brief Get a vector containing all CUDA-device names. The index of the device can be used for
 * initializing a connection.
 *
 */
std::vector<std::string> getCUDADevices();

/**
 * @brief Get a vector containing all CUDA-device names. The index of the device can be used for
 * initializing a connection.
 *
 */
int getDeviceMaxCUDA();

#ifdef __CUDACC__
/**
 * @brief Helper function to check for errors of CUDA runtime functions (like cudaMalloc, cudaMemcpy
 * and cudaFree).
 * @param err Return value of CUDA runtime functions.
 *
 */
void checkCUDAError(cudaError_t err);
#endif // __CUDACC__

} /* namespace GPUTools */

} // namespace Aux
#endif // NETWORKIT_AUXILIARY_GPU_TOOLS_HPP_
