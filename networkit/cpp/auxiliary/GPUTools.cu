/*
 * GPUTools.cu
 *
 *  Created on: 21.03.2022
 *      Author: Fabian Brandt-Tumescheit, fabratu
 *              Lucas Petersen
 */

#include <networkit/auxiliary/GPUTools.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace Aux {

namespace GPUTools {

bool initAndTestCUDA(int dev) {
    // Test result is only true if CUDA toolkit is found on the host and a device context could be established.
    bool result = false;

#ifdef __CUDACC__
    auto checkInitError = [](CUresult error, std::string msg) {
        if (error != CUDA_SUCCESS) {
            printf("%s: %d\n", msg.c_str(), error);
            return false;
        }
        return true;
    };

    CUdevice cuDevice;
    CUcontext cuContext;

    //initialize CUDA
    cuInit(0);
    
    if(checkInitError(cuDeviceGet(&cuDevice, dev), "cannot get device " + std::to_string(dev)) && checkInitError(cuCtxCreate(&cuContext, 0, cuDevice), "cannot create context")) {
        result = true;
    }
#endif
    (void)dev;
    return result;
}

void synchronizeCUDA() {
#ifdef __CUDACC__
    cudaDeviceSynchronize();
#else
    WARN("NetworKit core was built without GPU-support. Therefore synchronizeCUDA() provides no functionality.");
#endif
}

void resetCUDA() {
#ifdef __CUDACC__
    cudaDeviceReset();
#else
    WARN("NetworKit core was built without GPU-support. Therefore resetCUDA() provides no functionality.");
#endif
}

std::vector<std::string> getCUDADevices() {
    std::vector<std::string> devList(0);
    
#ifdef __CUDACC__  
    int devCount;
    cudaGetDeviceCount(&devCount);
    cudaDeviceProp props;
    for(int i = 0; i < devCount; i++) {
        cudaGetDeviceProperties(&props, i);
        devList.push_back(props.name);
    }
#endif
    return devList;
}

#ifdef __CUDACC__ 
void checkCUDAError(cudaError_t err) {

    if (cudaSuccess != err) {
        throw std::runtime_error("CUDA Error = " + std::to_string(err) + ": " + cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}
#endif

int getDeviceMaxCUDA() {
#ifdef __CUDACC__ 
	int devCount;
	checkCUDAError(cudaGetDeviceCount(&devCount));
	cudaDeviceProp prop;

	int maxcc=0, bestdev=0;
	for(int i=0; i<devCount; i++){
		checkCUDAError(cudaGetDeviceProperties(&prop,i));
		if((prop.major + 0.1*prop.minor) > maxcc){
				maxcc = prop.major + 0.1*prop.minor;
				bestdev = i;
		}	
	}
	checkCUDAError(cudaSetDevice(bestdev));
	checkCUDAError(cudaGetDeviceProperties(&prop,bestdev));
	return bestdev;
#else
    WARN("NetworKit core was built without GPU-support. Therefore getDeviceMaxCUDA() always returns an invalid device index.");
    return -1;
#endif
}



} /* namespace GPUTools */

} // namespace Aux
