// File: cuda_main.cu
// C/Fortran interface to CUDA C API.

// includes standard headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// includes cuda headers
#include <cuda_runtime.h>
// includes project headers
#include "cuda_globals.h"

#ifdef __PARA
//#undef SEEK_SET  // remove compilation errors
//#undef SEEK_CUR  // with C++ binding of MPI
//#undef SEEK_END
#include <mpi.h>
#endif



// global variables
int NUM_STREAMS=0;              // number of CUDA streams
cudaStream_t *stream;		// CUDA stream
double *d_reduce, *d_reduce1;	// arrays for parallel reduction
cuDoubleComplex *d_zreduce, *d_zreduce1;  // arrays for parallel reduction
devptr_t *d_ptrs, *d_ptrs1;	// arrays of device pointers
int nPE_, myPE_;


/******************************************************/
// CUDA C wrappers for init, used in VAMP

extern "C"
void cuda_init_(int *nstreams, int *nsim)
{
    int i;

    /* Get MPI Information */
#ifdef __PARA
    MPI_Comm_size(MPI_COMM_WORLD, &nPE_);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE_);
#else
    nPE_ = 1; myPE_ = 0;
#endif


    // check number of CUDA streams requested
    if(*nstreams<=0)
    {
        printf("Nstreams is %d\n", *nstreams);
        ERROR( "GPU Library", "Invalid number of CUDA streams:pick a number greater than zero!");
    }
    NUM_STREAMS=*nstreams;
    printf("creating %d CUDA streams...\n",NUM_STREAMS);

    // create CUDA streams
    stream=(cudaStream_t*)malloc(NUM_STREAMS*sizeof(cudaStream_t));
    for(i=0;i<NUM_STREAMS;i++)
	CUDA_ERROR( cudaStreamCreate(&stream[i]), "Failed to create CUDA stream!" );

    // allocate parallel reduction arrays
    CUDA_ERROR( cudaMalloc((void **)&d_zreduce,MAX_THREADS*sizeof(cuDoubleComplex)),
		"Failed to allocate device memory!" );
    CUDA_ERROR( cudaMalloc((void **)&d_zreduce1,MAX_THREADS*sizeof(cuDoubleComplex)),
                "Failed to allocate device memory!" );
    // set parallel reduction arrays
    d_reduce = (double *)d_zreduce;
    d_reduce1 = (double *)d_zreduce1;

    // allocate device pointer arryas
    CUDA_ERROR( cudaMalloc((void **)&d_ptrs,(*nsim)*sizeof(devptr_t)),
		"Failed to allocate device memory!" );
    CUDA_ERROR( cudaMalloc((void **)&d_ptrs1,(*nsim)*sizeof(devptr_t)),
                "Failed to allocate device memory!" );
}

extern "C"
void cuda_destroy_(void)
{
    int i;
    // destroy CUDA streams
    for(i=0;i<NUM_STREAMS;i++)
        CUDA_ERROR( cudaStreamDestroy(stream[i]), "Failed to destroy CUDA stream!" );
    free(stream);

    // free parallel reduction arrays
    CUDA_ERROR( cudaFree(d_reduce), "Failed to allocate device memory!" );
    CUDA_ERROR( cudaFree(d_reduce1), "Failed to allocate device memory!" );

    // free device pointer arrays
    CUDA_ERROR( cudaFree(d_ptrs), "Failed to allocate device memory!" );
    CUDA_ERROR( cudaFree(d_ptrs1), "Failed to allocate device memory!" );
}

// TODO: replace with init from exact exchange?
extern "C"
void cuda_mpi_init_(int *CudaDevice)
{
#ifndef EMULATION
    int deviceCount, gpu_rank;
    cudaDeviceProp deviceProp;

    /* Get MPI Information */
#ifdef __PARA
    MPI_Comm_size(MPI_COMM_WORLD, &nPE_);
    MPI_Comm_rank(MPI_COMM_WORLD, &myPE_);
#else
    nPE_ = 1; myPE_ = 0;
#endif



    CUDA_ERROR( cudaGetDeviceCount(&deviceCount), "No CUDA-supporting devices found!" );

    gpu_rank = (*CudaDevice)*deviceCount/nPE_;
    CUDA_ERROR( cudaGetDeviceProperties(&deviceProp, gpu_rank),
		"Device does not support CUDA!" );
    if(deviceProp.major < 1)
    {
        printf( "CUDA ERROR: Devices does not support CUDA!\n");
        cudaDeviceReset();
	exit(1);
    }
    printf("Using device %d (rank %d) : %s\n", gpu_rank,*CudaDevice,deviceProp.name);
    CUDA_ERROR(cudaSetDevice(gpu_rank), "Failed to set the device!" );
#endif
}

/******************************************************/
// CUDA C wrappers for thread sync, in VASP

extern "C"
void cuda_device_reset_(void)
{
    printf("Reseting the CUDA device...\n");
    CUDA_ERROR( cudaDeviceReset(), "Failed to reset the device!" );
}

// synchronze the device
extern "C"
void cuda_devicesynchronize_(char *msg)
{
    CUDA_ERROR( cudaDeviceSynchronize(), msg );
}

// in fortran source
extern "C"
void threadsynchronize_(void)
{
    CUDA_ERROR( cudaThreadSynchronize(), "Failed to synchronize the device!" );
}

extern "C"
void cuda_streamsynchronize_(int *sid)
{
    cudaStream_t st = CUDA_STREAM(*sid);  // CUDA stream
    CUDA_ERROR( cudaStreamSynchronize(st), "Failed to synchronize the CUDA stream!" );
}

extern "C"
void cuda_all_stream_synchronize_(void)
{
    CUDA_ERROR( cudaStreamSynchronize(0), "Failed to synchronize all CUDA streams!" );
}

/******************************************************/
