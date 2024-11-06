!
! Simple Fortan90 program that multiplies a matrix and a vector
! using CUBLAS or CUDA
!
! This example is using cudaMallocHost to allocate the matrices in pinned memory.
! Pinned memory enables fast PCI-e transfer.
!
! The code is using the iso_c_binding Fortran 2003 extension.
! It has been tested with the Intel Fortran compiler v10.1 and g95.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interface to cudaMallocHost and cudaFree
module cuda_alloc
 use iso_c_binding
 interface
! These are disabled and replaced by wrappers in fortran.c to have more control
! over the type of allocation (pinned,not pinned)


! cudaMallocHost
 integer (C_INT) function MycudaMallocHost(buffer, size)  bind(C,name="cudaMallocHost")
  use iso_c_binding
  implicit none
  type (C_PTR)  :: buffer
  integer (C_SIZE_T), value :: size
 end function MycudaMallocHost
! cudaFreeHost
 integer (C_INT) function cudaFreeHost(buffer)  bind(C,name="cudaFreeHost")
  use iso_c_binding
  implicit none
  type (C_PTR), value :: buffer
 end function cudaFreeHost
 end interface
end module cuda_alloc
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

