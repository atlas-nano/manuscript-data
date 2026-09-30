#include <mpi.h>
#include <stdio.h>
#if defined(OPEN_MPI) && OPEN_MPI
#include <mpi-ext.h>
#endif
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
#if defined(MPIX_CUDA_AWARE_SUPPORT)
  printf("compile-time MPIX_CUDA_AWARE_SUPPORT=%d runtime MPIX_Query_cuda_support()=%d\n", MPIX_CUDA_AWARE_SUPPORT, MPIX_Query_cuda_support());
#else
  printf("MPIX_CUDA_AWARE_SUPPORT not defined\n");
#endif
  MPI_Finalize(); return 0;
}
