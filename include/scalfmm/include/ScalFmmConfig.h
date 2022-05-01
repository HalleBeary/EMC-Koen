// See LICENCE file at project root
///
/// File generated by cmake

// Do not remove any line
///////////////////////////////////////////////////////
#ifndef SSCALFMMCONFIG_H
#define SSCALFMMCONFIG_H
///////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////
#define SCALFMM_API_2

///////////////////////////////////////////////////////
// Debug
///////////////////////////////////////////////////////

// Uncomment the next line to use debug mode
/* #undef SCALFMM_USE_LOG */

///////////////////////////////////////////////////////
// Blas
///////////////////////////////////////////////////////

#define SCALFMM_USE_BLAS
/* #undef SCALFMM_USE_MKL_AS_BLAS */
// Fortran Mangling
#define SCALFMM_BLAS_ADD_
/* #undef SCALFMM_BLAS_UPCASE */
/* #undef SCALFMM_BLAS_NOCHANGE */
////////////////////////////////////////////////////////
// FFT
///////////////////////////////////////////////////////

/* #undef SCALFMM_USE_FFT */
/* #undef SCALFMM_USE_MKL_AS_FFTW */
/* #undef SCALFMM_USE_ESSL_AS_FFTW */

//////////////////////////////////////////////////////
// MPI
///////////////////////////////////////////////////////

/* #undef SCALFMM_USE_MPI */

///////////////////////////////////////////////////////
// Memory trace
///////////////////////////////////////////////////////

/* #undef SCALFMM_USE_MEM_STATS */

///////////////////////////////////////////////////////
// CUDA
///////////////////////////////////////////////////////

/* #undef SCALFMM_USE_CUDA */

///////////////////////////////////////////////////////
// OPENCL
///////////////////////////////////////////////////////

/* #undef SCALFMM_USE_OPENCL */

///////////////////////////////////////////////////////
// STARPU
///////////////////////////////////////////////////////

/* #undef SCALFMM_USE_STARPU */
/* #undef SCALFMM_DISABLE_NATIVE_OMP4 */

///////////////////////////////////////////////////////
// EZTRACE
///////////////////////////////////////////////////////

/* #undef SCALFMM_USE_EZTRACE */
/* #undef SCALFMM_TRACE_ALGO */
/* #undef SCALFMM_TRACE_P2P */
/* #undef SCALFMM_TRACE_P2M */
/* #undef SCALFMM_TRACE_M2L */
/* #undef SCALFMM_TRACE_L2L */
/* #undef SCALFMM_TRACE_L2P */


///////////////////////////////////////////////////////
// Assert tests
///////////////////////////////////////////////////////

#define SCALFMM_USE_ASSERT


#ifdef __INTEL_COMPILER
#pragma warning (disable : 858 )
#pragma warning (disable : 2326 )
#endif

///////////////////////////////////////////////////////
// Path to the SCALFMM DATA (For UTests)
///////////////////////////////////////////////////////
#include <string>
const std::string SCALFMMDataPath("/home/p274598/EnsembleMC/ScalFMM/Data/");


///////////////////////////////////////////////////////
// Flags and libs used to compile
///////////////////////////////////////////////////////
const std::string SCALFMMCompileFlags("-fpic -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wconversion -Wcast-align -Woverloaded-virtual");
const std::string SCALFMMCompileLibs("OpenMP BLAS LAPACK");

///////////////////////////////////////////////////////
// To use commute for KSTAR OMP4
///////////////////////////////////////////////////////

/* #undef OPENMP_SUPPORT_COMMUTE */

///////////////////////////////////////////////////////
// To use priority for KSTAR OMP4
///////////////////////////////////////////////////////

/* #undef OPENMP_SUPPORT_PRIORITY */

///////////////////////////////////////////////////////
// To use a taskname clause for tasks with KSTAR OMP4
///////////////////////////////////////////////////////

/* #undef OPENMP_SUPPORT_TASK_NAME */

///////////////////////////////////////////////////////
// To record omp4 task times for statistics
///////////////////////////////////////////////////////

/* #undef SCALFMM_TIME_OMPTASKS */

///////////////////////////////////////////////////////
// To catch signals and print backtrace
///////////////////////////////////////////////////////

/* #undef SCALFMM_USE_SIGNALS */

///////////////////////////////////////////////////////
// To control starpu config
///////////////////////////////////////////////////////

/* #undef SCALFMM_STARPU_USE_COMMUTE */
/* #undef SCALFMM_STARPU_USE_REDUX */
/* #undef SCALFMM_STARPU_USE_PRIO */
/* #undef SCALFMM_STARPU_FORCE_NO_SCHEDULER */
/* #undef SCALFMM_USE_STARPU_EXTRACT */

///////////////////////////////////////////////////////
// To control simgrid config
///////////////////////////////////////////////////////

/* #undef SCALFMM_SIMGRID_NODATA */
/* #undef STARPU_SIMGRID_MLR_MODELS */

#endif // CONFIG_H