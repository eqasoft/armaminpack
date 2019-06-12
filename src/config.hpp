#ifndef EQA_CONFIG_HPP
#define EQA_CONFIG_HPP


// -------------------------- STL------------------------------
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdint>
#include <cstddef>      // size_t
#include <iomanip>
#include <limits>
#include <type_traits>  // for is_integral
#include <random>
#include <complex>
#include <algorithm>
#include <iterator>
#include <string>
#include <math.h>



using namespace std;

#ifdef USE_ARMA
#include "armadillo"
using namespace arma;
template<typename T>
using not_arma_mat = std::enable_if<!(std::is_same<T, arma::mat>::value)>;
#endif

#ifdef USE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif




#if defined(_OPENMP) && !defined(EQA_DONT_USE_OPENMP) && !defined(EQA_USE_OPENMP)
#define EQA_USE_OPENMP
#include <omp.h>
#endif

#if defined(_MSC_VER)
#ifdef  EQA_EXPORTS
#define EQA_API __declspec(dllexport)
#else
#define EQA_API __declspec(dllimport)
#endif
#else
#define EQA_API
#endif



#ifndef __stats_pointer_settings__
#if defined(__clang__) || defined(__GNUC__)
#define __stats_pointer_settings__ __restrict__
#elif defined(_MSC_VER)
#define __stats_pointer_settings__ __restrict
#else
#define __stats_pointer_settings__
#endif
#endif

#if defined(_WIN32) && !defined(WIN32)
#define WIN32
#endif


#endif // !EQA_CONFIG_HPP
