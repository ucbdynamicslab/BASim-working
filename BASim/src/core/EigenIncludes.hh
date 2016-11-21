/**
 * \file EigenIncludes.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/15/2009
 */

#ifndef EIGENINCLUDES_HH
#define EIGENINCLUDES_HH

//#define EIGEN_DONT_ALIGN

#define EIGEN_VECTOR_IO Eigen::IOFormat(8, Eigen::Raw, ",", ",", "", "", "{", "}")
#define EIGEN_MATRIX_IO Eigen::IOFormat(8, Eigen::Raw, ",", ",", "{", "}", "{", "}")
#define EIGEN_SPACES_ONLY_IO Eigen::IOFormat(8, Eigen::Raw, " ", " ", "", "", "", "")
#define EIGEN_DEFAULT_IO_FORMAT EIGEN_MATRIX_IO

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO 1

#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#endif // EIGENINCLUDES_HH
