// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2012 Desire Nuentsa Wakam <desire.nuentsa_wakam@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
#include "sparse.h"
#include <Eigen/SparseQR>


template<typename MatrixType,typename DenseMat>
int generate_sparse_rectangular_problem(MatrixType& A, DenseMat& dA, int maxRows = 300, int maxCols = 300)
{
  eigen_assert(maxRows >= maxCols);
  typedef typename MatrixType::Scalar Scalar;
  int rows = internal::random<int>(1,maxRows);
  int cols = internal::random<int>(1,rows);
  double density = (std::max)(8./(rows*cols), 0.01);
  
  A.resize(rows,rows);
  dA.resize(rows,rows);
  initSparse<Scalar>(density, dA, A,ForceNonZeroDiag);
  A.makeCompressed();
  return rows;
}

template<typename Scalar> void test_sparseqr_scalar()
{
  typedef SparseMatrix<Scalar,ColMajor> MatrixType; 
  MatrixType A;
  Matrix<Scalar,Dynamic,Dynamic> dA;
  typedef Matrix<Scalar,Dynamic,1> DenseVector;
  DenseVector refX,x,b; 
  SparseQR<MatrixType, AMDOrdering<int> > solver; 
  generate_sparse_rectangular_problem(A,dA);
  
  int n = A.cols();
  b = DenseVector::Random(n);
  solver.compute(A);
  if (solver.info() != Success)
  {
    std::cerr << "sparse QR factorization failed\n";
    exit(0);
    return;
  }
  x = solver.solve(b);
  if (solver.info() != Success)
  {
    std::cerr << "sparse QR factorization failed\n";
    exit(0);
    return;
  }  
  //Compare with a dense QR solver
  refX = dA.colPivHouseholderQr().solve(b);
  VERIFY(x.isApprox(refX,test_precision<Scalar>()));
}
void test_sparseqr()
{
  CALL_SUBTEST_1(test_sparseqr_scalar<double>());
  CALL_SUBTEST_2(test_sparseqr_scalar<std::complex<double> >());
}