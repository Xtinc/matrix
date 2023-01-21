# PPX quick start

A quick start for ppx users.

## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
- [License](#license)

## Background

ppx is a cluster of algorithms for numerical calculation. It consists of 3 parts: linear algebra, non-linear optimization, statistics.

This library contains:

1. Linear algebra operations. ppx provides a STL-like container to describe matrix. And related operations such as cofactor,determinant...(Be careful that lazy evaluation is always used in element-wised matrix operators)
2. Non-linear optimization algorithms. For 1D optimization, some line search methods are given.
For N-Dimensions optimization problems, Direct method and gradient method are given.  
3. Statics toolkits. to be developed.

## Install

This project requires CPP standard 14, and it's header-only.

## Usage

```cpp
#include "liegroup.hpp"
// declarations of  matrix
Matrix<4, 4> a = {1, 2, 3, 4, 5, 6, 7, 8};
PRINT_SINGLE_ELEMENTS(a);
PRINT_SINGLE_ELEMENTS(a(maxloc(a)), "max of a: ");
// slice of a matrix
a({2, 3}, {2, 3}) = Matrix<2, 2>{1, 1, 1, 1};
a(0, {0, -1}) = {89.0, 23.0, 44.0, 9.8};
PRINT_SINGLE_ELEMENTS(a, "a = ");
PRINT_LISTED_ELEMENTS(a, "a in list: ");
// functions for matrix operations
PRINT_SINGLE_ELEMENTS(determinant(a), "determinant(a) = ");
PRINT_SINGLE_ELEMENTS(inverse(a), "inverse(a) = ");
// expression templates for matrix operations
Matrix<4, 4> x = {1, 2, 3, 4, 5, 6, 7, 8};
Matrix<4, 4> y = x.T();
Matrix<4, 4> z = x * 2 + y * 2 + 3.3 + x * y;
PRINT_SINGLE_ELEMENTS(z, "2*(x + x^T) + 3.3 +x*x^T = ");
//  lie group operations
Matrix<3, 3> so3mat = {0, 3, -2, -3, 0, 1, 2, -1, 0};
PRINT_SINGLE_ELEMENTS(vee(so3mat).exp(), "exp(s) = ");
SO3 SO3mat = {0, 1, 0, 0, 0, 1, 1, 0, 0};
PRINT_SINGLE_ELEMENTS(SO3mat.log(), "log(S) = ");
PRINT_SINGLE_ELEMENTS(SE3(SO3mat, {1, 1, 1}), "T = ");
SE3 SE3mat = {1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 3, 1};
PRINT_SINGLE_ELEMENTS(SE3mat, "SE3mat = ");
PRINT_SINGLE_ELEMENTS(SE3mat.log(), "log(s) = ");
PRINT_SINGLE_ELEMENTS(hat(SE3mat.log()));
PRINT_SINGLE_ELEMENTS(se3{1.5708, 0.0, 0.0, 0.0, 2.3562, 2.3562}.exp());
// linear solver for equations
Matrix<4, 3> X{3.5, 1.6, 3.7, 4.3,
                   2.7, -5.7, -0.8, -9.8,
                   -3.1, -6.0, 1.9, 6.9};
Matrix<4, 1> Y{1, 1, 1, 1};
PRINT_SINGLE_ELEMENTS(linsolve<factorization::QR>(X, Y), "solved by QR = ");
PRINT_SINGLE_ELEMENTS(linsolve<factorization::SVD>(X, Y), "solved by SVD = ");
// factorization for matrix
Matrix<4, 3> u{1, 4, 7, 11,
                2, 5, 8, 1,
`               3, 6, 9, 5};
Matrix<3, 1> v{};
Matrix<3, 3> w{};
bool sing = false;
u = svdcmp(u, v, w, sing);
PRINT_SINGLE_ELEMENTS(u, "U = ");
PRINT_SINGLE_ELEMENTS(v, "S = ");
PRINT_SINGLE_ELEMENTS(w, "V = ");
PRINT_SINGLE_ELEMENTS(u * v.diag() * w.T(), "U*S*V = ");
// eigenvalue for symmetry matrix
Matrix<4, 4> g{2, -1, 1, 4,
                   -1, 2, -1, 1,
                   1, -1, 2, -2,
                   4, 1, -2, 3};
PRINT_SINGLE_ELEMENTS(g, "g = ");
PRINT_SINGLE_ELEMENTS(eig<eigensystem::SymValAndVecSorted>(g).vec, "eigen vector of g : ");
PRINT_SINGLE_ELEMENTS(eig<eigensystem::SymOnlyVal>(g), "eigen value of g : ");

```

## License

[MIT](LICENSE) © XTinc
