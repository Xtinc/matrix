#include "matrix.hpp"
#include <random>

int main(int, char **)
{
    Matrix<4, 4> a = {1, 2, 3, 4, 5, 6, 7, 8};
    PRINT_SINGLE_ELEMENTS(a);
    a({2, 3}, {2, 3}) = Matrix<2, 2>{1, 1, 1, 1};
    a(0, {0, -1}) = {89.0, 23.0, 44.0, 9.8};
    PRINT_SINGLE_ELEMENTS(a, "a = ");
    PRINT_SINGLE_ELEMENTS(determinant(a), "determinant(a) = ");
    PRINT_SINGLE_ELEMENTS(inverse(a), "inverse(a) = ");
    Matrix<4, 4> A(std::move(a));
    PRINT_LISTED_ELEMENTS(a.flat(), "Move Sementics, a = ");
    PRINT_SINGLE_ELEMENTS(A, "Move Sementics, A = ");
    Matrix<2, 2> b = {1, 2, 4, 3};
    PRINT_SINGLE_ELEMENTS(inverse(b), "inverse(b) = ");
    Matrix<1, 1> c = {3};
    PRINT_SINGLE_ELEMENTS(inverse(c), "inverse(c) = ");
    Matrix<10, 10> M{};
    std::default_random_engine e;
    std::uniform_real_distribution<> u(-5, 5);
    for (auto &i : M.flat())
    {
        i = u(e);
    }
    PRINT_SINGLE_ELEMENTS(M, "M = ");
    PRINT_SINGLE_ELEMENTS(determinant(M), "determinant(M) = ");
    PRINT_SINGLE_ELEMENTS(adjoint(M),"adjoint(M) = ");
}
