#include "LieGroup.hpp"
#include <random>

using namespace details;

int main(int, char **)
{
    Matrix<4, 4> a = {1, 2, 3, 4, 5, 6, 7, 8};
    PRINT_SINGLE_ELEMENTS(a);
    a({2, 3}, {2, 3}) = Matrix<2, 2>{1, 1, 1, 1};
    a(0, {0, -1}) = {89.0, 23.0, 44.0, 9.8};
    PRINT_SINGLE_ELEMENTS(a, "a = ");
    PRINT_LISTED_ELEMENTS(a, "a in list: ");
    PRINT_SINGLE_ELEMENTS(determinant(a), "determinant(a) = ");
    PRINT_SINGLE_ELEMENTS(inverse(a), "inverse(a) = ");
    Matrix<4, 4> A(std::move(a));
    //  PRINT_LISTED_ELEMENTS(a.flat(), "Move Sementics, a = ");
    PRINT_SINGLE_ELEMENTS(A, "Move Sementics, A = ");
    Matrix<2, 2> b = {1, 2, 4, 3};
    PRINT_SINGLE_ELEMENTS(inverse(b), "inverse(b) = ");
    Matrix<1, 1> c = {3};
    PRINT_SINGLE_ELEMENTS(c.I(), "inverse(c) = ");
    Matrix<30, 30> M{};
    std::default_random_engine eni;
    std::uniform_real_distribution<> uf(-5, 5);
    for (auto &i : M)
    {
        i = uf(eni);
    }
    PRINT_SINGLE_ELEMENTS(M, "M = ");
    PRINT_SINGLE_ELEMENTS(determinant(M), "determinant(M) = ");
    PRINT_SINGLE_ELEMENTS(adjugate(M), "adjoint(M) = ");
    Matrix<4, 4> x = {1, 2, 3, 4, 5, 6, 7, 8};
    Matrix<4, 4> y = x.T();
    Matrix<4, 4> z = x * 2 + y * 2 + 3.3 + x * y;
    PRINT_SINGLE_ELEMENTS(z, "2*(x + x^T) + 3.3 +x*x^T = ");
    Matrix<3, 3> so3mat = {0, 3, -2, -3, 0, 1, 2, -1, 0};
    PRINT_SINGLE_ELEMENTS(vee(so3mat).exp(), "exp(s) = ");
    SO3 SO3mat = {0, 1, 0, 0, 0, 1, 1, 0, 0};
    PRINT_SINGLE_ELEMENTS(SO3mat.log(), "log(S) = ");
}
