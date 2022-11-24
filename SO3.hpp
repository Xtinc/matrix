#ifndef VVERY_SIMPLE_SO3_HEADER
#define VVERY_SIMPLE_SO3_HEADER

#include "matrix.hpp"
#include <cassert>

// forward declare
class so3;

class SO3 : public Matrix<3, 3>
{
    using Rep = Matrix<3, 3>;
    using Tangent = so3;

public:
    using Matrix::Matrix;
    template <typename T>
    SO3(T &&other)
        : Matrix(std::forward<T>(other))
    {
    }
    // Member functions
    SO3 I() const
    {
        return this->T();
    }
    Tangent log() const
    {
        const auto &a = *this;
        auto phi = acos((a[0] + a[4] + a[8] - 1.0) / 2.0);
        
    }
    bool valid() const
    {
        auto AAT = (*this) * this->T();
        return (AAT == Rep::eye() && is_same(this->det(), 1.0));
    }
};

class so3 : public Matrix<3, 1>
{
    using Rep1 = Matrix<3, 1>;
    using Rep2 = Matrix<3, 3>;

public:
    using Matrix::Matrix;

    // member function
    Rep2 hat() const
    {
        const auto &a = *this;
        return {0.0, a[2], -a[1], -a[2], 0.0, a[0], a[1], -a[0], 0.0};
    }

    SO3 exp() const
    {
        auto a = *this;
        auto s = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        a = a / s;
        return cos(s) * Rep2::eye() + (1 - cos(s)) * (a * a.T()) + sin(s) * a.hat();
    }
};

#endif