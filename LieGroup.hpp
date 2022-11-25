#ifndef VVERY_SIMPLE_SO3_HEADER
#define VVERY_SIMPLE_SO3_HEADER

#include "matrix.hpp"

class SO3;

class so3 : public Matrix<3, 1>
{

public:
    using Matrix::Matrix;
    so3() = default;
    ~so3() = default;
    so3(const so3 &) = default;
    so3(so3 &&) = default;

    SO3 exp() const;
};

class SO3 : public Matrix<3, 3>
{
    using Rep = Matrix<3, 3>;

public:
    using Matrix::Matrix;
    SO3(Rep &&other) : Matrix(std::forward<Rep>(other)) {}
    SO3() = default;
    ~SO3() = default;
    SO3(const SO3 &) = default;
    SO3(SO3 &&) = default;

    SO3 Adt() const
    {
        return *this;
    }
    SO3 I() const
    {
        return this->T();
    }
    so3 log() const;
};

inline SO3 hat(const Matrix<3, 1> &vec)
{
    return {0.0, vec[2], -vec[1], -vec[2], 0.0, vec[0], vec[1], -vec[0], 0.0};
}

inline so3 vee(const Matrix<3, 3> &mat)
{
    return {mat[5], mat[6], mat[1]};
}

SO3 so3::exp() const
{
    const auto &s = *this;
    auto ret = SO3::eye();
    double theta = norm2(s);
    if (details::near_zero(theta))
    {
        return ret;
    }
    else
    {
        auto omg = hat(s / theta);
        return ret + sin(theta) * omg + (1 - cos(theta)) * (omg * omg);
    }
}

so3 SO3::log() const
{
    Matrix<3, 3> ret{};
    const auto &R = *this;
    auto acosinput = (this->trace() - 1) / 2.0;
    if (acosinput >= 1.0)
    {
        return {};
    }
    else if (acosinput <= -1.0)
    {
        Matrix<3, 1> omg{};
        if (!details::near_zero(1 + R(2, 2)))
        {
            omg = (1.0 / sqrt(2 + 2 * R(2, 2))) * Matrix<3, 1>{R(0, 2), R(1, 2), R(2, 2)};
        }
        else if (!details::near_zero(1 + R(1, 1)))
        {
            omg = (1.0 / sqrt(2 + 2 * R(1, 1))) * Matrix<3, 1>{R(0, 1), R(1, 1), +R(2, 1)};
        }
        else
        {
            omg = (1.0 / sqrt(2 + 2 * R(0, 0))) * Matrix<3, 1>{R(0, 0), R(1, 0), +R(2, 0)};
        }
        return gl_rep_pi * hat(omg);
    }
    else
    {
        double theta = acos(acosinput);
        ret = (R - R.T()) * theta / (2.0 * sin(theta));
        return vee(ret);
    }
}

#endif