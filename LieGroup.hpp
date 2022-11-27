#ifndef VVERY_SIMPLE_SO3_HEADER
#define VVERY_SIMPLE_SO3_HEADER

#include "matrix.hpp"

class SO3;

class so3 : public Matrix<3, 1>
{
    using Rep = Matrix<3, 1>;

public:
    using Matrix::Matrix;
    so3() = default;
    so3(Rep &&other) : Matrix(std::forward<Rep>(other))
    {
    }

    SO3 exp() const;
};

class SO3 : public Matrix<3, 3>
{
    using Rep = Matrix<3, 3>;

public:
    using Matrix::Matrix;
    SO3(Rep &&other) : Matrix(std::forward<Rep>(other)) {}
    SO3() : Matrix(Rep::eye()){};

    Matrix<3, 3> Adt() const
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
    double theta = norm2(s);
    if (details::near_zero(theta))
    {
        return {};
    }
    else
    {
        auto omg = hat(s);
        auto coff1 = sin(theta) / theta;
        auto coff2 = (1 - cos(theta)) / (theta * theta);
        return SO3::eye() + coff1 * omg + coff2 * (omg * omg);
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

class SE3;

class se3 : public Matrix<6, 1>
{
    using T3 = Matrix<3, 1>;

public:
    using Matrix::Matrix;
    se3(const T3 &w, const T3 &v)
    {
        m_data[0] = w[0];
        m_data[1] = w[1];
        m_data[2] = w[2];
        m_data[3] = v[0];
        m_data[4] = v[1];
        m_data[5] = v[2];
    }
    se3() = default;

    T3 w() const
    {
        return {m_data[0], m_data[1], m_data[2]};
    }
    T3 v() const
    {
        return {m_data[3], m_data[4], m_data[5]};
    }
    SE3 exp() const;
};

class SE3 : public Matrix<4, 4>
{
    using Rep = Matrix<4, 4>;
    using T3 = Matrix<3, 1>;

public:
    using Matrix::Matrix;
    SE3(Rep &&other) : Matrix(std::forward<Rep>(other)) {}
    SE3(const SO3 &Rot, const T3 &Pos) : Matrix(Rep::eye())
    {
        (*this)({0, 2}, {0, 2}) = Rot;
        (*this)({0, 2}, 3) = Pos;
    }
    SE3() : Matrix(Rep::eye()) {}

    SO3 Rot() const
    {
        return slice<3, 3>(*this, 0, 0);
    }
    T3 Pos() const
    {
        return slice<3, 1>(*this, 0, 3);
    }
    Matrix<6, 6> Adt() const
    {
        Matrix<6, 6> result{};
        const auto &R = Rot();
        const auto &p = Pos();
        result({0, 2}, {0, 2}) = R;
        result({3, 5}, {0, 2}) = hat(p) * R;
        result({3, 5}, {3, 5}) = R;
        return result;
    }
    SE3 I() const
    {
        return {Rot().T(), -1 * (Rot().T() * Pos())};
    }
    se3 log() const;
};

inline SE3 hat(const Matrix<6, 1> &vec)
{
    return {0.0, vec[2], -vec[1], 0.0, -vec[2], 0.0, vec[0], 0.0, vec[1], -vec[0], 0.0, 0.0, vec[3], vec[4], vec[5], 0.0};
}

inline se3 vee(const Matrix<4, 4> &mat)
{
    return {mat[6], mat[8], mat[1], mat[12], mat[13], mat[14]};
}

SE3 se3::exp() const
{
    auto phi = norm2(w());
    if (details::near_zero(phi))
    {
        return {SO3::zero(), v()};
    }
    else
    {
        const auto ksi = hat(*this);
        auto coff1 = (1 - cos(phi)) / (phi * phi);
        auto coff2 = (phi - sin(phi)) / (phi * phi * phi);
        return SE3::eye() + ksi + coff1 * (ksi * ksi) + coff2 * (ksi * ksi * ksi);
    }
}

se3 SE3::log() const
{
    const auto &R = Rot();
    const auto &p = Pos();
    if (R == SO3::zero())
    {
        return {T3{}, p};
    }
    else
    {
        auto theta = acos((R.trace() - 1) / 2.0);
        auto omg = R.log();
        auto omgmat = hat(omg);
        auto coff1 = (1.0 / theta - 1.0 / (tan(theta / 2.0) * 2.0)) / theta;
        Matrix<3, 3> logExpand = SO3::eye() - omgmat / 2.0 + coff1 * (omgmat * omgmat);
        return {omg, logExpand * p};
    }
}

#endif