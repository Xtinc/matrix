#ifndef VVERY_SIMPLE_LIEGROUP_HEADER
#define VVERY_SIMPLE_LIEGROUP_HEADER

#include "statistics.hpp"

namespace ppx
{
    using so3 = MatrixS<3, 1>;
    using se3 = MatrixS<6, 1>;
    using T3 = MatrixS<3, 1>;

    inline MatrixS<3, 3> hat(const so3 &vec)
    {
        return {0.0, vec[2], -vec[1], -vec[2], 0.0, vec[0], vec[1], -vec[0], 0.0};
    }

    inline so3 vee(const MatrixS<3, 3> &mat)
    {
        return {mat[5], mat[6], mat[1]};
    }

    inline MatrixS<4, 4> hat(const se3 &vec)
    {
        return {0.0, vec[2], -vec[1], 0.0, -vec[2], 0.0, vec[0], 0.0, vec[1], -vec[0], 0.0, 0.0, vec[3], vec[4], vec[5], 0.0};
    }

    inline se3 vee(const MatrixS<4, 4> &mat)
    {
        return {mat[6], mat[8], mat[1], mat[12], mat[13], mat[14]};
    }

    class SO3 : public MatrixS<3, 3>
    {
        using Rep = MatrixS<3, 3>;

    public:
        template <typename T>
        SO3(std::initializer_list<T> list) : MatrixS(list) {}
        SO3(const Rep &other) : MatrixS(other) {}
        SO3() : MatrixS(eye<3>()){};
        SO3(const MatrixS<3, 1> &xAxis, const MatrixS<3, 1> &yAxis, const MatrixS<3, 1> &zAxis)
            : MatrixS<3, 3>{xAxis[0], xAxis[1], xAxis[2], yAxis[0], yAxis[1], yAxis[2], zAxis[0], zAxis[1], zAxis[2]} {}
        SO3(double m[3][3]) : MatrixS({m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2]}) {}

        SO3 operator*(const SO3 &other) const
        {
            return MatrixS::operator*(other);
        }

        MatrixS<3, 3> Adt() const
        {
            return *this;
        }
        SO3 I() const
        {
            return this->T();
        }
        so3 log() const
        {
            const auto &R = *this;
            auto acosinput = (this->trace() - 1) / 2.0;
            if (acosinput >= 1.0)
            {
                return {};
            }
            else if (acosinput <= -1.0)
            {
                MatrixS<3, 1> omg{};
                if (!details::is_same(1 + R(2, 2), 0.0))
                {
                    omg = (1.0 / sqrt(2 + 2 * R(2, 2))) * T3{R(0, 2), R(1, 2), 1 + R(2, 2)};
                }
                else if (!details::is_same(1 + R(1, 1), 0.0))
                {
                    omg = (1.0 / sqrt(2 + 2 * R(1, 1))) * T3{R(0, 1), 1 + R(1, 1), +R(2, 1)};
                }
                else
                {
                    omg = (1.0 / sqrt(2 + 2 * R(0, 0))) * T3{1 + R(0, 0), R(1, 0), +R(2, 0)};
                }
                return PI * omg;
            }
            else
            {
                MatrixS<3, 3> ret{};
                double theta = acos(acosinput);
                ret = (R - R.T()) * theta / (2.0 * sin(theta));
                return vee(ret);
            }
        }

        static SO3 RotX(double theta)
        {
            return {1.0, 0.0, 0.0,
                    0.0, cos(theta), -sin(theta),
                    0.0, sin(theta), cos(theta)};
        }

        static SO3 RotY(double theta)
        {
            return {cos(theta), 0.0, -sin(theta),
                    0.0, 1.0, 0.0,
                    sin(theta), 0.0, cos(theta)};
        }

        static SO3 RotZ(double theta)
        {
            return {cos(theta), sin(theta), 0.0,
                    -sin(theta), cos(theta), 0.0,
                    0.0, 0.0, 1.0};
        }
    };

    template <size_t M, size_t N>
    template <size_t A, size_t B, std::enable_if_t<A == 3 && B == 1> *>
    SO3 MatrixS<M, N>::exp() const
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
            return eye<3>() + coff1 * omg + coff2 * (omg * omg);
        }
    }

    template <size_t M, size_t N>
    template <size_t A, size_t B, std::enable_if_t<A == 3 && B == 1> *>
    MatrixS<3, 3> MatrixS<M, N>::adt() const
    {
        return hat(*this);
    }

    template <size_t M, size_t N>
    template <size_t A, size_t B, std::enable_if_t<A == 3 && B == 1> *>
    MatrixS<3, 3> MatrixS<M, N>::ljac() const
    {
        auto result = eye<3>();
        auto x2 = this->m_data[0] * this->m_data[0] + this->m_data[1] * this->m_data[1] + this->m_data[2] * this->m_data[2];
        auto X = hat(*this);
        if (x2 < EPS_DP)
        {
            return result + 0.5 * X;
        }
        auto x = sqrt(x2);
        return result + (1 - cos(x)) / x2 * X + (x - sin(x)) / (x2 * x) * X * X;
    }

    template <size_t M, size_t N>
    template <size_t A, size_t B, std::enable_if_t<A == 3 && B == 1> *>
    MatrixS<3, 3> MatrixS<M, N>::ljacinv() const
    {
        auto result = eye<3>();
        auto x2 = this->m_data[0] * this->m_data[0] + this->m_data[1] * this->m_data[1] + this->m_data[2] * this->m_data[2];
        auto X = hat(*this);
        if (x2 < EPS_DP)
        {
            return result - 0.5 * X;
        }
        auto x = sqrt(x2);
        return result - 0.5 * X + (1 / x2 - (1 + cos(x)) / (2 * x * sin(x))) * X * X;
    }

    template <size_t M, size_t N>
    template <size_t A, size_t B, std::enable_if_t<A == 3 && B == 1> *>
    MatrixS<3, 3> MatrixS<M, N>::rjac() const
    {
        return ljac().T();
    }

    template <size_t M, size_t N>
    template <size_t A, size_t B, std::enable_if_t<A == 3 && B == 1> *>
    MatrixS<3, 3> MatrixS<M, N>::rjacinv() const
    {
        return ljacinv().T();
    }

    class SE3 : public MatrixS<4, 4>
    {
        using Rep = MatrixS<4, 4>;

    public:
        template <typename T>
        SE3(std::initializer_list<T> list) : MatrixS(list) {}
        SE3(const Rep &other) : MatrixS(other) {}
        SE3(const SO3 &Rot, const T3 &Pos) : MatrixS(eye<4>())
        {
            (*this).sub<3, 3, false>(0, 0) = Rot;
            (*this).sub<3, 1, false>(0, 3) = Pos;
        }
        SE3() : MatrixS(eye<4>()) {}
        SE3(double m[4][4]) : MatrixS({m[0][0], m[1][0], m[2][0], 0.0, m[0][1], m[1][1], m[2][1], 0.0, m[0][2], m[1][2], m[2][2], 0.0, m[0][3], m[1][3], m[2][3], 1.0}) {}

        SE3 operator*(const SE3 &other) const
        {
            return MatrixS::operator*(other);
        }

        SO3 Rot() const
        {
            return (*this).sub<3, 3>(0, 0);
        }
        T3 Pos() const
        {
            return (*this).sub<3, 1>(0, 3);
        }
        MatrixS<6, 6> Adt() const
        {
            MatrixS<6, 6> result{};
            const auto &R = Rot();
            const auto &p = Pos();
            result.sub<3, 3, false>(0, 0) = R;
            result.sub<3, 3, false>(3, 0) = hat(p) * R;
            result.sub<3, 3, false>(3, 3) = R;
            return result;
        }
        SE3 I() const
        {
            return {Rot().T(), -1 * (Rot().T() * Pos())};
        }
        se3 log() const
        {
            const auto &R = Rot();
            const auto &p = Pos();
            auto ctheta = (R.trace() - 1) / 2.0;
            assert(fabs(ctheta) < 1 + EPS_SP);
            if (fabs(ctheta) < 1)
            {
                auto theta = acos(ctheta);
                auto omg = R.log();
                auto omgmat = hat(omg);
                auto coff1 = (1.0 / theta - 1.0 / (tan(theta / 2.0) * 2.0)) / theta;
                MatrixS<3, 3> logExpand = eye<3>() - omgmat / 2.0 + coff1 * (omgmat * omgmat);
                return {omg, logExpand * p};
            }
            else
            {
                return {T3{}, p};
            }
        }
    };

    template <size_t M, size_t N>
    template <size_t A, size_t B, std::enable_if_t<A == 6 && B == 1> *>
    SE3 MatrixS<M, N>::exp() const
    {
        auto phi = norm2(_1());
        if (details::is_same(phi, 0))
        {
            return {SO3(), _2()};
        }
        else
        {
            const auto ksi = hat(*this);
            auto coff1 = (1 - cos(phi)) / (phi * phi);
            auto coff2 = (phi - sin(phi)) / (phi * phi * phi);
            return eye<4>() + ksi + coff1 * (ksi * ksi) + coff2 * (ksi * ksi * ksi);
        }
    }

    template <size_t M, size_t N>
    template <size_t A, size_t B, std::enable_if_t<A == 6 && B == 1> *>
    MatrixS<6, 6> MatrixS<M, N>::adt() const
    {
        MatrixS<6, 6> result;
        result.sub<3, 3, false>(0, 0) = hat(_1());
        result.sub<3, 3, false>(3, 0) = hat(_2());
        result.sub<3, 3, false>(3, 3) = hat(_1());
        return result;
    }

}
#endif