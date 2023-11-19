#ifndef VVERY_SIMPLE_LIEGROUP_HEADER
#define VVERY_SIMPLE_LIEGROUP_HEADER

#include "statistics.hpp"

namespace ppx
{
    class SO3;
    class SE3;
    using T3 = MatrixS<3, 1>;

    class so3 : public MatrixS<3, 1>
    {
    public:
        so3() = default;

        template <typename T, details::enable_arith_type_t<T> * = nullptr>
        so3(std::initializer_list<T> list) : MatrixS<3, 1>(list) {}

        so3(const MatrixS<3, 1> &mat) : MatrixS<3, 1>(mat) {}

        SO3 exp() const;

        MatrixS<3, 3> adt() const;

        MatrixS<3, 3> ljac() const;

        MatrixS<3, 3> ljacinv() const;

        MatrixS<3, 3> rjac() const;

        MatrixS<3, 3> rjacinv() const;
    };

    inline MatrixS<3, 3> hat(const so3 &vec)
    {
        return {0.0, vec[2], -vec[1], -vec[2], 0.0, vec[0], vec[1], -vec[0], 0.0};
    }

    inline so3 vee(const MatrixS<3, 3> &mat)
    {
        return {mat[5], mat[6], mat[1]};
    }

    class se3 : public MatrixS<6, 1>
    {
    public:
        se3() = default;

        template <typename T, details::enable_arith_type_t<T> * = nullptr>
        se3(std::initializer_list<T> list) : MatrixS<6, 1>(list) {}

        se3(const MatrixS<6, 1> &mat) : MatrixS(mat) {}

        se3(const MatrixS<3, 1> &elem1, const MatrixS<3, 1> &elem2)
        {
            m_data[0] = elem1[0];
            m_data[1] = elem1[1];
            m_data[2] = elem1[2];

            m_data[3] = elem2[0];
            m_data[4] = elem2[1];
            m_data[5] = elem2[2];
        }

        SE3 exp() const;

        MatrixS<6, 6> adt() const;

        MatrixS<3, 1> _1() const
        {
            return {(*this)[0], (*this)[1], (*this)[2]};
        }

        MatrixS<3, 1> _2() const
        {
            return {(*this)[3], (*this)[4], (*this)[5]};
        }
    };

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
        SO3() : MatrixS(Rep::eye()){};
        SO3(const MatrixS<3, 1> &xAxis, const MatrixS<3, 1> &yAxis, const MatrixS<3, 1> &zAxis)
            : MatrixS<3, 3>{xAxis[0], xAxis[1], xAxis[2], yAxis[0], yAxis[1], yAxis[2], zAxis[0], zAxis[1], zAxis[2]} {}
        SO3(double m[3][3]) : MatrixS({m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2]}) {}

        SO3 operator*(const SO3 &other) const
        {
            return MatrixS::operator*(other);
        }
        template <size_t L>
        MatrixS<3, L> operator*(const MatrixS<3, L> &other) const
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

    inline SO3 so3::exp() const
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

    inline MatrixS<3, 3> so3::adt() const
    {
        return hat(*this);
    }

    inline MatrixS<3, 3> so3::ljac() const
    {
        auto result = MatrixS<3, 3>::eye();
        auto x2 = m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2];
        auto X = hat(*this);
        if (x2 < EPS_DP)
        {
            return result + 0.5 * X;
        }
        auto x = sqrt(x2);
        return result + (1 - cos(x)) / x2 * X + (x - sin(x)) / (x2 * x) * X * X;
    }

    inline MatrixS<3, 3> so3::ljacinv() const
    {
        auto result = MatrixS<3, 3>::eye();
        auto x2 = m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2];
        auto X = hat(*this);
        if (x2 < EPS_DP)
        {
            return result - 0.5 * X;
        }
        auto x = sqrt(x2);
        return result - 0.5 * X + (1 / x2 - (1 + cos(x)) / (2 * x * sin(x))) * X * X;
    }

    inline MatrixS<3, 3> so3::rjac() const
    {
        return ljac().T();
    }

    inline MatrixS<3, 3> so3::rjacinv() const
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
        SE3(const SO3 &Rot, const T3 &Pos) : MatrixS(Rep::eye())
        {
            (*this).sub<3, 3, false>(0, 0) = Rot;
            (*this).sub<3, 1, false>(0, 3) = Pos;
        }
        SE3() : MatrixS(Rep::eye()) {}
        SE3(double m[4][4]) : MatrixS({m[0][0], m[1][0], m[2][0], 0.0, m[0][1], m[1][1], m[2][1], 0.0, m[0][2], m[1][2], m[2][2], 0.0, m[0][3], m[1][3], m[2][3], 1.0}) {}

        SE3 operator*(const SE3 &other) const
        {
            return MatrixS::operator*(other);
        }
        template <size_t L>
        MatrixS<4, L> operator*(const MatrixS<4, L> &other) const
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
                MatrixS<3, 3> logExpand = SO3::eye() - omgmat / 2.0 + coff1 * (omgmat * omgmat);
                return {omg, logExpand * p};
            }
            else
            {
                return {T3{}, p};
            }
        }
    };

    inline SE3 se3::exp() const
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
            return SE3::eye() + ksi + coff1 * (ksi * ksi) + coff2 * (ksi * ksi * ksi);
        }
    }

    inline MatrixS<6, 6> se3::adt() const
    {
        MatrixS<6, 6> result;
        result.sub<3, 3, false>(0, 0) = hat(_1());
        result.sub<3, 3, false>(3, 0) = hat(_2());
        result.sub<3, 3, false>(3, 3) = hat(_1());
        return result;
    }

}
#endif