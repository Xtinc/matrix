#ifndef VVERY_SIMPLE_MATRIXD_HEADER
#define VVERY_SIMPLE_MATRIXD_HEADER

#include "exprtmpl.hpp"
#include <iostream>
#include <iomanip>
#include <memory>
#include <array>
#include <vector>
#include <numeric>
#include <cassert>

namespace ppx
{

    class MatrixD
    {
        using value_type = double;
        using cast_type = MatrixD;
        using elem_tag = details::ElemTags::Matrix;
        using iterator = typename std::vector<double>::iterator;
        using const_iterator = typename std::vector<double>::const_iterator;

        iterator begin() noexcept
        {
            return this->m_data.begin();
        }
        iterator end() noexcept
        {
            return this->m_data.end();
        }
        const_iterator begin() const noexcept
        {
            return this->m_data.begin();
        }
        const_iterator end() const noexcept
        {
            return this->m_data.end();
        }
        const_iterator cbegin() const noexcept
        {
            return this->m_data.begin();
        }
        const_iterator cend() const noexcept
        {
            return this->m_data.end();
        }

    public:
        MatrixD() : m(0), n(0) {}

        MatrixD(size_t M, size_t N) : m(M), n(N), m_data(m * n, 0.0) {}

        template <typename T, size_t L, details::enable_arith_type_t<T> * = nullptr>
        MatrixD(size_t M, size_t N, const std::array<T, L> &list) : m(M), n(N), m_data(m * n, 0.0)
        {
            std::copy_n(list.begin(), std::min(L, m * n), m_data.begin());
        }

        template <typename T, details::enable_arith_type_t<T> * = nullptr>
        MatrixD(size_t M, size_t N, const std::initializer_list<T> &list) : m(M), n(N), m_data(m * n, 0.0)
        {
            std::copy_n(list.begin(), std::min(list.size(), m * n), m_data.begin());
        }

        template <typename T, details::enable_arith_type_t<T> * = nullptr>
        MatrixD(size_t M, size_t N, const std::vector<T> &list) : m(M), n(N), m_data(m * n, 0.0)
        {
            std::copy_n(list.begin(), std::min(list.size(), m * n), m_data.begin());
        }

        // Member functions
        double *data()
        {
            return this->m_data.data();
        }

        const double *data() const
        {
            return this->m_data.data();
        }

        MatrixD row(size_t idx) const
        {
            assert(idx < m);
            MatrixD result(n, 1);
            for (size_t i = 0; i < n; ++i)
            {
                result[i] = m_data.at(idx + i * m);
            }
            return result;
        }

        MatrixD col(size_t idx) const
        {
            assert(idx < n);
            MatrixD result(m, 1);
            for (size_t i = 0; i < m; ++i)
            {
                result[i] = m_data.at(idx * m + i);
            }
            return result;
        }

        void fill(double val)
        {
            std::fill(m_data.begin(), m_data.end(), val);
        }

        MatrixD T() const
        {
            MatrixD res(n, m);
            for (size_t i = 0; i < m; i++)
            {
                for (size_t j = 0; j < n; j++)
                {
                    res(j, i) = (*this)(i, j);
                }
            }
            return res;
        }

        MatrixD I() const
        {
            if (m == n)
            {
                std::vector<int> indx(m, 0);
                auto even = true;
                auto LU = ludcmp(*this, indx, even);
                if (LU.s == StatusCode::SINGULAR)
                {
                    return {};
                }
                auto result = MatrixD::eye(m, m);
                for (size_t j = 0; j < m; j++)
                {
                    ludbksb(LU.x, indx, result.data() + j * m);
                }
                return result;
            }
            else
            {
                return {};
            }
        }

        double det() const
        {
            assert(m == n);
            auto even = true;
            std::vector<int> indx(m, 0.0);
            auto LU = ludcmp(*this, indx, even);
            if (LU.s == StatusCode::SINGULAR)
            {
                return {};
            }
            auto D = even ? 1.0 : -1.0;
            for (size_t i = 0; i < m; i++)
            {
                D *= LU.x(i, i);
            }
            return D;
        }

        double trace() const
        {
            assert(m == n);
            double res = 0.0;
            for (size_t i = 0; i < m; i++)
            {
                res += (*this)(i, i);
            }
            return res;
        }

        // Overloaded Operators
        double &operator()(size_t row, size_t col)
        {
            assert(row < m && col < n);
            return m_data.at(row + col * m);
        }

        const double &operator()(size_t row, size_t col) const
        {
            assert(row < m && col < n);
            return m_data.at(row + col * m);
        }

        double &operator()(const std::pair<size_t, size_t> &idx)
        {
            assert(idx.first < m && idx.second < n);
            return m_data.at(idx.first + idx.second * m);
        }

        const double &operator()(const std::pair<size_t, size_t> &idx) const
        {
            assert(idx.first < m && idx.second < n);
            return m_data.at(idx.first + idx.second * m);
        }

        double &operator[](size_t idx)
        {
            assert(idx < m * n);
            return m_data.at(idx);
        }

        const double &operator[](size_t idx) const
        {
            assert(idx < m * n);
            return m_data.at(idx);
        }

        // Generate function.
        friend std::ostream &operator<<(std::ostream &os, const MatrixD &self)
        {
            auto m = self.m;
            auto n = self.n;
            os << "Matrix<" << m << "," << n << ">:\n";
            for (auto i = 0u; i < m; i++)
            {
                for (auto j = 0u; j < n; j++)
                {
                    os << std::setw(12) << self(i, j) << "\t";
                }
                os << std::endl;
            }
            return os;
        }

        size_t rows() const
        {
            return m;
        }

        size_t cols() const
        {
            return n;
        }

        size_t size()
        {
            return m * n;
        }

        static MatrixD eye(size_t M, size_t N)
        {
            MatrixD result(M, N);
            auto real_idx = M < N ? M : N;
            for (size_t i = 0; i < real_idx; i++)
            {
                result.data()[i * (M + 1)] = 1.0;
            }
            return result;
        }

    private:
        size_t m;
        size_t n;
        std::vector<double> m_data;
    };
} // namespace ppx

#endif