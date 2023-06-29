#ifndef VVERY_SIMPLE_MATRIXS_HEADER
#define VVERY_SIMPLE_MATRIXS_HEADER

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
    // forward declare
    template <size_t M, size_t N>
    class MatrixS;

    class SO3;
    class SE3;

    enum class StatusCode : char
    {
        NORMAL,
        CONVERGED,
        DIVERGED,
        SINGULAR
    };

    template <size_t M, size_t N>
    struct FacResult
    {
        MatrixS<M, N> x;
        StatusCode s = StatusCode::NORMAL;
    };

    inline std::ostream &operator<<(std::ostream &os, const StatusCode &self)
    {
        switch (self)
        {
        case StatusCode::NORMAL:
            os << "NORMAL";
            break;
        case StatusCode::CONVERGED:
            os << "CONVERGED";
            break;
        case StatusCode::DIVERGED:
            os << "DIVERGED";
            break;
        case StatusCode::SINGULAR:
            os << "SINGULAR";
        default:
            break;
        }
        return os;
    }

    template <size_t N>
    FacResult<N, N> ludcmp(MatrixS<N, N> A, MatrixS<N, 1> &indx, bool &even);

    template <size_t M, size_t N>
    FacResult<M, N> svdcmp(MatrixS<M, N> u, MatrixS<N, 1> &w, MatrixS<N, N> &v);

    template <size_t N>
    void ludbksb(const MatrixS<N, N> &A, const MatrixS<N, 1> &indx, double *b);

    template <size_t M, size_t N>
    void svbksb(const MatrixS<M, N> &u, const MatrixS<N, 1> &w, const MatrixS<N, N> &v, double *b);

    namespace details
    {
        constexpr size_t MAX_SIZE_LIMIT = 260;

        inline bool is_same(double a, double b)
        {
            return fabs(a - b) < EPS_SP;
        }

        inline bool near_zero(double a)
        {
            return fabs(a) < 1.0e-5;
        }

        constexpr size_t is_small_matrix_v(size_t A, size_t B)
        {
            return A * B < MAX_SIZE_LIMIT ? 1 : 0;
        }

        template <size_t M, size_t N, std::size_t A = is_small_matrix_v(M, N)>
        class MatrixBase;

        template <std::size_t M, std::size_t N>
        class MatrixBase<M, N, 1>
        {
        protected:
            using iterator = typename std::array<double, M * N>::iterator;
            using const_iterator = typename std::array<double, M * N>::const_iterator;

            MatrixBase() : m_data{} {}

            template <typename T, size_t L, details::enable_arith_type_t<T> * = nullptr>
            MatrixBase(const std::array<T, L> &list) : m_data{}
            {
                constexpr auto real_idx = std::min(L, M * N);
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

            template <typename T, details::enable_arith_type_t<T> * = nullptr>
            MatrixBase(const std::initializer_list<T> &list) : m_data{}
            {
                auto real_idx = list.size() < M * N ? list.size() : M * N;
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

            template <typename T, details::enable_arith_type_t<T> * = nullptr>
            MatrixBase(const std::vector<T> &list) : m_data{}
            {
                auto real_idx = list.size() < M * N ? list.size() : M * N;
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

        protected:
            std::array<double, M * N> m_data;
        };

        template <std::size_t M, std::size_t N>
        class MatrixBase<M, N, 0>
        {
        protected:
            using iterator = typename std::vector<double>::iterator;
            using const_iterator = typename std::vector<double>::const_iterator;

            MatrixBase() : m_data(M * N, 0.0) {}

            template <typename T, size_t L, details::enable_arith_type_t<T> * = nullptr>
            MatrixBase(const std::array<T, L> &list) : m_data(M * N, 0.0)
            {
                constexpr auto real_idx = std::min(L, M * N);
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

            MatrixBase(const std::initializer_list<int> &list) : m_data(M * N, 0.0)
            {
                auto real_idx = list.size() < M * N ? list.size() : M * N;
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

            template <typename T, details::enable_arith_type_t<T> * = nullptr>
            MatrixBase(const std::vector<T> &list) : m_data(M * N, 0.0)
            {
                auto real_idx = list.size() < M * N ? list.size() : M * N;
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

        protected:
            std::vector<double> m_data;
        };

    } // namespace details

    template <std::size_t M, std::size_t N>
    class MatrixS : public details::MatrixBase<M, N>
    {
        template <size_t A, typename RT = void>
        using enable_when_squre_t = std::enable_if_t<A == N, RT>;
        template <size_t A, typename RT = void>
        using disable_when_squre_t = std::enable_if_t<A != N, RT>;

        template <size_t A, size_t B>
        class SubMatrixBase
        {

        public:
            using elem_tag = details::ElemTags::Mblock;
            using cast_type = MatrixS<A, B>;

            SubMatrixBase(MatrixS<M, N> &self, size_t r, size_t c)
                : row_idx(r), col_idx(c), data(self)
            {
                assert(row_idx + A <= M && col_idx + B <= N);
            }

            const cast_type &snap() const
            {
                return *copy;
            }

        protected:
            size_t row_idx;
            size_t col_idx;
            MatrixS<M, N> &data;
            std::unique_ptr<MatrixS<A, B>> copy;

            void take_snap()
            {
                copy = std::make_unique<MatrixS<A, B>>();
                for (size_t i = 0; i < A; i++)
                {
                    for (size_t j = 0; j < B; j++)
                    {
                        (*copy)(i, j) = data(row_idx + i, col_idx + j);
                    }
                }
            }

            void copy_data(const SubMatrixBase &other)
            {
                for (size_t i = 0; i < A; i++)
                {
                    for (size_t j = 0; j < B; j++)
                    {
                        data(row_idx + i, col_idx + j) =
                            other.data(other.row_idx + i, other.col_idx + j);
                    }
                }
            }

            template <typename T>
            void ctor_by_list(const T &list)
            {
                auto iter = list.begin();
                for (auto j = 0u; j < B; j++)
                {
                    for (auto i = 0u; i < A; i++)
                    {
                        auto value = iter == list.end() ? typename T::value_type{} : *(iter++);
                        data(i + row_idx, j + col_idx) = value;
                        // which one is better?
                        // if (iter != list.end())
                        // {
                        //     data(i + row_idx, j + col_idx) = *(iter++);
                        // }
                    }
                }
            }
        };

        template <size_t A, size_t B, bool RefAndVal = true>
        class SubMatrix : public SubMatrixBase<A, B>
        {
            using base_type = SubMatrixBase<A, B>;

        public:
            using elem_tag = typename base_type::elem_tag;
            using cast_type = typename base_type::cast_type;

            SubMatrix(MatrixS<M, N> &self, size_t r, size_t c)
                : SubMatrixBase<A, B>(self, r, c)
            {
                SubMatrixBase<A, B>::take_snap();
            }

            SubMatrix(const SubMatrix &other)
            {
                base_type::copy_data(other);
                *base_type::copy = *other.copy;
            }

            SubMatrix(SubMatrix &&other) noexcept : base_type::copy(std::move(other.copy))
            {
                base_type::copy_data(other);
            }

            SubMatrix &operator=(const SubMatrix &other)
            {
                base_type::copy_data(other);
                *base_type::copy = *other.copy;
                return *this;
            }

            SubMatrix &operator=(SubMatrix &&other) noexcept
            {
                base_type::copy_data(other);
                base_type::copy = std::move(other.copy);
                return *this;
            }

            SubMatrix &operator=(const MatrixS<A, B> &other)
            {
                base_type::ctor_by_list(other);
                return *this;
            }

            template <typename T, size_t L, details::enable_arith_type_t<T> * = nullptr>
            SubMatrix &operator=(const std::array<T, L> &list)
            {
                base_type::ctor_by_list(list);
                return *this;
            }

            template <typename T, details::enable_arith_type_t<T> * = nullptr>
            SubMatrix &operator=(const std::initializer_list<T> &list)
            {
                base_type::ctor_by_list(list);
                return *this;
            }

            template <typename T, details::enable_arith_type_t<T> * = nullptr>
            SubMatrix &operator=(const std::vector<T> &list)
            {
                base_type::ctor_by_list(list);
                return *this;
            }

            template <typename T>
            details::enable_expr_type_t<T, SubMatrix &>
            operator=(const T &expr)
            {
                for (size_t j = 0; j < B; j++)
                {
                    for (size_t i = 0; i < A; i++)
                    {
                        base_type::data(base_type::row_idx + i, base_type::col_idx + j) = expr[i + j * A];
                    }
                }
                return *this;
            }

            operator cast_type()
            {
                MatrixS<A, B> result;
                for (size_t i = 0; i < A; i++)
                {
                    for (size_t j = 0; j < B; j++)
                    {
                        result(i, j) = base_type::data(base_type::row_idx + i, base_type::col_idx + j);
                    }
                }
                return result;
            }
        };

        template <size_t A, size_t B>
        class SubMatrix<A, B, false> : public SubMatrixBase<A, B>
        {
            using base_type = SubMatrixBase<A, B>;

        public:
            using elem_tag = typename base_type::elem_tag;
            using cast_type = typename base_type::cast_type;

            SubMatrix(MatrixS<M, N> &self, size_t r, size_t c)
                : SubMatrixBase<A, B>(self, r, c)
            {
            }

            SubMatrix(const SubMatrix &other)
            {
                base_type::copy_data(other);
            }

            SubMatrix(SubMatrix &&other) noexcept : base_type::copy(std::move(other.copy))
            {
                base_type::copy_data(other);
            }

            SubMatrix &operator=(const SubMatrix &other)
            {
                base_type::copy_data(other);
                return *this;
            }

            SubMatrix &operator=(SubMatrix &&other) noexcept
            {
                base_type::copy_data(other);
                return *this;
            }

            SubMatrix &operator=(const MatrixS<A, B> &other)
            {
                base_type::ctor_by_list(other);
                return *this;
            }

            template <typename T, size_t L, details::enable_arith_type_t<T> * = nullptr>
            SubMatrix &operator=(const std::array<T, L> &list)
            {
                base_type::ctor_by_list(list);
                return *this;
            }

            template <typename T, details::enable_arith_type_t<T> * = nullptr>
            SubMatrix &operator=(const std::initializer_list<T> &list)
            {
                base_type::ctor_by_list(list);
                return *this;
            }

            template <typename T, details::enable_arith_type_t<T> * = nullptr>
            SubMatrix &operator=(const std::vector<T> &list)
            {
                base_type::ctor_by_list(list);
                return *this;
            }

            template <typename T>
            details::enable_expr_type_t<T, SubMatrix &>
            operator=(const T &expr)
            {
                for (size_t j = 0; j < B; j++)
                {
                    for (size_t i = 0; i < A; i++)
                    {
                        base_type::data(base_type::row_idx + i, base_type::col_idx + j) = expr[i + j * A];
                    }
                }
                return *this;
            }

            operator cast_type()
            {
                MatrixS<A, B> result;
                for (size_t i = 0; i < A; i++)
                {
                    for (size_t j = 0; j < B; j++)
                    {
                        result(i, j) = base_type::data(base_type::row_idx + i, base_type::col_idx + j);
                    }
                }
                return result;
            }
        };

    public:
        using value_type = double;
        using cast_type = MatrixS<M, N>;
        using elem_tag = details::ElemTags::Matrix;
        using iterator = typename details::MatrixBase<M, N>::iterator;
        using const_iterator = typename details::MatrixBase<M, N>::const_iterator;

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
        MatrixS() = default;

        template <typename T, size_t L, details::enable_arith_type_t<T> * = nullptr>
        MatrixS(const std::array<T, L> &list) : details::MatrixBase<M, N>(list)
        {
        }

        template <typename T, details::enable_arith_type_t<T> * = nullptr>
        MatrixS(const std::initializer_list<T> &list) : details::MatrixBase<M, N>(list)
        {
        }

        template <typename T, details::enable_arith_type_t<T> * = nullptr>
        MatrixS(const std::vector<T> &list) : details::MatrixBase<M, N>(list)
        {
        }

        template <size_t A = M, size_t B = N, std::enable_if_t<A == 6 && B == 1> * = nullptr>
        MatrixS(const MatrixS<3, 1> &elem1, const MatrixS<3, 1> &elem2)
            : details::MatrixBase<M, N>({elem1[0], elem1[1], elem1[2], elem2[0], elem2[1], elem2[2]})
        {
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

        MatrixS<N, 1> row(size_t idx) const
        {
            assert(idx < M);
            MatrixS<N, 1> result;
            for (auto i = 0u; i < N; ++i)
            {
                result[i] = this->m_data.at(idx + i * M);
            }
            return result;
        }

        MatrixS<M, 1> col(size_t idx) const
        {
            assert(idx < N);
            MatrixS<M, 1> result;
            for (auto i = 0u; i < M; ++i)
            {
                result[i] = this->m_data.at(idx * M + i);
            }
            return result;
        }

        template <size_t A, size_t B, bool RefAndVal = true>
        SubMatrix<A, B, RefAndVal> sub(size_t row_start, size_t col_start)
        {
            return {*this, row_start, col_start};
        }

        template <size_t A, size_t B>
        MatrixS<A, B> sub(size_t row_start, size_t col_start) const
        {
            assert(row_start + A <= M && col_start + B <= N);
            MatrixS<A, B> result;
            for (size_t i = 0; i < A; i++)
            {
                for (size_t j = 0; j < B; j++)
                {
                    result(i, j) = (*this)(row_start + i, col_start + j);
                }
            }
            return result;
        }

        void fill(double val)
        {
            std::fill(this->m_data.begin(), this->m_data.end(), val);
        }

        MatrixS<N, M> T() const
        {
            MatrixS<N, M> res{};
            for (auto i = 0u; i < M; i++)
            {
                for (auto j = 0u; j < N; j++)
                {
                    res(j, i) = (*this)(i, j);
                }
            }
            return res;
        }

        template <size_t A = M>
        enable_when_squre_t<A, MatrixS> I() const
        {
            std::array<int, M> indx{};
            auto even = true;
            auto LU = ludcmp(*this, indx, even);
            if (LU.s == StatusCode::SINGULAR)
            {
                return {};
            }
            auto result = MatrixS<M, M>::eye();
            for (size_t j = 0; j < M; j++)
            {
                ludbksb(LU.x, indx, result.data() + j * M);
            }
            return result;
        }

        template <size_t A = M>
        disable_when_squre_t<A, MatrixS<N, M>> I() const
        {
            MatrixS<N, 1> w;
            MatrixS<N, N> V;
            auto U = svdcmp(*this, w, V);
            if (U.s == StatusCode::SINGULAR)
            {
                return {};
            }
            MatrixS<N, N> W;
            for (size_t i = 0; i < N; i++)
            {
                if (fabs(w[i]) > EPS_SP)
                {
                    W(i, i) = 1.0 / w[i];
                }
            }
            return V * W * U.x.T();
        }

        template <size_t A = M>
        enable_when_squre_t<A, double> det() const
        {
            auto even = true;
            std::array<int, M> indx{};
            auto LU = ludcmp(*this, indx, even);
            if (LU.s == StatusCode::SINGULAR)
            {
                return {};
            }
            auto D = even ? 1.0 : -1.0;
            for (size_t i = 0; i < M; i++)
            {
                D *= LU.x(i, i);
            }
            return D;
        }

        template <size_t A = M>
        enable_when_squre_t<A, double> trace() const
        {
            double res = 0.0;
            for (size_t i = 0; i < A; i++)
            {
                res += (*this)(i, i);
            }
            return res;
        }

        template <size_t A = std::min(M, N)>
        enable_when_matrix_t<M, N, MatrixS<A, 1>> diag() const
        {
            MatrixS<A, 1> result;
            for (size_t i = 0; i < A; i++)
            {
                result[i] = (*this)(i, i);
            }
            return result;
        }

        template <size_t A = std::max(M, N)>
        enable_when_array_t<M, N, MatrixS<A, A>> diag() const
        {
            MatrixS<A, A> result;
            for (size_t i = 0; i < A; i++)
            {
                result(i, i) = (*this)[i];
            }
            return result;
        }
        // Lie Group
        template <size_t A = M, size_t B = N>
        std::enable_if_t<A == 3 && B == 1, SO3> exp() const;

        template <size_t A = M, size_t B = N>
        std::enable_if_t<A == 3 && B == 1, MatrixS<3, 3>> adt() const;

        template <size_t A = M, size_t B = N>
        std::enable_if_t<A == 6 && B == 1, SE3> exp() const;

        template <size_t A = M, size_t B = N>
        std::enable_if_t<A == 6 && B == 1, MatrixS<6, 6>> adt() const;

        template <size_t A = M, size_t B = N>
        std::enable_if_t<A == 6 && B == 1, MatrixS<3, 1>> _1() const
        {
            return {(*this)[0], (*this)[1], (*this)[2]};
        }

        template <size_t A = M, size_t B = N>
        std::enable_if_t<A == 6 && B == 1, MatrixS<3, 1>> _2() const
        {
            return {(*this)[3], (*this)[4], (*this)[5]};
        }

        // Overloaded Operators
        double &operator()(size_t row, size_t col)
        {
            assert(row < M && col < N);
            return this->m_data.at(row + col * M);
        }

        const double &operator()(size_t row, size_t col) const
        {
            assert(row < M && col < N);
            return this->m_data.at(row + col * M);
        }

        double &operator()(const std::pair<size_t, size_t> &idx)
        {
            assert(idx.first < M && idx.second < N);
            return this->m_data.at(idx.first + idx.second * M);
        }

        const double &operator()(const std::pair<size_t, size_t> &idx) const
        {
            assert(idx.first < M && idx.second < N);
            return this->m_data.at(idx.first + idx.second * M);
        }

        double &operator[](size_t idx)
        {
            assert(idx < M * N);
            return this->m_data.at(idx);
        }

        const double &operator[](size_t idx) const
        {
            assert(idx < M * N);
            return this->m_data.at(idx);
        }

        template <size_t L>
        MatrixS<M, L> operator*(const MatrixS<N, L> &other) const
        {
            MatrixS<M, L> result;
            for (size_t k = 0; k < N; k++)
            {
                for (size_t j = 0; j < L; j++)
                {
                    for (size_t i = 0; i < M; i++)
                    {
                        result(i, j) += (*this)(i, k) * other(k, j);
                    }
                }
            }
            return result;
        }

        MatrixS &operator+=(const MatrixS &other)
        {
            for (size_t i = 0; i < M * N; i++)
            {
                this->m_data[i] += other.data()[i];
            }
            return *this;
        }

        MatrixS &operator-=(const MatrixS &other)
        {
            for (size_t i = 0; i < M * N; i++)
            {
                this->m_data[i] -= other.data()[i];
            }
            return *this;
        }

        bool operator==(const MatrixS &other) const
        {
            return std::equal(this->m_data.begin(), this->m_data.end(), other.data(),
                              [](double ele1, double ele2)
                              { return details::is_same(ele1, ele2); });
        }

        template <typename T>
        details::enable_arith_type_t<T, MatrixS &> operator+=(T ele)
        {
            for (auto &i : this->m_data)
            {
                i += ele;
            }
            return *this;
        }

        template <typename T>
        details::enable_arith_type_t<T, MatrixS &> operator-=(T ele)
        {
            for (auto &i : this->m_data)
            {
                i -= ele;
            }
            return *this;
        }

        template <typename T>
        details::enable_arith_type_t<T, MatrixS &> operator*=(T ele)
        {
            for (auto &i : this->m_data)
            {
                i *= ele;
            }
            return *this;
        }

        template <typename T>
        details::enable_arith_type_t<T, MatrixS &> operator/=(T ele)
        {
            for (auto &i : this->m_data)
            {
                i /= ele;
            }
            return *this;
        }

        // Generate function.
        template <size_t A = M>
        friend enable_when_matrix_t<A, N, std::ostream &> operator<<(std::ostream &os, const MatrixS &self)
        {
            os << "MAT<" << M << ',' << N << ">:\n";
            for (auto i = 0u; i < M; i++)
            {
                for (auto j = 0u; j < N; j++)
                {
                    os << std::setw(12) << self(i, j) << '\t';
                }
                os << (i == M - 1 ? ' ' : '\n');
            }
            return os;
        }

        template <size_t A = M>
        friend enable_when_array_t<A, N, std::ostream &> operator<<(std::ostream &os, const MatrixS &self)
        {
            os << "ARR<" << M * N << ">:";
            for (auto j = 0u; j < M * N; j++)
            {
                os << std::setw(12) << self[j] << "\t";
            }
            return os;
        }
        
        constexpr static auto LEN = M * N;

        // Static function.
        static MatrixS eye()
        {
            MatrixS<M, N> result{};
            auto real_idx = M < N ? M : N;
            for (size_t i = 0; i < real_idx; i++)
            {
                result.data()[i * (M + 1)] = 1.0;
            }
            return result;
        }

        static MatrixS zeros()
        {
            return {};
        }

        constexpr size_t rows() const
        {
            return M;
        }

        constexpr size_t cols() const
        {
            return N;
        }

        constexpr size_t size() const
        {
            return M * N;
        }
    };
}
#endif
