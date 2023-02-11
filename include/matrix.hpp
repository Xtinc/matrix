#ifndef VVERY_SIMPLE_MATRIX_HEADER
#define VVERY_SIMPLE_MATRIX_HEADER

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cassert>

namespace ppx
{
    // constexpr
    constexpr double PI = 3.141592653589793;
    constexpr double EPS_SP = std::numeric_limits<float>::epsilon();
    constexpr double EPS_DP = std::numeric_limits<double>::epsilon();
    constexpr double MAX_SP = std::numeric_limits<float>::max();
    constexpr double MAX_DP = std::numeric_limits<double>::max();

    constexpr double DEG_RAD(double deg)
    {
        return deg * PI / 180;
    }

    constexpr double RAD_DEG(double rad)
    {
        return 180 * rad / PI;
    }

    template <typename T, typename RT = void>
    using enable_arith_type_t = std::enable_if_t<std::is_arithmetic<T>::value, RT>;
    template <size_t A, size_t B, typename RT = void>
    using enable_when_array_t = std::enable_if_t<A == 1 || B == 1, RT>;
    template <size_t A, size_t B, typename RT = void>
    using enable_when_matrix_t = std::enable_if_t<A != 1 && B != 1, RT>;

    enum class StatusCode : char
    {
        NORMAL,
        CONVERGED,
        DIVERGED,
        SINGULAR
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

        template <typename T>
        inline void PRINT_SINGLE_ELEMENTS(const T &coll, const std::string &optcsrt = "")
        {
            std::cout << optcsrt << coll << std::endl;
        }

        template <typename T>
        inline void PRINT_LISTED_ELEMENTS(const T &coll, const std::string &optcsrt = "")
        {
            std::cout << optcsrt;
            for (const auto ele : coll)
            {
                std::cout << ele << ' ';
            }
            std::cout << std::endl;
        }

        // operators
        struct expr_plus_t
        {
            constexpr explicit expr_plus_t() = default;

            template <typename LType, typename RType>
            auto operator()(const LType &lhs, const RType &rhs) const
            {
                return lhs + rhs;
            }
        };

        struct expr_minus_t
        {
            constexpr explicit expr_minus_t() = default;

            template <typename LType, typename RType>
            auto operator()(const LType &lhs, const RType &rhs) const
            {
                return lhs - rhs;
            }
        };

        struct expr_mul_t
        {
            constexpr explicit expr_mul_t() = default;

            template <typename LType, typename RType>
            auto operator()(const LType &lhs, const RType &rhs) const
            {
                return lhs * rhs;
            }
        };

        struct expr_div_t
        {
            constexpr explicit expr_div_t() = default;

            template <typename LType, typename RType>
            auto operator()(const LType &lhs, const RType &rhs) const
            {
                return lhs / rhs;
            }
        };

        struct expr_abs_t
        {
            constexpr explicit expr_abs_t() = default;

            template <typename Type>
            auto operator()(const Type &ele) const
            {
                return std::fabs(ele);
            }
        };

        constexpr expr_plus_t expr_plus{};
        constexpr expr_minus_t expr_minus{};
        constexpr expr_mul_t expr_mul{};
        constexpr expr_div_t expr_div{};
        constexpr expr_abs_t expr_abs{};

        struct ElemTags
        {
            struct Scalar;
            struct Matrix;
            struct Mblock;
        };
        
        struct expr_mtx_process;
        struct expr_sca_process;

        // expr templates
        template <typename T>
        struct expr_scalar
        {
        private:
            const T &s;

        public:
            constexpr expr_scalar(const T &v) : s(v) {}

            constexpr T const &

            operator[](std::size_t) const
            {
                return s;
            }

            constexpr std::size_t

            size() const
            {
                return 0;
            }
        };

        template <typename T>
        struct expr_traits
        {
            using value_tag = expr_mtx_process;
            using value_type = T;
            using expr_ref = T const &;
        };

        template <typename T>
        struct expr_traits<expr_scalar<T>>
        {
            using value_tag = expr_sca_process;
            using value_type = T;
            using expr_ref = expr_scalar<T>;
        };

        template <typename T>
        class expr
        {
        public:
            const T &self() const
            {
                return static_cast<const T &>(*this);
            }

            T &self()
            {
                return static_cast<T &>(*this);
            }

            auto eval() const
            {
                return (typename T::value_type)self();
            }

        protected:
            explicit expr() = default;

            constexpr size_t size()
            {
                return self().size_impl();
            }

            auto operator[](size_t idx) const
            {
                return self().at_impl(idx);
            }

            auto operator()() const
            {
                return self()();
            }
        };

        template <typename T>
        class expr_elem : expr<expr_elem<T>>
        {
        public:
            using value_tag = typename expr_traits<T>::value_tag;
            using value_type = typename expr_traits<T>::value_type;
            using base_type = expr<expr_elem<T>>;
            using base_type::size;
            using base_type::operator[];
            friend base_type;

            explicit expr_elem(const T &val) : value(val) {}

            constexpr size_t size_impl()
            {
                return value.size();
            };

            auto at_impl(size_t idx) const
            {
                return value[idx];
            };

            decltype(auto) operator()() const
            {
                return (value);
            }

        private:
            typename expr_traits<T>::expr_ref value;
        };

        template <typename Ops, typename Expr>
        class unops : public expr<unops<Ops, Expr>>
        {
        public:
            using value_tag = typename Expr::value_tag;
            using value_type = typename Expr::value_type;

            using base_type = expr<unops<Ops, Expr>>;
            using base_type::size;
            using base_type::operator[];
            friend base_type;

            explicit unops(const Ops &ops, const Expr &expr)
                : m_ops(ops), m_expr(expr) {}

            constexpr size_t size_impl()
            {
                return m_expr.size();
            }

            auto at_impl(size_t idx) const
            {
                return m_ops(m_expr[idx]);
            }

            template <typename T>
            operator T() const
            {
                T res{};
                for (size_t idx = 0; idx < res.size(); ++idx)
                {
                    res[idx] = (*this)[idx];
                }
                return res;
            }

        private:
            Ops m_ops;
            Expr m_expr;
        };

        template <typename Ops, typename lExpr, typename rExpr>
        class biops : public expr<biops<Ops, lExpr, rExpr>>
        {
        public:
            using value_tag =
                std::conditional_t<std::is_same<typename lExpr::value_tag, expr_sca_process>::value,
                                   typename rExpr::value_tag, typename lExpr::value_tag>;
            using value_type =
                std::conditional_t<std::is_same<typename lExpr::value_tag, expr_sca_process>::value,
                                   typename rExpr::value_type, typename lExpr::value_type>;
            using base_type = expr<biops<Ops, lExpr, rExpr>>;
            using base_type::size;
            using base_type::operator[];
            friend base_type;

            explicit biops(const Ops &ops, const lExpr &lxpr, const rExpr &rxpr) : m_ops(ops), m_lxpr(lxpr), m_rxpr(rxpr) {}

            constexpr size_t size_impl()
            {
                return std::max(m_lxpr.size(), m_rxpr.size());
            }

            auto at_impl(size_t idx) const
            {
                return m_ops(m_lxpr[idx], m_rxpr[idx]);
            }
            // dangerous but powerful;
            template <typename T>
            operator T() const
            {
                T res{};
                for (size_t idx = 0; idx < res.size(); ++idx)
                {
                    res[idx] = (*this)[idx];
                }
                return res;
            }

        private:
            Ops m_ops;
            lExpr m_lxpr;
            rExpr m_rxpr;
        };

        template <typename T>
        constexpr bool is_expr_v()
        {
            return std::is_base_of<details::expr<T>, T>::value;
        }

        template <typename T1, typename T2>
        using enable_expr_expr_t = std::enable_if_t<is_expr_v<T1>() && is_expr_v<T2>()>;
        template <typename T1, typename T2>
        using enable_expr_num_t = std::enable_if_t<is_expr_v<T1>() && std::is_arithmetic<T2>::value>;
        template <typename T1, typename T2>
        using enable_num_expr_t = std::enable_if_t<std::is_arithmetic<T1>::value && is_expr_v<T2>()>;

    } // namespace details

    template <typename T, typename RT = void>
    using enable_expr_type_t = std::enable_if_t<details::is_expr_v<T>(), RT>;

    // forward declare
    template <size_t M, size_t N>
    class Matrix;

    template <size_t N>
    Matrix<N, N> ludcmp(Matrix<N, N> A, std::array<int, N> &indx, bool &even, bool &sing);

    template <size_t M, size_t N>
    Matrix<M, N> svdcmp(Matrix<M, N> u, Matrix<N, 1> &w, Matrix<N, N> &v, bool &sing);

    template <size_t N>
    void ludbksb(const Matrix<N, N> &A, const std::array<int, N> &indx, double *b);

    template <size_t M, size_t N>
    void svbksb(const Matrix<M, N> &u, const Matrix<N, 1> &w, const Matrix<N, N> &v, double *b);

    namespace details
    {
        constexpr size_t is_samll_matrix(size_t A, size_t B)
        {
            return A * B < MAX_SIZE_LIMIT ? 1 : 0;
        }

        template <size_t M, size_t N, std::size_t A = is_samll_matrix(M, N)>
        class MatrixBase;

        template <std::size_t M, std::size_t N>
        class MatrixBase<M, N, 1>
        {
        protected:
            MatrixBase() : m_data{} {}

            template <typename T, size_t L, enable_arith_type_t<T> * = nullptr>
            explicit MatrixBase(const std::array<T, L> &list) : m_data{}
            {
                constexpr auto real_idx = std::min(L, M * N);
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

            template <typename T, enable_arith_type_t<T> * = nullptr>
            MatrixBase(const std::initializer_list<T> &list) : m_data{}
            {
                auto real_idx = list.size() < M * N ? list.size() : M * N;
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

            template <typename T, enable_arith_type_t<T> * = nullptr>
            explicit MatrixBase(const std::vector<T> &list) : m_data{}
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
            MatrixBase() : m_data(M * N, 0.0) {}

            template <typename T, size_t L, enable_arith_type_t<T> * = nullptr>
            explicit MatrixBase(const std::array<T, L> &list) : m_data(M * N, 0.0)
            {
                constexpr auto real_idx = std::min(L, M * N);
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

            MatrixBase(const std::initializer_list<int> &list) : m_data(M * N, 0.0)
            {
                auto real_idx = list.size() < M * N ? list.size() : M * N;
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

            template <typename T, enable_arith_type_t<T> * = nullptr>
            explicit MatrixBase(const std::vector<T> &list) : m_data(M * N, 0.0)
            {
                auto real_idx = list.size() < M * N ? list.size() : M * N;
                std::copy_n(list.begin(), real_idx, m_data.begin());
            }

        protected:
            std::vector<double> m_data;
        };

    } // namespace details

    template <std::size_t M, std::size_t N>
    class Matrix : public details::MatrixBase<M, N>
    {
        template <size_t A, typename RT = void>
        using enable_when_squre_t = std::enable_if_t<A == N, RT>;
        template <size_t A, typename RT = void>
        using disable_when_squre_t = std::enable_if_t<A != N, RT>;
        using IndexRange = std::pair<int, int>;

        template <size_t A, size_t B>
        class SubMatrix
        {

        public:
            SubMatrix(Matrix<M, N> &self, size_t r, size_t c)
                : row_idx(r), col_idx(c), data(self)
            {
                assert(row_idx + A <= M && col_idx + B <= N);
            }

            SubMatrix(const SubMatrix &) = delete;

            SubMatrix(SubMatrix &&) = delete;

            SubMatrix &operator=(const SubMatrix &other) = delete;

            SubMatrix &operator=(SubMatrix &&other) = delete;

            SubMatrix &operator=(const Matrix<A, B> &other)
            {
                assert(row_idx <= row_end && col_idx <= col_end);
                for (size_t i = 0; i < A; i++)
                {
                    for (size_t j = 0; j < B; j++)
                    {
                        data(row_idx + i, col_idx + j) = other(i, j);
                    }
                }
                return *this;
            }

            template <typename T>
            enable_expr_type_t<T, SubMatrix &>
            operator=(const T &expr)
            {
                (*this) = expr.eval();
                return *this;
            }

        private:
            size_t row_idx;
            size_t col_idx;
            Matrix<M, N> &data;

            Matrix<A, B> val() const
            {
                Matrix<A, B> result;
                for (size_t i = 0; i < M; i++)
                {
                    for (size_t j = 0; j < N; j++)
                    {
                        result(i, j) = data(row_idx + i, col_idx + j);
                    }
                }
                return result;
            }
        };

        struct SubPart
        {
        public:
            SubPart(Matrix<M, N> &self, IndexRange r, IndexRange c)
                : row_idx(0), col_idx(0), row_end(0), col_end(0), data(self)
            {
                row_idx = r.first >= 0 ? r.first : 0;
                col_idx = c.first >= 0 ? c.first : 0;
                row_end = r.second >= 0 ? r.second : M - 1;
                col_end = c.second >= 0 ? c.second : N - 1;
            }

            SubPart(const SubPart &) = delete;

            SubPart(SubPart &&) = delete;

            SubPart &operator=(const SubPart &other) = delete;

            SubPart &operator=(SubPart &&other) = delete;

            template <size_t A, size_t B>
            SubPart &operator=(const Matrix<A, B> &other)
            {
                assert(row_idx <= row_end && col_idx <= col_end);
                auto real_row_counts = row_end - row_idx + 1;
                auto real_col_counts = col_end - col_idx + 1;
                real_row_counts = std::min(real_row_counts, A);
                real_col_counts = std::min(real_col_counts, B);
                for (size_t i = 0; i < real_row_counts; i++)
                {
                    for (size_t j = 0; j < real_col_counts; j++)
                    {
                        data(row_idx + i, col_idx + j) = other(i, j);
                    }
                }
                return *this;
            }

            template <typename T>
            std::enable_if_t<details::is_expr_v<T>(), SubPart &>
            operator=(const T &expr)
            {
                (*this) = expr.eval();
                return *this;
            }

            template <typename T, size_t A, enable_arith_type_t<T> * = nullptr>
            SubPart &operator=(const std::array<T, A> &list)
            {
                generator_by_list(list);
                return *this;
            }

            template <typename T, enable_arith_type_t<T> * = nullptr>
            SubPart &operator=(const std::initializer_list<T> &list)
            {
                generator_by_list(list);
                return *this;
            }

            template <typename T, enable_arith_type_t<T> * = nullptr>
            SubPart &operator=(const std::vector<T> &list)
            {
                generator_by_list(list);
                return *this;
            }

        private:
            size_t row_idx;
            size_t col_idx;
            size_t row_end;
            size_t col_end;
            Matrix<M, N> &data;

            template <typename T>
            void generator_by_list(const T &list)
            {
                auto A = list.size();
                auto iter = list.begin();
                auto real_row_counts = row_end - row_idx + 1;
                auto real_col_counts = col_end - col_idx + 1;
                if (real_col_counts == 1)
                {
                    real_row_counts = std::min(real_row_counts, A);
                    for (size_t i = 0; i < real_row_counts; i++)
                    {
                        data(row_idx + i, col_idx) = *(iter++);
                    }
                    return;
                }
                if (real_row_counts == 1)
                {
                    real_col_counts = std::min(real_col_counts, A);
                    for (size_t i = 0; i < real_col_counts; i++)
                    {
                        data(row_idx, col_idx + i) = *(iter++);
                    }
                    return;
                }
            }
        };

    public:
        using value_type = double;
        struct iterator : public std::iterator<std::random_access_iterator_tag, double>
        {
        public:
            using iterator_category = std::random_access_iterator_tag;
            iterator(value_type *ptr) noexcept : m_ptr(ptr) {}
            // iterator(const iterator &itr) noexcept : m_ptr(itr.m_ptr) {}

            pointer operator->() const noexcept
            {
                return m_ptr;
            }
            reference operator*() const noexcept
            {
                return *m_ptr;
            }

            iterator &operator=(const iterator &rhs) noexcept
            {
                m_ptr = rhs.m_ptr;
                return *this;
            }
            bool operator==(const iterator &rhs) const noexcept
            {
                return m_ptr == rhs.m_ptr;
            }
            bool operator!=(const iterator &rhs) const noexcept
            {
                return m_ptr != rhs.m_ptr;
            }
            bool operator>(const iterator &rhs) const noexcept
            {
                return m_ptr > rhs.m_ptr;
            }
            bool operator>=(const iterator &rhs) const noexcept
            {
                return m_ptr >= rhs.m_ptr;
            }
            bool operator<(const iterator &rhs) const noexcept
            {
                return m_ptr < rhs.m_ptr;
            }
            bool operator<=(const iterator &rhs) const noexcept
            {
                return m_ptr <= rhs.m_ptr;
            }

            iterator &operator++() noexcept
            {
                m_ptr += 1;
                return *this;
            }
            iterator operator++(int) noexcept
            {
                value_type *ret = m_ptr;
                m_ptr += 1;
                return ret;
            }
            iterator &operator--() noexcept
            {
                m_ptr -= 1;
                return *this;
            }
            iterator operator--(int) noexcept
            {
                value_type *ret = m_ptr;
                m_ptr -= 1;
                return ret;
            }
            iterator &operator+=(difference_type step) noexcept
            {
                m_ptr += step;
                return *this;
            }
            iterator &operator-=(difference_type step) noexcept
            {
                m_ptr -= step;
                return *this;
            }
            iterator operator+(difference_type step) noexcept
            {
                value_type *ret = m_ptr;
                ret += step;
                return ret;
            }
            iterator operator-(difference_type step) noexcept
            {
                value_type *ret = m_ptr;
                ret -= step;
                return ret;
            }
            iterator operator[](difference_type n) noexcept
            {
                m_ptr += n;
                return *this;
            }
            difference_type operator-(const iterator &rhs) const noexcept
            {
                return m_ptr - rhs.m_ptr;
            }

        private:
            value_type *m_ptr{};
        };
        struct const_iterator : public std::iterator<std::random_access_iterator_tag, const double>
        {
        public:
            using iterator_category = std::random_access_iterator_tag;
            const_iterator(value_type *const ptr) noexcept : m_ptr(ptr)
            {
            }

            pointer operator->() const noexcept
            {
                return m_ptr;
            }
            reference operator*() const noexcept
            {
                return *m_ptr;
            }

            const_iterator(const const_iterator &itr) noexcept : m_ptr(itr.m_ptr)
            {
            }
            const_iterator &operator=(const const_iterator &rhs) noexcept
            {
                m_ptr = rhs.m_ptr;
                return *this;
            }
            bool operator==(const const_iterator &rhs) const noexcept
            {
                return m_ptr == rhs.m_ptr;
            }
            bool operator!=(const const_iterator &rhs) const noexcept
            {
                return m_ptr != rhs.m_ptr;
            }
            bool operator>(const const_iterator &rhs) const noexcept
            {
                return m_ptr > rhs.m_ptr;
            }
            bool operator>=(const const_iterator &rhs) const noexcept
            {
                return m_ptr >= rhs.m_ptr;
            }
            bool operator<(const const_iterator &rhs) const noexcept
            {
                return m_ptr < rhs.m_ptr;
            }
            bool operator<=(const const_iterator &rhs) const noexcept
            {
                return m_ptr <= rhs.m_ptr;
            }

            const_iterator &operator++() noexcept
            {
                m_ptr++;
                return *this;
            }
            const_iterator operator++(int) noexcept
            {
                value_type const *ret = m_ptr;
                m_ptr++;
                return ret;
            }
            const_iterator &operator--() noexcept
            {
                m_ptr--;
                return *this;
            }
            const_iterator operator--(int) noexcept
            {
                value_type const *ret = m_ptr;
                m_ptr--;
                return ret;
            }
            const_iterator &operator+=(ptrdiff_t step) noexcept
            {
                m_ptr += step;
                return *this;
            }
            const_iterator &operator-=(ptrdiff_t step) noexcept
            {
                m_ptr -= step;
                return *this;
            }
            const_iterator operator+(ptrdiff_t step) noexcept
            {
                value_type const *ret = m_ptr;
                ret += step;
                return ret;
            }
            const_iterator operator-(ptrdiff_t step) noexcept
            {
                value_type const *ret = m_ptr;
                ret -= step;
                return ret;
            }
            const_iterator operator[](difference_type n) noexcept
            {
                m_ptr += n;
                return *this;
            }
            difference_type operator-(const const_iterator &rhs) const noexcept
            {
                return m_ptr - rhs.m_ptr;
            }

        private:
            value_type *m_ptr{};
        };
        iterator begin() noexcept
        {
            return this->m_data.data();
        }
        iterator end() noexcept
        {
            return this->m_data.data() + M * N;
        }
        const_iterator begin() const noexcept
        {
            return this->m_data.data();
        }
        const_iterator end() const noexcept
        {
            return this->m_data.data() + M * N;
        }
        const_iterator cbegin() const noexcept
        {
            return this->m_data.data();
        }
        const_iterator cend() const noexcept
        {
            return this->m_data.data() + M * N;
        }

    public:
        Matrix() = default;

        template <typename T, size_t L, enable_arith_type_t<T> * = nullptr>
        Matrix(const std::array<T, L> &list) : details::MatrixBase<M, N>(list)
        {
        }

        template <typename T, enable_arith_type_t<T> * = nullptr>
        Matrix(const std::initializer_list<T> &list) : details::MatrixBase<M, N>(list)
        {
        }

        template <typename T, enable_arith_type_t<T> * = nullptr>
        Matrix(const std::vector<T> &list) : details::MatrixBase<M, N>(list)
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

        Matrix<N, 1> row(size_t idx) const
        {
            assert(idx < M);
            Matrix<N, 1> result;
            for (auto i = 0u; i < N; ++i)
            {
                result[i] = this->m_data.at(idx + i * M);
            }
            return result;
        }

        Matrix<M, 1> col(size_t idx) const
        {
            assert(idx < N);
            Matrix<M, 1> result;
            for (auto i = 0u; i < M; ++i)
            {
                result[i] = this->m_data.at(idx * M + i);
            }
            return result;
        }

        SubPart operator()(const IndexRange &row_range, const IndexRange &col_range)
        {
            return {*this, row_range, col_range};
        }

        SubPart operator()(size_t row_idx, const IndexRange &col_range)
        {
            return {*this, {(int)row_idx, (int)row_idx}, col_range};
        }

        SubPart operator()(const IndexRange &row_range, size_t col_idx)
        {
            return {*this, row_range, {(int)col_idx, (int)col_idx}};
        }

        template <size_t A, size_t B>
        SubMatrix<A, B> sub(size_t row_start, size_t col_start)
        {
            return {*this, row_start, col_start};
        }

        void fill(double val)
        {
            std::fill(this->m_data.begin(), this->m_data.end(), val);
        }

        Matrix<N, M> T() const
        {
            Matrix<N, M> res{};
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
        enable_when_squre_t<A, Matrix> I() const
        {
            std::array<int, M> indx{};
            auto even = true;
            auto sing = false;
            auto LU = ludcmp(*this, indx, even, sing);
            if (sing)
            {
                return {};
            }
            auto result = Matrix<M, M>::eye();
            for (size_t j = 0; j < M; j++)
            {
                ludbksb(LU, indx, result.data() + j * M);
            }
            return result;
        }

        template <size_t A = M>
        disable_when_squre_t<A, Matrix<N, M>> I() const
        {
            Matrix<N, 1> w{};
            Matrix<N, N> V{};
            auto U = svdcmp(*this, w, V);
            Matrix<N, N> W{};
            for (size_t i = 0; i < N; i++)
            {
                if (fabs(w[i]) > EPS_SP)
                {
                    W(i, i) = 1.0 / w[i];
                }
            }
            return V * W * U.T();
        }

        template <size_t A = M>
        enable_when_squre_t<A, double> det() const
        {
            auto even = true;
            auto sing = false;
            std::array<int, M> indx{};
            auto LU = ludcmp(*this, indx, even, sing);
            if (sing)
            {
                return {};
            }
            auto D = even ? 1.0 : -1.0;
            for (size_t i = 0; i < M; i++)
            {
                D *= LU(i, i);
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
        enable_when_matrix_t<M, N, Matrix<A, 1>> diag() const
        {
            Matrix<A, 1> result;
            for (size_t i = 0; i < A; i++)
            {
                result[i] = (*this)(i, i);
            }
            return result;
        }

        template <size_t A = std::max(M, N)>
        enable_when_array_t<M, N, Matrix<A, A>> diag() const
        {
            Matrix<A, A> result;
            for (size_t i = 0; i < A; i++)
            {
                result(i, i) = (*this)[i];
            }
            return result;
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
        Matrix<M, L> operator*(const Matrix<N, L> &other) const
        {
            Matrix<M, L> result;
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
            // Matrix<N, M> AT = this->T();
            // auto iter_A = AT.data();
            // auto iter_B = other.data();
            // for (size_t i = 0; i < M; i++)
            // {
            //     iter_B = other.data();
            //     for (size_t j = 0; j < L; j++)
            //     {
            //         result(i, j) = std::inner_product(iter_A, iter_A + N, iter_B, 0.0);
            //         iter_B += N;
            //     }
            //     iter_A += N;
            // }
            return result;
        }

        Matrix &operator+=(const Matrix &other)
        {
            for (size_t i = 0; i < M * N; i++)
            {
                this->m_data[i] += other.data()[i];
            }
            return *this;
        }

        Matrix &operator-=(const Matrix &other)
        {
            for (size_t i = 0; i < M * N; i++)
            {
                this->m_data[i] -= other.data()[i];
            }
            return *this;
        }

        bool operator==(const Matrix &other) const
        {
            return std::equal(this->m_data.begin(), this->m_data.end(), other.data(),
                              [](double ele1, double ele2)
                              { return details::is_same(ele1, ele2); });
        }

        template <typename T>
        enable_arith_type_t<T, Matrix &> operator+=(T ele)
        {
            for (auto &i : this->m_data)
            {
                i += ele;
            }
            return *this;
        }

        template <typename T>
        enable_arith_type_t<T, Matrix &> operator-=(T ele)
        {
            for (auto &i : this->m_data)
            {
                i -= ele;
            }
            return *this;
        }

        template <typename T>
        enable_arith_type_t<T, Matrix &> operator*=(T ele)
        {
            for (auto &i : this->m_data)
            {
                i *= ele;
            }
            return *this;
        }

        template <typename T>
        enable_arith_type_t<T, Matrix &> operator/=(T ele)
        {
            for (auto &i : this->m_data)
            {
                i /= ele;
            }
            return *this;
        }

        // Generate function.
        friend std::ostream &operator<<(std::ostream &os, const Matrix &self)
        {
            os << "Matrix<" << M << "," << N << ">:\n";
            for (auto i = 0u; i < M; i++)
            {
                for (auto j = 0u; j < N; j++)
                {
                    os << std::setw(12) << self(i, j) << "\t";
                }
                os << std::endl;
            }
            return os;
        }

        // Static function.
        static Matrix eye()
        {
            Matrix<M, N> result{};
            auto real_idx = M < N ? M : N;
            for (size_t i = 0; i < real_idx; i++)
            {
                result.data()[i * (M + 1)] = 1.0;
            }
            return result;
        }

        static Matrix zero()
        {
            Matrix<M, N> result{};
            result.fill(0.0);
            return result;
        }

        constexpr size_t row_counts()
        {
            return M;
        }

        constexpr size_t col_counts()
        {
            return N;
        }

        constexpr size_t size()
        {
            return M * N;
        }
    };

    namespace details
    {
        template <size_t M, size_t N>
        using result_t = details::expr_elem<Matrix<M, N>>;

        template <typename T>
        using result_s = details::expr_elem<details::expr_scalar<T>>;
    } // namespace details

    // Ops abs
    template <typename T, std::enable_if_t<details::is_expr_v<T>()> * = nullptr>
    auto Abs(const T &t)
    {
        return details::unops<details::expr_abs_t, T>(details::expr_abs, t);
    }

    template <size_t M, size_t N>
    auto Abs(const Matrix<M, N> &t)
    {
        return details::unops<details::expr_abs_t,
                              details::result_t<M, N>>(details::expr_abs, details::result_t<M, N>(t));
    }

    // Ops +
    template <typename T1, typename T2, details::enable_expr_expr_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t, T1, T2>(details::expr_plus, t1, t2);
    }

    template <typename T1, typename T2, details::enable_expr_num_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t, T1,
                              details::result_s<T2>>(details::expr_plus, t1, details::result_s<T2>(t2));
    }

    template <typename T1, size_t M, size_t N, enable_expr_type_t<T1> * = nullptr>
    auto operator+(const T1 &t1, const Matrix<M, N> &t2)
    {
        return details::biops<details::expr_plus_t, T1,
                              details::result_t<M, N>>(details::expr_plus, t1, details::result_t<M, N>(t2));
    }

    template <typename T1, typename T2, details::enable_num_expr_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t,
                              details::result_s<T1>, T2>(details::expr_plus, details::result_s<T1>(t1), t2);
    }

    template <typename T1, size_t M, size_t N, enable_arith_type_t<T1> * = nullptr>
    auto operator+(const T1 &t1, const Matrix<M, N> &t2)
    {
        return details::biops<details::expr_plus_t, details::result_s<T1>,
                              details::result_t<M, N>>(details::expr_plus, details::result_s<T1>(t1), details::result_t<M, N>(t2));
    }

    template <size_t M, size_t N, typename T2, enable_expr_type_t<T2> * = nullptr>
    auto operator+(const Matrix<M, N> &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t,
                              details::result_t<M, N>, T2>(details::expr_plus, details::result_t<M, N>(t1), t2);
    }

    template <size_t M, size_t N, typename T2, enable_arith_type_t<T2> * = nullptr>
    auto operator+(const Matrix<M, N> &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t, details::result_t<M, N>,
                              details::result_s<T2>>(details::expr_plus, details::result_t<M, N>(t1), details::result_s<T2>(t2));
    }

    template <size_t M, size_t N>
    auto operator+(const Matrix<M, N> &t1, const Matrix<M, N> &t2)
    {
        return details::biops<details::expr_plus_t, details::result_t<M, N>,
                              details::result_t<M, N>>(details::expr_plus, details::result_t<M, N>(t1), details::result_t<M, N>(t2));
    }

    // Ops -
    template <typename T1, typename T2, details::enable_expr_expr_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t, T1, T2>(details::expr_minus, t1, t2);
    }

    template <typename T1, typename T2, details::enable_expr_num_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t, T1,
                              details::result_s<T2>>(details::expr_minus, t1, details::result_s<T2>(t2));
    }

    template <typename T1, size_t M, size_t N, enable_expr_type_t<T1> * = nullptr>
    auto operator-(const T1 &t1, const Matrix<M, N> &t2)
    {
        return details::biops<details::expr_minus_t, T1,
                              details::result_t<M, N>>(details::expr_minus, t1, details::result_t<M, N>(t2));
    }

    template <typename T1, typename T2, details::enable_num_expr_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t,
                              details::result_s<T1>, T2>(details::expr_minus, details::result_s<T1>(t1), t2);
    }

    template <typename T1, size_t M, size_t N, enable_arith_type_t<T1> * = nullptr>
    auto operator-(const T1 &t1, const Matrix<M, N> &t2)
    {
        return details::biops<details::expr_minus_t, details::result_s<T1>,
                              details::result_t<M, N>>(details::expr_minus, details::result_s<T1>(t1), details::result_t<M, N>(t2));
    }

    template <size_t M, size_t N, typename T2, enable_expr_type_t<T2> * = nullptr>
    auto operator-(const Matrix<M, N> &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t,
                              details::result_t<M, N>, T2>(details::expr_minus, details::result_t<M, N>(t1), t2);
    }

    template <size_t M, size_t N, typename T2, enable_arith_type_t<T2> * = nullptr>
    auto operator-(const Matrix<M, N> &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t, details::result_t<M, N>,
                              details::result_s<T2>>(details::expr_minus, details::result_t<M, N>(t1), details::result_s<T2>(t2));
    }

    template <size_t M, size_t N>
    auto operator-(const Matrix<M, N> &t1, const Matrix<M, N> &t2)
    {
        return details::biops<details::expr_minus_t, details::result_t<M, N>,
                              details::result_t<M, N>>(details::expr_minus, details::result_t<M, N>(t1), details::result_t<M, N>(t2));
    }

    // Ops *
    template <typename T1, typename T2, details::enable_expr_expr_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_mul_t, T1, T2>(details::expr_mul, t1, t2);
    }

    template <typename T1, typename T2, details::enable_expr_num_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_mul_t, T1,
                              details::result_s<T2>>(details::expr_mul, t1, details::result_s<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_num_expr_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_mul_t,
                              details::result_s<T1>, T2>(details::expr_mul, details::result_s<T1>(t1), t2);
    }

    template <typename T1, size_t M, size_t N, enable_arith_type_t<T1> * = nullptr>
    auto operator*(const T1 &t1, const Matrix<M, N> &t2)
    {
        return details::biops<details::expr_mul_t, details::result_s<T1>,
                              details::result_t<M, N>>(details::expr_mul, details::result_s<T1>(t1), details::result_t<M, N>(t2));
    }

    template <size_t M, size_t N, typename T2, enable_arith_type_t<T2> * = nullptr>
    auto operator*(const Matrix<M, N> &t1, const T2 &t2)
    {
        return details::biops<details::expr_mul_t, details::result_t<M, N>,
                              details::result_s<T2>>(details::expr_mul, details::result_t<M, N>(t1), details::result_s<T2>(t2));
    }

    // Ops /
    template <typename T1, typename T2, details::enable_expr_expr_t<T1, T2> * = nullptr>
    auto operator/(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_div_t, T1, T2>(details::expr_div, t1, t2);
    }

    template <typename T1, typename T2, details::enable_expr_num_t<T1, T2> * = nullptr>
    auto operator/(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_div_t, T1,
                              details::result_s<T2>>(details::expr_div, t1, details::result_s<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_num_expr_t<T1, T2> * = nullptr>
    auto operator/(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_div_t,
                              details::result_s<T1>, T2>(details::expr_div, details::result_s<T1>(t1), t2);
    }

    template <typename T1, size_t M, size_t N, enable_arith_type_t<T1> * = nullptr>
    auto operator/(const T1 &t1, const Matrix<M, N> &t2)
    {
        return details::biops<details::expr_div_t, details::result_s<T1>,
                              details::result_t<M, N>>(details::expr_div, details::result_s<T1>(t1), details::result_t<M, N>(t2));
    }

    template <size_t M, size_t N, typename T2, enable_arith_type_t<T2> * = nullptr>
    auto operator/(const Matrix<M, N> &t1, const T2 &t2)
    {
        return details::biops<details::expr_div_t, details::result_t<M, N>,
                              details::result_s<T2>>(details::expr_div, details::result_t<M, N>(t1), details::result_s<T2>(t2));
    }
}
#endif
