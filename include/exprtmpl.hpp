#ifndef VVERY_SIMPLE_EXPR_TMPL_HEADER
#define VVERY_SIMPLE_EXPR_TMPL_HEADER

#include <limits>
#include <cmath>
#include <algorithm>

namespace ppx
{
    // constexpr
    constexpr double PI = 3.141592653589793;
    constexpr double EPS_SP = std::numeric_limits<float>::epsilon();
    constexpr double EPS_DP = std::numeric_limits<double>::epsilon();
    constexpr double MAX_SP = std::numeric_limits<float>::max();
    constexpr double MAX_DP = std::numeric_limits<double>::max();
    constexpr int MAJOR_VERSION = 1;
    constexpr int MINOR_VERSION = 1;

    constexpr double DEG_RAD(double deg)
    {
        return deg * PI / 180.0;
    }

    constexpr double RAD_DEG(double rad)
    {
        return 180.0 * rad / PI;
    }

    template <size_t A, size_t B, typename RT = void>
    using enable_when_array_t = std::enable_if_t<A == 1 || B == 1, RT>;
    template <size_t A, size_t B, typename RT = void>
    using enable_when_matrix_t = std::enable_if_t<A != 1 && B != 1, RT>;

    namespace details
    {
        // void_t
        template <typename... Ts>
        struct make_void
        {
            typedef void type;
        };

        template <typename... Ts>
        using void_t = typename make_void<Ts...>::type;

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

        // expr templates
        template <typename T>
        struct scalar
        {
        private:
            const T &s;

        public:
            using elem_tag = ElemTags::Scalar;

            constexpr scalar(const T &v) : s(v) {}

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

        template <typename T, typename Tag = typename T::elem_tag>
        struct elem_traits
        {
            using elem_tag = Tag;
            using elem_type = typename T::cast_type;
            using elem_ref = T const &;
        };

        template <typename T>
        struct elem_traits<scalar<T>, ElemTags::Scalar>
        {
            using elem_tag = ElemTags::Scalar;
            using elem_type = T;
            using elem_ref = scalar<T>;
        };

        template <typename T>
        struct elem_traits<T, ElemTags::Mblock>
        {
            using elem_tag = ElemTags::Mblock;
            using elem_type = typename T::cast_type;
            using elem_ref = typename T::cast_type const &;
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
                return (typename T::elem_type)self();
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
            using elem_tag = typename elem_traits<T>::elem_tag;
            using elem_type = typename elem_traits<T>::elem_type;
            using base_type = expr<expr_elem<T>>;
            using base_type::size;
            using base_type::operator[];
            friend base_type;

            template <typename U = elem_tag, std::enable_if_t<std::is_same<U, ElemTags::Scalar>::value> * = nullptr>
            explicit expr_elem(const T &val) : value(val)
            {
            }
            template <typename U = elem_tag, std::enable_if_t<std::is_same<U, ElemTags::Matrix>::value> * = nullptr>
            explicit expr_elem(const T &val) : value(val)
            {
            }
            template <typename U = elem_tag, std::enable_if_t<std::is_same<U, ElemTags::Mblock>::value> * = nullptr>
            explicit expr_elem(const T &val) : value(val.snap())
            {
            }

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
            typename elem_traits<T>::elem_ref value;
        };

        template <typename Ops, typename Expr>
        class unops : public expr<unops<Ops, Expr>>
        {
        public:
            using elem_tag = typename Expr::elem_tag;
            using elem_type = typename Expr::elem_type;

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
            using elem_tag =
                std::conditional_t<std::is_same<typename lExpr::elem_tag, ElemTags::Scalar>::value,
                                   typename rExpr::elem_tag, typename lExpr::elem_tag>;
            using elem_type =
                std::conditional_t<std::is_same<typename lExpr::elem_tag, ElemTags::Scalar>::value,
                                   typename rExpr::elem_type, typename lExpr::elem_type>;
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

        template <typename T, typename Tag = void>
        struct is_matrix : public std::false_type
        {
        };

        template <typename T>
        struct is_matrix<T, void_t<typename T::cast_type>> : public std::true_type
        {
        };

        template <typename T>
        constexpr bool is_expr_v()
        {
            return std::is_base_of<details::expr<T>, T>::value;
        }

        template <typename T>
        constexpr bool is_matrix_v()
        {
            return is_matrix<T>::value;
        }

        template <typename T, typename RT = void>
        using enable_arith_type_t = std::enable_if_t<std::is_arithmetic<T>::value, RT>;

        template <typename T, typename RT = void>
        using enable_expr_type_t = std::enable_if_t<details::is_expr_v<T>(), RT>;

        template <typename T, typename RT = void>
        using enable_matrix_type_t = std::enable_if_t<is_matrix_v<T>(), RT>;

        template <typename T>
        using result_t = details::expr_elem<T>;

        template <typename T>
        using result_s = details::expr_elem<details::scalar<T>>;

        template <typename T1, typename T2>
        using enable_expr_expr_t = std::enable_if_t<is_expr_v<T1>() && is_expr_v<T2>()>;
        template <typename T1, typename T2>
        using enable_expr_num_t = std::enable_if_t<is_expr_v<T1>() && std::is_arithmetic<T2>::value>;
        template <typename T1, typename T2>
        using enable_expr_mat_t = std::enable_if_t<is_expr_v<T1>() && is_matrix_v<T2>()>;

        template <typename T1, typename T2>
        using enable_num_expr_t = std::enable_if_t<std::is_arithmetic<T1>::value && is_expr_v<T2>()>;
        template <typename T1, typename T2>
        using enable_num_mat_t = std::enable_if_t<std::is_arithmetic<T1>::value && is_matrix_v<T2>()>;

        template <typename T1, typename T2>
        using enable_mat_mat_t = std::enable_if_t<is_matrix_v<T1>() && is_matrix_v<T2>()>;
        template <typename T1, typename T2>
        using enable_mat_expr_t = std::enable_if_t<is_matrix_v<T1>() && is_expr_v<T2>()>;
        template <typename T1, typename T2>
        using enable_mat_num_t = std::enable_if_t<is_matrix_v<T1>() && std::is_arithmetic<T2>::value>;
        template <typename T1, typename T2>
        using enable_mat_block_t = std::enable_if_t<std::is_same<typename T2::elem_tag, ElemTags::Mblock>::value && std::is_same<typename T2::elem_tag, ElemTags::Matrix>::value>;
        template <typename T1, typename T2>
        using enable_block_mat_t = std::enable_if_t<std::is_same<typename T1::elem_tag, ElemTags::Matrix>::value && std::is_same<typename T2::elem_tag, ElemTags::Mblock>::value>;
        template <typename T1, typename T2>
        using enable_block_block_t = std::enable_if_t<std::is_same<typename T1::elem_tag, ElemTags::Mblock>::value && std::is_same<typename T2::elem_tag, ElemTags::Mblock>::value>;
    } // namespace details

    // Ops abs
    template <typename T, details::enable_expr_type_t<T> * = nullptr>
    auto Abs(const T &t)
    {
        return details::unops<details::expr_abs_t, T>(details::expr_abs, t);
    }

    template <typename T, details::enable_matrix_type_t<T> * = nullptr>
    auto Abs(const T &t)
    {
        return details::unops<details::expr_abs_t,
                              details::result_t<T>>(details::expr_abs, details::result_t<T>(t));
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

    template <typename T1, typename T2, details::enable_expr_mat_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t, T1,
                              details::result_t<T2>>(details::expr_plus, t1, details::result_t<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_num_expr_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t,
                              details::result_s<T1>, T2>(details::expr_plus, details::result_s<T1>(t1), t2);
    }

    template <typename T1, typename T2, details::enable_num_mat_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t, details::result_s<T1>,
                              details::result_t<T2>>(details::expr_plus, details::result_s<T1>(t1), details::result_t<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_mat_expr_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t,
                              details::result_t<T1>, T2>(details::expr_plus, details::result_t<T1>(t1), t2);
    }

    template <typename T1, typename T2, details::enable_mat_num_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t, details::result_t<T1>,
                              details::result_s<T2>>(details::expr_plus, details::result_t<T1>(t1), details::result_s<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_mat_mat_t<T1, T2> * = nullptr>
    auto operator+(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_plus_t, details::result_t<T1>,
                              details::result_t<T2>>(details::expr_plus, details::result_t<T1>(t1), details::result_t<T2>(t2));
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

    template <typename T1, typename T2, details::enable_expr_mat_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t, T1,
                              details::result_t<T2>>(details::expr_minus, t1, details::result_t<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_num_expr_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t,
                              details::result_s<T1>, T2>(details::expr_minus, details::result_s<T1>(t1), t2);
    }

    template <typename T1, typename T2, details::enable_num_mat_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t, details::result_s<T1>,
                              details::result_t<T2>>(details::expr_minus, details::result_s<T1>(t1), details::result_t<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_mat_expr_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t,
                              details::result_t<T1>, T2>(details::expr_minus, details::result_t<T1>(t1), t2);
    }

    template <typename T1, typename T2, details::enable_mat_num_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t, details::result_t<T1>,
                              details::result_s<T2>>(details::expr_minus, details::result_t<T1>(t1), details::result_s<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_mat_mat_t<T1, T2> * = nullptr>
    auto operator-(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_minus_t, details::result_t<T1>,
                              details::result_t<T2>>(details::expr_minus, details::result_t<T1>(t1), details::result_t<T2>(t2));
    }

    // Ops *
    template <typename T1, typename T2, details::enable_expr_expr_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return t1.eval() * t2.eval();
    }

    template <typename T1, typename T2, details::enable_expr_num_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_mul_t, T1,
                              details::result_s<T2>>(details::expr_mul, t1, details::result_s<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_expr_mat_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return t1.eval() * t2;
    }

    template <typename T1, typename T2, details::enable_num_expr_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_mul_t,
                              details::result_s<T1>, T2>(details::expr_mul, details::result_s<T1>(t1), t2);
    }

    template <typename T1, typename T2, details::enable_num_mat_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_mul_t, details::result_s<T1>,
                              details::result_t<T2>>(details::expr_mul, details::result_s<T1>(t1), details::result_t<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_mat_num_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_mul_t, details::result_t<T1>,
                              details::result_s<T2>>(details::expr_mul, details::result_t<T1>(t1), details::result_s<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_mat_expr_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return t1 * t2.eval();
    }

    template <typename T1, typename T2, details::enable_mat_block_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return t1 * t2.snap();
    }

    template <typename T1, typename T2, details::enable_block_mat_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return t1.snap() * t2;
    }

    template <typename T1, typename T2, details::enable_block_block_t<T1, T2> * = nullptr>
    auto operator*(const T1 &t1, const T2 &t2)
    {
        return t1.snap() * t2.snap();
    }

    // Ops /
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

    template <typename T1, typename T2, details::enable_num_mat_t<T1, T2> * = nullptr>
    auto operator/(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_div_t, details::result_s<T1>,
                              details::result_t<T2>>(details::expr_div, details::result_s<T1>(t1), details::result_t<T2>(t2));
    }

    template <typename T1, typename T2, details::enable_mat_num_t<T1, T2> * = nullptr>
    auto operator/(const T1 &t1, const T2 &t2)
    {
        return details::biops<details::expr_div_t, details::result_t<T1>,
                              details::result_s<T2>>(details::expr_div, details::result_t<T1>(t1), details::result_s<T2>(t2));
    }

}

#endif