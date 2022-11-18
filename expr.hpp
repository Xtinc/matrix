#ifndef VVERY_SIMPLE_EXPRESSION_HEADER
#define VVERY_SIMPLE_EXPRESSION_HEADER

#include <iostream>
#include <iomanip>

constexpr double gl_rep_eps = std::numeric_limits<float>::epsilon();
constexpr size_t gl_get_less(size_t A, size_t B)
{
    return A < B ? A : B;
}
constexpr size_t gl_get_more(size_t A, size_t B)
{
    return A < B ? B : A;
}

template <typename T, typename RT = void>
using enable_arith_type_t = typename std::enable_if<std::is_arithmetic<T>::value, RT>::type;
template <typename T, typename RT = void>
using disable_arith_type_t = typename std::enable_if<!std::is_arithmetic<T>::value, RT>::type;

namespace details
{
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
    template <typename T>
    inline void PRINT_MAPPED_ELEMENTS(const T &coll, const std::string &optcsrt = "")
    {
        std::cout << optcsrt;
        for (auto ele : coll)
        {
            std::cout << '[' << ele.first << ',' << ele.second << "] ";
        }
        std::cout << std::endl;
    }

    // expr templates
    template <typename T>
    struct expr_scalar
    {
    private:
        T const &s;

    public:
        constexpr expr_scalar(T const &v)
            : s(v)
        {
        }
        constexpr T const &operator[](std::size_t) const
        {
            return s;
        }
        constexpr std::size_t size() const
        {
            return 0;
        }
    };

    template <typename T>
    class expr
    {
    public:
        using expr_type = expr<T>;
        const T &self() const { return static_cast<const T &>(*this); }
        T &self() { return static_cast<T &>(*this); }

    protected:
        explicit expr(){};
        constexpr size_t size() { return self().size_impl(); }
        auto operator[](size_t idx) const { return self().at_impl(idx); }
        auto operator()() const { return self()(); };
    };

    template <typename T>
    class expr_result : expr<expr_result<T>>
    {
    public:
        using base_type = expr<expr_result<T>>;
        using base_type::size;
        using base_type::operator[];
        friend base_type;

        explicit expr_result(const T &val) : value(val) {}
        size_t size_impl() const { return value.size(); };
        auto at_impl(size_t idx) const { return value[idx]; };
        decltype(auto) operator()() const { return (value); }

        // template <typename T, disable_arith_type_t<T> * = nullptr>
        // auto operator+(const T &rhs)
        // {
        //     return biops<expr_plus_t, biops<Ops, expr_result<T>, rExpr>, T>(expr_plus, *this, rhs);
        // }
        // template <typename T, enable_arith_type_t<T> * = nullptr>
        // auto operator+(const T &rhs)
        // {
        //     return biops<expr_plus_t, biops<Ops, expr_result<T>, rExpr>, expr_scalar<T>>(expr_plus, *this, rhs);
        // }
        // template <typename T, disable_arith_type_t<T> * = nullptr>
        // auto operator-(const T &rhs)
        // {
        //     return biops<expr_minus_t, biops<Ops, expr_result<T>, rExpr>, T>(expr_minus, *this, rhs);
        // }
        // template <typename T, enable_arith_type_t<T> * = nullptr>
        // auto operator-(const T &rhs)
        // {
        //     return biops<expr_minus_t, biops<Ops, expr_result<T>, rExpr>, expr_scalar<T>>(expr_minus, *this, rhs);
        // }
        // template <typename T, enable_arith_type_t<T> * = nullptr>
        // auto operator*(const T &rhs)
        // {
        //     return biops<expr_mul_t, biops<Ops, expr_result<T>, rExpr>, expr_scalar<T>>(expr_mul, *this, rhs);
        // }
        // template <typename T, enable_arith_type_t<T> * = nullptr>
        // auto operator/(const T &rhs)
        // {
        //     return biops<expr_div_t, biops<Ops, expr_result<T>, rExpr>, expr_scalar<T>>(expr_div, *this, rhs);
        // }

    private:
        const T &value;
    };

    template <typename Ops, typename lExpr, typename rExpr>
    class biops : public expr<biops<Ops, lExpr, rExpr>>
    {
    public:
        using base_type = expr<biops<Ops, lExpr, rExpr>>;
        using base_type::size;
        using base_type::operator[];
        friend base_type;

        explicit biops(const Ops &ops, const lExpr &lxpr, const rExpr &rxpr)
            : m_ops(ops), m_lxpr(lxpr), m_rxpr(rxpr){};
        constexpr size_t size_impl() { return gl_get_more(m_lxpr.size(), m_rxpr.size()); };
        auto at_impl(size_t idx) const { return m_ops(m_lxpr[idx], m_rxpr[idx]); };
        template <typename T>
        operator T()
        {
            T res{};
            for (size_t idx = 0; idx < res.size(); ++idx)
            {
                res[idx] = (*this)[idx];
            }
            return res;
        }
        template <typename T, disable_arith_type_t<T> * = nullptr, typename T::expr_type>
        auto operator+(const T &rhs)
        {
            return biops<expr_plus_t, biops<Ops, lExpr, rExpr>, T>(expr_plus, *this, rhs);
        }
        template <typename T, disable_arith_type_t<T> * = nullptr>
        auto operator+(const T &rhs)
        {
            using result_t = expr_result<T>;
            return biops<expr_plus_t, biops<Ops, lExpr, rExpr>, result_t>(expr_plus, *this, result_t(rhs));
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        auto operator+(const T &rhs)
        {
            using result_s = expr_result<expr_scalar<T>>;
            return biops<expr_plus_t, biops<Ops, lExpr, rExpr>, result_s>(expr_plus, *this, result_s(rhs));
        }
        template <typename T, disable_arith_type_t<T> * = nullptr, typename T::expr_type>
        auto operator-(const T &rhs)
        {
            return biops<expr_minus_t, biops<Ops, lExpr, rExpr>, T>(expr_minus, *this, rhs);
        }
        template <typename T, disable_arith_type_t<T> * = nullptr>
        auto operator-(const T &rhs)
        {
            using result_t = expr_result<T>;
            return biops<expr_minus_t, biops<Ops, lExpr, rExpr>, result_t>(expr_minus, *this, result_t(rhs));
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        auto operator-(const T &rhs)
        {
            using result_s = expr_result<expr_scalar<T>>;
            return biops<expr_minus_t, biops<Ops, lExpr, rExpr>, result_s>(expr_minus, *this, result_s(rhs));
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        auto operator*(const T &rhs)
        {
            using result_s = expr_result<expr_scalar<T>>;
            return biops<expr_mul_t, biops<Ops, lExpr, rExpr>, result_s>(expr_mul, *this, result_s(rhs));
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        auto operator/(const T &rhs)
        {
            using result_s = expr_result<expr_scalar<T>>;
            return biops<expr_div_t, biops<Ops, lExpr, rExpr>, result_s>(expr_div, *this, result_s(rhs));
        }

    private:
        Ops m_ops;
        lExpr m_lxpr;
        rExpr m_rxpr;
    };

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

    constexpr expr_plus_t expr_plus{};
    constexpr expr_minus_t expr_minus{};
    constexpr expr_mul_t expr_mul{};
    constexpr expr_div_t expr_div{};
}
#endif