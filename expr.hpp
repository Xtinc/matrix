#ifndef VVERY_SIMPLE_EXPRESSION_HEADER
#define VVERY_SIMPLE_EXPRESSION_HEADER

namespace details
{
    template <typename LHS, typename RHS>
    class MatrixSum
    {
    public:
        using value_type = typename LHS::value_type;
        MatrixSum(const LHS &lhs, const RHS &rhs) : rhs(rhs), lhs(lhs)
        {
        }
        value_type operator[](size_t idx) const
        {
            return lhs[idx] + rhs[idx];
        }
        template <typename T>
        MatrixSum<MatrixSum<LHS, RHS>, T> operator+(const T &rhs)
        {
            return MatrixSum<MatrixSum<LHS, RHS>, T>(*this, rhs);
        }

    private:
        const LHS &lhs;
        const RHS &rhs;
    };
}
#endif