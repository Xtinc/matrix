#ifndef VERY_SIMPLE_MATRIX_HEADER
#define VERY_SIMPLE_MATRIX_HEADER

#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>

constexpr double gl_rep_eps = std::numeric_limits<float>::epsilon();
constexpr size_t gl_get_less(size_t A, size_t B)
{
    return A < B ? A : B;
}
constexpr size_t gl_get_greater(size_t A, size_t B)
{
    return A < B ? B : A;
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
    for (auto ele : coll)
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

template <std::size_t M, std::size_t N>
class Matrix
{
    template <typename T, typename RT = void>
    using enable_arith_type_t = typename std::enable_if<std::is_arithmetic<T>::value, RT>::type;
    using IndexRange = std::pair<int, int>;

    struct SubPart
    {
    public:
        SubPart(Matrix<M, N> &self, IndexRange r, IndexRange c)
            : row_idx(0), row_end(0), col_idx(0), col_end(0), data(self)
        {
            row_idx = r.first >= 0 ? r.first : 0;
            col_idx = c.first >= 0 ? c.first : 0;
            row_end = r.second >= 0 ? r.second : M - 1;
            col_end = c.second >= 0 ? c.second : N - 1;
        }
        SubPart(const SubPart &) = delete;
        SubPart(SubPart &&) = delete;

        SubPart &operator=(const SubPart &) = delete;
        SubPart &operator=(SubPart &&) = delete;

        template <size_t A, size_t B>
        void operator=(const Matrix<A, B> &other)
        {
            if (row_idx > row_end || col_idx > col_end)
            {
                return;
            }
            auto real_row_counts = row_end - row_idx + 1;
            auto real_col_counts = col_end - col_idx + 1;
            real_row_counts = real_row_counts > A ? A : real_row_counts;
            real_col_counts = real_col_counts > B ? B : real_col_counts;
            for (size_t i = 0; i < real_row_counts; i++)
            {
                for (size_t j = 0; j < real_col_counts; j++)
                {
                    data(row_idx + i, col_idx + j) = other(i, j);
                }
            }
            // for (size_t j = 0; j < real_col_counts; j++)
            // {
            //     for (size_t i = 0; i < real_row_counts; i++)
            //     {
            //         data(row_idx + i, col_idx + j) = other(i, j);
            //     }
            // }
            return;
        }
        template <typename T, size_t A, enable_arith_type_t<T> * = nullptr>
        void operator=(const std::array<T, A> &list)
        {
            generator_by_list(list);
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        void operator=(const std::initializer_list<T> &list)
        {
            generator_by_list(list);
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        void operator=(const std::vector<T> &list)
        {
            generator_by_list(list);
        }
        std::array<double, M> operator*(const std::array<double, N> &x) const
        {
            return data * x;
        }
        template <size_t L>
        Matrix<M, L> operator*(const Matrix<N, L> &other) const
        {
            return data * other;
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
                for (size_t i = 0; i < real_row_counts; i++)
                {
                    data(row_idx + i, col_idx) = *(iter++);
                }
                return;
            }
            if (real_row_counts == 1)
            {
                real_col_counts = real_col_counts > A ? A : real_col_counts;
                for (size_t i = 0; i < real_col_counts; i++)
                {
                    data(row_idx, col_idx + i) = *(iter++);
                }
                return;
            }
        }
    };

public:
    // Constructors
    Matrix() : m_data(M * N, 0.0){};
    ~Matrix() = default;
    Matrix(const Matrix &) = default;
    Matrix(Matrix &&) = default;
    Matrix &operator=(const Matrix &) = default;
    Matrix &operator=(Matrix &&) = default;

    template <typename T, size_t L, enable_arith_type_t<T> * = nullptr>
    Matrix(const std::array<T, L> &list) : m_data(M * N, 0.0)
    {
        auto real_idx = list.size() < M * N ? list.size() : M * N;
        std::copy_n(list.begin(), real_idx, m_data.begin());
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    Matrix(const std::initializer_list<T> &list) : m_data(M * N, 0.0)
    {
        auto real_idx = list.size() < M * N ? list.size() : M * N;
        std::copy_n(list.begin(), real_idx, m_data.begin());
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    Matrix(const std::vector<T> &list) : m_data(M * N, 0.0)
    {
        auto real_idx = list.size() < M * N ? list.size() : M * N;
        std::copy_n(list.begin(), real_idx, m_data.begin());
    }

    // Member functions
    double *data()
    {
        return m_data.data();
    }
    const double *data() const
    {
        return m_data.data();
    }
    std::vector<double> &flat()
    {
        return m_data;
    }
    const std::vector<double> &flat() const
    {
        return m_data;
    }
    std::array<double, N> row(size_t idx) const
    {
        auto real_idx = idx < M ? idx : M;
        std::array<double, N> result;
        for (auto i = 0u; i < N; ++i)
        {
            result[i++] = m_data.at(real_idx + i * M);
        }
        return result;
    }
    std::array<double, M> col(size_t idx) const
    {
        auto real_idx = idx < N ? idx : N;
        std::array<double, M> result;
        for (auto i = 0u; i < M; ++i)
        {
            result[i++] = m_data.at(real_idx * M + i);
        }
        return result;
    }
    double &operator()(size_t row, size_t col)
    {
        return m_data.at(row + col * M);
    }
    const double &operator()(size_t row, size_t col) const
    {
        return m_data.at(row + col * M);
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
    void fill(double val)
    {
        std::fill(m_data.begin(), m_data.end(), val);
    }
    Matrix<N, M> T() const
    {
        std::vector<double> res(M * N, 0.0);
        for (auto i = 0u; i < M; i++)
        {
            for (auto j = 0u; j < N; j++)
            {
                res.at(i + j) = (*this)(i, j);
            }
        }
        return res;
    }

    // Overloaded Operators
    Matrix operator+(const Matrix &other) const
    {
        std::vector<double> res(M * N, 0.0);
        for (auto i = 0u; i < M * N; i++)
        {
            result[i] = m_data[i] + other.data()[i];
        }
        return result;
    }
    Matrix operator-(const Matrix &other) const
    {
        std::vector<double> res(M * N, 0.0);
        for (auto i = 0u; i < M * N; i++)
        {
            result[i] = m_data[i] - other.data()[i];
        }
        return result;
    }
    template <size_t L>
    Matrix<M, L> operator*(const Matrix<N, L> &other) const
    {
        Matrix<M, L> result{};
        Matrix<N, M> AT = this->T();
        auto iter_A = AT.data();
        auto iter_B = other.data();
        for (size_t i = 0; i < M; i++)
        {
            iter_B = other.data();
            for (size_t j = 0; j < L; j++)
            {
                result(i, j) = std::inner_product(iter_A, iter_A + N, iter_B, 0.0);
                iter_B += N;
            }
            iter_A += N;
        }
        return result;
        // Matrix<M, L> result{};
        // for (size_t i = 0; i < M; i++)
        // {
        //     for (size_t j = 0; j < L; j++)
        //     {
        //         for (size_t k = 0; k < N; k++)
        //         {
        //             result(i, j) += (*this)(i, k) * other(k, j);
        //         }
        //     }
        // }
        // return result;
    }
    std::array<double, M> operator*(const std::array<double, N> &x) const
    {
        std::vector<double> res(M * N, 0.0);
        for (size_t i = 0; i < M; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                result[i] = (*this)(i, j) * x[j];
            }
        }
        return result;
    }
    // template <size_t L>
    // Matrix<M, L> Mul2(const Matrix<N, L> &other) const
    // {
    //     Matrix<M, L> result{};
    //     Matrix<N, M> AT = this->T();
    //     auto iter_A = AT.data();
    //     auto iter_B = other.data();
    //     for (size_t i = 0; i < M; i++)
    //     {
    //         iter_B = other.data();
    //         for (size_t j = 0; j < L; j++)
    //         {
    //             double d = 0.0;
    //             for (size_t k = 0; k < N; k++)
    //             {
    //                 d += (*(iter_A + k)) * (*(iter_B + k));
    //             }
    //             result(i, j) = d;
    //             iter_B += N;
    //         }
    //         iter_A += N;
    //     }
    //     return result;
    // }
    // template <size_t L>
    // Matrix<M, L> Mul3(const Matrix<N, L> &other) const
    // {
    //     Matrix<M, L> result{};
    //     for (size_t i = 0; i < M; i++)
    //     {
    //         for (size_t j = 0; j < L; j++)
    //         {
    //             for (size_t k = 0; k < N; k++)
    //             {
    //                 result(i, j) += (*this)(i, k) * other(k, j);
    //             }
    //         }
    //     }
    //     return result;
    // }
    Matrix operator+=(const Matrix &other) const
    {
        for (size_t i = 0; i < M * N; i++)
        {
            m_data[i] += other.data()[i];
        }
    }
    Matrix operator-=(const Matrix &other) const
    {
        for (size_t i = 0; i < M * N; i++)
        {
            m_data[i] -= other.data()[i];
        }
    }
    bool operator==(const Matrix &other) const
    {
        return std::equal(m_data.begin(), m_data.end(), other.data(), [](double ele1, double ele2)
                          { return fabs(ele1 - ele2) < gl_rep_eps; });
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator+(T ele)
    {
        std::vector<double> res(M * N, 0.0);
        for (auto i = 0u; i < M * N; i++)
        {
            result[i] = m_data[i] + ele;
        }
        return result;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator-(T ele)
    {
        std::vector<double> res(M * N, 0.0);
        for (auto i = 0u; i < M * N; i++)
        {
            result[i] = m_data[i] - ele;
        }
        return result;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator*(T ele)
    {
        std::vector<double> res(M * N, 0.0);
        for (auto i = 0u; i < M * N; i++)
        {
            result[i] = m_data[i] * ele;
        }
        return result;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator/(T ele)
    {
        std::vector<double> res(M * N, 0.0);
        for (auto i = 0u; i < M * N; i++)
        {
            result[i] = m_data[i] / ele;
        }
        return result;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator+=(T ele)
    {
        for (auto &i : m_data)
        {
            i += ele;
        }
        return *this;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator-=(T ele)
    {
        for (auto &i : m_data)
        {
            i -= ele;
        }
        return *this;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator*=(T ele)
    {
        for (auto &i : m_data)
        {
            i *= ele;
        }
        return *this;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator/=(T ele)
    {
        for (auto &i : m_data)
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
    static size_t row_counts()
    {
        return M;
    }
    static size_t col_counts()
    {
        return N;
    }

private:
    // data member
    std::vector<double> m_data;
};

template <size_t M, size_t N>
void zeros(Matrix<M, N> &m)
{
    m.flat() = {};
}

template <size_t M, size_t N>
void ones(Matrix<M, N> &m)
{
    m.fill(1);
}

template <size_t M, size_t N, size_t A, size_t B>
Matrix<gl_get_greater(M, A), N + B> catcol(const Matrix<M, N> &m1, const Matrix<A, B> &m2)
{
    constexpr size_t N_M = gl_get_greater(M, A);
    Matrix<N_M, N + B> result{};
    for (size_t j = 0; j < N; j++)
    {
        for (size_t i = 0; i < M; i++)
        {
            result(i, j) = m1(i, j);
        }
    }
    for (size_t j = N, idx = 0; j < N + B; ++j, ++idx)
    {
        for (size_t i = 0; i < A; ++i)
        {
            result(i, j) = m2(i, idx);
        }
    }
    return result;
}

template <size_t M, size_t N, size_t A, size_t B>
Matrix<M + A, gl_get_greater(N, B)> catrow(const Matrix<M, N> &m1, const Matrix<A, B> &m2)
{
    constexpr size_t N_N = gl_get_greater(N, B);
    Matrix<M + A, N_N> result{};
    for (size_t j = 0; j < N; j++)
    {
        for (size_t i = 0; i < M; i++)
        {
            result(i, j) = m1(i, j);
        }
    }
    for (size_t j = 0; j < B; ++j)
    {
        for (size_t i = A, idx = 0; i < M + A; ++i, ++idx)
        {
            result(i, j) = m2(idx, j);
        }
    }
    return result;
}

template <size_t M, size_t N>
Matrix<M - 1u, N - 1u> cofactor(const Matrix<M, N> &mat, size_t p, size_t q)
{
    Matrix<M - 1u, N - 1u> result{};
    size_t i = 0, j = 0;
    for (size_t row = 0; row < M; row++)
    {
        for (size_t col = 0; col < N; col++)
        {
            if (row != p && col != q)
            {
                result(i, j++) = mat(row, col);
                if (j == N - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return result;
}

// template <>
// double determinant(const Matrix<1, 1, nullptr> &mat)
// {
//     return mat.flat()[0];
// }

// template <size_t M>
// Matrix<M, M> adjoint(const Matrix<M, M> &mat)
// {
//     Matrix<M, M> result{};
//     for (int i = 0; i < M; i++)
//     {
//         for (int j = 0; j < M; j++)
//         {
//             auto sign = (i + j) % 2 == 0 ? 1 : -1;
//             result(j, i) = sign * (determinant(cofactor(mat, i, j)));
//         }
//     }
//     return result;
// }

// template <>
// Matrix<1, 1> adjoint(const Matrix<1, 1, nullptr> &mat)
// {
//     return {1};
// }

template <size_t N>
Matrix<N, N> ludcmp(Matrix<N, N> A, std::array<int, N> &indx, bool &even)
{
    even = true;
    for (int i = 0; i < N; i++)
    {
        indx[i] = i;
    }
    for (int k = 0; k < N - 1; k++)
    {
        auto valmax = fabs(A(k, k));
        auto ip = k;
        for (int row = k + 1; row < N; row++)
        {
            double tmp = fabs(A(row, k));
            if (valmax < tmp)
            {
                valmax = tmp;
                ip = row;
            }
        }
        if (valmax < gl_rep_eps)
        {
            return {};
        }
        if (ip != k)
        {
            for (int col = k; col < N; col++)
            {
                std::swap(A(ip, col), A(k, col));
            }
            std::swap(indx[ip], indx[k]);
            for (int col = 0; col < k; col++)
            {
                std::swap(A(k, col), A(ip, col));
            }
            even = !even;
        }
        for (int row = k + 1; row < N; row++)
        {
            double weight = A(row, k) / A(k, k);
            A(row, k) = weight;
            for (int col = k + 1; col < N; col++)
            {
                A(row, col) -= weight * A(k, col);
            }
        }
    }
    return A;
}

template <size_t N>
void ludbksb(const Matrix<N, N> &A, const std::array<int, N> &indx, double *b)
{
    std::array<double, N> y{};
    y[0] = b[indx[0]];
    for (int row = 1; row < N; row++)
    {
        double sum = 0.0;
        for (int col = 0; col < row; col++)
        {
            sum += A(row, col) * y[col];
        }
        y[row] = b[indx[row]] - sum;
    }

    int n = N;

    b[n - 1] = y[n - 1] / A(n - 1, n - 1);
    for (int row = n - 2; row >= 0; row--)
    {
        auto id = row + 1;
        double sum = 0.0;
        for (int col = id; col < n; col++)
        {
            sum += A(row, col) * b[col];
        }
        b[row] = (y[row] - sum) / A(row, row);
    }
}

template <size_t M>
double determinant(const Matrix<M, M> &mat)
{
    auto even = true;
    std::array<int, M> indx{};
    auto LU = ludcmp(mat, indx, even);
    auto D = even ? 1.0 : -1.0;
    for (size_t i = 0; i < M; i++)
    {
        D *= LU(i, i);
    }
    return D;
}

template <size_t M>
Matrix<M, M> inverse(const Matrix<M, M> &mat)
{
    std::array<int, M> indx{};
    auto even = true;
    auto result = Matrix<M, M>::eye();
    auto LU = ludcmp(mat, indx, even);
    for (int j = 0; j < M; j++)
    {
        ludbksb(LU, indx, result.data() + j * M);
    }
    return result;
}

template <size_t A, size_t B, size_t M, size_t N>
Matrix<A, B> slice(const Matrix<M, N> &m, size_t row_start, size_t col_start)
{
    if (row_start + A > M || col_start + B > N)
    {
        return {};
    }
    Matrix<A, B> res{};
    for (size_t i = row_start, ix = 0u; i < row_start + A; ++i, ++ix)
    {
        for (size_t j = col_start, iy = 0u; j < col_start + B; ++j, ++iy)
        {
            res(ix, iy) = m(i, j);
        }
    }
    return res;
}

template <size_t M, size_t N>
Matrix<N, M> transpose(const Matrix<M, N> &m)
{
    return m.T();
}

#endif
