#ifndef VVERY_SIMPLE_MATRIX_HEADER
#define VVERY_SIMPLE_MATRIX_HEADER

#include "expr.hpp"
#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>

// forward declare
template <size_t M, size_t N>
class Matrix;

template <size_t N>
Matrix<N, N> ludcmp(Matrix<N, N> A, std::array<int, N> &indx, bool &even);

template <size_t N>
void ludbksb(const Matrix<N, N> &A, const std::array<int, N> &indx, double *b);

namespace details
{
    constexpr size_t gl_sm(size_t A, size_t B)
    {
        return A * B < 260 ? 1 : 0;
    }

    template <size_t M, size_t N, std::size_t A = gl_sm(M, N)>
    class MatrixBase;

    template <std::size_t M, std::size_t N>
    class MatrixBase<M, N, 1>
    {
    protected:
        MatrixBase() : m_data{} {}
        ~MatrixBase() = default;
        MatrixBase(const MatrixBase &) = default;
        MatrixBase(MatrixBase &&) = default;
        MatrixBase &operator=(const MatrixBase &) = default;
        MatrixBase &operator=(MatrixBase &&) = default;

        template <typename T, size_t L, enable_arith_type_t<T> * = nullptr>
        MatrixBase(const std::array<T, L> &list) : m_data{}
        {
            auto real_idx = list.size() < M * N ? list.size() : M * N;
            std::copy_n(list.begin(), real_idx, m_data.begin());
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        MatrixBase(const std::initializer_list<T> &list) : m_data{}
        {
            auto real_idx = list.size() < M * N ? list.size() : M * N;
            std::copy_n(list.begin(), real_idx, m_data.begin());
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
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
        MatrixBase() : m_data(M * N, 0.0) {}
        ~MatrixBase() = default;
        MatrixBase(const MatrixBase &) = default;
        MatrixBase(MatrixBase &&) = default;
        MatrixBase &operator=(const MatrixBase &) = default;
        MatrixBase &operator=(MatrixBase &&) = default;

        template <typename T, size_t L, enable_arith_type_t<T> * = nullptr>
        MatrixBase(const std::array<T, L> &list) : m_data(M * N, 0.0)
        {
            auto real_idx = list.size() < M * N ? list.size() : M * N;
            std::copy_n(list.begin(), real_idx, m_data.begin());
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        MatrixBase(const std::initializer_list<T> &list) : m_data(M * N, 0.0)
        {
            auto real_idx = list.size() < M * N ? list.size() : M * N;
            std::copy_n(list.begin(), real_idx, m_data.begin());
        }
        template <typename T, enable_arith_type_t<T> * = nullptr>
        MatrixBase(const std::vector<T> &list) : m_data(M * N, 0.0)
        {
            auto real_idx = list.size() < M * N ? list.size() : M * N;
            std::copy_n(list.begin(), real_idx, m_data.begin());
        }

    protected:
        std::vector<double> m_data;
    };
}

template <std::size_t M, std::size_t N>
class Matrix : public details::MatrixBase<M, N>
{
    template <size_t A, typename RT = void>
    using enable_when_squre_t = typename std::enable_if<A == N, RT>::type;
    using IndexRange = std::pair<int, int>;
    using result_t = details::expr_result<Matrix>;
    template <typename T>
    using result_s = details::expr_result<details::expr_scalar<T>>;

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
    struct iterator : public std::iterator<std::random_access_iterator_tag, double>
    {
    public:
        using iterator_category = std::random_access_iterator_tag;
        inline iterator(value_type *ptr) noexcept : m_ptr(ptr) {}
        inline iterator(const iterator &itr) noexcept : m_ptr(itr.m_ptr) {}
        inline iterator &operator=(const iterator &rhs) noexcept
        {
            m_ptr = rhs.m_ptr;
            return *this;
        }
        inline bool operator==(const iterator &rhs) const noexcept
        {
            return m_ptr == rhs.m_ptr;
        }
        inline bool operator!=(const iterator &rhs) const noexcept
        {
            return m_ptr != rhs.m_ptr;
        }
        inline bool operator>(const iterator &rhs) const noexcept
        {
            return m_ptr > rhs.m_ptr;
        }
        inline bool operator<=(const iterator &rhs) const noexcept
        {
            return m_ptr <= rhs.m_ptr;
        }
        inline pointer operator->() const noexcept
        {
            return m_ptr;
        }
        inline reference operator*() const noexcept
        {
            return *m_ptr;
        }
        inline iterator &operator++() noexcept
        {
            m_ptr += 1;
            return *this;
        }
        inline iterator operator++(int) noexcept
        {
            value_type *ret = m_ptr;
            m_ptr += 1;
            return ret;
        }
        inline iterator &operator--() noexcept
        {
            m_ptr -= 1;
            return *this;
        }
        inline iterator operator--(int) noexcept
        {
            value_type *ret = m_ptr;
            m_ptr -= 1;
            return ret;
        }
        inline iterator &operator+=(difference_type step) noexcept
        {
            m_ptr += step;
            return *this;
        }
        inline iterator &operator-=(difference_type step) noexcept
        {
            m_ptr -= step;
            return *this;
        }
        inline iterator operator+(difference_type step) noexcept
        {
            value_type *ret = m_ptr;
            ret += step;
            return ret;
        }
        inline iterator operator-(difference_type step) noexcept
        {
            value_type *ret = m_ptr;
            ret -= step;
            return ret;
        }
        inline iterator operator[](difference_type n) noexcept
        {
            m_ptr += n;
            return *this;
        }
        inline difference_type operator-(const iterator &rhs) const noexcept
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
        inline const_iterator(value_type *const ptr) noexcept : m_ptr(ptr)
        {
        }
        inline const_iterator(const const_iterator &itr) noexcept : m_ptr(itr.m_ptr)
        {
        }
        inline const_iterator &operator=(const const_iterator &rhs) noexcept
        {
            m_ptr = rhs.m_ptr;
            return *this;
        }
        inline bool operator==(const const_iterator &rhs) const noexcept
        {
            return m_ptr == rhs.m_ptr;
        }
        inline bool operator!=(const const_iterator &rhs) const noexcept
        {
            return m_ptr != rhs.m_ptr;
        }
        inline bool operator>(const const_iterator &rhs) const noexcept
        {
            return m_ptr > rhs.m_ptr;
        }
        inline bool operator<=(const const_iterator &rhs) const noexcept
        {
            return m_ptr <= rhs.m_ptr;
        }
        inline const pointer operator->() const noexcept
        {
            return m_ptr;
        }
        inline const reference operator*() const noexcept
        {
            return *m_ptr;
        }
        inline const_iterator &operator++() noexcept
        {
            m_ptr++;
            return *this;
        }
        inline const_iterator operator++(int) noexcept
        {
            value_type const *ret = m_ptr;
            m_ptr++;
            return ret;
        }
        inline const_iterator &operator--() noexcept
        {
            m_ptr--;
            return *this;
        }
        inline const_iterator operator--(int) noexcept
        {
            value_type const *ret = m_ptr;
            m_ptr--;
            return ret;
        }
        inline const_iterator &operator+=(ptrdiff_t step) noexcept
        {
            m_ptr += step;
            return *this;
        }
        inline const_iterator &operator-=(ptrdiff_t step) noexcept
        {
            m_ptr -= step;
            return *this;
        }
        inline const_iterator operator+(ptrdiff_t step) noexcept
        {
            value_type const *ret = m_ptr;
            ret += step;
            return ret;
        }
        inline const_iterator operator-(ptrdiff_t step) noexcept
        {
            value_type const *ret = m_ptr;
            ret -= step;
            return ret;
        }
        inline const_iterator operator[](difference_type n) noexcept
        {
            m_ptr += n;
            return *this;
        }
        inline difference_type operator-(const const_iterator &rhs) const noexcept
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
    ~Matrix() = default;
    Matrix(const Matrix &) = default;
    Matrix(Matrix &&) = default;
    Matrix &operator=(const Matrix &) = default;
    Matrix &operator=(Matrix &&) = default;

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
    std::array<double, N> row(size_t idx) const
    {
        auto real_idx = idx < M ? idx : M;
        std::array<double, N> result;
        for (auto i = 0u; i < N; ++i)
        {
            result[i++] = this->m_data.at(real_idx + i * M);
        }
        return result;
    }
    std::array<double, M> col(size_t idx) const
    {
        auto real_idx = idx < N ? idx : N;
        std::array<double, M> result;
        for (auto i = 0u; i < M; ++i)
        {
            result[i++] = this->m_data.at(real_idx * M + i);
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
        auto result = Matrix<M, M>::eye();
        auto LU = ludcmp(*this, indx, even);
        for (int j = 0; j < M; j++)
        {
            ludbksb(LU, indx, result.data() + j * M);
        }
        return result;
    }
    template <size_t A = M>
    enable_when_squre_t<A, double> det() const
    {
        auto even = true;
        std::array<int, M> indx{};
        auto LU = ludcmp(*this, indx, even);
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

    // Overloaded Operators
    double &operator()(size_t row, size_t col)
    {
        return this->m_data.at(row + col * M);
    }
    const double &operator()(size_t row, size_t col) const
    {
        return this->m_data.at(row + col * M);
    }
    double &operator[](size_t idx)
    {
        return this->m_data.at(idx);
    }
    const double &operator[](size_t idx) const
    {
        return this->m_data.at(idx);
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    auto operator+(const T &other) const
    {
        return details::biops<details::expr_plus_t, result_t, result_s<T>>(details::expr_plus, result_t(*this), result_s<T>(other));
    }
    template <typename T, disable_arith_type_t<T> * = nullptr>
    auto operator+(const T &other) const
    {
        return details::biops<details::expr_plus_t, result_t, T>(details::expr_plus, result_t(*this), other);
    }
    auto operator+(const Matrix &other) const
    {
        return details::biops<details::expr_plus_t, result_t, result_t>(details::expr_plus, result_t(*this), result_t(other));
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    auto operator-(const T &other) const
    {
        return details::biops<details::expr_minus_t, result_t, result_s<T>>(details::expr_minus, result_t(*this), result_s<T>(other));
    }
    template <typename T, disable_arith_type_t<T> * = nullptr>
    auto operator-(const T &other) const
    {
        return details::biops<details::expr_minus_t, result_t, T>(details::expr_minus, result_t(*this), other);
    }
    auto operator-(const Matrix &other) const
    {
        return details::biops<details::expr_minus_t, result_t, result_t>(details::expr_minus, result_t(*this), result_t(other));
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    auto operator*(const T &t) const
    {
        return details::biops<details::expr_mul_t, result_t, result_s<T>>(details::expr_mul, result_t(*this), result_s<T>(t));
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    auto operator/(const T &t) const
    {
        return details::biops<details::expr_div_t, result_t, result_s<T>>(details::expr_div, result_t(*this), result_s<T>(t));
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
    }
    std::array<double, M> operator*(const std::array<double, N> &x) const
    {
        std::vector<double> result(M * N, 0.0);
        for (size_t i = 0; i < M; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                result[i] = (*this)(i, j) * x[j];
            }
        }
        return result;
    }
    Matrix operator+=(const Matrix &other) const
    {
        for (size_t i = 0; i < M * N; i++)
        {
            this->m_data[i] += other.data()[i];
        }
    }
    Matrix operator-=(const Matrix &other) const
    {
        for (size_t i = 0; i < M * N; i++)
        {
            this->m_data[i] -= other.data()[i];
        }
    }
    bool operator==(const Matrix &other) const
    {
        return std::equal(this->m_data.begin(), this->m_data.end(), other.data(), [](double ele1, double ele2)
                          { return details::is_same(ele1, ele2); });
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator+=(T ele)
    {
        for (auto &i : this->m_data)
        {
            i += ele;
        }
        return *this;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator-=(T ele)
    {
        for (auto &i : this->m_data)
        {
            i -= ele;
        }
        return *this;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator*=(T ele)
    {
        for (auto &i : this->m_data)
        {
            i *= ele;
        }
        return *this;
    }
    template <typename T>
    enable_arith_type_t<T, Matrix<M, N>> operator/=(T ele)
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

    template <typename T, enable_arith_type_t<T> * = nullptr>
    friend auto operator+(const T &t, const Matrix<M, N> &self)
    {
        return details::biops<details::expr_plus_t, result_t, result_s<T>>(details::expr_plus, result_t(self), result_s<T>(t));
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    friend auto operator-(const T &t, const Matrix<M, N> &self)
    {
        return details::biops<details::expr_minus_t, result_t, result_s<T>>(details::expr_minus, result_t(self), result_s<T>(t));
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    friend auto operator*(const T &t, const Matrix<M, N> &self)
    {
        return details::biops<details::expr_mul_t, result_t, result_s<T>>(details::expr_mul, result_t(self), result_s<T>(t));
    }
    template <typename T, enable_arith_type_t<T> * = nullptr>
    friend auto operator/(const T &t, const Matrix<M, N> &self)
    {
        return details::biops<details::expr_div_t, result_t, result_s<T>>(details::expr_div, result_t(self), result_s<T>(t));
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

template <size_t M, size_t N>
void zeros(Matrix<M, N> &m)
{
    m.fill(0.0);
}

template <size_t M, size_t N>
void ones(Matrix<M, N> &m)
{
    m.fill(1.0);
}

template <size_t M, size_t N, size_t A, size_t B>
Matrix<gl_get_more(M, A), N + B> catcol(const Matrix<M, N> &m1, const Matrix<A, B> &m2)
{
    constexpr size_t N_M = gl_get_more(M, A);
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
Matrix<M + A, gl_get_more(N, B)> catrow(const Matrix<M, N> &m1, const Matrix<A, B> &m2)
{
    constexpr size_t N_N = gl_get_more(N, B);
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

template <size_t M>
Matrix<M, M> adjugate(const Matrix<M, M> &mat)
{
    Matrix<M, M> result{};
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            auto sign = (i + j) % 2 == 0 ? 1 : -1;
            result(j, i) = sign * (determinant(cofactor(mat, i, j)));
        }
    }
    return result;
}

template <>
Matrix<1, 1> adjugate(const Matrix<1, 1> &mat)
{
    return {1};
}

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
    return mat.det();
}

template <size_t M>
Matrix<M, M> inverse(const Matrix<M, M> &mat)
{
    return mat.I();
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

template <size_t M, size_t N>
enable_when_array_t<M, N, double> norm2(const Matrix<M, N> &mat)
{
    double res = 0.0;
    for (auto ele : mat)
    {
        res += ele * ele;
    }
    return sqrt(res);
}

template <size_t M, size_t N>
double trace(const Matrix<M, N> &mat)
{
    return mat.trace();
}

#endif