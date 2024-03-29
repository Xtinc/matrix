#ifndef VVERY_SIMPLE_LINEAR_ALGORITHM_HEADER
#define VVERY_SIMPLE_LINEAR_ALGORITHM_HEADER

#include "matrixs.hpp"
#include <complex>

namespace ppx
{
    template <typename T>
    constexpr int SIGN(T a)
    {
        if (a > T{})
        {
            return 1;
        }
        else if (a < T{})
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }

    template <typename T>
    constexpr T SIGN(T a, T b)
    {
        return b >= T{} ? std::abs(a) : -std::abs(a);
    }

    template <typename T>
    constexpr T SQR(T a)
    {
        return a * a;
    }

    template <typename T>
    constexpr T CUB(T a)
    {
        return a * a * a;
    }

    template <class T, class Compare>
    constexpr const T &CLAMP(const T &v, const T &lo, const T &hi, Compare comp)
    {
        return comp(v, lo) ? lo : comp(hi, v) ? hi
                                              : v;
    }

    template <class T>
    constexpr const T &CLAMP(const T &v, const T &lo, const T &hi)
    {
        return CLAMP(v, lo, hi, std::less<T>{});
    }

    template <class T>
    constexpr const T &CLAMP(const T &v, const std::pair<T, T> &ra)
    {
        return CLAMP(v, ra.first, ra.second, std::less<T>{});
    }

    template <typename T, std::enable_if_t<std::is_floating_point<T>::value, void> * = nullptr>
    std::vector<T> LINSPACE(T lo, T hi, size_t num)
    {
        if (num < 2)
        {
            // throw new exception ?;
            return {};
        }
        auto partitions = num - 1;
        std::vector<T> pts(num, 0.0);
        T length = (hi - lo) / partitions;
        pts[0] = lo;
        for (size_t i = 1; i < partitions; i++)
        {
            pts[i] = lo + i * length;
        }
        pts[num - 1] = hi;
        return pts;
    }

    template <typename T>
    std::vector<T> LINSPACE(std::pair<T, T> range, size_t num)
    {
        return LINSPACE(range.first, range.second, num);
    }

    // matrix related

    inline double sum(const double *cbegin, size_t len)
    {
#if defined(PPX_USE_AVX)
        return avxt::sum(cbegin, len);
#else
        return std::accumulate(cbegin, cbegin + len, 0.0);
#endif
    }

    template <size_t M, size_t N>
    void zeros(MatrixS<M, N> &m)
    {
        m.fill(0.0);
    }

    template <size_t M, size_t N>
    void ones(MatrixS<M, N> &m)
    {
        m.fill(1.0);
    }

    template <size_t N>
    MatrixS<N, N> eye()
    {
        MatrixS<N, N> result;
        for (size_t i = 0; i < N; i++)
        {
            result(i, i) = 1;
        }
        return result;
    }

    template <size_t M, size_t N>
    MatrixS<M - 1u, N - 1u> cofactor(const MatrixS<M, N> &mat, size_t p, size_t q)
    {
        MatrixS<M - 1u, N - 1u> result{};
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
    MatrixS<M, M> adjugate(const MatrixS<M, M> &mat)
    {
        MatrixS<M, M> result{};
        for (size_t i = 0; i < M; i++)
        {
            for (size_t j = 0; j < M; j++)
            {
                auto sign = (i + j) % 2 == 0 ? 1 : -1;
                result(j, i) = sign * (determinant(cofactor(mat, i, j)));
            }
        }
        return result;
    }

    template <>
    inline MatrixS<1, 1> adjugate(const MatrixS<1, 1> &)
    {
        return {1};
    }

    template <size_t M>
    double determinant(const MatrixS<M, M> &mat)
    {
        return mat.det();
    }

    template <size_t M>
    MatrixS<M, M> inverse(const MatrixS<M, M> &mat)
    {
        return mat.I();
    }

    template <size_t M, size_t N>
    MatrixS<N, M> transpose(const MatrixS<M, N> &m)
    {
        return m.T();
    }

    template <size_t M, size_t N>
    double norm1(const MatrixS<M, N> &mat)
    {
        auto res = 0.0;
        for (size_t j = 0; j < N; j++)
        {
            auto norm1_colj = 0.0;
            for (size_t i = 0; i < M; i++)
            {
                norm1_colj += std::abs(mat(i, j));
            }
            res = std::max(res, norm1_colj);
        }
        return res;
    }

    template <size_t M, size_t N>
    double norminf(const MatrixS<M, N> &mat)
    {
        auto res = 0.0;
        for (size_t i = 0; i < M; i++)
        {
            auto norm_inf_rowi = 0.0;
            for (size_t j = 0; j < N; j++)
            {
                norm_inf_rowi += std::abs(mat(i, j));
            }
            res = std::max(res, norm_inf_rowi);
        }
        return res;
    }

    template <size_t M, size_t N, enable_when_array_t<M, N> * = nullptr>
    double norm2(const MatrixS<M, N> &mat)
    {
#if defined(PPX_USE_AVX)
        return std::sqrt(avxt::inrpdt(mat.data(), mat.data(), M * N));
#else
        return std::sqrt(inner_product(mat, mat));
#endif
    }

    template <size_t M, size_t N>
    double trace(const MatrixS<M, N> &mat)
    {
        return mat.trace();
    }

    template <size_t N>
    double inner_product(const MatrixS<N, 1> &a, const MatrixS<N, 1> &b)
    {
#if defined(PPX_USE_AVX)
        return avxt::inrpdt(a.data(), b.data(), N);
#else
        return std::inner_product(a.cbegin(), a.cend(), b.cbegin(), 0.0);
#endif
    }

    template <size_t M, size_t N>
    double inner_product(const MatrixS<M, 1> &a, const MatrixS<N, 1> &b, const MatrixS<M, N> &Q)
    {
        auto result = 0.0;
        for (size_t i = 0; i < N; i++)
        {
            auto sum = 0.0;
            for (size_t j = 0; j < M; j++)
            {
                sum += a[j] * Q(i, j);
            }
            result += b[i] * sum;
        }
        return result;
    }

    template <size_t M, size_t N>
    std::pair<size_t, size_t> maxloc(const MatrixS<M, N> &mat)
    {
        auto max_pos = std::max_element(mat.begin(), mat.end());
        auto max_dis = std::div(std::distance(mat.begin(), max_pos), M);
        return {max_dis.rem, max_dis.quot};
    }

    template <Ori ori, size_t M, size_t N>
    std::enable_if_t<ori == Ori::Row, MatrixS<2 * M, N>> concat(const MatrixS<M, N> &a,
                                                                const MatrixS<M, N> &b)
    {
        MatrixS<2 * M, N> result;
        for (size_t i = 0; i < 2 * M; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                result(i, j) = i < M ? a(i, j) : b(i - M, j);
            }
        }
        return result;
    }

    template <Ori ori, size_t M, size_t N>
    std::enable_if_t<ori == Ori::Col, MatrixS<M, 2 * N>> concat(const MatrixS<M, N> &a,
                                                                const MatrixS<M, N> &b)
    {
        MatrixS<M, 2 * N> result;
        for (size_t i = 0; i < M; ++i)
        {
            for (size_t j = 0; j < 2 * N; ++j)
            {
                result(i, j) = j < N ? a(i, j) : b(i, j - N);
            }
        }
        return result;
    }

    template <size_t M, size_t N, size_t L>
    MatrixS<M + N, L> concat(const MatrixS<M, L> &a, const MatrixS<N, L> &b)
    {
        MatrixS<M + N, L> result;
        for (size_t i = 0; i < M + N; ++i)
        {
            for (size_t j = 0; j < L; ++j)
            {
                result(i, j) = i < M ? a(i, j) : b(i - M, j);
            }
        }
        return result;
    }

    template <size_t M, size_t N, size_t L>
    MatrixS<L, M + N> concat(const MatrixS<L, M> &a, const MatrixS<L, N> &b)
    {
        MatrixS<L, M + N> result;
        for (size_t i = 0; i < L; ++i)
        {
            for (size_t j = 0; j < M + N; ++j)
            {
                result(i, j) = j < M ? a(i, j) : b(i, j - M);
            }
        }
        return result;
    }

    template <size_t M>
    MatrixS<M + 1, 1> concat(const MatrixS<M, 1> &a, double b)
    {
        MatrixS<M + 1, 1> result;
        std::copy_n(a.begin(), M, result.begin());
        result[M] = b;
        return result;
    }

    template <size_t M>
    MatrixS<1, M + 1> concat(const MatrixS<1, M> &a, double b)
    {
        MatrixS<1, M + 1> result;
        std::copy_n(a.begin(), M, result.begin());
        result[M] = b;
        return result;
    }

    template <size_t M>
    MatrixS<M + 1, 1> concat(double a, const MatrixS<M, 1> &b)
    {
        MatrixS<M + 1, 1> result;
        std::copy_n(b.begin(), M, result.begin() + 1);
        result[0] = a;
        return result;
    }

    template <size_t M>
    MatrixS<1, M + 1> concat(double a, const MatrixS<1, M> &b)
    {
        MatrixS<1, M + 1> result;
        std::copy_n(b.begin(), M, result.begin() + 1);
        result[0] = a;
        return result;
    }

    // solver linear system
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

    template <size_t N>
    struct EqnResult
    {
        MatrixS<N, 1> x;
        StatusCode s = StatusCode::NORMAL;

        friend std::ostream &operator<<(std::ostream &os, const EqnResult &self)
        {
            os << "EqnResult<" << N << ">:\n"
               << "Status:\t" << self.s << "\n"
               << "x     =\t" << self.x << std::endl;
            return os;
        }
    };

    template <size_t M, size_t N>
    MatrixS<M, N> MGS(MatrixS<M, N> A)
    {
        for (int j = 0; j < (int)N; j++)
        {
            auto pivot = 0.0;
            for (size_t i = 0; i < M; i++)
            {
                pivot = std::max(pivot, std::abs(A(i, j)));
            }
            if (pivot < EPS_DP)
            {
                continue;
            }

            for (int k = 0; k < j; k++)
            {
                auto sum = 0.0;
                for (size_t i = 0; i < N; i++)
                {
                    sum += A(i, j) * A(i, k);
                }

                for (size_t i = 0; i < M; i++)
                {
                    A(i, j) -= sum * A(i, k);
                }
            }
            auto norm_j = 0.0;
            for (size_t i = 0; i < N; i++)
            {
                norm_j += A(i, j) * A(i, j);
            }
            norm_j = sqrt(norm_j);
            for (size_t i = 0; i < N; i++)
            {
                A(i, j) /= norm_j;
            }
        }
        return A;
    }

    template <int n>
    class LU
    {
    public:
        LU(const MatrixS<n, n> &mat)
            : A(mat), even(true), s(StatusCode::NORMAL)
        {
            static_assert(n > 0, "dimension must be positive!");
            for (int i = 0; i < n; i++)
            {
                indx[i] = i;
            }
            for (int k = 0; k < n; k++)
            {
                auto valmax = fabs(A(k, k));
                auto ip = k;
                for (int row = k + 1; row < n; row++)
                {
                    double tmp = fabs(A(row, k));
                    if (valmax < tmp)
                    {
                        valmax = tmp;
                        ip = row;
                    }
                }
                if (valmax < EPS_DP)
                {
                    s = StatusCode::SINGULAR;
                    break;
                }
                if (ip != k)
                {
                    for (int col = k; col < n; col++)
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
                for (int row = k + 1; row < n; row++)
                {
                    double weight = A(row, k) / A(k, k);
                    A(row, k) = weight;
                    for (int col = k + 1; col < n; col++)
                    {
                        A(row, col) -= weight * A(k, col);
                    }
                }
            }
        }

        void solve(double *b) const
        {
            std::array<double, n> y{};
            y[0] = b[indx[0]];
            for (int row = 1; row < n; row++)
            {
                double sum = 0.0;
                for (int col = 0; col < row; col++)
                {
                    sum += A(row, col) * y[col];
                }
                y[row] = b[indx[row]] - sum;
            }
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

        MatrixS<n, n> A;
        std::array<int, n> indx{};
        bool even;
        StatusCode s;
    };

    template <int m, int n>
    class QR
    {
    public:
        QR(const MatrixS<m, n> &mat) : A(mat), s(StatusCode::NORMAL)
        {
            static_assert(m > 0 && n > 0, "dimension must be a positive number!");
            for (int k = 0; k < std::min(m, n); k++)
            {
                auto sum = 0.0;
                for (int i = k; i < m; i++)
                {
                    sum += A(i, k) * A(i, k);
                }
                if (fabs(sum) < EPS_DP)
                {
                    s = StatusCode::SINGULAR;
                }
                auto sigma = SIGN(sqrt(sum), A(k, k));
                A(k, k) += sigma;
                c[k] = sigma * A(k, k);
                d[k] = -sigma;
                for (int j = k + 1; j < n; j++)
                {
                    auto tau = 0.0;
                    for (int i = k; i < m; i++)
                    {
                        tau += A(i, k) * A(i, j);
                    }
                    tau = fabs(c[k]) > EPS_DP ? tau / c[k] : 0.0;
                    for (int i = k; i < m; i++)
                    {
                        A(i, j) -= tau * A(i, k);
                    }
                }
            }
        }

        void solve(double *b) const
        {
            for (int j = 0; j < n; j++)
            {
                auto sum = 0.0;
                for (int i = j; i < m; i++)
                {
                    sum += A(i, j) * b[i];
                }
                sum /= c[j];
                for (int i = j; i < m; i++)
                {
                    b[i] -= sum * A(i, j);
                }
            }
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++)
                {
                    sum += A(i, j) * b[j];
                }
                b[i] = (b[i] - sum) / d[i];
            }
        }

        MatrixS<m, n> A;
        MatrixS<n, 1> c;
        MatrixS<n, 1> d;
        StatusCode s;
    };

    template <int m, int n>
    class SVD
    {
    public:
        SVD(const MatrixS<m, n> &mat) : u(mat)
        {
            static_assert(m > 0 && n > 0, "dimension must be a positive number!");
            int i, its, j, jj, k, nm;
            double c, f, h, s, x, y, z;
            std::array<double, n> rv1{};
            double g = 0.0;
            double scale = 0.0;
            double anorm = 0.0;
            int l = 0;
            for (i = 0; i < n; i++)
            {
                l = i + 2;
                rv1[i] = scale * g;
                g = 0.0;
                s = 0.0;
                scale = 0.0;
                if (i < m)
                {
                    for (k = i; k < m; k++)
                    {
                        scale += fabs(u(k, i));
                    }
                    if (scale != 0.0)
                    {
                        for (k = i; k < m; k++)
                        {
                            u(k, i) /= scale;
                            s += u(k, i) * u(k, i);
                        }
                        f = u(i, i);
                        g = -SIGN(sqrt(s), f);
                        h = f * g - s;
                        u(i, i) = f - g;
                        for (j = l - 1; j < n; j++)
                        {
                            for (s = 0.0, k = i; k < m; k++)
                            {
                                s += u(k, i) * u(k, j);
                            }
                            f = s / h;
                            for (k = i; k < m; k++)
                            {
                                u(k, j) += f * u(k, i);
                            }
                        }
                        for (k = i; k < m; k++)
                        {
                            u(k, i) *= scale;
                        }
                    }
                }
                w[i] = scale * g;
                g = s = scale = 0.0;
                if (i + 1 <= m && i + 1 != n)
                {
                    for (k = l - 1; k < n; k++)
                    {
                        scale += fabs(u(i, k));
                    }
                    if (scale != 0.0)
                    {
                        for (k = l - 1; k < n; k++)
                        {
                            u(i, k) /= scale;
                            s += u(i, k) * u(i, k);
                        }
                        f = u(i, l - 1);
                        g = -SIGN(sqrt(s), f);
                        h = f * g - s;
                        u(i, l - 1) = f - g;
                        for (k = l - 1; k < n; k++)
                        {
                            rv1[k] = u(i, k) / h;
                        }
                        for (j = l - 1; j < m; j++)
                        {
                            for (s = 0.0, k = l - 1; k < n; k++)
                            {
                                s += u(j, k) * u(i, k);
                            }
                            for (k = l - 1; k < n; k++)
                            {
                                u(j, k) += s * rv1[k];
                            }
                        }
                        for (k = l - 1; k < n; k++)
                        {
                            u(i, k) *= scale;
                        }
                    }
                }
                anorm = std::max(anorm, fabs(w[i]) + fabs(rv1[i]));
            }
            for (i = n - 1; i >= 0; i--)
            {
                if (i < n - 1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < n; j++)
                        {
                            v(j, i) = (u(i, j) / u(i, l)) / g;
                        }
                        for (j = l; j < n; j++)
                        {
                            for (s = 0.0, k = l; k < n; k++)
                            {
                                s += u(i, k) * v(k, j);
                            }
                            for (k = l; k < n; k++)
                            {
                                v(k, j) += s * v(k, i);
                            }
                        }
                    }
                    for (j = l; j < n; j++)
                    {
                        v(i, j) = 0.0;
                        v(j, i) = 0.0;
                    }
                }
                v(i, i) = 1.0;
                g = rv1[i];
                l = i;
            }
            for (i = std::min(m, n) - 1; i >= 0; i--)
            {
                l = i + 1;
                g = w[i];
                for (j = l; j < n; j++)
                {
                    u(i, j) = 0.0;
                }
                if (g != 0.0)
                {
                    g = 1.0 / g;
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++)
                        {
                            s += u(k, i) * u(k, j);
                        }
                        f = (s / u(i, i)) * g;
                        for (k = i; k < m; k++)
                        {
                            u(k, j) += f * u(k, i);
                        }
                    }
                    for (j = i; j < m; j++)
                    {
                        u(j, i) *= g;
                    }
                }
                else
                {
                    for (j = i; j < m; j++)
                    {
                        u(j, i) = 0.0;
                    }
                }
                u(i, i) += 1.0;
            }
            for (k = n - 1; k >= 0; k--)
            {
                for (its = 0; its < 30; its++)
                {
                    auto flag = true;
                    for (l = k; l >= 0; l--)
                    {
                        nm = l - 1;
                        if (l == 0 || fabs(rv1[l]) < EPS_DP * anorm)
                        {
                            flag = false;
                            break;
                        }
                        if (fabs(w[nm]) < EPS_DP * anorm)
                        {
                            break;
                        }
                    }
                    if (flag)
                    {
                        c = 0.0;
                        s = 1.0;
                        for (i = l; i < k + 1; i++)
                        {
                            f = s * rv1[i];
                            rv1[i] = c * rv1[i];
                            if (fabs(f) < EPS_DP * anorm)
                            {
                                break;
                            }
                            g = w[i];
                            h = sqrt(f * f + g * g);
                            w[i] = h;
                            h = 1.0 / h;
                            c = g * h;
                            s = -f * h;
                            for (j = 0; j < m; j++)
                            {
                                y = u(j, nm);
                                z = u(j, i);
                                u(j, nm) = y * c + z * s;
                                u(j, i) = z * c - y * s;
                            }
                        }
                    }
                    z = w[k];
                    if (l == k)
                    {
                        if (z < 0.0)
                        {
                            w[k] = -z;
                            for (j = 0; j < n; j++)
                            {
                                v(j, k) *= -1;
                            }
                        }
                        break;
                    }
                    x = w[l];
                    nm = k - 1;
                    y = w[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g = sqrt(f * f + 1.0);
                    f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
                    c = s = 1.0;
                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = w[i];
                        h = s * g;
                        g = c * g;
                        z = sqrt(f * f + h * h);
                        rv1[j] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = g * c - x * s;
                        h = y * s;
                        y *= c;
                        for (jj = 0; jj < n; jj++)
                        {
                            x = v(jj, j);
                            z = v(jj, i);
                            v(jj, j) = x * c + z * s;
                            v(jj, i) = z * c - x * s;
                        }
                        z = sqrt(f * f + h * h);
                        w[j] = z;
                        if (fabs(z) > EPS_DP)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            s = h * z;
                        }
                        f = c * g + s * y;
                        x = c * y - s * g;
                        for (jj = 0; jj < m; jj++)
                        {
                            y = u(jj, j);
                            z = u(jj, i);
                            u(jj, j) = y * c + z * s;
                            u(jj, i) = z * c - y * s;
                        }
                    }
                    rv1[l] = 0.0;
                    rv1[k] = f;
                    w[k] = x;
                }
            }
        }

        void solve(double *b) const
        {
            std::array<double, n> tmp{};
            auto eigen_max = *std::max_element(w.cbegin(), w.cend());
            auto tsh = 0.5 * sqrt(m + n + 1) * eigen_max * EPS_DP;
            for (int j = 0; j < n; j++)
            {
                auto s = 0.0;
                if (w[j] > tsh)
                {
                    for (int i = 0; i < m; i++)
                    {
                        s += u(i, j) * b[i];
                    }
                    s /= w[j];
                }
                tmp[j] = s;
            }
            for (int j = 0; j < n; j++)
            {
                auto s = 0.0;
                for (int jj = 0; jj < n; jj++)
                {
                    s += v(j, jj) * tmp[jj];
                }
                b[j] = s;
            }
        }

        double norm() const
        {
            return *std::max_element(w.cbegin(), w.cend());
        }

        MatrixS<m, n> u;
        MatrixS<n, 1> w;
        MatrixS<n, n> v;
        // StatusCode state;
    };

    template <size_t M, size_t N, enable_when_matrix_t<M, N> * = nullptr>
    double norm2(const MatrixS<M, N> &mat)
    {
        SVD<M, N> svd(mat);
        return svd.norm();
    }

    template <size_t M, size_t N>
    MatrixS<M, N> MatrixS<M, N>::I() const
    {
        static_assert(M == N, "only square matrix has an inverse. For rectangular matrix ,query for pinv.");
        LU<M> lu(*this);
        if (lu.s == StatusCode::SINGULAR)
        {
            return {};
        }
        auto result = eye<M>();
        for (size_t j = 0; j < M; j++)
        {
            lu.solve(result.data() + j * M);
        }
        return result;
    }

    template <size_t M, size_t N>
    double MatrixS<M, N>::det() const
    {
        static_assert(M == N, "only square matrix has a determinant.");
        LU<M> lu(*this);
        if (lu.s == StatusCode::SINGULAR)
        {
            return {};
        }
        auto D = lu.even ? 1.0 : -1.0;
        for (size_t i = 0; i < M; i++)
        {
            D *= lu.A(i, i);
        }
        return D;
    }

    template <int n>
    class EigenValue
    {
        bool m_srt;
        MatrixS<n, 1> e;

        void tred2()
        {
            for (int i = n - 1; i > 0; i--)
            {
                int l = i - 1;
                double h = 0.0;
                double scale = 0.0;
                if (l > 0)
                {
                    for (int k = 0; k < i; k++)
                    {
                        scale += fabs(z(i, k));
                    }
                    if (fabs(scale) < EPS_DP)
                    {
                        e[i] = z(i, l);
                    }
                    else
                    {
                        for (int k = 0; k < i; k++)
                        {
                            z(i, k) /= scale;
                            h += z(i, k) * z(i, k);
                        }
                        double f = z(i, l);
                        double g = f > 0.0 ? -sqrt(h) : sqrt(h);
                        e[i] = scale * g;
                        h -= f * g;
                        z(i, l) = f - g;
                        f = 0.0;
                        for (int j = 0; j < i; j++)
                        {
                            z(j, i) = z(i, j) / h;
                            g = 0.0;
                            for (int k = 0; k < j + 1; k++)
                            {
                                g += z(j, k) * z(i, k);
                            }
                            for (int k = j + 1; k < i; k++)
                            {
                                g += z(k, j) * z(i, k);
                            }
                            e[j] = g / h;
                            f += e[j] * z(i, j);
                        }
                        double hh = f / (2 * h);
                        for (int j = 0; j < i; j++)
                        {
                            f = z(i, j);
                            e[j] = g = e[j] - hh * f;
                            for (int k = 0; k < j + 1; k++)
                            {
                                z(j, k) -= f * e[k] + g * z(i, k);
                            }
                        }
                    }
                }
                else
                {
                    e[i] = z(i, l);
                }
                d[i] = h;
            }
            d[0] = 0.0;
            e[0] = 0.0;
            for (int i = 0; i < n; i++)
            {
                if (fabs(d[i]) > EPS_DP)
                {
                    for (int j = 0; j < i; j++)
                    {
                        double g = 0.0;
                        for (int k = 0; k < i; k++)
                        {
                            g += z(i, k) * z(k, j);
                        }
                        for (int k = 0; k < i; k++)
                        {
                            z(k, j) -= g * z(k, i);
                        }
                    }
                }
                d[i] = z(i, i);
                z(i, i) = 1.0;
                for (int j = 0; j < i; j++)
                {
                    z(j, i) = 0.0;
                    z(i, j) = 0.0;
                }
            }
        }

        void tliq()
        {
            for (int i = 1; i < n; i++)
            {
                e[i - 1] = e[i];
            }
            e[n - 1] = 0.0;
            for (int l = 0; l < n; l++)
            {
                int iter = 0;
                int m;
                do
                {
                    for (m = l; m < n - 1; m++)
                    {
                        if (fabs(e[m]) < EPS_DP * (fabs(d[m]) + fabs(d[m + 1])))
                        {
                            break;
                        }
                    }
                    if (m != l)
                    {
                        if (iter++ == 30)
                        // how ? throw("Too many iterations in tqli");
                        {
                            return;
                        }
                        double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                        double r = sqrt(g * g + 1.0);
                        g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                        double s = 1.0;
                        double c = 1.0;
                        double p = 0.0;
                        int i = 0;
                        for (i = m - 1; i >= l; i--)
                        {
                            double f = s * e[i];
                            double b = c * e[i];
                            r = sqrt(f * f + g * g);
                            e[i + 1] = r;
                            if (fabs(r) < EPS_DP)
                            {
                                d[i + 1] -= p;
                                e[m] = 0.0;
                                break;
                            }
                            s = f / r;
                            c = g / r;
                            g = d[i + 1] - p;
                            r = (d[i] - g) * s + 2.0 * c * b;
                            p = s * r;
                            d[i + 1] = g + p;
                            g = c * r - b;
                            for (int k = 0; k < n; k++)
                            {
                                f = z(k, i + 1);
                                z(k, i + 1) = s * z(k, i) + c * f;
                                z(k, i) = c * z(k, i) - s * f;
                            }
                        }
                        if (fabs(r) < EPS_DP && i >= l)
                        {
                            continue;
                        }
                        d[l] -= p;
                        e[l] = g;
                        e[m] = 0.0;
                    }
                } while (m != l);
            }
        }

        void eigsrt()
        {
            for (auto i = 0; i < n - 1; i++)
            {
                auto j = std::distance(d.begin() + i, std::min_element(d.begin() + i, d.end())) + i;
                if (j != i)
                {
                    std::swap(d[i], d[j]);
                    for (auto k = 0; k < n; k++)
                    {
                        std::swap(z(k, i), z(k, j));
                    }
                }
            }
        }

    public:
        EigenValue(const MatrixS<n, n> &mat, bool sorted = false)
            : m_srt(sorted), z(mat)
        {
            tred2();
            tliq();
            if (m_srt)
            {
                eigsrt();
            }
        }

        MatrixS<n, n> z;
        MatrixS<n, 1> d;
    };

    enum class Factorization : char
    {
        LU,
        QR,
        SVD
    };

    template <Factorization type, size_t M, size_t N>
    std::enable_if_t<type == Factorization::LU, EqnResult<N>>
    linsolve(const MatrixS<M, N> &A, MatrixS<M, 1> b)
    {
        static_assert(M == N, "LU decomposation only supports square matrix.");
        LU<M> lu(A);
        if (lu.s != StatusCode::SINGULAR)
        {
            lu.solve(b.data());
        }
        return {b, lu.s};
    }

    template <Factorization type, size_t M, size_t N>
    std::enable_if_t<type == Factorization::QR, EqnResult<N>>
    linsolve(const MatrixS<M, N> &A, MatrixS<M, 1> b)
    {
        static_assert(M >= N, "QR decomposation only supports column full-rank matrix.");
        QR<M, N> qr(A);
        if (qr.s != StatusCode::SINGULAR)
        {
            qr.solve(b.data());
        }
        return {b.template sub<N, 1>(0, 0), qr.s};
    }

    template <Factorization type, size_t M, size_t N>
    std::enable_if_t<type == Factorization::SVD, EqnResult<N>>
    linsolve(const MatrixS<M, N> &A, MatrixS<M, 1> b)
    {
        SVD<M, N> svd(A);
        svd.solve(b.data());
        return {b.template sub<N, 1>(0, 0), StatusCode::NORMAL};
    }

    template <size_t M, size_t N>
    MatrixS<N, M> pinv(const MatrixS<M, N> &mat)
    {
        SVD<M, N> svd(mat);
        const auto &w = svd.w;
        MatrixS<N, N> W{};
        auto eigen_max = *std::max_element(w.cbegin(), w.cend());
        auto tsh = 0.5 * sqrt(M + N + 1) * eigen_max * EPS_DP;
        for (size_t i = 0; i < N; i++)
        {
            if (fabs(w[i]) > tsh)
            {
                W(i, i) = 1.0 / w[i];
            }
        }
        return svd.v * W * svd.u.T();
    }

    inline EqnResult<2> quadsolve(double a, double b, double c)
    {
        auto det = b * b - 4 * a * c;
        if (det < 0.0 || fabs(a) < EPS_DP)
        {
            return {{}, StatusCode::SINGULAR};
        }
        auto q = -0.5 * (b + SIGN(sqrt(det), b));
        return {{q / a, c / q}, StatusCode::NORMAL};
    }
}
#endif