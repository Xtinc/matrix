#ifndef VVERY_SIMPLE_LINEAR_ALGORITHM_HEADER
#define VVERY_SIMPLE_LINEAR_ALGORITHM_HEADER

#include "matrixs.hpp"
#include <complex>

namespace ppx
{
    template <typename T>
    int SIGN(T a)
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
    T SIGN(T a, T b)
    {
        return b >= T{} ? std::abs(a) : -std::abs(a);
    }

    template <typename T>
    T SQR(T a)
    {
        return a * a;
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

    // matrix related
    enum class Factorization : char
    {
        LU,
        QR,
        SVD
    };

    enum class EigenSystem : char
    {
        SymOnlyVal,
        SymValAndVec,
        SymValAndVecSorted,
        // GemOnlyVal,
        // GemValAndVec,
        // GemValAndVecSorted
        // todo
    };

    template <size_t N>
    struct EigResult
    {
        MatrixS<N, N> vec;
        MatrixS<N, 1> val;
    };

    template <size_t N>
    struct EqnResult
    {
        MatrixS<N, 1> x;
        StatusCode s = StatusCode::NORMAL;
    };

    template <size_t N>
    std::ostream &operator<<(std::ostream &os, const EqnResult<N> &self)
    {
        os << "EqnResult<" << N << ">:\n"
           << "Status:\t" << self.s << "\n"
           << "x     =\t" << self.x << std::endl;
        return os;
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
    enable_when_array_t<M, N, double> norm2(const MatrixS<M, N> &mat)
    {
        double res = 0.0;
        for (auto ele : mat)
        {
            res += ele * ele;
        }
        return sqrt(res);
    }

    template <size_t M, size_t N>
    enable_when_array_t<M, N, double> norminf(const MatrixS<M, N> &mat)
    {
        double max = -std::numeric_limits<double>::max();
        for (auto i : mat)
        {
            auto ti = fabs(i);
            if (ti > max)
            {
                max = ti;
            }
        }
        return max;
    }

    template <size_t M, size_t N>
    double trace(const MatrixS<M, N> &mat)
    {
        return mat.trace();
    }

    template <size_t N>
    double inner_product(const MatrixS<N, 1> &a, const MatrixS<N, 1> &b)
    {
        return std::inner_product(a.cbegin(), a.cend(), b.cbegin(), 0.0);
    }

    template <size_t M, size_t N>
    std::pair<size_t, size_t> maxloc(const MatrixS<M, N> &mat)
    {
        auto max_pos = std::max_element(mat.begin(), mat.end());
        auto max_dis = std::div(std::distance(mat.begin(), max_pos), M);
        return {max_dis.rem, max_dis.quot};
    }

    enum class Orientation
    {
        Horizontal,
        Vertical
    };

    template <Orientation ori, size_t M, size_t N>
    std::enable_if_t<ori == Orientation::Vertical, MatrixS<2 * M, N>> concat(const MatrixS<M, N> &a, const MatrixS<M, N> &b)
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

    template <Orientation ori, size_t M, size_t N>
    std::enable_if_t<ori == Orientation::Horizontal, MatrixS<M, 2 * N>> concat(const MatrixS<M, N> &a,
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

    template <size_t N>
    FacResult<N, N> ludcmp(MatrixS<N, N> A, std::array<int, N> &indx, bool &even)
    {
        constexpr int n = N;
        even = true;
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
            if (valmax < EPS_SP)
            {
                return {{}, StatusCode::SINGULAR};
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
        return {A, StatusCode::CONVERGED};
    }

    template <size_t N>
    void ludbksb(const MatrixS<N, N> &A, const std::array<int, N> &indx, double *b)
    {
        constexpr int n = N;
        std::array<double, N> y{};
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

    template <size_t M, size_t N>
    FacResult<M, N> qrdcmp(MatrixS<M, N> a, MatrixS<N, 1> &c, MatrixS<N, 1> &d)
    {
        constexpr int m = M;
        constexpr int n = N;
        StatusCode sing = StatusCode::NORMAL;
        double scale, sigma, sum, tau;
        for (int k = 0; k < n; k++)
        {
            scale = 0.0;
            for (int i = k; i < n; i++)
            {
                scale = std::max(scale, fabs(a(i, k)));
            }
            if (scale < EPS_SP)
            {
                sing = StatusCode::SINGULAR;
                c[k] = 0.0;
                d[k] = 0.0;
            }
            else
            {
                sum = 0.0;
                for (int i = k; i < m; i++)
                {
                    sum += a(i, k) * a(i, k);
                }
                sigma = fabs(a(k, k)) < EPS_SP ? sqrt(sum) : SIGN(sqrt(sum), a(k, k));
                a(k, k) += sigma;
                c[k] = sigma * a(k, k);
                // d[k] = -scale * sigma;
                d[k] = -sigma;
                for (int j = k + 1; j < n; j++)
                {
                    sum = 0.0;
                    for (int i = k; i < m; i++)
                    {
                        sum += a(i, k) * a(i, j);
                    }
                    tau = sum / c[k];
                    for (int i = k; i < m; i++)
                    {
                        a(i, j) -= tau * a(i, k);
                    }
                    // details::PRINT_SINGLE_ELEMENTS(a);
                }
            }
        }
        // d[N - 1] = a(N - 1, N - 1);
        if (fabs(d[n - 1]) < EPS_SP)
        {
            sing = StatusCode::SINGULAR;
        }
        return {a, sing};
    }

    template <size_t M, size_t N>
    void qrsolv(const MatrixS<M, N> &a, const MatrixS<N, 1> &c, const MatrixS<N, 1> &d, double *b)
    {
        constexpr int m = M;
        constexpr int n = N;
        for (int j = 0; j < n; j++)
        {
            double sum = 0.0;
            for (int i = j; i < m; i++)
            {
                sum += a(i, j) * b[i];
            }
            double tau = sum / c[j];
            for (int i = j; i < m; i++)
            {
                b[i] -= tau * a(i, j);
            }
        }
        // b[N - 1] /= d[N - 1];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++)
            {
                sum += a(i, j) * b[j];
            }
            b[i] = (b[i] - sum) / d[i];
        }
    }

    template <size_t M, size_t N>
    FacResult<M, N> svdcmp(MatrixS<M, N> u, MatrixS<N, 1> &w, MatrixS<N, N> &v)
    {
        constexpr int m = M;
        constexpr int n = N;
        int i, its, j, jj, k, nm;
        double c, f, h, s, x, y, z;
        std::array<double, N> rv1{};
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
                if (scale > EPS_SP)
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
                if (scale > EPS_SP)
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
                if (fabs(g) > EPS_SP)
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
            if (fabs(g) > EPS_SP)
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
                bool flag = true;
                nm = k - 1;
                for (l = k; l >= 0; l--)
                {
                    nm = l - 1;
                    if (l == 0 || fabs(rv1[l]) <= EPS_SP * anorm)
                    {
                        flag = false;
                        break;
                    }
                    if (fabs(w[nm]) <= EPS_SP * anorm)
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
                        if (fabs(f) <= EPS_SP * anorm)
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
                            v(j, k) = -v(j, k);
                        }
                    }
                    break;
                }
                if (its == 29)
                {
                    return {MatrixS<M, N>(), StatusCode::SINGULAR};
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
                    if (fabs(z) > EPS_SP)
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
        return {u, StatusCode::CONVERGED};
    }

    template <size_t M, size_t N>
    void svbksb(const MatrixS<M, N> &u, const MatrixS<N, 1> &w, const MatrixS<N, N> &v, double *b)
    {
        const int m = M;
        const int n = N;
        std::array<double, N> tmp{};
        auto eigen_max = *std::max_element(w.cbegin(), w.cend());
        auto tsh = 0.5 * sqrt(m + n + 1) * eigen_max * EPS_SP;
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

    template <Factorization type, size_t M, size_t N>
    std::enable_if_t<type == Factorization::LU, EqnResult<N>>
    linsolve(const MatrixS<M, N> &A, MatrixS<M, 1> b)
    {
        std::array<int, M> indx{};
        auto even = true;
        auto LU = ludcmp(A, indx, even);
        if (LU.s == StatusCode::CONVERGED)
        {
            ludbksb(LU.x, indx, b.data());
        }
        return {b, LU.s};
    }

    template <Factorization type, size_t M, size_t N>
    std::enable_if_t<type == Factorization::QR, EqnResult<N>>
    linsolve(const MatrixS<M, N> &A, MatrixS<M, 1> b)
    {
        MatrixS<N, 1> c, d;
        auto R = qrdcmp(A, c, d);
        if (R.s != StatusCode::SINGULAR)
        {
            qrsolv(R.x, c, d, b.data());
        }
        return {b.template sub<N, 1>(0, 0), R.s};
    }

    template <Factorization type, size_t M, size_t N>
    std::enable_if_t<type == Factorization::SVD, EqnResult<N>>
    linsolve(const MatrixS<M, N> &A, MatrixS<M, 1> b)
    {
        MatrixS<N, 1> w{};
        MatrixS<N, N> V{};
        auto U = svdcmp(A, w, V);
        if (U.s == StatusCode::CONVERGED)
        {
            svbksb(U.x, w, V, b.data());
        }
        return {b.template sub<N, 1>(0, 0), U.s};
    }

    template <size_t M, size_t N>
    MatrixS<N, M> pinv(const MatrixS<M, N> &mat)
    {
        MatrixS<N, 1> w{};
        MatrixS<N, N> V{};
        bool sing = false;
        auto U = svdcmp(mat, w, V, sing);
        if (sing)
        {
            return {};
        }
        MatrixS<N, N> W{};
        auto eigen_max = *std::max_element(w.cbegin(), w.cend());
        auto tsh = 0.5 * sqrt(M + N + 1) * eigen_max * EPS_SP;
        for (size_t i = 0; i < N; i++)
        {
            if (fabs(w[i]) > tsh)
            {
                W(i, i) = 1.0 / w[i];
            }
        }
        return V * W * U.T();
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

    // eigen system
    namespace details
    {
        template <size_t N>
        MatrixS<N, N> tred2(MatrixS<N, N> z, MatrixS<N, 1> &d, MatrixS<N, 1> &e)
        {
            constexpr int n = N;
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
                    if (fabs(scale) < EPS_SP)
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
                if (fabs(d[i]) > EPS_SP)
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
            return z;
        }

        template <size_t N>
        void tliq(MatrixS<N, N> &z, MatrixS<N, 1> &d, MatrixS<N, 1> &e)
        {
            constexpr int n = N;
            for (int i = 1; i < n; i++)
            {
                e[i - 1] = e[i];
            }
            e[n - 1] = 0.0;
            for (int l = 0; l < n; l++)
            {
                int iter = 0;
                int m = 0;
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
                            if (fabs(r) < EPS_SP)
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
                        if (fabs(r) < EPS_SP && i >= l)
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

        template <size_t N>
        void eigsrt(MatrixS<N, N> &mat, MatrixS<N, 1> &vec)
        {
            constexpr int n = N;
            for (auto i = 0; i < n - 1; i++)
            {
                auto j = std::distance(vec.begin() + i, std::min_element(vec.begin() + i, vec.end())) + i;
                if (j != i)
                {
                    std::swap(vec[i], vec[j]);
                    for (auto k = 0; k < n; k++)
                    {
                        std::swap(mat(k, i), mat(k, j));
                    }
                }
            }
        }
    }

    template <EigenSystem type, size_t N>
    std::enable_if_t<type == EigenSystem::SymValAndVec, EigResult<N>>
    eig(const MatrixS<N, N> &mat)
    {
        MatrixS<N, 1> e, eig_value;
        auto result = details::tred2(mat, eig_value, e);
        details::tliq(result, eig_value, e);
        return {result, eig_value};
    }

    template <EigenSystem type, size_t N>
    std::enable_if_t<type == EigenSystem::SymValAndVecSorted, EigResult<N>>
    eig(const MatrixS<N, N> &mat)
    {
        MatrixS<N, 1> e, eig_value;
        auto result = details::tred2(mat, eig_value, e);
        details::tliq(result, eig_value, e);
        details::eigsrt(result, eig_value);
        return {result, eig_value};
    }

    template <EigenSystem type, size_t N>
    std::enable_if_t<type == EigenSystem::SymOnlyVal, MatrixS<N, 1>>
    eig(const MatrixS<N, N> &mat)
    {
        MatrixS<N, 1> e, eig_value;
        auto result = details::tred2(mat, eig_value, e);
        details::tliq(result, eig_value, e);
        std::sort(eig_value.begin(), eig_value.end());
        return eig_value;
    }

}
#endif