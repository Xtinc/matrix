#ifndef VVERY_SIMPLE_MATRIX_HEADER
#define VVERY_SIMPLE_MATRIX_HEADER

#include "matrix.hpp"

namespace ppx
{
    // matrix related

    enum class factorization : char
    {
        LU,
        QR,
        SVD
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
    inline Matrix<1, 1> adjugate(const Matrix<1, 1> &mat)
    {
        return {1};
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

    // solver linear system

    inline double SIGN(double a, double b)
    {
        return b > 0.0 ? fabs(a) : -fabs(a);
    };

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

    template <size_t M, size_t N>
    Matrix<M, N> qrdcmp(Matrix<M, N> a, Matrix<N, 1> &c, Matrix<N, 1> &d, bool &sing)
    {
        sing = false;
        double scale, sigma, sum, tau;
        for (int k = 0; k < N - 1; k++)
        {
            scale = 0.0;
            for (int i = k; i < N; i++)
            {
                scale = gl_get_more_dynamic(scale, fabs(a(i, k)));
            }
            if (scale < gl_rep_eps)
            {
                sing = true;
                c[k] = 0.0;
                d[k] = 0.0;
            }
            else
            {
                for (int i = k; i < M; i++)
                {
                    a(i, k) /= scale;
                }
                sum = 0.0;
                for (int i = k; i < M; i++)
                {
                    sum += a(i, k) * a(i, k);
                }
                sigma = fabs(a(k, k)) < gl_rep_eps ? sqrt(sum) : SIGN(sqrt(sum), a(k, k));
                a(k, k) += sigma;
                c[k] = sigma * a(k, k);
                d[k] = -scale * sigma;
                for (int j = k + 1; j < N; j++)
                {
                    sum = 0.0;
                    for (int i = k; i < M; i++)
                    {
                        sum += a(i, k) * a(i, j);
                    }
                    tau = sum / c[k];
                    for (int i = k; i < M; i++)
                    {
                        a(i, j) -= tau * a(i, k);
                        auto ccc = a(i, j);
                        std::cout << i << " " << j << " " << a(i, j) << std::endl;
                    }
                }
            }
        }
        d[N - 1] = a(N - 1, N - 1);
        if (fabs(d[N - 1]) < gl_rep_eps)
        {
            sing = true;
        }
        return a;
    }

    template <size_t M, size_t N>
    void qrsolv(const Matrix<M, N> &a, const Matrix<N, 1> &c, const Matrix<N, 1> &d, double *b)
    {
        for (int j = 0; j < N - 1; j++)
        {
            double sum = 0.0;
            for (int i = j; i < N; i++)
            {
                sum += a(i, j) * b[i];
            }
            double tau = sum / c[j];
            for (int i = j; i < N; i++)
            {
                b[i] -= tau * a(i, j);
            }
        }
        b[N - 1] /= d[N - 1];
        for (int i = N - 2; i >= 0; i--)
        {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++)
            {
                sum += a(i, j) * b[j];
            }
            b[i] = (b[i] - sum) / d[i];
        }
    }

    template <size_t M, size_t N>
    Matrix<M, N> svdcmp(Matrix<M, N> u, Matrix<N, 1> &w, Matrix<N, N> &v)
    {
        auto PYTHAG = [](double a, double b)
        {
            double absa = fabs(a);
            double absb = fabs(b);
            if (absa > absb)
            {
                auto tmp = absb / absa;
                return absa * sqrt(1.0 + tmp * tmp);
            }
            else if (fabs(absb) < gl_rep_eps)
            {
                return 0.0;
            }
            else
            {
                auto tmp = absa / absb;
                return absb * sqrt(1.0 + tmp * tmp);
            }
        };
        bool flag;
        int i, its, j, jj, k, l, nm;
        double anorm, c, f, g, h, s, scale, x, y, z;
        double rv1[N];
        g = 0.0;
        scale = 0.0;
        anorm = 0.0;
        for (i = 0; i < N; i++)
        {
            l = i + 2;
            rv1[i] = scale * g;
            g = 0.0;
            s = 0.0;
            scale = 0.0;
            if (i < M)
            {
                for (k = i; k < M; k++)
                {
                    scale += fabs(u(k, i));
                }
                if (scale > gl_rep_eps)
                {
                    for (k = i; k < M; k++)
                    {
                        u(k, i) /= scale;
                        s += u(k, i) * u(k, i);
                    }
                    f = u(i, i);
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    u(i, i) = f - g;
                    for (j = l - 1; j < N; j++)
                    {
                        for (s = 0.0, k = i; k < M; k++)
                        {
                            s += u(k, i) * u(k, j);
                        }
                        f = s / h;
                        for (k = i; k < M; k++)
                        {
                            u(k, j) += f * u(k, i);
                        }
                    }
                    for (k = i; k < M; k++)
                    {
                        u(k, i) *= scale;
                    }
                }
            }
            w[i] = scale * g;
            g = s = scale = 0.0;
            if (i + 1 <= M && i + 1 != N)
            {
                for (k = l - 1; k < N; k++)
                {
                    scale += fabs(u(i, k));
                }
                if (scale > gl_rep_eps)
                {
                    for (k = l - 1; k < N; k++)
                    {
                        u(i, k) /= scale;
                        s += u(i, k) * u(i, k);
                    }
                    f = u(i, l - 1);
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    u(i, l - 1) = f - g;
                    for (k = l - 1; k < N; k++)
                    {
                        rv1[k] = u(i, k) / h;
                    }
                    for (j = l - 1; j < M; j++)
                    {
                        for (s = 0.0, k = l - 1; k < N; k++)
                        {
                            s += u(j, k) * u(i, k);
                        }
                        for (k = l - 1; k < N; k++)
                        {
                            u(j, k) += s * rv1[k];
                        }
                    }
                    for (k = l - 1; k < N; k++)
                    {
                        u(i, k) *= scale;
                    }
                }
            }
            anorm = gl_get_more_dynamic(anorm, fabs(w[i]) + fabs(rv1[i]));
        }
        for (i = N - 1; i >= 0; i--)
        {
            if (i < N - 1)
            {
                if (g != 0.0)
                {
                    for (j = l; j < N; j++)
                    {
                        v(j, i) = (u(i, j) / u(i, l)) / g;
                    }
                    for (j = l; j < N; j++)
                    {
                        for (s = 0.0, k = l; k < N; k++)
                        {
                            s += u(i, k) * v(k, j);
                        }
                        for (k = l; k < N; k++)
                        {
                            v(k, j) += s * v(k, i);
                        }
                    }
                }
                for (j = l; j < N; j++)
                {
                    v(i, j) = 0.0;
                    v(j, i) = 0.0;
                }
            }
            v(i, i) = 1.0;
            g = rv1[i];
            l = i;
        }
        for (i = (int)gl_get_less(M, N) - 1; i >= 0; i--)
        {
            l = i + 1;
            g = w[i];
            for (j = l; j < N; j++)
            {
                u(i, j) = 0.0;
            }
            if (fabs(g) > gl_rep_eps)
            {
                g = 1.0 / g;
                for (j = l; j < N; j++)
                {
                    for (s = 0.0, k = l; k < M; k++)
                    {
                        s += u(k, i) * u(k, j);
                    }
                    f = (s / u(i, i)) * g;
                    for (k = i; k < M; k++)
                    {
                        u(k, j) += f * u(k, i);
                    }
                }
                for (j = i; j < M; j++)
                {
                    u(j, i) *= g;
                }
            }
            else
            {
                for (j = i; j < M; j++)
                {
                    u(j, i) = 0.0;
                }
            }
            u(i, i) += 1.0;
        }
        for (k = N - 1; k >= 0; k--)
        {
            for (its = 0; its < 30; its++)
            {
                flag = true;
                for (l = k; l >= 0; l--)
                {
                    nm = l - 1;
                    if (l == 0 || fabs(rv1[l]) <= gl_rep_eps * anorm)
                    {
                        flag = false;
                        break;
                    }
                    if (fabs(w[nm]) <= gl_rep_eps * anorm)
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
                        if (fabs(f) <= gl_rep_eps * anorm)
                        {
                            break;
                        }
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.0 / h;
                        c = g * h;
                        s = -f * h;
                        for (j = 0; j < M; j++)
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
                        for (j = 0; j < N; j++)
                        {
                            v(j, k) = -v(j, k);
                        }
                    }
                    break;
                }
                if (its == 29)
                {
                    return {};
                }
                x = w[l];
                nm = k - 1;
                y = w[nm];
                g = rv1[nm];
                h = rv1[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = PYTHAG(f, 1.0);
                f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
                c = s = 1.0;
                for (j = l; j <= nm; j++)
                {
                    i = j + 1;
                    g = rv1[i];
                    y = w[i];
                    h = s * g;
                    g = c * g;
                    z = PYTHAG(f, h);
                    rv1[j] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y *= c;
                    for (jj = 0; jj < N; jj++)
                    {
                        x = v(jj, j);
                        z = v(jj, i);
                        v(jj, j) = x * c + z * s;
                        v(jj, i) = z * c - x * s;
                    }
                    z = PYTHAG(f, h);
                    w[j] = z;
                    if (z)
                    {
                        z = 1.0 / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for (jj = 0; jj < M; jj++)
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
        return u;
    }

    template <size_t M, size_t N>
    void svbksb(const Matrix<M, N> &u, const Matrix<N, 1> &w, const Matrix<N, N> &v, double *b)
    {
        double tmp[N];
        auto eigen_max = *std::max_element(w.cbegin(), w.cend());
        auto tsh = 0.5 * sqrt(M + N + 1) * eigen_max * gl_rep_eps;
        for (int j = 0; j < N; j++)
        {
            auto s = 0.0;
            if (w[j] > tsh)
            {
                for (int i = 0; i < M; i++)
                {
                    s += u(i, j) * b[i];
                }
                s /= w[j];
            }
            tmp[j] = s;
        }
        for (int j = 0; j < N; j++)
        {
            auto s = 0.0;
            for (int jj = 0; jj < N; jj++)
            {
                s += v(j, jj) * tmp[jj];
            }
            b[j] = s;
        }
    }

    template <factorization type, size_t M, size_t N>
    typename std::enable_if<type == factorization::LU, Matrix<N, 1>>::type
    solve(const Matrix<M, N> &A, Matrix<M, 1> b)
    {
        std::array<int, M> indx{};
        auto even = true;
        auto LU = ludcmp(A, indx, even);
        ludbksb(LU, indx, b.data());
        return b;
    }

    template <factorization type, size_t M, size_t N>
    typename std::enable_if<type == factorization::QR, Matrix<N, 1>>::type
    solve(const Matrix<M, N> &A, Matrix<M, 1> b)
    {
        bool singular = false;
        Matrix<N, 1> c;
        Matrix<N, 1> d;
        auto R = qrdcmp(A, c, d, singular);
        if (singular)
        {
            return {};
        }
        qrsolv(R, c, d, b.data());
        return slice<N, 1>(b, 0, 0);
    }

    template <factorization type, size_t M, size_t N>
    typename std::enable_if<type == factorization::SVD, Matrix<N, 1>>::type
    solve(const Matrix<M, N> &A, Matrix<M, 1> b)
    {
        Matrix<N, 1> w{};
        Matrix<N, N> V{};
        auto U = svdcmp(A, w, V);
        svbksb(U, w, V, b.data());
        return slice<N, 1>(b, 0, 0);
    }

    template <size_t M, size_t N>
    Matrix<N, M> pinv(const Matrix<M, N> &mat)
    {
        Matrix<N, 1> w{};
        Matrix<N, N> V{};
        auto U = svdcmp(mat, w, V);
        Matrix<N, N> W{};
        for (size_t i = 0; i < N; i++)
        {
            if (fabs(w[i]) > gl_rep_eps)
            {
                W(i, i) = 1.0 / w[i];
            }
        }
        return V * W * U.T();
    }
}
#endif