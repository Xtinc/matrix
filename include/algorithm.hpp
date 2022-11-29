#ifndef VVERY_SIMPLE_MATRIX_HEADER
#define VVERY_SIMPLE_MATRIX_HEADER

#include "matrix.hpp"

namespace ppx
{
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

    inline double pythag(double a, double b)
    {
        double at = fabs(a);
        double bt = fabs(b);
        double ct{};
        double result{};

        if (at > bt)
        {
            ct = bt / at;
            result = at * sqrt(1.0 + ct * ct);
        }
        else if (bt > 0.0)
        {
            ct = at / bt;
            result = bt * sqrt(1.0 + ct * ct);
        }
        else
        {
            result = 0.0;
        }
        return result;
    }

    template <size_t M, size_t N>
    Matrix<M, N> svdcmp(Matrix<M, N> a, Matrix<N, 1> &w, Matrix<N, N> &v)
    {
        int flag, i, its, j, jj, k, l, nm;
        double c, f, h, s, x, y, z;
        double anorm = 0.0;
        double g = 0.0;
        double scale = 0.0;
        double rv1[N];

        auto SIGN = [](double a, double b)
        {
            return b > 0.0 ? fabs(a) : -fabs(a);
        };

        if (M < N)
        {
            return {};
        }

        /* Householder reduction to bidiagonal form */
        for (i = 0; i < N; i++)
        {
            /* left-hand reduction */
            l = i + 1;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i < M)
            {
                for (k = i; k < M; k++)
                {
                    scale += fabs(a(k, i));
                }
                if (scale > gl_rep_eps)
                {
                    for (k = i; k < M; k++)
                    {
                        a(k, i) = a(k, i) / scale;
                        s += a(k, i) * a(k, i);
                    }
                    f = a(i, i);
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    a(i, i) = f - g;
                    if (i != N - 1)
                    {
                        for (j = l; j < N; j++)
                        {
                            for (s = 0.0, k = i; k < M; k++)
                            {
                                s += a(k, i) * a(k, j);
                            }
                            f = s / h;
                            for (k = i; k < M; k++)
                            {
                                a(k, j) += f * a(k, i);
                            }
                        }
                    }
                    for (k = i; k < M; k++)
                    {
                        a(k, i) = a(k, i) * scale;
                    }
                }
            }
            w[i] = scale * g;

            /* right-hand reduction */
            g = 0.0;
            s = 0.0;
            scale = 0.0;
            if (i < M && i != N - 1)
            {
                for (k = l; k < N; k++)
                {
                    scale += fabs(a(i, k));
                }
                if (scale > gl_rep_eps)
                {
                    for (k = l; k < N; k++)
                    {
                        a(i, k) = a(i, k) / scale;
                        s += a(i, k) * a(i, k);
                    }
                    f = a(i, l);
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    a(i, l) = f - g;
                    for (k = l; k < N; k++)
                    {
                        rv1[k] = a(i, k) / h;
                    }
                    if (i != M - 1)
                    {
                        for (j = l; j < M; j++)
                        {
                            for (s = 0.0, k = l; k < N; k++)
                            {
                                s += a(j, k) * a(i, k);
                            }
                            for (k = l; k < N; k++)
                            {
                                a(j, k) += s * rv1[k];
                            }
                        }
                    }
                    for (k = l; k < N; k++)
                    {
                        a(i, k) = a(i, k) * scale;
                    }
                }
            }
            anorm = gl_get_more_dynamic(anorm, (fabs(w[i]) + fabs(rv1[i])));
        }

        /* accumulate the right-hand transformation */
        for (i = N - 1; i >= 0; i--)
        {
            if (i < N - 1)
            {
                if (g > gl_rep_eps)
                {
                    for (j = l; j < N; j++)
                    {
                        v(j, i) = (a(i, j) / a(i, l)) / g;
                    }
                    /* double division to avoid underflow */
                    for (j = l; j < N; j++)
                    {
                        for (s = 0.0, k = l; k < N; k++)
                        {
                            s += a(i, k) * v(k, j);
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

        /* accumulate the left-hand transformation */
        for (i = N - 1; i >= 0; i--)
        {
            l = i + 1;
            g = w[i];
            if (i < N - 1)
            {
                for (j = l; j < N; j++)
                {
                    a(i, j) = 0.0;
                }
            }
            if (g > gl_rep_eps)
            {
                g = 1.0 / g;
                if (i != N - 1)
                {
                    for (j = l; j < N; j++)
                    {
                        for (s = 0.0, k = l; k < M; k++)
                        {
                            s += a(k, i) * a(k, j);
                        }
                        f = (s / a(i, i)) * g;
                        for (k = i; k < M; k++)
                        {
                            a(k, j) += f * a(k, i);
                        }
                    }
                }
                for (j = i; j < M; j++)
                {
                    a(j, i) = a(j, i) * g;
                }
            }
            else
            {
                for (j = i; j < M; j++)
                {
                    a(j, i) = 0.0;
                }
            }
            a(i, i) += 1.0;
        }

        /* diagonalize the bidiagonal form */
        for (k = N - 1; k >= 0; k--)
        { /* loop over singular values */
            for (its = 0; its < 30; its++)
            { /* loop over allowed iterations */
                flag = 1;
                for (l = k; l >= 0; l--)
                { /* test for splitting */
                    nm = l - 1;
                    if (details::is_same((rv1[l]) + anorm, anorm))
                    {
                        flag = 0;
                        break;
                    }
                    if (details::is_same((w[nm]) + anorm, anorm))
                    {
                        break;
                    }
                }
                if (flag > gl_rep_eps)
                {
                    c = 0.0;
                    s = 1.0;
                    for (i = l; i <= k; i++)
                    {
                        f = s * rv1[i];
                        if (fabs(f) + anorm != anorm)
                        {
                            g = w[i];
                            h = pythag(f, g);
                            w[i] = h;
                            h = 1.0 / h;
                            c = g * h;
                            s = -f * h;
                            for (j = 0; j < M; j++)
                            {
                                y = a(j, nm);
                                z = a(j, i);
                                a(j, nm) = y * c + z * s;
                                a(j, i) = z * c - y * s;
                            }
                        }
                    }
                }
                z = w[k];
                if (l == k)
                { /* convergence */
                    if (z < 0.0)
                    { /* make singular value nonnegative */
                        w[k] = -z;
                        for (j = 0; j < N; j++)
                        {
                            v(j, k) = -v(j, k);
                        }
                    }
                    break;
                }
                if (its >= 30)
                {
                    // fprintf(stderr, "No convergence after 30,000! iterations \N");
                    return {};
                }

                /* shift from bottom 2 x 2 minor */
                x = w[l];
                nm = k - 1;
                y = w[nm];
                g = rv1[nm];
                h = rv1[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = pythag(f, 1.0);
                f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

                /* next QR transformation */
                c = 1.0;
                s = 1.0;
                for (j = l; j <= nm; j++)
                {
                    i = j + 1;
                    g = rv1[i];
                    y = w[i];
                    h = s * g;
                    g = c * g;
                    z = pythag(f, h);
                    rv1[j] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y = y * c;
                    for (jj = 0; jj < N; jj++)
                    {
                        x = v(jj, j);
                        z = v(jj, i);
                        v(jj, j) = x * c + z * s;
                        v(jj, i) = z * c - x * s;
                    }
                    z = pythag(f, h);
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
                        y = a(jj, j);
                        z = a(jj, i);
                        a(jj, j) = y * c + z * s;
                        a(jj, i) = z * c - y * s;
                    }
                }
                rv1[l] = 0.0;
                rv1[k] = f;
                w[k] = x;
            }
        }
        return a;
    }
}
#endif