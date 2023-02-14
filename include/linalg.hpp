#ifndef VVERY_SIMPLE_ALGORITHM1_HEADER
#define VVERY_SIMPLE_ALGORITHM1_HEADER

#include "matrixs.hpp"
#include <complex>

namespace ppx
{
    inline double SIGN(double a, double b)
    {
        return b > 0.0 ? fabs(a) : -fabs(a);
    }

    inline double SQR(double a)
    {
        return a * a;
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

    template <size_t A, size_t B, size_t M, size_t N>
    MatrixS<A, B> slice(const MatrixS<M, N> &m, size_t row_start, size_t col_start)
    {
        if (row_start + A > M || col_start + B > N)
        {
            return {};
        }
        MatrixS<A, B> res{};
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

    // solver linear system

    template <typename MatNN, typename VecN1>
    FacResult<MatNN> ludcmp(MatNN A, VecN1 &indx, bool &even)
    {
        // static_assert(std::is_same<typename VecN1::value_type, int>::value);
        assert(A.rows() == A.cols());
        assert(indx.rows() == A.cols() && indx.cols() == 1);
        const int n = A.rows();
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
                return {MatNN(), StatusCode::SINGULAR};
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

    template <typename MatNN, typename VecN1>
    void ludbksb(const MatNN &A, const VecN1 &indx, double *b)
    {
        // static_assert(std::is_same<typename VecN1::value_type, int>::value);
        assert(A.rows() == A.cols());
        assert(indx.rows() == A.cols() && indx.cols() == 1);
        const int n = A.rows();
        std::vector<double> y(n, 0.0);
        y[0] = b[indx[0]];
        for (int row = 1; row < n; row++)
        {
            auto sum = 0.0;
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
            auto sum = 0.0;
            for (int col = id; col < n; col++)
            {
                sum += A(row, col) * b[col];
            }
            b[row] = (y[row] - sum) / A(row, row);
        }
    }

    template <size_t M, size_t N>
    MatrixS<M, N> qrdcmp(MatrixS<M, N> a, MatrixS<N, 1> &c, MatrixS<N, 1> &d, bool &sing)
    {
        constexpr int IN = N;
        constexpr int IM = M;
        sing = false;
        double scale, sigma, sum, tau;
        for (int k = 0; k < IN; k++)
        {
            scale = 0.0;
            for (int i = k; i < IN; i++)
            {
                scale = std::max(scale, fabs(a(i, k)));
            }
            if (scale < EPS_SP)
            {
                sing = true;
                c[k] = 0.0;
                d[k] = 0.0;
            }
            else
            {
                // for (int i = k; i < M; i++)
                // {
                //     a(i, k) /= scale;
                // }
                sum = 0.0;
                for (int i = k; i < IM; i++)
                {
                    sum += a(i, k) * a(i, k);
                }
                sigma = fabs(a(k, k)) < EPS_SP ? sqrt(sum) : SIGN(sqrt(sum), a(k, k));
                a(k, k) += sigma;
                c[k] = sigma * a(k, k);
                // d[k] = -scale * sigma;
                d[k] = -sigma;
                for (int j = k + 1; j < IN; j++)
                {
                    sum = 0.0;
                    for (int i = k; i < IM; i++)
                    {
                        sum += a(i, k) * a(i, j);
                    }
                    tau = sum / c[k];
                    for (int i = k; i < IM; i++)
                    {
                        a(i, j) -= tau * a(i, k);
                    }
                    // details::PRINT_SINGLE_ELEMENTS(a);
                }
            }
        }
        // d[N - 1] = a(N - 1, N - 1);
        if (fabs(d[IN - 1]) < EPS_SP)
        {
            sing = true;
        }
        return a;
    }

    template <size_t M, size_t N>
    void qrsolv(const MatrixS<M, N> &a, const MatrixS<N, 1> &c, const MatrixS<N, 1> &d, double *b)
    {
        constexpr int IN = N;
        constexpr int IM = M;
        for (int j = 0; j < IN; j++)
        {
            double sum = 0.0;
            for (int i = j; i < IM; i++)
            {
                sum += a(i, j) * b[i];
            }
            double tau = sum / c[j];
            for (int i = j; i < IM; i++)
            {
                b[i] -= tau * a(i, j);
            }
        }
        // b[N - 1] /= d[N - 1];
        for (int i = IN - 1; i >= 0; i--)
        {
            double sum = 0.0;
            for (int j = i + 1; j < IN; j++)
            {
                sum += a(i, j) * b[j];
            }
            b[i] = (b[i] - sum) / d[i];
        }
    }

    template <size_t M, size_t N>
    MatrixS<M, N> svdcmp(MatrixS<M, N> u, MatrixS<N, 1> &w, MatrixS<N, N> &v, bool &sing)
    {
        constexpr int IN = N;
        constexpr int IM = M;
        sing = false;
        bool flag;
        int i, its, j, jj, k, l, nm;
        double anorm, c, f, g, h, s, scale, x, y, z;
        double rv1[N];
        g = 0.0;
        scale = 0.0;
        anorm = 0.0;
        l = 0;
        for (i = 0; i < IN; i++)
        {
            l = i + 2;
            rv1[i] = scale * g;
            g = 0.0;
            s = 0.0;
            scale = 0.0;
            if (i < IM)
            {
                for (k = i; k < IM; k++)
                {
                    scale += fabs(u(k, i));
                }
                if (scale > EPS_SP)
                {
                    for (k = i; k < IM; k++)
                    {
                        u(k, i) /= scale;
                        s += u(k, i) * u(k, i);
                    }
                    f = u(i, i);
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    u(i, i) = f - g;
                    for (j = l - 1; j < IN; j++)
                    {
                        for (s = 0.0, k = i; k < IM; k++)
                        {
                            s += u(k, i) * u(k, j);
                        }
                        f = s / h;
                        for (k = i; k < IM; k++)
                        {
                            u(k, j) += f * u(k, i);
                        }
                    }
                    for (k = i; k < IM; k++)
                    {
                        u(k, i) *= scale;
                    }
                }
            }
            w[i] = scale * g;
            g = s = scale = 0.0;
            if (i + 1 <= IM && i + 1 != N)
            {
                for (k = l - 1; k < IN; k++)
                {
                    scale += fabs(u(i, k));
                }
                if (scale > EPS_SP)
                {
                    for (k = l - 1; k < IN; k++)
                    {
                        u(i, k) /= scale;
                        s += u(i, k) * u(i, k);
                    }
                    f = u(i, l - 1);
                    g = -SIGN(sqrt(s), f);
                    h = f * g - s;
                    u(i, l - 1) = f - g;
                    for (k = l - 1; k < IN; k++)
                    {
                        rv1[k] = u(i, k) / h;
                    }
                    for (j = l - 1; j < IM; j++)
                    {
                        for (s = 0.0, k = l - 1; k < IN; k++)
                        {
                            s += u(j, k) * u(i, k);
                        }
                        for (k = l - 1; k < IN; k++)
                        {
                            u(j, k) += s * rv1[k];
                        }
                    }
                    for (k = l - 1; k < IN; k++)
                    {
                        u(i, k) *= scale;
                    }
                }
            }
            anorm = std::max(anorm, fabs(w[i]) + fabs(rv1[i]));
        }
        for (i = IN - 1; i >= 0; i--)
        {
            if (i < IN - 1)
            {
                if (fabs(g) > EPS_SP)
                {
                    for (j = l; j < IN; j++)
                    {
                        v(j, i) = (u(i, j) / u(i, l)) / g;
                    }
                    for (j = l; j < IN; j++)
                    {
                        for (s = 0.0, k = l; k < IN; k++)
                        {
                            s += u(i, k) * v(k, j);
                        }
                        for (k = l; k < IN; k++)
                        {
                            v(k, j) += s * v(k, i);
                        }
                    }
                }
                for (j = l; j < IN; j++)
                {
                    v(i, j) = 0.0;
                    v(j, i) = 0.0;
                }
            }
            v(i, i) = 1.0;
            g = rv1[i];
            l = i;
        }
        for (i = std::min(IM, IN) - 1; i >= 0; i--)
        {
            l = i + 1;
            g = w[i];
            for (j = l; j < IN; j++)
            {
                u(i, j) = 0.0;
            }
            if (fabs(g) > EPS_SP)
            {
                g = 1.0 / g;
                for (j = l; j < IN; j++)
                {
                    for (s = 0.0, k = l; k < IM; k++)
                    {
                        s += u(k, i) * u(k, j);
                    }
                    f = (s / u(i, i)) * g;
                    for (k = i; k < IM; k++)
                    {
                        u(k, j) += f * u(k, i);
                    }
                }
                for (j = i; j < IM; j++)
                {
                    u(j, i) *= g;
                }
            }
            else
            {
                for (j = i; j < IM; j++)
                {
                    u(j, i) = 0.0;
                }
            }
            u(i, i) += 1.0;
        }
        for (k = IN - 1; k >= 0; k--)
        {
            for (its = 0; its < 30; its++)
            {
                flag = true;
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
                        for (j = 0; j < IM; j++)
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
                        for (j = 0; j < IN; j++)
                        {
                            v(j, k) = -v(j, k);
                        }
                    }
                    break;
                }
                if (its == 29)
                {
                    sing = true;
                    return {};
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
                    for (jj = 0; jj < IN; jj++)
                    {
                        x = v(jj, j);
                        z = v(jj, i);
                        v(jj, j) = x * c + z * s;
                        v(jj, i) = z * c - x * s;
                    }
                    z = sqrt(f * f + h * h);
                    w[j] = z;
                    if (z)
                    {
                        z = 1.0 / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for (jj = 0; jj < IM; jj++)
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
    void svbksb(const MatrixS<M, N> &u, const MatrixS<N, 1> &w, const MatrixS<N, N> &v, double *b)
    {
        constexpr int IN = N;
        constexpr int IM = M;
        double tmp[N];
        auto eigen_max = *std::max_element(w.cbegin(), w.cend());
        auto tsh = 0.5 * sqrt(M + N + 1) * eigen_max * EPS_SP;
        for (int j = 0; j < IN; j++)
        {
            auto s = 0.0;
            if (w[j] > tsh)
            {
                for (int i = 0; i < IM; i++)
                {
                    s += u(i, j) * b[i];
                }
                s /= w[j];
            }
            tmp[j] = s;
        }
        for (int j = 0; j < IN; j++)
        {
            auto s = 0.0;
            for (int jj = 0; jj < IN; jj++)
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
        MatrixS<N, 1> c;
        MatrixS<N, 1> d;
        auto sing = false;
        auto R = qrdcmp(A, c, d, sing);
        if (!sing)
        {
            qrsolv(R, c, d, b.data());
        }
        StatusCode s = sing ? StatusCode::DIVERGED : StatusCode::NORMAL;
        return {slice<N, 1>(b, 0, 0), s};
    }

    template <Factorization type, size_t M, size_t N>
    std::enable_if_t<type == Factorization::SVD, EqnResult<N>>
    linsolve(const MatrixS<M, N> &A, MatrixS<M, 1> b)
    {
        MatrixS<N, 1> w{};
        MatrixS<N, N> V{};
        auto sing = false;
        auto U = svdcmp(A, w, V, sing);
        if (!sing)
        {
            svbksb(U, w, V, b.data());
        }
        StatusCode s = sing ? StatusCode::DIVERGED : StatusCode::NORMAL;
        return {slice<N, 1>(b, 0, 0), s};
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

    // eigen system
    namespace details
    {
        template <size_t N>
        void tred2_no_vec(MatrixS<N, N> z, MatrixS<N, 1> &d, MatrixS<N, 1> &e)
        {
            constexpr int IN = N;
            for (int i = IN - 1; i > 0; i--)
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
            e[0] = 0.0;
            for (int i = 0; i < IN; i++)
            {
                d[i] = z(i, i);
            }
        }

        template <size_t N>
        void tliq_no_vec(MatrixS<N, 1> &d, MatrixS<N, 1> &e)
        {
            constexpr int IN = N;
            // constexpr double EPS_SP = std::numeric_limits<double>::epsilon();
            for (int i = 1; i < IN; i++)
            {
                e[i - 1] = e[i];
            }
            e[N - 1] = 0.0;
            for (int l = 0; l < IN; l++)
            {
                int iter = 0;
                int m = 0;
                do
                {
                    for (m = l; m < IN - 1; m++)
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
        MatrixS<N, N> tred2_with_vec(MatrixS<N, N> z, MatrixS<N, 1> &d, MatrixS<N, 1> &e)
        {
            constexpr int IN = N;
            for (int i = IN - 1; i > 0; i--)
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
                        double g = (f > 0.0 ? -sqrt(h) : sqrt(h));
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
            for (int i = 0; i < IN; i++)
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
        void tliq_with_vec(MatrixS<N, N> &z, MatrixS<N, 1> &d, MatrixS<N, 1> &e)
        {
            constexpr int IN = N;
            // constexpr double EPS_SP = std::numeric_limits<double>::epsilon();
            for (int i = 1; i < IN; i++)
            {
                e[i - 1] = e[i];
            }
            e[IN - 1] = 0.0;
            for (int l = 0; l < IN; l++)
            {
                int iter = 0;
                int m = 0;
                do
                {
                    for (m = l; m < IN - 1; m++)
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
                            for (int k = 0; k < IN; k++)
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
            constexpr int IN = N;
            for (size_t i = 0; i < IN - 1; i++)
            {
                auto j = std::distance(vec.begin() + i, std::min_element(vec.begin() + i, vec.end())) + i;
                if (j != i)
                {
                    std::swap(vec[i], vec[j]);
                    for (size_t k = 0; k < N; k++)
                    {
                        std::swap(mat(k, i), mat(k, j));
                    }
                }
            }
        }
        /*
        template <size_t N>
        struct Unsymmeig
        {
            using complexd = std::complex<double>;
            constexpr static int n = N;
            Matrix<N, N> a, zz;
            Matrix<N, 1> scale;
            std::array<int, N> perm;
            std::array<complexd, N> wri;
            bool yesvecs, hessen;

            void balance()
            {
                constexpr double RADIX = std::numeric_limits<double>::radix;
                bool done = false;
                double sqrdx = RADIX * RADIX;
                while (!done)
                {
                    done = true;
                    for (int i = 0; i < n; i++)
                    {
                        double r = 0.0, c = 0.0;
                        for (int j = 0; j < n; j++)
                            if (j != i)
                            {
                                c += fabs(a(j, i));
                                r += fabs(a(i, j));
                            }
                        if (c != 0.0 && r != 0.0)
                        {
                            double g = r / RADIX;
                            double f = 1.0;
                            double s = c + r;
                            while (c < g)
                            {
                                f *= RADIX;
                                c *= sqrdx;
                            }
                            g = r * RADIX;
                            while (c > g)
                            {
                                f /= RADIX;
                                c /= sqrdx;
                            }
                            if ((c + r) / f < 0.95 * s)
                            {
                                done = false;
                                g = 1.0 / f;
                                scale[i] *= f;
                                for (int j = 0; j < n; j++)
                                {
                                    a(i, j) *= g;
                                }
                                for (int j = 0; j < n; j++)
                                {
                                    a(j, i) *= f;
                                }
                            }
                        }
                    }
                }
            }

            void elmhes()
            {
                for (int m = 1; m < n - 1; m++)
                {
                    double x = 0.0;
                    int i = m;
                    for (int j = m; j < n; j++)
                    {
                        if (fabs(a(j, m - 1)) > fabs(x))
                        {
                            x = a(j, m - 1);
                            i = j;
                        }
                    }
                    perm[m] = i;
                    if (i != m)
                    {
                        for (int j = m - 1; j < n; j++)
                            std::swap(a(i, j), a(m, j));
                        for (int j = 0; j < n; j++)
                            std::swap(a(j, i), a(j, m));
                    }
                    if (x != 0.0)
                    {
                        for (i = m + 1; i < n; i++)
                        {
                            double y = a(i, m - 1);
                            if (y != 0.0)
                            {
                                y /= x;
                                a(i, m - 1) = y;
                                for (int j = m; j < n; j++)
                                    a(i, j) -= y * a(m, j);
                                for (int j = 0; j < n; j++)
                                    a(j, m) += y * a(j, i);
                            }
                        }
                    }
                }
            }

            void eltran()
            {
                for (int mp = n - 2; mp > 0; mp--)
                {
                    for (int k = mp + 1; k < n; k++)
                        zz(k, mp) = a(k, mp - 1);
                    int i = perm[mp];
                    if (i != mp)
                    {
                        for (int j = mp; j < n; j++)
                        {
                            zz(mp, j) = zz(i, j);
                            zz(i, j) = 0.0;
                        }
                        zz(i, mp) = 1.0;
                    }
                }
            }

            void hqr()
            {
                int nn, m, l, k, j, its, i, mmin;
                double z, y, x, w, v, u, t, s, r, q, p, anorm = 0.0;
                for (i = 0; i < n; i++)
                    for (j = std::max(i - 1, 0); j < n; j++)
                        anorm += fabs(a(i, j));
                nn = n - 1;
                t = 0.0;
                while (nn >= 0)
                {
                    its = 0;
                    do
                    {
                        for (l = nn; l > 0; l--)
                        {
                            s = fabs(a(l - 1, l - 1)) + fabs(a(l, l));
                            if (s == 0.0)
                                s = anorm;
                            if (fabs(a(l, l - 1)) <= EPS_DP * s)
                            {
                                a(l, l - 1) = 0.0;
                                break;
                            }
                        }
                        x = a(nn, nn);
                        if (l == nn)
                        {
                            wri[nn--] = x + t;
                        }
                        else
                        {
                            y = a(nn - 1, nn - 1);
                            w = a(nn, nn - 1) * a(nn - 1, nn);
                            if (l == nn - 1)
                            {
                                p = 0.5 * (y - x);
                                q = p * p + w;
                                z = sqrt(fabs(q));
                                x += t;
                                if (q >= 0.0)
                                {
                                    z = p + SIGN(z, p);
                                    wri[nn - 1] = wri[nn] = x + z;
                                    if (z != 0.0)
                                        wri[nn] = x - w / z;
                                }
                                else
                                {
                                    wri[nn] = complexd(x + p, -z);
                                    wri[nn - 1] = conj(wri[nn]);
                                }
                                nn -= 2;
                            }
                            else
                            {
                                if (its == 30)
                                    throw("Too many iterations in hqr");
                                if (its == 10 || its == 20)
                                {
                                    t += x;
                                    for (i = 0; i < nn + 1; i++)
                                        a(i, i) -= x;
                                    s = fabs(a(nn, nn - 1)) + fabs(a(nn - 1, nn - 2));
                                    y = x = 0.75 * s;
                                    w = -0.4375 * s * s;
                                }
                                ++its;
                                for (m = nn - 2; m >= l; m--)
                                {
                                    z = a(m, m);
                                    r = x - z;
                                    s = y - z;
                                    p = (r * s - w) / a(m + 1, m) + a(m, m + 1);
                                    q = a(m + 1, m + 1) - z - r - s;
                                    r = a(m + 2, m + 1);
                                    s = fabs(p) + fabs(q) + fabs(r);
                                    p /= s;
                                    q /= s;
                                    r /= s;
                                    if (m == l)
                                        break;
                                    u = fabs(a(m, m - 1)) * (fabs(q) + fabs(r));
                                    v = fabs(p) * (fabs(a(m - 1, m - 1)) + fabs(z) + fabs(a(m + 1, m + 1)));
                                    if (u <= EPS_DP * v)
                                        break;
                                }
                                for (i = m; i < nn - 1; i++)
                                {
                                    a(i + 2, i) = 0.0;
                                    if (i != m)
                                        a(i + 2, i - 1) = 0.0;
                                }
                                for (k = m; k < nn; k++)
                                {
                                    if (k != m)
                                    {
                                        p = a(k, k - 1);
                                        q = a(k + 1, k - 1);
                                        r = 0.0;
                                        if (k + 1 != nn)
                                            r = a(k + 2, k - 1);
                                        if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0)
                                        {
                                            p /= x;
                                            q /= x;
                                            r /= x;
                                        }
                                    }
                                    if ((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0)
                                    {
                                        if (k == m)
                                        {
                                            if (l != m)
                                                a(k, k - 1) *= -1.0;
                                        }
                                        else
                                            a(k, k - 1) = -s * x;
                                        p += s;
                                        x = p / s;
                                        y = q / s;
                                        z = r / s;
                                        q /= p;
                                        r /= p;
                                        for (j = k; j < nn + 1; j++)
                                        {
                                            p = a(k, j) + q * a(k + 1, j);
                                            if (k + 1 != nn)
                                            {
                                                p += r * a(k + 2, j);
                                                a(k + 2, j) -= p * z;
                                            }
                                            a(k + 1, j) -= p * y;
                                            a(k, j) -= p * x;
                                        }
                                        mmin = nn < k + 3 ? nn : k + 3;
                                        for (i = l; i < mmin + 1; i++)
                                        {
                                            p = x * a(i, k) + y * a(i, k + 1);
                                            if (k + 1 != nn)
                                            {
                                                p += z * a(i, k + 2);
                                                a(i, k + 2) -= p * r;
                                            }
                                            a(i, k + 1) -= p * q;
                                            a(i, k) -= p;
                                        }
                                    }
                                }
                            }
                        }
                    } while (l + 1 < nn);
                }
            }

            void hqr2()
            {
                int nn, m, l, k, j, its, i, mmin, na;
                double z, y, x, w, v, u, t, s, r, q, p, anorm = 0.0, ra, sa, vr, vi;
                for (i = 0; i < n; i++)
                    for (j = std::max(i - 1, 0); j < n; j++)
                        anorm += fabs(a(i, j));
                nn = n - 1;
                t = 0.0;
                while (nn >= 0)
                {
                    its = 0;
                    do
                    {
                        for (l = nn; l > 0; l--)
                        {
                            s = fabs(a(l - 1, l - 1)) + fabs(a(l, l));
                            if (s == 0.0)
                                s = anorm;
                            if (fabs(a(l, l - 1)) <= EPS_DP * s)
                            {
                                a(l, l - 1) = 0.0;
                                break;
                            }
                        }
                        x = a(nn, nn);
                        if (l == nn)
                        {
                            wri[nn] = a(nn, nn) = x + t;
                            nn--;
                        }
                        else
                        {
                            y = a(nn - 1, nn - 1);
                            w = a(nn, nn - 1) * a(nn - 1, nn);
                            if (l == nn - 1)
                            {
                                p = 0.5 * (y - x);
                                q = p * p + w;
                                z = sqrt(fabs(q));
                                x += t;
                                a(nn, nn) = x;
                                a(nn - 1, nn - 1) = y + t;
                                if (q >= 0.0)
                                {
                                    z = p + SIGN(z, p);
                                    wri[nn - 1] = wri[nn] = x + z;
                                    if (z != 0.0)
                                        wri[nn] = x - w / z;
                                    x = a(nn, nn - 1);
                                    s = fabs(x) + fabs(z);
                                    p = x / s;
                                    q = z / s;
                                    r = sqrt(p * p + q * q);
                                    p /= r;
                                    q /= r;
                                    for (j = nn - 1; j < n; j++)
                                    {
                                        z = a(nn - 1, j);
                                        a(nn - 1, j) = q * z + p * a(nn, j);
                                        a(nn, j) = q * a(nn, j) - p * z;
                                    }
                                    for (i = 0; i <= nn; i++)
                                    {
                                        z = a(i, nn - 1);
                                        a(i, nn - 1) = q * z + p * a(i, nn);
                                        a(i, nn) = q * a(i, nn) - p * z;
                                    }
                                    for (i = 0; i < n; i++)
                                    {
                                        z = zz(i, nn - 1);
                                        zz(i, nn - 1) = q * z + p * zz(i, nn);
                                        zz(i, nn) = q * zz(i, nn) - p * z;
                                    }
                                }
                                else
                                {
                                    wri[nn] = complexd(x + p, -z);
                                    wri[nn - 1] = conj(wri[nn]);
                                }
                                nn -= 2;
                            }
                            else
                            {
                                if (its == 30)
                                    throw("Too many iterations in hqr");
                                if (its == 10 || its == 20)
                                {
                                    t += x;
                                    for (i = 0; i < nn + 1; i++)
                                        a(i, i) -= x;
                                    s = fabs(a(nn, nn - 1)) + fabs(a(nn - 1, nn - 2));
                                    y = x = 0.75 * s;
                                    w = -0.4375 * s * s;
                                }
                                ++its;
                                for (m = nn - 2; m >= l; m--)
                                {
                                    z = a(m, m);
                                    r = x - z;
                                    s = y - z;
                                    p = (r * s - w) / a(m + 1, m) + a(m, m + 1);
                                    q = a(m + 1, m + 1) - z - r - s;
                                    r = a(m + 2, m + 1);
                                    s = fabs(p) + fabs(q) + fabs(r);
                                    p /= s;
                                    q /= s;
                                    r /= s;
                                    if (m == l)
                                        break;
                                    u = fabs(a(m, m - 1)) * (fabs(q) + fabs(r));
                                    v = fabs(p) * (fabs(a(m - 1, m - 1)) + fabs(z) + fabs(a(m + 1, m + 1)));
                                    if (u <= EPS_DP * v)
                                        break;
                                }
                                for (i = m; i < nn - 1; i++)
                                {
                                    a(i + 2, i) = 0.0;
                                    if (i != m)
                                        a(i + 2, i - 1) = 0.0;
                                }
                                for (k = m; k < nn; k++)
                                {
                                    if (k != m)
                                    {
                                        p = a(k, k - 1);
                                        q = a(k + 1, k - 1);
                                        r = 0.0;
                                        if (k + 1 != nn)
                                            r = a(k + 2, k - 1);
                                        if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0)
                                        {
                                            p /= x;
                                            q /= x;
                                            r /= x;
                                        }
                                    }
                                    if ((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0)
                                    {
                                        if (k == m)
                                        {
                                            if (l != m)
                                                a(k, k - 1) = -a(k, k - 1);
                                        }
                                        else
                                            a(k, k - 1) = -s * x;
                                        p += s;
                                        x = p / s;
                                        y = q / s;
                                        z = r / s;
                                        q /= p;
                                        r /= p;
                                        for (j = k; j < n; j++)
                                        {
                                            p = a(k, j) + q * a(k + 1, j);
                                            if (k + 1 != nn)
                                            {
                                                p += r * a(k + 2, j);
                                                a(k + 2, j) -= p * z;
                                            }
                                            a(k + 1, j) -= p * y;
                                            a(k, j) -= p * x;
                                        }
                                        mmin = nn < k + 3 ? nn : k + 3;
                                        for (i = 0; i < mmin + 1; i++)
                                        {
                                            p = x * a(i, k) + y * a(i, k + 1);
                                            if (k + 1 != nn)
                                            {
                                                p += z * a(i, k + 2);
                                                a(i, k + 2) -= p * r;
                                            }
                                            a(i, k + 1) -= p * q;
                                            a(i, k) -= p;
                                        }
                                        for (i = 0; i < n; i++)
                                        {
                                            p = x * zz(i, k) + y * zz(i, k + 1);
                                            if (k + 1 != nn)
                                            {
                                                p += z * zz(i, k + 2);
                                                zz(i, k + 2) -= p * r;
                                            }
                                            zz(i, k + 1) -= p * q;
                                            zz(i, k) -= p;
                                        }
                                    }
                                }
                            }
                        }
                    } while (l + 1 < nn);
                }
                if (anorm != 0.0)
                {
                    for (nn = n - 1; nn >= 0; nn--)
                    {
                        p = real(wri[nn]);
                        q = imag(wri[nn]);
                        na = nn - 1;
                        if (q == 0.0)
                        {
                            m = nn;
                            a(nn, nn) = 1.0;
                            for (i = nn - 1; i >= 0; i--)
                            {
                                w = a(i, i) - p;
                                r = 0.0;
                                for (j = m; j <= nn; j++)
                                    r += a(i, j) * a(j, nn);
                                if (imag(wri[i]) < 0.0)
                                {
                                    z = w;
                                    s = r;
                                }
                                else
                                {
                                    m = i;

                                    if (imag(wri[i]) == 0.0)
                                    {
                                        t = w;
                                        if (t == 0.0)
                                            t = EPS_DP * anorm;
                                        a(i, nn) = -r / t;
                                    }
                                    else
                                    {
                                        x = a(i, i + 1);
                                        y = a(i + 1, i);
                                        q = SQR(real(wri[i]) - p) + SQR(imag(wri[i]));
                                        t = (x * s - z * r) / q;
                                        a(i, nn) = t;
                                        if (fabs(x) > fabs(z))
                                            a(i + 1, nn) = (-r - w * t) / x;
                                        else
                                            a(i + 1, nn) = (-s - y * t) / z;
                                    }
                                    t = fabs(a(i, nn));
                                    if (EPS_DP * t * t > 1)
                                        for (j = i; j <= nn; j++)
                                            a(j, nn) /= t;
                                }
                            }
                        }
                        else if (q < 0.0)
                        {
                            m = na;
                            if (fabs(a(nn, na)) > fabs(a(na, nn)))
                            {
                                a(na, na) = q / a(nn, na);
                                a(na, nn) = -(a(nn, nn) - p) / a(nn, na);
                            }
                            else
                            {
                                auto temp = complexd(0.0, -a(na, nn)) / complexd(a(na, na) - p, q);
                                a(na, na) = temp.real();
                                a(na, nn) = temp.imag();
                            }
                            a(nn, na) = 0.0;
                            a(nn, nn) = 1.0;
                            for (i = nn - 2; i >= 0; i--)
                            {
                                w = a(i, i) - p;
                                ra = sa = 0.0;
                                for (j = m; j <= nn; j++)
                                {
                                    ra += a(i, j) * a(j, na);
                                    sa += a(i, j) * a(j, nn);
                                }
                                if (imag(wri[i]) < 0.0)
                                {
                                    z = w;
                                    r = ra;
                                    s = sa;
                                }
                                else
                                {
                                    m = i;
                                    if (imag(wri[i]) == 0.0)
                                    {
                                        auto temp = complexd(-ra, -sa) / complexd(w, q);
                                        a(i, na) = temp.real();
                                        a(i, nn) = temp.imag();
                                    }
                                    else
                                    {
                                        x = a(i, i + 1);
                                        y = a(i + 1, i);
                                        vr = SQR(real(wri[i]) - p) + SQR(imag(wri[i])) - q * q;
                                        vi = 2.0 * q * (real(wri[i]) - p);
                                        if (vr == 0.0 && vi == 0.0)
                                            vr = EPS_DP * anorm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(z));
                                        auto temp = complexd(x * r - z * ra + q * sa, x * s - z * sa - q * ra) /
                                                    complexd(vr, vi);
                                        a(i, na) = temp.real();
                                        a(i, nn) = temp.imag();
                                        if (fabs(x) > fabs(z) + fabs(q))
                                        {
                                            a(i + 1, na) = (-ra - w * a(i, na) + q * a(i, nn)) / x;
                                            a(i + 1, nn) = (-sa - w * a(i, nn) - q * a(i, na)) / x;
                                        }
                                        else
                                        {
                                            auto temp = complexd(-r - y * a(i, na), -s - y * a(i, nn)) /
                                                        complexd(z, q);
                                            a(i + 1, na) = temp.real();
                                            a(i + 1, nn) = temp.imag();
                                        }
                                    }
                                }
                                t = std::max(fabs(a(i, na)), fabs(a(i, nn)));
                                if (EPS_DP * t * t > 1)
                                    for (j = i; j <= nn; j++)
                                    {
                                        a(j, na) /= t;
                                        a(j, nn) /= t;
                                    }
                            }
                        }
                    }
                    for (j = n - 1; j >= 0; j--)
                        for (i = 0; i < n; i++)
                        {
                            z = 0.0;
                            for (k = 0; k <= j; k++)
                                z += zz(i, k) * a(k, j);
                            zz(i, j) = z;
                        }
                }
            }

            void balbak()
            {
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        zz(i, j) *= scale[i];
            }

            void sort()
            {
                int i;
                for (int j = 1; j < n; j++)
                {
                    auto x = wri[j];
                    for (i = j - 1; i >= 0; i--)
                    {
                        if (real(wri[i]) >= real(x))
                            break;
                        wri[i + 1] = wri[i];
                    }
                    wri[i + 1] = x;
                }
            }

            void sortvecs()
            {
                int i;
                Matrix<N, 1> temp;
                for (int j = 1; j < n; j++)
                {
                    complexd x = wri[j];
                    for (int k = 0; k < n; k++)
                        temp[k] = zz(k, j);
                    for (i = j - 1; i >= 0; i--)
                    {
                        if (real(wri[i]) >= real(x))
                            break;
                        wri[i + 1] = wri[i];
                        for (int k = 0; k < n; k++)
                            zz(k, i + 1) = zz(k, i);
                    }
                    wri[i + 1] = x;
                    for (int k = 0; k < n; k++)
                        zz(k, i + 1) = temp[k];
                }
            }

            void operator()(const Matrix<N, N> &aa, bool yesvec = true, bool hessenb = false)
            {
                a = aa;
                yesvecs = yesvec;
                hessen = hessenb;
                ones(scale);
                balance();
                if (!hessen)
                    elmhes();
                if (yesvecs)
                {
                    for (size_t i = 0; i < n; i++)
                    {
                        zz(i, i) = 1.0;
                    }
                    if (!hessen)
                        eltran();
                    hqr2();
                    balbak();
                    sortvecs();
                }
                else
                {
                    hqr();
                    sort();
                }
            }
        };*/
    }

    template <EigenSystem type, size_t N>
    std::enable_if_t<type == EigenSystem::SymValAndVec, EigResult<N>>
    eig(const MatrixS<N, N> &mat)
    {
        MatrixS<N, 1> e, eig_value;
        auto result = details::tred2_with_vec(mat, eig_value, e);
        details::tliq_with_vec(result, eig_value, e);
        return {result, eig_value};
    }

    template <EigenSystem type, size_t N>
    std::enable_if_t<type == EigenSystem::SymValAndVecSorted, EigResult<N>>
    eig(const MatrixS<N, N> &mat)
    {
        MatrixS<N, 1> e, eig_value;
        auto result = details::tred2_with_vec(mat, eig_value, e);
        details::tliq_with_vec(result, eig_value, e);
        details::eigsrt(result, eig_value);
        return {result, eig_value};
    }

    template <EigenSystem type, size_t N>
    std::enable_if_t<type == EigenSystem::SymOnlyVal, MatrixS<N, 1>>
    eig(const MatrixS<N, N> &mat)
    {
        MatrixS<N, 1> e, eig_value;
        details::tred2_no_vec(mat, eig_value, e);
        details::tliq_no_vec(eig_value, e);
        std::sort(eig_value.begin(), eig_value.end());
        return eig_value;
    }

}
#endif