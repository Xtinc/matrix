#ifndef VVERY_SIMPLE_OPTIMIZATION_HEADER
#define VVERY_SIMPLE_OPTIMIZATION_HEADER

#include "linalg.hpp"

namespace ppx
{
    enum class Optimization : char
    {
        GoldenSearch,
        QuadraticSearch,
        BrentSearch,
        Powell,
        GradientDescent,
        ConjugateGradient,
        BGFS
    };

    template <size_t N>
    struct OptResult
    {
        MatrixS<N, 1> x;
        double y = 0.0;
        StatusCode s = StatusCode::NORMAL;
    };

    template <>
    struct OptResult<1>
    {
        double x = 0.0;
        double y = 0.0;
        StatusCode s = StatusCode::NORMAL;
    };

    template <size_t N>
    std::ostream &operator<<(std::ostream &os, const OptResult<N> &self)
    {
        os << "OptResult<" << N << ">:\n"
           << "Status:\t" << self.s << "\n"
           << "xmin =\t" << self.x << "\n"
           << "ymin =\t" << self.y << std::endl;
        return os;
    }

    namespace details
    {
        // 1D
        struct Univariate
        {
        protected:
            double a;
            double b;
            double c;
            double ya;
            double yb;
            double yc;

        public:
            void already_bounded(double x0, double x1)
            {
                a = std::min(x0, x1);
                c = std::max(x0, x1);
                b = 0.5 * (a + c);
            }

            template <typename T>
            void bracket_minimum(const T &func, double x0)
            {
                // deal with monotony?
                constexpr double k = 2.0;
                constexpr size_t itmax = 1000;
                double s = 0.05;
                a = x0;
                ya = func(a);
                b = a + s;
                yb = func(b);
                if (yb > ya)
                {
                    std::swap(a, b);
                    std::swap(ya, yb);
                    s = -s;
                }
                for (size_t its = 0; its < itmax; ++its)
                {
                    c = b + s;
                    yc = func(c);
                    if (yc > yb)
                    {
                        if (a > c)
                        {
                            std::swap(a, c);
                            std::swap(ya, yc);
                        }
                        return;
                    }
                    a = b;
                    ya = yb;
                    b = c;
                    yb = yc;
                    s *= k;
                }
            }
        };

        struct GoldenSearch : Univariate
        {
            OptResult<1> R;
            size_t ITMAX = 200;

            template <typename T>
            OptResult<1> operator()(const T &func)
            {
                constexpr double rho = 0.6180339887;
                constexpr double phi = 1.0 - rho;
                auto d = rho * c + phi * a;
                auto yd = func(d);
                auto its = 0u;
                while (fabs(a - c) > EPS_SP && its < ITMAX)
                {
                    R.x = rho * a + phi * c;
                    R.y = func(R.x);
                    if (R.y < yd)
                    {
                        c = d;
                        d = R.x;
                        yd = R.y;
                    }
                    else
                    {
                        a = c;
                        c = R.x;
                    }
                    ++its;
                    // printf("iter = %d, xc = %g, residual = %g\n", its, xmin, fabs(a - b));
                }
                R.s = its == ITMAX ? StatusCode::DIVERGED : StatusCode::CONVERGED;
                return R;
            }
        };

        struct QuadraticSearch : Univariate
        {
            OptResult<1> R;
            size_t ITMAX = 200;

            template <typename T>
            OptResult<1> operator()(const T &func)
            {
                ya = func(a);
                yb = func(b);
                yc = func(c);
                auto its = 0u;
                while (fabs(a - c) > EPS_SP && fabs(b - c) > EPS_SP && fabs(a - b) > EPS_SP && its < ITMAX)
                {
                    R.x = 0.5 * (ya * (b * b - c * c) + yb * (c * c - a * a) + yc * (a * a - b * b)) /
                          (ya * (b - c) + yb * (c - a) + yc * (a - b));
                    R.y = func(R.x);
                    if (R.x > b)
                    {
                        if (R.y > yb)
                        {
                            c = R.x;
                            yc = R.y;
                        }
                        else
                        {
                            a = b;
                            ya = yb;
                            b = R.x;
                            yb = R.y;
                        }
                    }
                    else
                    {
                        if (R.y > yb)
                        {
                            a = R.x;
                            ya = R.y;
                        }
                        else
                        {
                            c = b;
                            yc = yb;
                            b = R.x;
                            yb = R.y;
                        }
                    }
                    // printf("%g %g\n", xmin, ymin);
                    ++its;
                }
                R.s = its == ITMAX ? StatusCode::DIVERGED : StatusCode::CONVERGED;
                return R;
            }
        };

        struct BrentSearch : Univariate
        {
            OptResult<1> R;
            size_t ITMAX = 200;

            template <typename T>
            OptResult<1> operator()(const T &func)
            {
                constexpr double rho = 0.6180339887;
                constexpr double phi = 1.0 - rho;
                auto x = c;
                auto w = c;
                auto v = c;
                auto fx = func(x);
                auto fw = fx;
                auto fv = fx;
                for (auto its = 0u; its < ITMAX; its++)
                {
                    double d = 0.0;
                    double e = 0.0;
                    double u = 0.0;
                    auto xm = 0.5 * (a + c);
                    auto tol1 = EPS_SP * fabs(x);
                    auto tol2 = 2.0 * tol1;
                    if (fabs(x - xm) < tol2 - 0.5 * (c - a))
                    {
                        R.y = fx;
                        R.x = x;
                        R.s = StatusCode::CONVERGED;
                        return R;
                    }
                    if (fabs(e) > tol1)
                    {
                        auto r = (x - w) * (fx - fv);
                        auto q = (x - v) * (fx - fw);
                        auto p = (x - v) * q - (x - w) * r;
                        q = 2.0 * (q - r);
                        if (q > 0.0)
                        {
                            p = -p;
                        }
                        q = fabs(q);
                        r = e;
                        e = d;
                        if (fabs(p) > fabs(0.5 * q * r) ||
                            p < q * (a - x) || p > q * (c - x))
                        {
                            e = x > xm ? a - x : c - x;
                            d = phi * e;
                        }
                        else
                        {
                            d = p / q;
                            u = x + d;
                            if (u - a < tol2 || c - u < tol2)
                            {
                                d = SIGN(tol1, xm - x);
                            }
                        }
                    }
                    else
                    {
                        e = x > xm ? a - x : c - x;
                        d = phi * e;
                    }
                    u = fabs(d) > tol1 ? x + d : x + SIGN(tol1, d);
                    auto fu = func(u);
                    if (fu < fx)
                    {
                        if (u > x)
                        {
                            a = x;
                        }
                        else
                        {
                            c = x;
                        }
                        v = w;
                        w = x;
                        x = u;
                        fv = fw;
                        fw = fx;
                        fx = fu;
                    }
                    else
                    {
                        if (u < x)
                        {
                            a = u;
                        }
                        else
                        {
                            c = u;
                        }
                        if (fu < fw || details::is_same(w, x))
                        {
                            v = w;
                            w = u;
                            fv = fw;
                            fw = fu;
                        }
                        else if (fu < fv || details::is_same(v, x) || details::is_same(v, w))
                        {
                            v = u;
                            fv = fu;
                        }
                    }
                }
                R.y = fx;
                R.x = x;
                R.s = StatusCode::DIVERGED;
                return R;
            }
        };

        template <typename T, size_t N>
        OptResult<1> lnsrch(const T &f, const MatrixS<N, 1> &x, const MatrixS<N, 1> &d)
        {
            auto func = [&f, &x, &d](double a)
            {
                return f(x + a * d);
            };
            BrentSearch brent;
            brent.bracket_minimum(func, 0.0);
            brent(func);
            return brent.R;
        }

        // nD
        template <size_t N>
        struct GradientDescent
        {
            OptResult<N> R;
            size_t ITMAX = 200;
            double FTOLA = EPS_SP;

            template <typename T1, typename T2>
            OptResult<N> operator()(const T1 &f, const T2 &df, const MatrixS<N, 1> &x)
            {
                R.x = x;
                MatrixS<N, 1> dx;
                dx.fill(1.0);
                auto its = 0u;
                while (norminf(dx) > FTOLA && its < ITMAX)
                {
                    dx = -1 * df(R.x);
                    auto lr = lnsrch(f, R.x, dx);
                    if (lr.s == StatusCode::DIVERGED)
                    {
                        break;
                    }
                    dx *= lr.x;
                    R.x += dx;
                    R.y = lr.y;
                    // printf("iter = %d, xerr = %g\n", its++, norminf(dx));
                }
                R.s = its == ITMAX ? StatusCode::DIVERGED : StatusCode::CONVERGED;
                return R;
            }
        };

        template <size_t N>
        struct ConjuateGradient
        {
            OptResult<N> R;
            size_t ITMAX = 200;
            double FTOLA = EPS_SP;

            template <typename T1, typename T2>
            OptResult<N> operator()(const T1 &f, const T2 &df, const MatrixS<N, 1> &x)
            {
                R.x = x;
                auto its = 0u;
                MatrixS<N, 1> g = df(R.x);
                MatrixS<N, 1> d = -1.0 * g;
                MatrixS<N, 1> dh;
                dh.fill(1.0);
                while (norminf(dh) > FTOLA && its < ITMAX)
                {
                    MatrixS<N, 1> gh = df(R.x);
                    double beta =
                        std::max(0.0, inner_product<N>(gh, gh - g) / inner_product(g, g));
                    dh = -1 * gh + beta * d;
                    auto lr = lnsrch(f, R.x, dh);
                    if (lr.s == StatusCode::DIVERGED)
                    {
                        break;
                    }
                    d = dh;
                    g = gh;
                    dh *= lr.x;
                    R.x += dh;
                    R.y = lr.y;
                    // printf("iter = %d, xerr = %g\n", its++, xerr);
                }
                R.s = its == ITMAX ? StatusCode::DIVERGED : StatusCode::CONVERGED;
                return R;
            }
        };

        template <size_t N>
        struct BFGS
        {
            OptResult<N> R;
            size_t ITMAX = 20;
            double FTOLA = EPS_SP;

            template <typename T1, typename T2>
            OptResult<N> operator()(const T1 &f, const T2 &df, const MatrixS<N, 1> &x)
            {
                R.x = x;
                auto Q = MatrixS<N, N>::eye();
                MatrixS<N, 1> g = df(x);
                MatrixS<N, 1> dg;
                MatrixS<N, 1> dx;
                dg.fill(1.0);
                dx.fill(1.0);
                auto its = 0u;
                while (norm2(dx) > FTOLA && norm2(dg) > FTOLA && its < ITMAX)
                {
                    dx = Q * g * (-1.0);
                    auto lr = lnsrch(f, R.x, dx);
                    if (lr.s == StatusCode::DIVERGED)
                    {
                        break;
                    }
                    dx *= lr.x;
                    R.x += dx;
                    R.y = lr.y;
                    auto gh = df(R.x);
                    dg = gh - g;
                    double nominator = inner_product(dx, dg);
                    auto dgT = dg.T();
                    auto dxT = dx.T();
                    double coff = (1 + dgT * Q * dg / nominator)[0] / nominator;
                    Q -= (dx * dgT * Q + Q * dg * dxT) / nominator + coff * (dx * dxT);
                    // Matrix<N, 1> tp = dx - Q * dg;
                    // Q += tp * tp.T() / inner_product(dg, tp);
                    g = gh;
                    // printf("iter = %d, xerr = %g\n", its++, norm2(dx));
                }
                R.s = its == ITMAX ? StatusCode::DIVERGED : StatusCode::CONVERGED;
                return R;
            }
        };

        template <size_t N>
        struct Powell
        {
            OptResult<N> R;
            size_t ITMAX = 200;
            double FTOLA = EPS_SP;

            template <typename T>
            OptResult<N> operator()(const T &f, const MatrixS<N, 1> &x0)
            {
                R.x = x0;
                auto U = MatrixS<N, N>::eye();
                auto dx = 1.0;
                auto its = 0u;
                while (dx > FTOLA && its < ITMAX)
                {
                    auto x = R.x;
                    MatrixS<N, 1> d;
                    for (size_t i = 0; i < N; i++)
                    {
                        d = U.col(i);
                        auto lr = lnsrch(f, x, d);
                        // if (lr.s == StatusCode::DIVERGED)
                        // {
                        //     R.s = StatusCode::DIVERGED;
                        //     return;
                        // }
                        x += lr.x * d;
                    }
                    for (size_t i = 0; i < N - 1; i++)
                    {

                        U.template sub<N, 1>(0, i) = U.template sub<N, 1>(0, i + 1);
                        // U({-1, -1}, i) = U.col(i + 1);
                    }
                    d = x - R.x;
                    U.template sub<N, 1>(0, N - 1) = d;
                    // U({-1, -1}, N - 1) = d;
                    auto lr = lnsrch(f, x, d);
                    if (lr.s == StatusCode::DIVERGED)
                    {
                        break;
                    }
                    x += lr.x * d;
                    dx = norm2<N, 1>(x - R.x);
                    R.x = x;
                    R.y = lr.y;
                    // printf("iters = %d, residual = %g\n", its++, dx);
                }
                R.s = its == ITMAX ? StatusCode::DIVERGED : StatusCode::CONVERGED;
                return R;
            }
        };
    }

    template <Optimization etype, typename T>
    std::enable_if_t<etype == Optimization::GoldenSearch, OptResult<1>>
    fminbnd(const T &fn, double a, double b)
    {
        details::GoldenSearch gs;
        gs.already_bounded(a, b);
        return gs(fn);
    }

    template <Optimization etype, typename T>
    std::enable_if_t<etype == Optimization::QuadraticSearch, OptResult<1>>
    fminbnd(const T &fn, double a, double b)
    {
        details::QuadraticSearch qs;
        qs.already_bounded(a, b);
        return qs(fn);
    }

    template <Optimization etype, typename T>
    std::enable_if_t<etype == Optimization::BrentSearch, OptResult<1>>
    fminbnd(const T &fn, double a, double b)
    {
        details::BrentSearch brent;
        brent.already_bounded(a, b);
        return brent(fn);
    }

    template <Optimization etype, typename T>
    std::enable_if_t<etype == Optimization::GoldenSearch, OptResult<1>>
    fminunc(const T &fn, double a)
    {
        details::GoldenSearch gs;
        gs.bracket_minimum(fn, a);
        return gs(fn);
    }

    template <Optimization etype, typename T>
    std::enable_if_t<etype == Optimization::QuadraticSearch, OptResult<1>>
    fminunc(const T &fn, double a)
    {
        details::QuadraticSearch qs;
        qs.bracket_minimum(fn, a);
        return qs(fn);
    }

    template <Optimization etype, typename T>
    std::enable_if_t<etype == Optimization::BrentSearch, OptResult<1>>
    fminunc(const T &fn, double a)
    {
        details::BrentSearch brent;
        brent.bracket_minimum(fn, a);
        return brent(fn);
    }

    template <Optimization etype, size_t N, typename T>
    std::enable_if_t<etype == Optimization::Powell, OptResult<N>>
    fminunc(const T &func, const MatrixS<N, 1> &x0)
    {
        details::Powell<N> pw;
        return pw(func, x0);
    }

    template <Optimization etype, size_t N, typename T, typename T2>
    std::enable_if_t<etype == Optimization::GradientDescent, OptResult<N>>
    fminunc(const T &func, const T2 &dfunc, const MatrixS<N, 1> &x0)
    {
        details::GradientDescent<N> gd;
        return gd(func, dfunc, x0);
    }

    template <Optimization etype, size_t N, typename T, typename T2>
    std::enable_if_t<etype == Optimization::ConjugateGradient, OptResult<N>>
    fminunc(const T &func, const T2 &dfunc, const MatrixS<N, 1> &x0)
    {
        details::ConjuateGradient<N> cg;
        return cg(func, dfunc, x0);
    }

    template <Optimization etype, size_t N, typename T, typename T2>
    std::enable_if_t<etype == Optimization::BGFS, OptResult<N>>
    fminunc(const T &func, const T2 &dfunc, const MatrixS<N, 1> &x0)
    {
        details::BFGS<N> bfgs;
        return bfgs(func, dfunc, x0);
    }

}

#endif