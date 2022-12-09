#ifndef VVERY_SIMPLE_ALGORITHM2_HEADER
#define VVERY_SIMPLE_ALGORITHM2_HEADER

#include "algorithm1.hpp"

namespace ppx
{
    enum class optimization : char
    {
        GoldenSection,
        Brent,
        Powell,
        FRPRCG,
        BGFS
    };

    namespace details
    {
        using range = std::pair<double, double>;
        constexpr double gl_rep_ftol = 3.0e-8;

        template <typename T>
        range bracket_minimum(const T &func, double x0, double x1)
        {
            // deal with monotony?
            constexpr double k = 2.0;
            auto s = fabs(x0 - x1);
            auto a = x1;
            auto ya = func(x1);
            auto b = x2;
            auto yb = func(x2);
            auto c = 0.0;
            if (yb > ya)
            {
                std::swap(a, b);
                std::swap(ya, yb);
                s = -s;
            }
            while (true)
            {
                c = b + s;
                auto yc = func(c);
                if (yc > yb)
                {
                    break;
                }
                a = b;
                ya = yb;
                b = c;
                yb = yc;
                s *= k;
            }
            return a < c ? {a, c} : {c, a};
        }

        template <typename T>
        double golden_section(const T &func, double a, double b)
        {
            constexpr double R = 0.6180339887;
            constexpr double C = 1.0 - R;
            auto d = R * b + C * a;
            auto yd = func(d);
            auto its = 0u;
            while (fabs(a - b) > gl_rep_ftol && its < 200)
            {
                auto c = R * a + C * b;
                auto yc = func(c);
                if (yc < yd)
                {
                    b = d;
                    d = c;
                    yd = yc;
                }
                else
                {
                    a = b;
                    b = c;
                }
                ++its;
                // printf("iter = %d, xc = %g, residual = %g\n", its, c, fabs(a - b));
            }
            return 0.5 * (a + b);
        }

        template <typename T>
        double brent(const T &func, double x0, double x1)
        {
            constexpr int ITMAX = 200;
            constexpr double C = 0.3819660;

            auto a = gl_get_less_dynamic(x0, x1);
            auto b = gl_get_more_dynamic(x0, x1);
            auto x = b;
            auto w = b;
            auto v = b;
            auto fx = func(x);
            auto fw = fx;
            auto fv = fx;
            for (int its = 0; its < ITMAX; its++)
            {
                double d = 0.0;
                double e = 0.0;
                double u = 0.0;
                auto xm = 0.5 * (a + b);
                auto tol1 = gl_rep_ftol * fabs(x);
                auto tol2 = 2.0 * tol1;
                // printf("its = %d, xc = %g, residual = %g\n", its, x, fabs(x - xm));
                // printf("%g %g\n", xm, func(xm));
                if (fabs(x - xm) < tol2 - 0.5 * (b - a))
                {
                    return x;
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
                    auto etemp = e;
                    e = d;
                    if (fabs(p) > fabs(0.5 * q * etemp) || p < q * (a - x) || p > q * (b - x))
                    {
                        e = x > xm ? a - x : b - x;
                        d = C * e;
                    }
                    else
                    {
                        d = p / q;
                        u = x + d;
                        if (u - a < tol2 || b - u < tol2)
                        {
                            d = fabs(xm - x);
                        }
                    }
                }
                else
                {
                    e = x > xm ? a - x : b - x;
                    d = C * e;
                }
                // u = fabs(d) > tol1 ? x + d : x + SIGN(tol1, d);
                u = x + d;
                auto fu = func(u);
                if (fu < fx)
                {
                    if (u > x)
                    {
                        a = x;
                    }
                    else
                    {
                        b = x;
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
                        b = u;
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
            return x;
        }

        template <size_t N, typename T>
        double testconv(const Matrix<N, 1> &a, const Matrix<N, 1> &b, const T &func)
        {
            double test = 0.0;
            double temp = 0.0;
            for (int i = 0; i < N; i++)
            {
                temp = func(a[i], b[i]);
                if (temp > test)
                {
                    test = temp;
                }
            }
            return test;
        }
        template <typename Callable, size_t N>
        void lnsrch(const Matrix<N, 1> &xold, double fold, const Matrix<N, 1> &g, Matrix<N, 1> &p,
                    Matrix<N, 1> &x, double &f, double stpmax, bool &check, Callable &func)
        {
            constexpr double ALF = 1.0e-4;
            double alam2 = 0.0;
            double f2 = 0.0;
            check = false;
            auto sum = norm2(p);
            if (sum > stpmax)
            {
                p *= stpmax / sum;
            }
            auto slope = inner_product(g, p);
            if (slope > 0.0)
            {
                return;
            }
            // throw("Roundoff problem in lnsrch.");
            auto test = testconv(p, xold, [](double a, double b)
                                 { return fabs(a) / gl_get_more_dynamic(fabs(b), 1.0); });
            auto alamin = gl_rep_eps / test;
            auto alam = 1.0;
            for (;;)
            {
                auto tmplam = 0.0;
                x = xold + alam * p;
                f = func(x);
                if (alam < alamin)
                {
                    x = xold;
                    check = true;
                    return;
                }
                else if (f < fold + ALF * alam * slope)
                {
                    return;
                }
                else
                {
                    if (alam == 1.0)
                    {
                        tmplam = -slope / (2.0 * (f - fold - slope));
                    }
                    else
                    {
                        auto rhs1 = f - fold - alam * slope;
                        auto rhs2 = f2 - fold - alam2 * slope;
                        auto a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
                        auto b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
                        if (fabs(a) < gl_rep_eps)
                        {
                            tmplam = -slope / (2.0 * b);
                        }
                        else
                        {
                            auto disc = b * b - 3.0 * a * slope;
                            if (disc < 0.0)
                            {
                                tmplam = 0.5 * alam;
                            }
                            else if (b < 0.0)
                            {
                                tmplam = (-b + sqrt(disc)) / (3.0 * a);
                            }
                            else
                            {
                                tmplam = -slope / (b + sqrt(disc));
                            }
                        }
                        if (tmplam > 0.5 * alam)
                        {
                            tmplam = 0.5 * alam;
                        }
                    }
                }
                alam2 = alam;
                f2 = f;
                alam = gl_get_more_dynamic(tmplam, 0.1 * alam);
            }
        }

        struct Bracket
        {
            double ax;
            double bx;
            double cx;
            double fa;
            double fb;
            double fc;

            template <typename T>
            void bounded(double a, double b, const T &func)
            {
                ax = a;
                bx = b;
                fa = func(ax);
                fb = func(bx);
                if (fb > fa)
                {
                    std::swap(ax, bx);
                    std::swap(fb, fa);
                }
                cx = (ax + bx) / 2.0;
                fc = func(cx);
            }

            template <typename T>
            void bracket(double a, double b, const T &func)
            {
                constexpr double GOLD = 1.618034;
                // how long should be ?
                double GLIMIT = 5 * fabs(b - a);
                ax = a;
                bx = b;
                fa = func(ax);
                fb = func(bx);
                if (fb > fa)
                {
                    std::swap(ax, bx);
                    std::swap(fb, fa);
                }
                cx = bx + GOLD * (bx - ax);
                fc = func(cx);
                while (fb > fc)
                {
                    double fu = 0.0;
                    double r = (bx - ax) * (fb - fc);
                    double q = (bx - cx) * (fb - fa);
                    double u = bx - ((bx - cx) * q - (bx - ax) * r) /
                                        (2.0 * SIGN(gl_get_more_dynamic(fabs(q - r), gl_rep_eps), q - r));
                    double ulim = bx + GLIMIT * (cx - bx);
                    if ((bx - u) * (u - cx) > 0.0)
                    {
                        fu = func(u);
                        if (fu < fc)
                        {
                            ax = bx;
                            bx = u;
                            fa = fb;
                            fb = fu;
                            return;
                        }
                        else if (fu > fb)
                        {
                            cx = u;
                            fc = fu;
                            return;
                        }
                        u = cx + GOLD * (cx - bx);
                        fu = func(u);
                    }
                    else if ((cx - u) * (u - ulim) > 0.0)
                    {
                        fu = func(u);
                        if (fu < fc)
                        {
                            shft3(bx, cx, u, u + GOLD * (u - cx));
                            shft3(fb, fc, fu, func(u));
                        }
                    }
                    else if ((u - ulim) * (ulim - cx) >= 0.0)
                    {
                        u = ulim;
                        fu = func(u);
                    }
                    else
                    {
                        u = cx + GOLD * (cx - bx);
                        fu = func(u);
                    }
                    shft3(ax, bx, cx, u);
                    shft3(fa, fb, fc, fu);
                }
            }

            void shft2(double &a, double &b, const double c)
            {
                a = b;
                b = c;
            }

            void shft3(double &a, double &b, double &c, const double d)
            {
                a = b;
                b = c;
                c = d;
            }

            void mov3(double &a, double &b, double &c, const double d, const double e,
                      const double f)
            {
                a = d;
                b = e;
                c = f;
            }
        };

        template <size_t N>
        struct LineMethod
        {
            Matrix<N, 1> p;
            Matrix<N, 1> xi;
            template <typename T>
            double minimize(const T &func)
            {
                auto F1dim = [&](double x)
                {
                    return func(p + x * xi);
                };
                auto ax = 0.0;
                auto xx = 1.0;
                Brent brent;
                brent.bracket(ax, xx, F1dim);
                auto xmin = brent.minimize(F1dim);
                xi *= xmin;
                p += xi;
                return brent.fmin;
            }
        };
    }

    struct GoldenSection : details::Bracket
    {
        double xmin;
        double fmin;
        double tol;

        GoldenSection(double tol_ = 3.0e-8) : tol(tol_) {}

        template <typename T>
        double minimize(const T &func)
        {
            constexpr double R = 0.61803399;
            constexpr double C = 1.0 - R;
            double x1, x2;
            double x0 = ax;
            double x3 = cx;
            if (fabs(cx - bx) > fabs(bx - ax))
            {
                x1 = bx;
                x2 = bx + C * (cx - bx);
            }
            else
            {
                x2 = bx;
                x1 = bx - C * (bx - ax);
            }
            double f1 = func(x1);
            double f2 = func(x2);
            while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2)))
            {
                if (f2 < f1)
                {
                    shft3(x0, x1, x2, R * x2 + C * x3);
                    shft2(f1, f2, func(x2));
                }
                else
                {
                    shft3(x3, x2, x1, R * x1 + C * x0);
                    shft2(f2, f1, func(x1));
                }
            }
            if (f1 < f2)
            {
                xmin = x1;
                fmin = f1;
            }
            else
            {
                xmin = x2;
                fmin = f2;
            }
            return xmin;
        }
    };

    struct Brent : details::Bracket
    {
        double xmin;
        double fmin;
        double tol;

        Brent(double tol_ = 3.0e-8) : tol(tol_) {}

        template <typename T>
        double minimize(const T &func)
        {
            constexpr int ITMAX = 100;
            constexpr double CGOLD = 0.3819660;
            double d = 0.0;

            auto a = (ax < cx ? ax : cx);
            auto b = (ax > cx ? ax : cx);
            auto x = bx;
            auto w = bx;
            auto v = bx;
            auto fx = func(x);
            auto fw = fx;
            auto fv = fx;
            for (int iter = 0; iter < ITMAX; iter++)
            {
                double e = 0.0;
                double u = 0.0;
                auto xm = 0.5 * (a + b);
                auto tol1 = tol * fabs(x);
                auto tol2 = 2.0 * tol1;
                if (fabs(x - xm) < tol2 - 0.5 * (b - a))
                {
                    fmin = fx;
                    xmin = x;
                    return xmin;
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
                    auto etemp = e;
                    e = d;
                    if (fabs(p) > fabs(0.5 * q * etemp) || p < q * (a - x) || p > q * (b - x))
                    {
                        e = x > xm ? a - x : b - x;
                        d = CGOLD * e;
                    }
                    else
                    {
                        d = p / q;
                        u = x + d;
                        if (u - a < tol2 || b - u < tol2)
                        {
                            d = SIGN(tol1, xm - x);
                        }
                    }
                }
                else
                {
                    e = x > xm ? a - x : b - x;
                    d = CGOLD * e;
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
                        b = x;
                    }
                    shft3(v, w, x, u);
                    shft3(fv, fw, fx, fu);
                }
                else
                {
                    if (u < x)
                    {
                        a = u;
                    }
                    else
                    {
                        b = u;
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
            return xmin;
            // Too many iterations in brent?
        }
    };

    template <optimization etype, typename T>
    typename std::enable_if<etype == optimization::GoldenSection, double>::type
    fminbnd(T fn, double a, double b)
    {
        GoldenSection gs;
        gs.bounded(a, b, fn);
        return gs.minimize(fn);
    }

    template <optimization etype, typename T>
    typename std::enable_if<etype == optimization::Brent, double>::type
    fminbnd(T fn, double a, double b)
    {
        Brent brent;
        brent.bounded(a, b, fn);
        return brent.minimize(fn);
    }

    template <size_t N>
    struct Powell : details::LineMethod<N>
    {
        double fret;
        double ftol;

        using LineMethod<N>::p;
        using LineMethod<N>::xi;

        Powell(double ftoll = 3.0e-8) : ftol(ftoll) {}

        template <typename T>
        Matrix<N, 1> minimize(const Matrix<N, 1> &x0, const T &func)
        {
            constexpr int ITMAX = 20000;
            double fptt;
            p = x0;
            Matrix<N, 1> pt;
            Matrix<N, 1> ptt;
            auto ximat = Matrix<N, N>::eye();
            fret = func(p);
            pt = p;
            for (int iter = 0; iter < ITMAX; ++iter)
            {
                double fp = fret;
                int ibig = 0;
                double del = 0.0;
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        xi[j] = ximat(j, i);
                    }
                    fptt = fret;
                    fret = LineMethod<N>::minimize(func);
                    if (fptt - fret > del)
                    {
                        del = fptt - fret;
                        ibig = i + 1;
                    }
                }
                if (2.0 * (fp - fret) < ftol * (fabs(fp) + fabs(fret)))
                {
                    return p;
                }
                else
                {
                    printf("iter = %d, residual = %g\n", iter, fp - fret);
                }
                ptt = 2.0 * p - pt;
                xi = p - pt;
                pt = p;
                fptt = func(ptt);
                if (fptt < fp)
                {
                    double t = 2.0 * (fp - 2.0 * fret + fptt);
                    t *= (fp - fret - del) * (fp - fret - del) - del * (fp - fptt) * (fp - fptt);
                    if (t < 0.0)
                    {
                        fret = LineMethod<N>::minimize(func);
                        for (int j = 0; j < N; j++)
                        {
                            ximat(j, ibig - 1) = ximat(j, N - 1);
                            ximat(j, N - 1) = xi[j];
                        }
                    }
                }
            }
            // should what?
            return {};
        }
    };

    template <size_t N>
    struct FRPRCG : details::LineMethod<N>
    {
        double fret;
        const double ftol;

        using LineMethod<N>::p;
        using LineMethod<N>::xi;

        FRPRCG(double ftoll = 3.0e-8) : ftol(ftoll) {}

        template <typename T, typename T2>
        Matrix<N, 1> minimize(const Matrix<N, 1> &x0, const T &func, const T2 &dfunc)
        {
            constexpr int ITMAX = 200;
            constexpr double GTOL = 1.0e-8;
            p = x0;
            double fp = func(p);
            Matrix<N, 1> g = -1 * dfunc(p);
            Matrix<N, 1> h(g);
            xi = g;
            for (int its = 0; its < ITMAX; its++)
            {
                fret = LineMethod<N>::minimize(func);
                if (2.0 * fabs(fret - fp) < ftol * (fabs(fret) + fabs(fp)))
                {
                    return p;
                }
                fp = fret;
                xi = dfunc(p);
                double den = gl_get_more_dynamic(fabs(fp), 1.0);
                if (testconv(xi, p, [den](double a, double b)
                             { return fabs(a) * gl_get_more_dynamic(fabs(b), 1.0) / den; }) < GTOL)
                {
                    return p;
                }
                double dgg = 0.0;
                double gg = 0.0;
                for (int j = 0; j < N; j++)
                {
                    gg += g[j] * g[j];
                    dgg += (xi[j] + g[j]) * xi[j];
                }
                if (fabs(gg) < gl_rep_eps)
                {
                    return p;
                }
                else
                {
                    printf("iter = %d, residual = %g\n", its, gg);
                }
                double gam = dgg / gg;
                g = -1 * xi;
                h = g + gam * h;
                xi = h;
            }
            return {};
        }
    };

    template <size_t N>
    struct BGFS : details::LineMethod<N>
    {
        double fret;
        const double gtol;

        using LineMethod<N>::p;
        using LineMethod<N>::xi;

        BGFS(double gtoll = 3.0e-8) : gtol(gtoll) {}

        template <typename T, typename T2>
        Matrix<N, 1> minimize(const Matrix<N, 1> &x0, const T &func, const T2 &dfunc)
        {
            constexpr int ITMAX = 200;
            constexpr double TOLX = 4 * gl_rep_eps;
            constexpr double STPMX = 100.0;
            p = x0;
            bool check = false;
            Matrix<N, 1> dg, hdg, pnew;
            auto fp = func(p);
            auto g = dfunc(p);
            xi = -1.0 * g;
            auto hessin = Matrix<N, N>::eye();
            auto sum = inner_product(p, p);
            auto stpmax = STPMX * gl_get_more_dynamic(sqrt(sum), double(N));
            for (int its = 0; its < ITMAX; its++)
            {
                lnsrch(p, fp, g, xi, pnew, fret, stpmax, check, func);
                fp = fret;
                xi = pnew - p;
                p = pnew;
                if (testconv(xi, p, [](double a, double b)
                             { return fabs(a) / gl_get_more_dynamic(fabs(b), 1.0); }) < TOLX)
                {
                    return p;
                }
                dg = g;
                g = dfunc(p);
                auto den = gl_get_more_dynamic(fabs(fret), 1.0);
                auto test = testconv(g, p, [den](double a, double b)
                                     { return fabs(a) * gl_get_more_dynamic(fabs(b), 1.0) / den; });
                if (test < gtol)
                {
                    return p;
                }
                else
                {
                    printf("iter = %d, residual = %g\n", its, test);
                }
                dg = g - dg;
                hdg = hessin * dg;
                auto fac = 0.0;
                auto fae = 0.0;
                auto fad = 0.0;
                auto sumdg = 0.0;
                auto sumxi = 0.0;
                for (int i = 0; i < N; i++)
                {
                    fac += dg[i] * xi[i];
                    fae += dg[i] * hdg[i];
                    sumdg += dg[i] * dg[i];
                    sumxi += xi[i] * xi[i];
                }
                if (fac > sqrt(gl_rep_eps * sumdg * sumxi))
                {
                    fac = 1.0 / fac;
                    fad = 1.0 / fae;
                    dg = fac * xi - fad * hdg;
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = i; j < N; j++)
                        {
                            hessin(i, j) += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
                            hessin(j, i) = hessin(i, j);
                        }
                    }
                }
                zeros(xi);
                xi -= hessin * g;
            }
            return {};
        }
    };

    template <optimization etype, size_t N, typename T>
    typename std::enable_if<etype == optimization::Powell, Matrix<N, 1>>::type
    fminunc(const Matrix<N, 1> &x0, T func)
    {
        Powell<N> pw;
        return pw.minimize(x0, func);
    }

    template <optimization etype, size_t N, typename T, typename T2>
    typename std::enable_if<etype == optimization::FRPRCG, Matrix<N, 1>>::type
    fminunc(const Matrix<N, 1> &x0, T func, T2 dfunc)
    {
        FRPRCG<N> cg;
        return cg.minimize(x0, func, dfunc);
    }

    template <optimization etype, size_t N, typename T, typename T2>
    typename std::enable_if<etype == optimization::BGFS, Matrix<N, 1>>::type
    fminunc(const Matrix<N, 1> &x0, T func, T2 dfunc)
    {
        BGFS<N> bgfs;
        return bgfs.minimize(x0, func, dfunc);
    }

}

#endif