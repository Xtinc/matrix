#ifndef VVERY_SIMPLE_ALGORITHM2_HEADER
#define VVERY_SIMPLE_ALGORITHM2_HEADER

#include "algorithm1.hpp"

namespace ppx
{
    enum class optimization : char
    {
        GoldenSection,
        Brent
    };

    namespace details
    {
        struct Bracket
        {
            double ax;
            double bx;
            double cx;
            double fa;
            double fb;
            double fc;

            template <typename T>
            void bracket(double a, double b, const T &func)
            {
                constexpr double GOLD = 1.618034;
                // how long should be ?
                double GLIMIT = 5 * fabs(b - a);
                ax = a;
                bx = b;
                double fu;
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
                auto tol1 = tol * fabs(x) + gl_rep_eps;
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
            printf("Too many iterations in brent\n");
            return xmin;
            // Too many iterations in brent?
        }
    };

    template <optimization etype, typename T>
    typename std::enable_if<etype == optimization::GoldenSection, double>::type
    fminbnd(T fn, double a, double b)
    {
        GoldenSection gs;
        gs.bracket(a, b, fn);
        return gs.minimize(fn);
    }

    template <optimization etype, typename T>
    typename std::enable_if<etype == optimization::Brent, double>::type
    fminbnd(T fn, double a, double b)
    {
        Brent brent;
        brent.bracket(a, b, fn);
        return brent.minimize(fn);
    }
}

#endif