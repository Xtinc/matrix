#ifndef VVERY_SIMPLE_OPTIMIZATION_HEADER
#define VVERY_SIMPLE_OPTIMIZATION_HEADER

#include "linalg.hpp"
#include <functional>

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
                auto Q = eye<N>();
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
                auto U = eye<N>();
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

        template <size_t M, size_t N>
        struct CoDo
        {
            size_t ITMAX = 10000;
            double FTOLA = 1e-6;

            using vecx = MatrixS<M, 1>;
            using vecy = MatrixS<N, 1>;
            using matyx = MatrixS<N, M>;
            using func = std::function<vecy(const vecx &)>;

            CoDo(const func &_f) : f(_f)
            {
                lo.fill(-MAX_SP);
                hi.fill(MAX_SP);
            }

            CoDo(const func &_f, const vecx &lower, const vecx &upper)
                : f(_f), lo(lower), hi(upper)
            {
                // Check up > lo.
            }

            OptResult<M> operator()(vecx x)
            {
                auto sigma = 0.99995;
                auto thetal = 0.99995;

                auto fx = f(x);
                auto fnrm = norm2(fx);

                // iteration
                auto itc = 0;
                auto lambda = 0.0;
                auto delta = 1.0;
                vecx grad, p;
                while (fnrm > FTOLA && itc++ < ITMAX)
                {
                    std::cout << "itc:\t" << itc << "\t" << fnrm << "\n";
                    auto fnrm0 = fnrm;
                    auto jac = NumericaJacobi(x, fx);
                    auto grad_old = grad;
                    grad = jac.T() * fx;

                    // calculation of the scaling matrices d (d), d^1/2 (dsqr), and inv(d^1/2) (dmsqr)

                    if (itc == 1)
                    {
                        lambda = std::max(1e-2, norm2(grad));
                    }
                    else
                    {
                        lambda = std::max(1e-2, (p.T() * (grad - grad_old))[0] / inner_product(p, p));
                    }

                    auto d = DMatrixCi(x, grad, lambda);
                    vecx dsqrt = Sqrt(d);
                    vecx dsqrtgrad = pwmul(dsqrt, grad);
                    vecx G = 1.0 / dsqrt;
                    vecx dgrad = pwmul(d, grad);
                    vecx Gdgrad = pwmul(G, dgrad);
                    vecx Gmgrad = pwdiv(grad, G);

                    if (norm2(dgrad) < FTOLA)
                    {
                        // the sequence has approached a minumum of f in the box
                        break;
                    }

                    // gradient descent Step
                    auto vert = SQR(norm2(dsqrtgrad) / norm2(jac * dgrad));
                    vecx pc = -vert * dgrad;
                    // if (itc != 1)
                    // {
                    //     delta = norm2(Gmgrad);
                    // }

                    // Projected Newton Step
                    auto result = linsolve<Factorization::LU>(jac, fx);
                    vecx sn = -1 * result.x;
                    if (result.s != StatusCode::SINGULAR)
                    {
                        sn = Max(x + sn, lo) - x;
                        sn = Min(x + sn, hi) - x;
                        sn = std::max(sigma, 1 - norm2(sn)) * sn;
                    }

                    // trust-region strategy
                    auto rhof = 0.0;
                    auto delta_tmp = delta;
                    vecx xp;
                    vecy fxp;
                    double fxpnrm;
                    while (rhof < 0.25 && delta > sqrt(EPS_DP))
                    {
                        auto pcv = pc;
                        // Cauchy Point
                        if (vert * norm2(Gdgrad) > delta)
                        {
                            pcv = -delta / norm2(Gdgrad) * dgrad;
                        }

                        auto npcv = norm2(pcv);
                        vecx alp;
                        for (size_t i = 0; i < vecx::LEN; i++)
                        {
                            if (std::abs(pcv[i]) < EPS_DP)
                            {
                                alp[i] = MAX_SP;
                            }
                            else
                            {
                                alp[i] = std::max((lo[i] - x[i]) / pcv[i], (hi[i] - x[i]) / pcv[i]);
                            }
                        }
                        auto alpha1 = *std::min_element(alp.cbegin(), alp.cend());
                        auto pciv = pcv;
                        if (alpha1 <= 1)
                        {
                            pciv = std::max(thetal, 1 - npcv) * alpha1 * pcv;
                        }

                        // path computation
                        if (result.s != StatusCode::SINGULAR)
                        {
                            vecy aa = fx + jac * pciv;
                            vecx seg = sn - pciv;
                            vecy bb = jac * seg;

                            auto gamma = 0.0;
                            if (std::abs(norm2(bb)) > 0)
                            {
                                auto gammamin = -inner_product(aa, bb) / inner_product(bb, bb);
                                vecx Gseg = pwmul(G, seg);
                                vecx Gpciv = pwmul(G, pciv);
                                auto a = SQR(norm2(Gseg));
                                auto b = inner_product(Gpciv, Gseg);
                                auto c = SQR(norm2(Gpciv)) - SQR(delta);
                                auto l1 = (-b + SIGN(std::sqrt(b * b - a * c), -b)) / a;
                                auto l2 = c / (l1 * a);
                                if (gammamin > 0)
                                {
                                    auto gammaend = std::max(l1, l2);
                                    gamma = std::min(gammamin, gammaend);
                                    if (gamma > 1)
                                    {
                                        for (size_t i = 0; i < vecx::LEN; i++)
                                        {
                                            if (fabs(seg[i]) > EPS_DP)
                                            {
                                                alp[i] = std::max((lo[i] - x[i] - pciv[i]) / seg[i],
                                                                  (hi[i] - x[i] - pciv[i]) / seg[i]);
                                            }
                                            else
                                            {
                                                alp[i] = MAX_SP;
                                            }
                                        }
                                        alpha1 = *std::min_element(alp.cbegin(), alp.cend());
                                        gamma = std::min(thetal * alpha1, gamma);
                                    }
                                }
                                else
                                {
                                    auto gammaend = std::min(l1, l2);
                                    gamma = std::max(gammamin, gammaend);
                                    if (gamma < 0)
                                    {
                                        for (size_t i = 0; i < vecx::LEN; i++)
                                        {
                                            if (fabs(seg[i]) > EPS_DP)
                                            {
                                                alp[i] = std::max((lo[i] - x[i] - pciv[i]) / seg[i],
                                                                  (hi[i] - x[i] - pciv[i]) / seg[i]);
                                            }
                                            else
                                            {
                                                alp[i] = MAX_SP;
                                            }
                                        }
                                        alpha1 = *std::min_element(alp.cbegin(), alp.cend());
                                        gamma = std::max(-thetal * alpha1, gamma);
                                    }
                                }
                            }
                            p = (1 - gamma) * pciv + gamma * sn;
                        }
                        else
                        {
                            p = pciv;
                        }

                        // accuracy requirements
                        xp = x + p;
                        fxp = f(xp);
                        fxpnrm = norm2(fxp);
                        rhof = (fnrm - fxpnrm) / (fnrm - norm2<N, 1>(fx + jac * p));
                        delta_tmp = delta;
                        delta = std::min(0.25 * delta, 0.5 * norm2<N, 1>(pwmul(G, p)));
                    }

                    if (delta < sqrt(EPS_DP) && rhof < 0.25)
                    {
                        // the trust region radius Delta has become too small
                        break;
                    }
                    delta = delta_tmp;

                    // update
                    x = xp;
                    bool inddel{false};
                    for (size_t i = 0; i < vecx::LEN; i++)
                    {
                        if (x[i] < lo[i] + 1e2 * EPS_DP)
                        {
                            x[i] = lo[i] + 1e2 * EPS_DP;
                            inddel = true;
                        }
                        if (x[i] > hi[i] - 1e2 * EPS_DP)
                        {
                            x[i] = hi[i] - 1e2 * EPS_DP;
                            inddel = true;
                        }
                    }
                    if (inddel)
                    {
                        fx = f(x);
                    }
                    else
                    {
                        fx = fxp;
                    }
                    fnrm = fxpnrm;
                    if (std::abs(fnrm - fnrm0) < 1e2 * EPS_DP * fnrm && fnrm > FTOLA)
                    {
                        // no improvement for the nonlinear reasidual could be obtained
                        break;
                    }
                    if (rhof > 0.75)
                    {
                        delta = std::max(delta, 2 * norm2<M, 1>(pwmul(G, p)));
                    }
                }
                return {x, 0.5 * norm2(fx), StatusCode::CONVERGED};
            }

        private:
            func f;
            vecx lo;
            vecx hi;

            auto DMatrixCi(const vecx &x, const vecx &grad, double lambda)
            {
                constexpr auto gamma = 1;
                vecx d;
                // Hager Scaling.
                for (size_t i = 0; i < vecx::LEN; i++)
                {
                    if (grad[i] < 0.0 && hi[i] <= MAX_SP)
                    {
                        auto diff = hi[i] - x[i];
                        d[i] = 1.0 / (lambda - grad[i] / diff);
                    }
                    else if (grad[i] > 0 && lo[i] >= -MAX_SP)
                    {
                        auto diff = x[i] - lo[i];
                        d[i] = 1.0 / (lambda + grad[i] / diff);
                    }
                    else
                    {
                        d[i] = 1.0 / lambda;
                    }
                }
                // Check d[i] gt 0.0
                return d;
            }

            auto NumericaJacobi(const vecx &x, const vecy &fx)
            {
                auto EPS_DIFF = std::sqrt(EPS_DP);
                matyx result;
                auto nrmx = norm1(x) / (double)(vecx::LEN);
                for (size_t j = 0; j < vecx::LEN; j++)
                {
                    auto xh = x;
                    auto h = std::abs(x[j]) > EPS_DP ? EPS_DIFF * SIGN(std::max(x[j], nrmx), x[j]) : EPS_DIFF;
                    xh[j] = x[j] + h;
                    if (xh[j] < lo[j] || xh[j] > hi[j])
                    {
                        h *= -1;
                        xh[j] = x[j] + h;
                    }
                    result.sub<vecy::LEN, 1, false>(0, j) = (f(xh) - fx) / h;
                }
                return result;
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