#ifndef VVERY_SIMPLE_SIGNALS_HEADER
#define VVERY_SIMPLE_SIGNALS_HEADER

#include "statistics.hpp"
#include <queue>

namespace ppx
{
    // signals
    inline double sinc(double x)
    {
        return fabs(x) < EPS_DP ? 1.0 : sin(x) / x;
    }

    template <size_t M, size_t N>
    MatrixS<N + M - 1, 1>
    convovle(const MatrixS<M, 1> &f, const MatrixS<N, 1> &g)
    {
        MatrixS<M + N - 1, 1> out;
        for (auto i = 0; i < M + N - 1; ++i)
        {
            int jmn = (i >= N - 1) ? i - (N - 1) : 0;
            int jmx = (i < M - 1) ? i : M - 1;
            for (auto j = jmn; j <= jmx; ++j)
            {
                out[i] += (f[j] * g[i - j]);
            }
        }
        return out;
    }

    namespace details
    {
        inline double BartlettWindow(int N, int i)
        {
            double period = N - 1;
            return 1 - 2.0 * fabs(i - period / 2.0) / period;
        }

        inline double BlackmanWindow(int N, int i)
        {
            double period = N - 1;
            return 0.42 - 0.5 * cos(2.0 * PI * i / period) +
                   0.08 * cos(4.0 * PI * i / period);
        }

        inline double HammingWindow(int N, int i)
        {
            double period = N - 1;
            return 0.54 - 0.46 * cos(2.0 * PI * i / period);
        }

        inline double HanningWindow(int N, int i)
        {
            double period = N - 1;
            return 0.5 * (1.0 - cos(2.0 * PI * i / period));
        }

        // The multiplication has the following form:
        // (x+p[0])*(x+p[1])*...*(x+p[n-1])
        // The p[i] coefficients are assumed to be complex and are passed to the
        // function as an array of doubles of length 2n. where p[2i] (i=0...n-1)
        // is assumed to be the real part of the coefficient of the ith binomial
        // and p[2i+1] is assumed to be the imaginary part. The resulting polynomial
        // has the following form:
        // x^n + a[0]*x^n-1 + a[1]*x^n-2 + ... +a[n-2]*x + a[n-1]
        template <size_t N>
        void binomial_mult(const MatrixS<N, 1> &p, MatrixS<N, 1> &a)
        {
            for (int i = 0; i < (int)N / 2; ++i)
            {
                for (int j = i; j > 0; --j)
                {
                    a[2 * j] += p[2 * i] * a[2 * (j - 1)] - p[2 * i + 1] * a[2 * (j - 1) + 1];
                    a[2 * j + 1] += p[2 * i] * a[2 * (j - 1) + 1] + p[2 * i + 1] * a[2 * (j - 1)];
                }
                a[0] += p[2 * i];
                a[1] += p[2 * i + 1];
            }
        }

        // trinomial_mult - multiplies a series of trinomials together and returns
        // the coefficients of the resulting polynomial. The multiplication has
        // the following form:
        //(x^2 + b[0]x + c[0])*(x^2 + b[1]x + c[1])*...*(x^2 + b[n-1]x + c[n-1])
        // The b[i] and c[i] coefficients are assumed to be complex and are passed
        // to the function as arrays of doubles of length 2n. The real
        // part of the coefficients are stored in the even numbered elements of the
        // array and the imaginary parts are stored in the odd numbered elements.
        // The resulting polynomial has the following form:
        // x^2n + a[0]*x^2n-1 + a[1]*x^2n-2 + ... +a[2n-2]*x + a[2n-1]
        template <size_t N>
        void trinomial_mult(const MatrixS<N, 1> &b, const MatrixS<N, 1> &c, MatrixS<2 * N, 1> &a)
        {
            // todo: static_assert of N > 4
            a[2] = c[0];
            a[3] = c[1];
            a[0] = b[0];
            a[1] = b[1];

            for (int i = 1; i < (int)N / 2; ++i)
            {
                a[2 * (2 * i + 1)] += c[2 * i] * a[2 * (2 * i - 1)] - c[2 * i + 1] * a[2 * (2 * i - 1) + 1];
                a[2 * (2 * i + 1) + 1] += c[2 * i] * a[2 * (2 * i - 1) + 1] + c[2 * i + 1] * a[2 * (2 * i - 1)];

                for (int j = 2 * i; j > 1; --j)
                {
                    a[2 * j] += b[2 * i] * a[2 * (j - 1)] - b[2 * i + 1] * a[2 * (j - 1) + 1] +
                                c[2 * i] * a[2 * (j - 2)] - c[2 * i + 1] * a[2 * (j - 2) + 1];
                    a[2 * j + 1] += b[2 * i] * a[2 * (j - 1) + 1] + b[2 * i + 1] * a[2 * (j - 1)] +
                                    c[2 * i] * a[2 * (j - 2) + 1] + c[2 * i + 1] * a[2 * (j - 2)];
                }

                a[2] += b[2 * i] * a[0] - b[2 * i + 1] * a[1] + c[2 * i];
                a[3] += b[2 * i] * a[1] + b[2 * i + 1] * a[0] + c[2 * i + 1];
                a[0] += b[2 * i];
                a[1] += b[2 * i + 1];
            }
        }
    }

    enum class FreqProperty : char
    {
        LowPass,
        HighPass,
        BandPass,
        BandStop
    };

    enum class FIRType : char
    {
        Boxcar,
        Bartlett,
        Blackman,
        Hamming,
        Hanning
    };

    enum class IIRType : char
    {
        ButterWorth
        // Chebyshev
    };

    class Filter
    {
    public:
        Filter(FreqProperty tp = FreqProperty::LowPass, bool isdeferred = true) : freqtype(tp), deferred(isdeferred), last_result(0.0), a{1.0} {}

        std::vector<double> &coff_a()
        {
            return a;
        }

        const std::vector<double> &coff_a() const
        {
            return a;
        }

        std::vector<double> &coff_b()
        {
            return b;
        }

        const std::vector<double> &coff_b() const
        {
            return b;
        }

        FreqProperty type() const
        {
            return freqtype;
        }

        void reset()
        {
            samples.clear();
            results.clear();
        }

        double operator()(double new_sample) const
        {
            auto s1 = 0.0;
            auto s2 = 0.0;
            auto n = b.size();
            auto m = a.size();

            samples.push_front(new_sample);

            auto samples_full = samples.size() > n;
            auto results_full = results.size() + 1 > m;

            if (!results.empty())
            {
                last_result = results.front();
            }

            if (samples_full)
            {
                samples.pop_back();
            }

            if (results_full)
            {
                results.pop_back();
            }

            for (size_t i = 0; i < samples.size(); i++)
            {
                s2 += samples.at(i) * b[i];
            }

            for (size_t i = 0; i < results.size(); i++)
            {
                s1 += results.at(i) * a[i + 1];
            }

            auto res = deferred && (!samples_full) ? new_sample : s2 - s1;
            results.push_front(res);
            return res;
        }

        double diff() const
        {
            if (!results.empty())
            {
                return results.front() - last_result;
            }
            return 0.0;
        }

    protected:
        const FreqProperty freqtype;
        const bool deferred;
        std::vector<double> a;
        std::vector<double> b;

        mutable std::deque<double> samples;
        mutable std::deque<double> results;
        mutable double last_result;
    };

    template <size_t N, FreqProperty U>
    class IIRFilter : public Filter
    {
    public:
        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::LowPass> * = nullptr>
        IIRFilter(double fcf, bool isdeferred = true) : Filter(U, isdeferred)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(a), N, 0.0);
            std::fill_n(std::back_inserter(b), N + 1, 0.0);
            cal_lphp_coffb(true);
            cal_lphp_coffa(fcf);
            sca_lphp_coffb(fcf, true);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::HighPass> * = nullptr>
        IIRFilter(double fcf, bool isdeferred = true) : Filter(U, isdeferred)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(a), N, 0.0);
            std::fill_n(std::back_inserter(b), N + 1, 0.0);
            cal_lphp_coffb(false);
            cal_lphp_coffa(fcf);
            sca_lphp_coffb(fcf, false);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::BandPass> * = nullptr>
        IIRFilter(double f1f, double f2f, bool isdeferred = true) : Filter(U, isdeferred)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(a), 2 * N, 0.0);
            std::fill_n(std::back_inserter(b), 2 * N + 1, 0.0);
            cal_bpbs_coffb(f1f, f2f, true);
            cal_bpbs_coffa(f1f, f2f, true);
            sca_bpbs_coffb(f1f, f2f, true);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::BandStop> * = nullptr>
        IIRFilter(double f1f, double f2f, bool isdeferred = true) : Filter(U, isdeferred)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(a), 2 * N, 0.0);
            std::fill_n(std::back_inserter(b), 2 * N + 1, 0.0);
            cal_bpbs_coffb(f1f, f2f, false);
            cal_bpbs_coffa(f1f, f2f, false);
            sca_bpbs_coffb(f1f, f2f, false);
        }

    private:
        void sca_lphp_coffb(double fcf, bool lp)
        {
            auto theta = PI * fcf;
            auto st = sin(theta);

            auto sf = 1.0;
            for (int k = 0; k < (int)N / 2; ++k)
            {
                sf *= 1.0 + st * sin((2 * k + 1) * PI / (2 * N));
            }
            auto fomega = lp ? sin(theta / 2.0) : cos(theta / 2.0);
            if (N % 2)
            {
                sf *= fomega + (lp ? cos(theta / 2.0) : sin(theta / 2.0));
            }
            sf = std::pow(fomega, N) / sf;
            for (auto &i : b)
            {
                i *= sf;
            }
        }

        void sca_bpbs_coffb(double f1f, double f2f, bool bp)
        {
            auto ctt = bp ? 1.0 / tan(PI * (f2f - f1f) / 2.0) : tan(PI * (f2f - f1f) / 2.0);
            auto sfr = 1.0;
            auto sfi = 0.0;

            for (int k = 0; k < N; ++k)
            {
                auto parg = PI * (2 * k + 1) / (double)(2 * N);
                auto sparg = ctt + sin(parg);
                auto cparg = cos(parg);
                auto sa = (sfr + sfi) * (sparg - cparg);
                auto sb = sfr * sparg;
                auto sc = -sfi * cparg;
                sfr = sb - sc;
                sfi = sa - sb - sc;
            }
            for (auto &i : b)
            {
                i /= sfr;
            }
        }

        void cal_lphp_coffb(bool lp)
        {
            b[0] = 1.0;
            b[1] = N;
            for (int i = 2; i <= (int)N / 2; ++i)
            {
                b[i] = (double)(N - i + 1) * b[i - 1] / i;
                b[N - i] = b[i];
            }
            b[N - 1] = N;
            b[N] = 1;
            if (!lp)
            {
                for (size_t i = 1; i < N + 1; i = i + 2)
                {
                    b[i] *= -1.0;
                }
            }
        }

        void cal_bpbs_coffb(double f1f, double f2f, bool bp)
        {
            if (bp)
            {
                cal_lphp_coffb(false);
                auto tcof = b;
                for (int i = 0; i < N; ++i)
                {
                    b[2 * i] = tcof[i];
                    b[2 * i + 1] = 0.0;
                }
                b[2 * N] = tcof[N];
            }
            else
            {
                auto alpha = -2.0 * cos(PI * (f2f + f1f) / 2.0) / cos(PI * (f2f - f1f) / 2.0);
                b[0] = 1.0;
                b[1] = alpha;
                b[2] = 1.0;

                for (int i = 1; i < N; ++i)
                {
                    b[2 * i + 2] += b[2 * i];
                    for (int j = 2 * i; j > 1; --j)
                    {
                        b[j + 1] += alpha * b[j] + b[j - 1];
                    }
                    b[2] += alpha * b[1] + 1.0;
                    b[1] += alpha;
                }
            }
        }

        void cal_lphp_coffa(double fcf)
        {
            auto theta = PI * fcf;
            auto st = sin(theta);
            auto ct = cos(theta);
            MatrixS<2 * N, 1> rcof, dcof;

            for (int k = 0; k < (int)N; ++k)
            {
                auto parg = PI * (2 * k + 1) / (double)(2 * N);
                auto dena = 1.0 + st * sin(parg);
                rcof[2 * k] = -ct / dena;
                rcof[2 * k + 1] = -st * cos(parg) / dena;
            }
            details::binomial_mult(rcof, dcof);

            dcof[1] = dcof[0];
            dcof[0] = 1.0;
            for (size_t i = 3; i <= N; i++)
            {
                dcof[i] = dcof[2 * i - 2];
            }
            std::copy_n(dcof.begin(), N + 1, a.begin());
        }

        void cal_bpbs_coffa(double f1f, double f2f, bool bp)
        {
            auto cp = cos(PI * (f2f + f1f) / 2.0);
            auto theta = PI * (f2f - f1f) / 2.0;
            auto st = sin(theta);
            auto ct = cos(theta);
            auto s2t = 2.0 * st * ct;
            auto c2t = 2.0 * ct * ct - 1.0;

            MatrixS<2 * N, 1> rcof, tcof;
            MatrixS<4 * N, 1> dcof;

            for (int k = 0; k < N; ++k)
            {
                auto parg = PI * (2 * k + 1) / (double)(2 * N);
                auto sa = 1.0 + s2t * sin(parg);
                rcof[2 * k] = c2t / sa;
                rcof[2 * k + 1] = (bp ? 1.0 : -1.0) * s2t * cos(parg) / sa;
                tcof[2 * k] = -2.0 * cp * (ct + st * sin(parg)) / sa;
                tcof[2 * k + 1] = (bp ? -2.0 : 2.0) * cp * st * cos(parg) / sa;
            }

            details::trinomial_mult(tcof, rcof, dcof);
            dcof[1] = dcof[0];
            dcof[0] = 1.0;
            for (size_t k = 3; k <= 2 * N; ++k)
            {
                dcof[k] = dcof[2 * k - 2];
            }
            std::copy_n(dcof.begin(), 2 * N + 1, a.begin());
        }
    };

    template <size_t N, FreqProperty U>
    class FIRFilter : public Filter
    {
    public:
        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::LowPass> * = nullptr>
        FIRFilter(double fcf, FIRType k_win = FIRType::Hamming, bool isdeferred = true) : Filter(U, isdeferred)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            static_assert(N % 2 != 0, "lowpass FIR filter degree must be odd");
            std::fill_n(std::back_inserter(b), N, 0.0);
            set_window_function(k_win);
            cal_lphp_coffb(fcf / 2.0, true);
            sca_lphp_coffb(true);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::HighPass> * = nullptr>
        FIRFilter(double fcf, FIRType k_win = FIRType::Hamming, bool isdeferred = true) : Filter(U, isdeferred)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            static_assert(N % 2 != 0, "lowpass FIR filter degree must be odd");
            std::fill_n(std::back_inserter(b), N, 0.0);
            set_window_function(k_win);
            cal_lphp_coffb(fcf / 2.0, false);
            sca_lphp_coffb(false);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::BandPass> * = nullptr>
        FIRFilter(double f1f, double f2f, FIRType k_win = FIRType::Hamming, bool isdeferred = true) : Filter(U, isdeferred)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(b), N, 0.0);
            set_window_function(k_win);
            cal_bpbs_coffb(f1f / 2.0, f2f / 2.0, true);
            sca_lphp_coffb(false);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::BandStop> * = nullptr>
        FIRFilter(double f1f, double f2f, FIRType k_win = FIRType::Hamming, bool isdeferred = true) : Filter(U, isdeferred)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(b), N, 0.0);
            set_window_function(k_win);
            cal_bpbs_coffb(f1f / 2.0, f2f / 2.0, false);
            sca_lphp_coffb(true);
        }

    private:
        void set_window_function(FIRType win_type)
        {
            switch (win_type)
            {
            case FIRType::Bartlett:
                winfunc = details::BartlettWindow;
                break;
            case FIRType::Blackman:
                winfunc = details::BlackmanWindow;
                break;
            case FIRType::Hanning:
                winfunc = details::HanningWindow;
                break;
            case FIRType::Hamming:
                winfunc = details::HammingWindow;
                break;
            case FIRType::Boxcar:
            default:
                winfunc = [](int, int)
                { return 1.0; };
                break;
            }
        }

        void sca_lphp_coffb(bool lp)
        {
            auto sum = lp ? std::accumulate(b.begin(), b.end(), 0.0, [](double y, double x)
                                            { return x + y; })
                          : std::accumulate(b.begin(), b.end(), 0.0, [](double y, double x)
                                            { return fabs(x) + y; });
            for (auto &i : b)
            {
                i /= sum;
            }
        }

        void cal_lphp_coffb(double fcf, bool lp)
        {
            auto omega = 2 * PI * fcf;
            auto bsign = lp ? 2.0 : -2.0;
            for (int i = 0; i < (int)N; i++)
            {
                auto ni = i - int(N / 2);
                b[i] = bsign * fcf * sinc(omega * ni) * winfunc(N, i);
            }
            if (!lp)
            {
                b[(N - 1) / 2] += 1;
            }
        }

        void cal_bpbs_coffb(double f1f, double f2f, bool bp)
        {
            auto bsign = bp ? 1.0 : -1.0;
            cal_lphp_coffb(f1f, true);
            auto h = b;
            cal_lphp_coffb(f2f, true);
            for (size_t i = 0; i < b.size(); i++)
            {
                b[i] = bsign * (b[i] - h[i]);
            }
            if (!bp)
            {
                b[(N - 1) / 2] += 1;
            }
        }

    private:
        std::function<double(int, int)> winfunc;
    };

    template <size_t N>
    class MovAvgFilter : public Filter
    {
    public:
        MovAvgFilter(bool isdeferred = true) : Filter(FreqProperty::LowPass, isdeferred)
        {
            std::fill_n(std::back_inserter(b), N, 1.0 / N);
        }
    };

    template <size_t M, size_t NL, size_t NR, size_t NP = NL + NR + 1>
    class SGFilter : public Filter
    {
    public:
        SGFilter(int derivative_order = 0) : Filter(FreqProperty::LowPass)
        {
            std::fill_n(std::back_inserter(b), NP, 0.0);
            savgo(derivative_order);
        }

    private:
        void savgo(int ld)
        {
            int nl = NL;
            int nr = NR;
            int np = NP;

            MatrixS<M + 1, M + 1> ta;
            MatrixS<M + 1, 1> tb;
            for (int ipj = 0; ipj <= (M << 1); ipj++)
            {
                auto sum = ipj ? 0.0 : 1.0;
                for (auto k = 1; k <= nr; k++)
                {
                    sum += std::pow(k, ipj);
                }
                for (int k = 1; k <= nl; k++)
                {
                    sum += std::pow(-k, ipj);
                }
                auto mm = std::min(ipj, 2 * int(M) - ipj);
                for (auto imj = -mm; imj <= mm; imj += 2)
                {
                    ta((ipj + imj) / 2, (ipj - imj) / 2) = sum;
                }
            }
            tb[ld] = 1.0;
            tb = linsolve<Factorization::LU>(ta, tb).x;
            for (auto k = -nl; k <= nr; k++)
            {
                auto sum = tb[0];
                auto fac = 1.0;
                for (auto mm = 1; mm <= M; mm++)
                {
                    sum += tb[mm] * (fac *= k);
                }
                // auto kk = (np - k) % np;
                b[k + nl] = sum;
            }
        }
    };

} // namespace ppx

#endif