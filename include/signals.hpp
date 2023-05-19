#ifndef VVERY_SIMPLE_ALGORITHM3_HEADER
#define VVERY_SIMPLE_ALGORITHM3_HEADER

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
            for (int i = 0; i < N / 2; ++i)
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

            for (int i = 1; i < N / 2; ++i)
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

    enum FIRType : char
    {
        Boxcar,
        Bartlett,
        Blackman,
        Hamming,
        Hanning
    };

    enum IIRType : char
    {
        ButterWorth
        // Chebyshev
    };

    class Filter
    {
    public:
        Filter(FreqProperty tp) : freqtype(tp) {}

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

            if (samples.size() > n)
            {
                samples.pop_back();
            }

            for (int i = 0; i < samples.size(); i++)
            {
                s2 += samples.at(i) * b[i];
            }

            for (int i = 0; i < results.size(); i++)
            {
                s1 += results.at(i) * a[i + 1];
            }

            results.push_front(s2 - s1);

            if (results.size() > (int)m - 1)
            {
                results.pop_back();
            }
            return results.front();
        }

    protected:
        const FreqProperty freqtype;
        std::vector<double> a;
        std::vector<double> b;

        mutable std::deque<double> samples;
        mutable std::deque<double> results;
    };

    template <size_t N, FreqProperty U>
    class IIRFilter : public Filter
    {
    public:
        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::LowPass> * = nullptr>
        IIRFilter(double fcf) : Filter(U)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(a), N + 1, 0.0);
            std::fill_n(std::back_inserter(b), N + 1, 0.0);
            cal_lphp_coffb(fcf, true);
            cal_lphp_coffa(fcf, true);
            sca_lphp_coffb(fcf, true);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::HighPass> * = nullptr>
        IIRFilter(double fcf) : Filter(U)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(a), N + 1, 0.0);
            std::fill_n(std::back_inserter(b), N + 1, 0.0);
            cal_lphp_coffb(fcf, false);
            cal_lphp_coffa(fcf, false);
            sca_lphp_coffb(fcf, false);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::BandPass> * = nullptr>
        IIRFilter(double f1f, double f2f) : Filter(U)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(a), 2 * N + 1, 0.0);
            std::fill_n(std::back_inserter(b), 2 * N + 1, 0.0);
            cal_bpbs_coffb(f1f, f2f, true);
            cal_bpbs_coffa(f1f, f2f, true);
            sca_bpbs_coffb(f1f, f2f, true);
        }

        template <FreqProperty T = U, std::enable_if_t<T == FreqProperty::BandStop> * = nullptr>
        IIRFilter(double f1f, double f2f) : Filter(U)
        {
            static_assert(N > 1, "filter order must greater than 1!");
            std::fill_n(std::back_inserter(a), 2 * N + 1, 0.0);
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
            auto ct = cos(theta);

            auto sf = 1.0;
            for (int k = 0; k < N / 2; ++k)
            {
                sf *= 1.0 + st * sin((2 * k + 1) * PI / (2 * N));
            }
            auto fomega = lp ? sin(theta / 2.0) : cos(theta / 2.0);
            if (N % 2)
            {
                sf *= fomega + lp ? cos(theta / 2.0) : sin(theta / 2.0);
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
                auto a = (sfr + sfi) * (sparg - cparg);
                auto b = sfr * sparg;
                auto c = -sfi * cparg;
                sfr = b - c;
                sfi = a - b - c;
            }
            for (auto &i : b)
            {
                i /= sfr;
            }
        }

        void cal_lphp_coffb(double fcf, bool lp)
        {
            auto theta = PI * fcf;
            auto st = sin(theta);
            auto ct = cos(theta);

            b[0] = 1.0;
            b[1] = N;
            for (int i = 2; i <= N / 2; ++i)
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
                cal_lphp_coffb(f2f, false);
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

        void cal_lphp_coffa(double fcf, bool lp)
        {
            auto theta = PI * fcf;
            auto st = sin(theta);
            auto ct = cos(theta);
            MatrixS<2 * N, 1> rcof, dcof;

            for (int k = 0; k < N; ++k)
            {
                auto parg = PI * (2 * k + 1) / (double)(2 * N);
                auto a = 1.0 + st * sin(parg);
                rcof[2 * k] = -ct / a;
                rcof[2 * k + 1] = -st * cos(parg) / a;
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
                auto a = 1.0 + s2t * sin(parg);
                rcof[2 * k] = c2t / a;
                rcof[2 * k + 1] = (bp ? 1.0 : -1.0) * s2t * cos(parg) / a;
                tcof[2 * k] = -2.0 * cp * (ct + st * sin(parg)) / a;
                tcof[2 * k + 1] = (bp ? -2.0 : 2.0) * cp * st * cos(parg) / a;
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

    template <size_t N>
    class FIRFilter
    {
    public:
        FIRFilter(
            double f1, double f2, FreqProperty k_type = FreqProperty::LowPass, FIRType k_win = FIRType::Hamming)
            : freqtype(k_type), idx(0)
        {
            f1 /= 2.0;
            f2 /= 2.0;
            // Calculate the coefficient corresponding to the filter type
            switch (k_win)
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

            switch (k_type)
            {
            case FreqProperty::LowPass:
                assert(N % 2 != 0);
                h = cal_lp_coff(f1);
                break;
            case FreqProperty::HighPass:
                h = cal_hp_coff(f2);
                break;
            case FreqProperty::BandPass:
                h = convovle(cal_hp_coff(f1), cal_lp_coff(f2)).template sub<N, 1>((N - 1) / 2, 0);
                break;
            case FreqProperty::BandStop:
                h = cal_lp_coff(f1) + cal_hp_coff(f2);
                break;
            default:
                break;
            }
        }

        MatrixS<N, 1> coff() const
        {
            return h;
        }

        FreqProperty type() const
        {
            return freqtype;
        }

        double operator()(double new_sample) const
        {
            double result = 0.0;
            samples[idx] = new_sample;
            // Calculate the output
            for (int n = 0; n < N; n++)
            {
                result += samples[(idx + n) % N] * h[n];
            }
            // Increase the round robin index
            idx = (idx + 1) % N;
            return result;
        }

        void reset()
        {
            idx = 0;
            samples.fill(0.0);
        }

    protected:
        MatrixS<N, 1> cal_lp_coff(double fc) const
        {
            MatrixS<N, 1> result;
            double omega = 2 * PI * fc;
            for (int i = 0; i < N; i++)
            {
                double ni = i - int(N / 2);
                result[i] = 2.0 * fc * sinc(omega * ni) * winfunc(N, i);
            }
            result /= std::accumulate(result.begin(), result.end(), 0.0);
            return result;
        }

        MatrixS<N, 1> cal_hp_coff(double fc) const
        {
            MatrixS<N, 1> result;
            double omega = 2 * PI * fc;
            for (int i = 0; i < N; i++)
            {
                double ni = i - int(N / 2);
                result[i] = 2.0 * fc * sinc(omega * ni) * winfunc(N, i);
            }
            result /= -std::accumulate(result.begin(), result.end(), 0.0);
            result[(N - 1) / 2] += 1;
            return result;
        }

    protected:
        const FreqProperty freqtype;
        std::function<double(int, int)> winfunc;
        MatrixS<N, 1> h;
        mutable size_t idx;
        mutable MatrixS<N, 1> samples;
    };

    template <size_t N>
    class IIRFilter2
    {
    public:
        IIRFilter2(
            double f1, double f2, FreqProperty k_type = FreqProperty::LowPass, IIRType func = IIRType::ButterWorth)
            : freqtype(k_type)
        {
            static_assert(N > 1, "order must greater than 1!");
            switch (k_type)
            {
            case FreqProperty::LowPass:
                std::fill_n(std::back_inserter(c), N + 1, 0.0);
                std::fill_n(std::back_inserter(d), N + 1, 0.0);
                cal_lp_coffc(f1);
                cal_lp_coffd(f1);
                sca_lp_coffc(f1);
                break;
            case FreqProperty::HighPass:
                std::fill_n(std::back_inserter(c), N + 1, 0.0);
                std::fill_n(std::back_inserter(d), N + 1, 0.0);
                cal_hp_coffc(f2);
                cal_hp_coffd(f2);
                sca_hp_coffc(f2);
                break;
            case FreqProperty::BandPass:
                std::fill_n(std::back_inserter(c), 2 * N + 1, 0.0);
                std::fill_n(std::back_inserter(d), 2 * N + 1, 0.0);
                cal_bp_coffc(f1, f2);
                cal_bp_coffd(f1, f2);
                sca_bp_coffc(f1, f2);
                break;
            case FreqProperty::BandStop:
                std::fill_n(std::back_inserter(c), 2 * N + 1, 0.0);
                std::fill_n(std::back_inserter(d), 2 * N + 1, 0.0);
                cal_bs_coffc(f1, f2);
                cal_bs_coffd(f1, f2);
                sca_bs_coffc(f1, f2);
                break;
            default:
                break;
            }
            std::cout << "c: " << std::endl;
            for (auto &&i : c)
            {
                std::cout << i << " ";
            }
            std::cout << std::endl;

            std::cout << "d: " << std::endl;
            for (auto &&i : d)
            {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }

        double operator()(double new_sample) const
        {
            double s1 = 0.0;
            double s2 = 0.0;
            samples.push_front(new_sample);
            if (samples.size() > N + 1)
            {
                samples.pop_back();
            }

            for (int i = 0; i < results.size(); i++)
            {
                s1 += results.at(i) * d[i + 1];
            }

            for (int i = 0; i < samples.size(); i++)
            {
                s2 += samples.at(i) * c[i];
            }

            results.push_front(s2 - s1);
            if (results.size() > N)
            {
                results.pop_back();
            }

            return results.front();
        }

    protected:
        void sca_lp_coffc(double fcf)
        {
            auto theta = PI * fcf;
            auto st = sin(theta);
            auto ct = cos(theta);

            auto sf = 1.0;
            for (int k = 0; k < N / 2; ++k)
            {
                sf *= 1.0 + st * sin((2 * k + 1) * PI / (2 * N));
            }
            auto fomega = sin(theta / 2.0);
            if (N % 2)
            {
                sf *= fomega + cos(theta / 2.0);
            }
            sf = pow(fomega, N) / sf;
            for (auto &i : c)
            {
                i *= sf;
            }
        }

        void sca_hp_coffc(double fcf)
        {
            auto theta = PI * fcf;
            auto st = sin(theta);
            auto ct = cos(theta);

            auto sf = 1.0;
            for (int k = 0; k < N / 2; ++k)
            {
                sf *= 1.0 + st * sin((2 * k + 1) * PI / (2 * N));
            }
            auto fomega = cos(theta / 2.0);
            if (N % 2)
            {
                sf *= fomega + sin(theta / 2.0);
            }
            sf = pow(fomega, N) / sf;
            for (auto &i : c)
            {
                i *= sf;
            }
        }

        void sca_bp_coffc(double f1f, double f2f)
        {
            auto ctt = 1.0 / tan(PI * (f2f - f1f) / 2.0);
            auto sfr = 1.0;
            auto sfi = 0.0;

            for (int k = 0; k < N; ++k)
            {
                auto parg = PI * (2 * k + 1) / (double)(2 * N);
                auto sparg = ctt + sin(parg);
                auto cparg = cos(parg);
                auto a = (sfr + sfi) * (sparg - cparg);
                auto b = sfr * sparg;
                auto c = -sfi * cparg;
                sfr = b - c;
                sfi = a - b - c;
            }
            for (auto &i : c)
            {
                i /= sfr;
            }
        }

        void sca_bs_coffc(double f1f, double f2f)
        {
            auto tt = tan(PI * (f2f - f1f) / 2.0);
            auto sfr = 1.0;
            auto sfi = 0.0;

            for (int k = 0; k < N; ++k)
            {
                auto parg = PI * (double)(2 * k + 1) / (double)(2 * N);
                auto sparg = tt + sin(parg);
                auto cparg = cos(parg);
                auto a = (sfr + sfi) * (sparg - cparg);
                auto b = sfr * sparg;
                auto c = -sfi * cparg;
                sfr = b - c;
                sfi = a - b - c;
            }
            for (auto &i : c)
            {
                i /= sfr;
            }
        }

        void cal_lp_coffc(double fcf)
        {
            auto theta = PI * fcf;
            auto st = sin(theta);
            auto ct = cos(theta);

            c[0] = 1.0;
            c[1] = N;
            for (int i = 2; i <= N / 2; ++i)
            {
                c[i] = (double)(N - i + 1) * c[i - 1] / i;
                c[N - i] = c[i];
            }
            c[N - 1] = N;
            c[N] = 1;
        }

        void cal_hp_coffc(double fcf)
        {
            cal_lp_coffc(fcf);
            for (size_t i = 1; i < N + 1; i = i + 2)
            {
                c[i] *= -1.0;
            }
        }

        void cal_bp_coffc(double f1f, double f2f)
        {
            cal_hp_coffc(f2f);

            auto tcof = c;

            for (int i = 0; i < N; ++i)
            {
                c[2 * i] = tcof[i];
                c[2 * i + 1] = 0.0;
            }
            c[2 * N] = tcof[N];
        }

        void cal_bs_coffc(double f1f, double f2f)
        {
            auto alpha = -2.0 * cos(PI * (f2f + f1f) / 2.0) / cos(PI * (f2f - f1f) / 2.0);
            c[0] = 1.0;
            c[1] = alpha;
            c[2] = 1.0;

            for (int i = 1; i < N; ++i)
            {
                c[2 * i + 2] += c[2 * i];
                for (int j = 2 * i; j > 1; --j)
                {
                    c[j + 1] += alpha * c[j] + c[j - 1];
                }
                c[2] += alpha * c[1] + 1.0;
                c[1] += alpha;
            }
        }

        void cal_lp_coffd(double fcf)
        {
            auto theta = PI * fcf;
            auto st = sin(theta);
            auto ct = cos(theta);
            MatrixS<2 * N, 1> rcof, dcof;

            for (int k = 0; k < N; ++k)
            {
                auto parg = PI * (2 * k + 1) / (double)(2 * N);
                auto a = 1.0 + st * sin(parg);
                rcof[2 * k] = -ct / a;
                rcof[2 * k + 1] = -st * cos(parg) / a;
            }
            details::binomial_mult(rcof, dcof);

            dcof[1] = dcof[0];
            dcof[0] = 1.0;
            for (size_t i = 3; i <= N; i++)
            {
                dcof[i] = dcof[2 * i - 2];
            }
            std::copy_n(dcof.begin(), N + 1, d.begin());
        }

        void cal_hp_coffd(double fcf)
        {
            cal_lp_coffd(fcf);
        }

        void cal_bp_coffd(double f1f, double f2f)
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
                auto a = 1.0 + s2t * sin(parg);
                rcof[2 * k] = c2t / a;
                rcof[2 * k + 1] = s2t * cos(parg) / a;
                tcof[2 * k] = -2.0 * cp * (ct + st * sin(parg)) / a;
                tcof[2 * k + 1] = -2.0 * cp * st * cos(parg) / a;
            }

            details::trinomial_mult(tcof, rcof, dcof);
            dcof[1] = dcof[0];
            dcof[0] = 1.0;
            for (size_t k = 3; k <= 2 * N; ++k)
            {
                dcof[k] = dcof[2 * k - 2];
            }
            std::copy_n(dcof.begin(), 2 * N + 1, d.begin());
        }

        void cal_bs_coffd(double f1f, double f2f)
        {
            auto cp = cos(PI * (f2f + f1f) / 2.0);
            auto theta = PI * (f2f - f1f) / 2.0;
            auto st = sin(theta);
            auto ct = cos(theta);
            auto s2t = 2.0 * st * ct;       // sine of 2*theta
            auto c2t = 2.0 * ct * ct - 1.0; // cosine 0f 2*theta

            MatrixS<2 * N, 1> rcof, tcof;
            MatrixS<4 * N, 1> dcof;

            for (int k = 0; k < N; ++k)
            {
                auto parg = PI * (2 * k + 1) / (double)(2 * N);
                auto sparg = sin(parg);
                auto cparg = cos(parg);
                auto a = 1.0 + s2t * sparg;
                rcof[2 * k] = c2t / a;
                rcof[2 * k + 1] = -s2t * cparg / a;
                tcof[2 * k] = -2.0 * cp * (ct + st * sparg) / a;
                tcof[2 * k + 1] = 2.0 * cp * st * cparg / a;
            }

            details::trinomial_mult(tcof, rcof, dcof);

            dcof[1] = dcof[0];
            dcof[0] = 1.0;
            for (int k = 3; k <= 2 * N; ++k)
            {
                dcof[k] = dcof[2 * k - 2];
            }
            std::copy_n(dcof.begin(), 2 * N + 1, d.begin());
        }

    protected:
        const FreqProperty freqtype;
        std::vector<double> c;
        std::vector<double> d;

        mutable std::deque<double> samples;
        mutable std::deque<double> results;
        // what is best way to deal with interv variables ?
    };

    template <size_t N>
    class MovAvgFilter : public FIRFilter<N>
    {
    public:
        MovAvgFilter() : FIRFilter<N>(0.0, 0.0, FreqProperty::LowPass, FIRType::Boxcar)
        {
            FIRFilter<N>::h.fill(1.0 / N);
        }
    };

    template <size_t M, size_t NL, size_t NR, size_t NP = NL + NR + 1>
    class SGFilter : public FIRFilter<NP>
    {
    public:
        SGFilter(int derivative_order = 0) : FIRFilter<NP>(0.0, 0.0, FreqProperty::LowPass, FIRType::Boxcar)
        {
            savgo(FIRFilter<NP>::h, derivative_order);
        }

    private:
        void savgo(MatrixS<NP, 1> &c, int ld)
        {
            int nl = NL;
            int nr = NR;
            int np = NP;

            MatrixS<M + 1, M + 1> a;
            MatrixS<M + 1, 1> b;
            for (int ipj = 0; ipj <= (M << 1); ipj++)
            {
                auto sum = (ipj ? 0.0 : 1.0);
                for (auto k = 1; k <= nr; k++)
                {
                    sum += pow(double(k), double(ipj));
                }
                for (int k = 1; k <= nl; k++)
                {
                    sum += pow(double(-k), double(ipj));
                }
                auto mm = std::min(ipj, 2 * int(M) - ipj);
                for (auto imj = -mm; imj <= mm; imj += 2)
                {
                    a((ipj + imj) / 2, (ipj - imj) / 2) = sum;
                }
            }
            b[ld] = 1.0;
            b = linsolve<Factorization::LU>(a, b).x;
            for (auto k = -nl; k <= nr; k++)
            {
                auto sum = b[0];
                auto fac = 1.0;
                for (auto mm = 1; mm <= M; mm++)
                {
                    sum += b[mm] * (fac *= k);
                }
                // auto kk = (np - k) % np;
                c[k + nl] = sum;
            }
        }
    };

} // namespace ppx

#endif