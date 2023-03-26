#ifndef VVERY_SIMPLE_ALGORITHM3_HEADER
#define VVERY_SIMPLE_ALGORITHM3_HEADER

#include "statistics.hpp"

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
    }

    enum class FreqProperty : char
    {
        LowPass,
        HighPass,
        BandPass,
        BandStop
    };

    enum WindowType : char
    {
        Boxcar,
        Bartlett,
        Blackman,
        Hamming,
        Hanning
    };

    template <size_t N>
    class FIRFilter
    {
    public:
        FIRFilter(
            double f1, double f2, FreqProperty k_type = FreqProperty::LowPass, WindowType k_win = WindowType::Hamming)
            : freqtype(k_type), idx(0)
        {
            f1 /= 2.0;
            f2 /= 2.0;
            // Calculate the coefficient corresponding to the filter type
            switch (k_win)
            {
            case WindowType::Bartlett:
                winfunc = details::BartlettWindow;
                break;
            case WindowType::Blackman:
                winfunc = details::BlackmanWindow;
                break;
            case WindowType::Hanning:
                winfunc = details::HanningWindow;
                break;
            case WindowType::Hamming:
                winfunc = details::HammingWindow;
                break;
            case WindowType::Boxcar:
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
    class MovAvgFilter : public FIRFilter<N>
    {
    public:
        MovAvgFilter() : FIRFilter<N>(0.0, 0.0, FreqProperty::LowPass, WindowType::Boxcar)
        {
            FIRFilter<N>::h.fill(1.0 / N);
        }
    };

    template <size_t M, size_t NL, size_t NR, size_t NP = NL + NR + 1>
    class SGFilter : public FIRFilter<NP>
    {
    public:
        SGFilter(int derivative_order = 0) : FIRFilter<NP>(0.0, 0.0, FreqProperty::LowPass, WindowType::Boxcar)
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