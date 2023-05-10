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

    template <size_t N>
    class ButterWorthFilter
    {
        struct SOS
        {
            double b0;
            double b1;
            double b2;
            double a1;
            double a2;
            double z1 = 0.0;
            double z2 = 0.0;

            double operator()(double input)
            {
                double out = input * b0 + z1;
                z1 = input * b1 - a1 * out + z2;
                z2 = input * b2 - a2 * out;
                return out;
            }

            void reset()
            {
                z1 = 0.0;
                z2 = 0.0;
            }
        };

    public:
        ButterWorthFilter(double f1, double f2, FreqProperty k_type = FreqProperty::LowPass)
        {
            // Need to pre-warp the frequency for the analog -> digital conversion
            double freq1 = 2.0 * std::tan(PI * f1 / 2.0);
            double freq2 = 2.0 * std::tan(PI * f2 / 2.0);

            double centerFreq = freq1;
            double bandwidth = 0.0;

            // Butterworth filter has N poles and no zeros. The poles all live right on the unit circle.
            constexpr std::complex<double> k_i(0, 1);
            constexpr std::complex<double> k_iPi(0, PI);
            std::vector<std::complex<double>> zeros;
            std::vector<std::complex<double>> polePairs; // a pole and its implied conjugate
            // Odd-order filters (3, 5, ...)
            bool hasOddPole = (N & 1) != 0;
            std::complex<double> oddPole = -1;
            // The poles (rounded down) get calculated evenly spaced in the left side of the space (negative real side).
            for (int i = 1; i < N; i += 2)
            {
                polePairs.push_back(std::exp(k_iPi * (double(i) / double(2 * N) + 0.5)));
            }
            double gain = 1.0;
            switch (k_type)
            {
            case FreqProperty::LowPass:
                // multiply the poles by the cutoff frequency
                for (auto &polePair : polePairs)
                {
                    polePair *= centerFreq;
                }
                if (hasOddPole)
                {
                    oddPole *= centerFreq;
                }
                // The gain simply adjusts by the center frequency power.
                gain = std::pow(centerFreq, N);
                break;
            case FreqProperty::HighPass:
                // flip the pole's value (via a complex divide) while multiplying in the center frequency
                for (auto &polePair : polePairs)
                {
                    polePair = std::complex<double>(centerFreq) / polePair;

                    // High pass has zeros at 0 (one per pole, which means two per pole pair)
                    zeros.push_back(0.0);
                    zeros.push_back(0.0);
                }
                if (hasOddPole)
                {
                    oddPole = std::complex<double>(centerFreq) / oddPole;
                    zeros.push_back(0.0);
                }
                break;
            case FreqProperty::BandPass:
            {
                // Band passes will end up with 2*order poles (so "order" pole pairs)
                // Calculate the center frequency and the bandwidth
                centerFreq = std::sqrt(freq1 * freq2);
                bandwidth = freq2 - freq1;
                std::vector<std::complex<double>> origPolePairs{polePairs};
                polePairs = {};
                for (auto &polePair : origPolePairs)
                {
                    // each pole gets modified by the following formula (leaving us still with these values and their conjugates, so the conjugate pairing here still works)
                    auto a = 0.5 * polePair * bandwidth;
                    auto b = 0.5 * std::sqrt(SQR(bandwidth) * SQR(polePair) - 4 * SQR(centerFreq));
                    polePairs.push_back(a + b);
                    polePairs.push_back(a - b);
                }
                // The odd pole itself becomes now a pair (new pole and conjugate), so stop factoring the odd pole in.
                if (hasOddPole)
                {
                    polePairs.push_back(0.5 * bandwidth * oddPole + 0.5 * std::sqrt(SQR(bandwidth) * SQR(oddPole) - 4 * SQR(centerFreq)));
                    hasOddPole = false;
                }

                for (int i = 0; i < N; i++)
                {
                    zeros.push_back(0.0);
                }
                gain = std::pow(bandwidth, N);
            }
            break;
            case FreqProperty::BandStop:
            {
                // it's a divide by the pole instead of a multiply.
                centerFreq = std::sqrt(freq1 * freq2);
                bandwidth = freq2 - freq1;
                std::vector<std::complex<double>> origPolePairs{polePairs};
                polePairs = {};
                for (auto &polePair : origPolePairs)
                {
                    auto a = 0.5 * bandwidth / polePair;
                    auto b = 0.5 * std::sqrt(SQR(bandwidth) / SQR(polePair) - 4 * SQR(centerFreq));
                    polePairs.push_back(a + b);
                    polePairs.push_back(a - b);
                }
                if (hasOddPole)
                {
                    polePairs.push_back(0.5 * bandwidth / oddPole + 0.5 * std::sqrt(SQR(bandwidth) / SQR(oddPole) - 4 * SQR(centerFreq)));
                    hasOddPole = false;
                }
                for (int i = 0; i < N; i++)
                {
                    zeros.push_back({0.0, centerFreq});
                    zeros.push_back({0.0, -centerFreq});
                }
            }
            break;
            }

            // map from the s-plane to the z-plane by bilinear transform
            for (auto &zero : zeros)
            {
                // Adjust the gain before modifying the zero
                gain *= std::abs(2.0 - zero);
                zero = (2.0 + zero) / (2.0 - zero);
            }
            for (auto &polePair : polePairs)
            {
                // Each polePair is a conjugate pair so factor each length into the gain twice
                double gainMod = std::abs(2.0 - polePair);
                gain /= gainMod * gainMod;
                polePair = (2.0 + polePair) / (2.0 - polePair);
            }
            if (hasOddPole)
            {
                gain /= std::abs(2.0 - oddPole);
                oddPole = (2.0 + oddPole) / (2.0 - oddPole);
            }
            // Pad out the zero array with -1*s
            auto zeroCount = polePairs.size() * 2 + (hasOddPole ? 1 : 0);
            for (auto i = zeros.size(); i < zeroCount; i++)
            {
                zeros.push_back(-1.0);
            }
            // Finally, generate second order sections from the poles
            if (hasOddPole)
            {
                // Odd filters have a special second order section that's really a first-order section.
                coff.push_back(SOS{
                    gain,
                    -zeros[zeros.size() - 1].real() * gain,
                    0.0,
                    -oddPole.real(),
                    0.0});
                gain = 1.0; // Only the first SOS in the set needs to have the gain factored in, so set the gain to 1 so it doesn't affect any more
            }

            for (unsigned int i = 0; i < polePairs.size(); i++)
            {
                auto poleA = polePairs[i];
                auto zeroA = zeros[i * 2 + 0];
                auto zeroB = zeros[i * 2 + 1];
                coff.push_back(SOS{
                    gain,
                    -(zeroA + zeroB).real() * gain,
                    (zeroA * zeroB).real() * gain,
                    -2.0 * poleA.real(),                                         // This is technically -(pole + conj(pole)), which simplifies down to -2*real
                    poleA.real() * poleA.real() + poleA.imag() * poleA.imag()}); // This is (pole * conj(pole))

                gain = 1.0; // Only the first SOS in the set needs to have the gain factored in, so set the gain to 1 so it doesn't affect any more
            }
        }

        double operator()(double new_sample) const
        {
            for (auto &i : coff)
            {
                new_sample = i(new_sample);
            }
            return new_sample;
        }

    private:
        mutable std::vector<SOS> coff;
    };

} // namespace ppx

#endif