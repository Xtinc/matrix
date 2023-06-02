#ifndef VVERY_SIMPLE_STATISTICS_HEADER
#define VVERY_SIMPLE_STATISTICS_HEADER

#include "optimization.hpp"
#include <random>
#include <functional>

namespace ppx
{
    // Moments of a Distribution: Mean, Variance, Skewness
    template <typename T>
    auto mean(const T &vec)
    {
        using value_t = typename T::value_type;
        return static_cast<value_t>(std::accumulate(vec.begin(), vec.end(), value_t()) / vec.size());
    }

    template <typename T>
    double var(const T &vec)
    {
        using value_t = typename T::value_type;
        value_t e = mean(vec);
        double sum = 0.0;
        for (size_t i = 0; i < vec.size(); i++)
        {
            value_t tmp = vec[i] - e;
            sum += tmp * tmp;
        }
        return sum / (vec.size() - 1);
    }

    template <size_t N>
    class MultiNormalDistribution
    {
    public:
        using samples = std::vector<MatrixS<N, 1>>;

        MultiNormalDistribution() : m_cov(MatrixS<N, N>::eye()), m_gen(std::random_device{}()) {}

        MultiNormalDistribution(const MatrixS<N, 1> &mu, const MatrixS<N, N> &sigma)
            : m_mean(mu), m_cov(sigma) {}

        MatrixS<N, 1> operator()() const
        {
            MatrixS<N, 1> x;
            std::normal_distribution<> d{0, 1};
            auto eigsys = eig<EigenSystem::SymValAndVec>(m_cov);
            MatrixS<N, N> diag;
            for (size_t i = 0; i < N; i++)
            {
                x[i] = d(m_gen);
                diag(i, i) = sqrt(eigsys.val[i]);
            }
            return eigsys.vec * diag * x + m_mean;
        }

        double pdf(const MatrixS<N, 1> &x) const
        {
            int n = static_cast<int>(N);
            MatrixS<2, 1> normalized_mu = x - m_mean;
            double quadform = (normalized_mu.T() * m_cov.I() * normalized_mu)[0];
            double norm = pow(sqrt(2 * PI), -n) * pow(m_cov.det(), -0.5);
            return norm * exp(-0.5 * quadform);
        }

        const MatrixS<N, 1> &mean() const
        {
            return m_mean;
        }

        MatrixS<N, 1> &mean()
        {
            return m_mean;
        }

        const MatrixS<N, N> &covariance() const
        {
            return m_cov;
        }

        MatrixS<N, N> &covariance()
        {
            return m_cov;
        }

        void loglikehood(const samples &data)
        {
            auto n = data.size();
            auto sum_m = std::accumulate(data.begin(), data.end(), MatrixS<N, 1>());
            m_mean = sum_m / n;
            MatrixS<N, N> sum_s;
            for (size_t i = 0; i < n; i++)
            {
                MatrixS<N, N> tpx = data.at(i) - m_mean;
                sum_s += tpx * tpx.T();
            }
            m_cov = sum_s / n;
        }

    private:
        MatrixS<N, 1> m_mean;
        MatrixS<N, N> m_cov;
        mutable std::mt19937 m_gen;
    };

    template <>
    class MultiNormalDistribution<1>
    {
    public:
        using samples = std::vector<double>;

        MultiNormalDistribution() : m_mean(0.0), m_cov(1.0), m_gen(std::random_device{}()) {}

        MultiNormalDistribution(double mu, double sigma)
            : m_mean(mu), m_cov(sigma) {}

        double operator()() const
        {
            std::normal_distribution<> d{m_mean, m_cov};
            return d(m_gen);
        }

        double pdf(double x) const
        {
            auto normalized_mu = (x - m_mean) / m_cov;
            double quadform = normalized_mu * normalized_mu;
            double norm = 1.0 / sqrt(2 * PI * m_cov * m_cov);
            return norm * exp(-0.5 * quadform);
        }

        const double &mean() const
        {
            return m_mean;
        }

        double &mean()
        {
            return m_mean;
        }

        const double &covariance() const
        {
            return m_cov;
        }

        double &covariance()
        {
            return m_cov;
        }

        void loglikehood(const samples &data)
        {
            auto n = data.size();
            auto sum_m = std::accumulate(data.begin(), data.end(), 0.0);
            m_mean = sum_m / n;
            double sum_s = 0.0;
            for (size_t i = 0; i < n; i++)
            {
                auto tpx = data.at(i) - m_mean;
                sum_s += tpx * tpx;
            }
            m_cov = sum_s / n;
        }

    private:
        double m_mean;
        double m_cov;
        mutable std::mt19937 m_gen;
    };

    template <size_t N>
    std::ostream &operator<<(std::ostream &os, const MultiNormalDistribution<N> &self)
    {
        os << "MultiNormalDistribution<" << N << ">:\n"
           << "mean = \t" << self.mean()
           << "cov  = \t" << self.covariance();
        return os;
    }

    template <size_t N, size_t K>
    class MixedNormalDistribution
    {
    public:
        using dist = MultiNormalDistribution<N>;
        using samples = std::vector<MatrixS<N, 1>>;

        MixedNormalDistribution() : m_prior(), m_gen(std::random_device{}()) {}

        double pdf(const MatrixS<N, 1> &x)
        {
            double sum = 0.0;
            for (size_t i = 0; i < K; i++)
            {
                sum += m_guassian[i].pdf(x) * m_prior[i];
            }
            return sum;
        }

        MatrixS<N, 1> operator()() const
        {
            std::uniform_real_distribution<> dis(0.0, 1.0);
            double sample = dis(m_gen);
            double sum = m_prior.front();
            size_t idx = 0;
            while (sample > sum && idx < K)
            {
                sum += m_prior[++idx];
            }
            return m_guassian[idx]();
        }

        const dist &distribution(size_t idx) const
        {
            return m_guassian.at(idx);
        }

        dist &distribution(size_t idx)
        {
            return m_guassian.at(idx);
        }

        const double &prior(size_t idx) const
        {
            return m_prior.at(idx);
        }

        double &prior(size_t idx)
        {
            return m_prior.at(idx);
        }

        void setcomp(size_t idx, const dist &d, double p)
        {
            m_guassian[idx] = d;
            m_prior[idx] = p;
        }

        void loglikehood(const samples &data)
        {
            constexpr auto c = K;
            auto n = data.size();
            double residual = 1.0;
            double last_p = 0.0;
            auto its = 0u;

            while (residual > EPS_SP && its < ITMAX)
            {
                std::vector<std::vector<double>> a(c, std::vector<double>(n));
                for (size_t k = 0; k < c; k++)
                {
                    std::transform(data.begin(), data.end(), a[k].begin(),
                                   [&](const MatrixS<N, 1> &x)
                                   { return m_prior[k] * m_guassian[k].pdf(x) / pdf(x); });
                }
                for (size_t k = 0; k < c; k++)
                {
                    auto sum_g = std::accumulate(a[k].begin(), a[k].end(), 0.0);
                    m_prior[k] = sum_g / n;
                    MatrixS<N, 1> sum_m = std::inner_product(a[k].begin(), a[k].end(), data.begin(), MatrixS<N, 1>());
                    sum_m /= sum_g;
                    MatrixS<N, N> sum_s;
                    for (size_t i = 0; i < n; i++)
                    {
                        MatrixS<N, 1> tpx = data.at(i) - sum_m;
                        sum_s += tpx * tpx.T() * a[k][i];
                    }
                    m_guassian[k].mean() = sum_m;
                    m_guassian[k].covariance() = sum_s / sum_g;
                }
                double tpp = std::accumulate(data.begin(), data.end(), 0.0,
                                             [&](double y0, const MatrixS<N, 1> &x)
                                             { return y0 + log(pdf(x)); });
                residual = fabs(last_p - tpp);
                last_p = tpp;
                ++its;
            }
        }

    private:
        std::array<dist, K> m_guassian;
        std::array<double, K> m_prior;
        mutable std::mt19937 m_gen;
        size_t ITMAX = 200;
    };

    template <size_t N, size_t K>
    std::ostream &operator<<(std::ostream &os, const MixedNormalDistribution<N, K> &self)
    {
        os << "MixedNormalDistribution<" << N << ',' << K << ">:\n";
        for (size_t i = 0; i < K; i++)
        {
            os << self.prior(i) << " of all, component " << i << ":\n";
            os << self.distribution(i) << "\n";
        }
        return os;
    }

    template <size_t N, size_t K>
    class Kmeans
    {
    public:
        using samples = std::vector<MatrixS<N, 1>>;

        explicit Kmeans(const samples &data)
            : m_data(data), m_assign(m_data.size(), 0) {}

    private:
        const samples &m_data;
        std::vector<int> m_assign;
        std::array<int, K> m_count;
        std::array<MatrixS<N, 1>, K> m_mean;

    private:
        int estep()
        {
            int nchg = 0;
            int kmin = 0;
            m_count.fill(0);
            for (auto i = 0; i < m_data.size(); i++)
            {
                auto dmin = std::numeric_limits<double>::max();
                for (auto j = 0; j < K; j++)
                {
                    auto d = norm2(m_data[i] - m_mean[j]);
                    if (d < dmin)
                    {
                        dmin = d;
                        kmin = j;
                    }
                }
                if (kmin != m_assign[i])
                {
                    nchg++;
                }
                m_assign[i] = kmin;
                m_count[kmin]++;
            }
            return nchg;
        }

        void mstep()
        {
            for (const auto &elem : m_mean)
            {
                elem.fill(0);
            }

            for (auto i = 0; i < m_data.size(); i++)
            {
                m_mean[m_assign[i]] += m_data[i];
            }

            for (auto i = 0; i < K; i++)
            {
                if (m_count[i] > 0)
                {
                    m_mean[i] /= m_count[i];
                }
            }
        }
    };
}
#endif