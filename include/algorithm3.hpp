#ifndef VVERY_SIMPLE_ALGORITHM3_HEADER
#define VVERY_SIMPLE_ALGORITHM3_HEADER

#include "algorithm2.hpp"
#include <random>

namespace ppx
{
    template <typename T>
    double mean(const T &vec)
    {
        return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    }

    template <typename T>
    double var(const T &vec)
    {
        auto e = mean(vec);
        double sum = 0.0;
        for (size_t i = 0; i < vec.size(); i++)
        {
            auto tmp = vec[i] - e;
            sum += tmp * tmp;
        }
        return sum / (vec.size() - 1);
    }

    template <size_t N>
    class MultiNormalDistribution
    {
    public:
        using samples = std::vector<Matrix<N, 1>>;
        MultiNormalDistribution() : m_cov(Matrix<N, N>::eye()) {}
        MultiNormalDistribution(const Matrix<N, 1> &mu, const Matrix<N, N> &sigma)
            : m_mean(mu), m_cov(sigma) {}
        Matrix<N, 1> operator()() const
        {
            Matrix<N, 1> x;
            std::random_device rd{};
            std::mt19937 gen{rd()};
            std::normal_distribution<> d{0, 1};
            auto eigsys = eig<eigensystem::SymValAndVec>(m_cov);
            Matrix<N, N> diag;
            for (size_t i = 0; i < N; i++)
            {
                x[i] = d(gen);
                diag(i, i) = sqrt(eigsys.val[i]);
            }
            return eigsys.vec * diag * x + m_mean;
        }

        double pdf(const Matrix<N, 1> &x) const
        {
            int n = static_cast<int>(N);
            Matrix<2, 1> normalized_mu = x - m_mean;
            double quadform = (normalized_mu.T() * m_cov.I() * normalized_mu)[0];
            double norm = pow(sqrt(2 * gl_rep_pi), -n) * pow(m_cov.det(), -0.5);
            return norm * exp(-0.5 * quadform);
        }
        const Matrix<N, 1> &mean() const
        {
            return m_mean;
        }
        Matrix<N, 1> &mean()
        {
            return m_mean;
        }
        const Matrix<N, N> &covariance() const
        {
            return m_cov;
        }
        Matrix<N, N> &covariance()
        {
            return m_cov;
        }

        void loglikehood(const samples &data)
        {
            auto n = data.size();
            auto sum_m = std::accumulate(data.begin(), data.end(), Matrix<N, 1>());
            m_mean = sum_m / n;
            Matrix<N, N> sum_s;
            for (size_t i = 0; i < n; i++)
            {
                Matrix<N, N> tpx = data.at(i) - m_mean;
                sum_s += tpx * tpx.T();
            }
            m_cov = sum_s / n;
        }

    private:
        Matrix<N, 1> m_mean;
        Matrix<N, N> m_cov;
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
        using samples = std::vector<Matrix<N, 1>>;

        double pdf(const Matrix<N, 1> &x)
        {
            double sum = 0.0;
            for (size_t i = 0; i < K; i++)
            {
                sum += m_guassian[i].pdf(x) * m_prior[i];
            }
            return sum;
        }

        Matrix<N, 1> operator()() const
        {
            // auto total_p = std::accumulate(m_prior.begin(), m_prior.end(), 0.0);
            // how?
            // if (!details::is_same(total_p, 1.0))
            // {
            // return {};
            // }
            std::random_device rd{};
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);
            double sample = dis(gen);
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

            while (residual > gl_rep_eps && its < ITMAX)
            {
                std::vector<std::vector<double>> a(c, std::vector<double>(n));
                for (size_t k = 0; k < c; k++)
                {
                    std::transform(data.begin(), data.end(), a[k].begin(), [&](const Matrix<N, 1> &x)
                                   { return m_prior[k] * m_guassian[k].pdf(x) / pdf(x); });
                }
                for (size_t k = 0; k < c; k++)
                {
                    auto sum_g = std::accumulate(a[k].begin(), a[k].end(), 0.0);
                    m_prior[k] = sum_g / n;
                    Matrix<N, 1> sum_m = std::inner_product(a[k].begin(), a[k].end(), data.begin(), Matrix<N, 1>());
                    sum_m /= sum_g;
                    Matrix<N, N> sum_s;
                    for (size_t i = 0; i < n; i++)
                    {
                        Matrix<N, 1> tpx = data.at(i) - sum_m;
                        sum_s += tpx * tpx.T() * a[k][i];
                    }
                    m_guassian[k].mean() = sum_m;
                    m_guassian[k].covariance() = sum_s / sum_g;
                }
                double tpp = std::accumulate(data.begin(), data.end(), 0.0, [&](double y0, const Matrix<N, 1> &x)
                                             { return y0 + log(pdf(x)); });
                residual = fabs(last_p - tpp);
                last_p = tpp;
                ++its;
            }
        }

    private:
        std::array<dist, K> m_guassian;
        std::array<double, K> m_prior;
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
}

#endif