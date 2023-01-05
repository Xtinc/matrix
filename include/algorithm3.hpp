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
    class multi_normal_distribution
    {
    public:
        using samples = std::vector<Matrix<N, 1>>;
        multi_normal_distribution() : m_cov(Matrix<N, N>::eye()) {}
        multi_normal_distribution(const Matrix<N, 1> &mu, const Matrix<N, N> &sigma)
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
    class mixed_normal_distribution
    {
    public:
        using samples = std::vector<Matrix<N, 1>>;

        template <typename T>
        void push_back(T &&elem, double prior)
        {
            m_guassian.push_back(std::forward<T>(elem));
            m_prior.push_back(prior);
        }
        double pdf(const Matrix<N, 1> &x)
        {
            double sum = 0.0;
            for (size_t i = 0; i < m_prior.size(); i++)
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
            //     return {};
            // }
            std::random_device rd{};
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(0.0, 1.0);
            double sample = dis(gen);
            double sum = m_prior.front();
            size_t idx = 0;
            size_t idx_max = m_prior.size();
            while (sample > sum && idx < idx_max)
            {
                sum += m_prior[++idx];
            }
            return m_guassian[idx]();
        }

        const multi_normal_distribution<N> &distribution(size_t idx) const
        {
            return m_guassian.at(idx);
        }
        multi_normal_distribution<N> &distribution(size_t idx)
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

        void loglikehood(const samples &data)
        {
            auto n = data.size();
            auto c = m_prior.size();

            for (size_t iter = 0; iter < 100; ++iter)
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
                    Matrix<N, 1> sum_m = std::inner_product(a[k].begin(), a[k].end(), data.begin(), Matrix<N, 1>()) / n;
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
            }
        }

    private:
        std::vector<multi_normal_distribution<N>> m_guassian;
        std::vector<double> m_prior;
    };
}

#endif