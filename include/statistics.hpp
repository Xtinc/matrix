#ifndef VVERY_SIMPLE_STATISTICS_HEADER
#define VVERY_SIMPLE_STATISTICS_HEADER

#include "optimization.hpp"
#include <random>
#include <functional>

namespace ppx
{
    // Moments of a Distribution: Mean, Variance, Skewness
    template <typename T>
    auto avg(const T &vec)
    {
        using value_t = typename T::value_type;
        return static_cast<value_t>(std::accumulate(vec.begin(), vec.end(), value_t{}) / vec.size());
    }

    template <typename T>
    double var(const T &vec)
    {
        assert(vec.size() >= 1);
        using value_t = typename T::value_type;
        value_t e = avg(vec);
        double sum = 0.0;
        for (size_t i = 0; i < vec.size(); i++)
        {
            value_t tmp = vec[i] - e;
            sum += tmp * tmp;
        }
        return sum / (vec.size() - 1);
    }

    template <size_t M, size_t N>
    void random(MatrixS<M, N> &mat, double lo = -MAX_SP, double hi = MAX_SP)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> distr(lo, hi);
        for (auto &i : mat)
        {
            i = distr(gen);
        }
    }

    template <size_t N>
    class MultiGaussianDistribution
    {
    public:
        using samples = std::vector<MatrixS<N, 1>>;

        MultiGaussianDistribution() : m_cov(eye<N>()), m_gen(std::random_device{}()) {}

        MultiGaussianDistribution(const MatrixS<N, 1> &mu, const MatrixS<N, N> &sigma)
            : m_mean(mu), m_cov(sigma), m_gen(std::random_device{}()) {}

        MatrixS<N, 1> operator()() const
        {
            MatrixS<N, 1> x;
            std::normal_distribution<> d{0, 1};
            EigenValue<N> eigsys(m_cov);
            MatrixS<N, N> diag;
            for (size_t i = 0; i < N; i++)
            {
                x[i] = d(m_gen);
                diag(i, i) = sqrt(eigsys.d[i]);
            }
            return eigsys.z * diag * x + m_mean;
        }

        double pdf(const MatrixS<N, 1> &x) const
        {
            int n = static_cast<int>(N);
            MatrixS<N, 1> normalized_mu = x - m_mean;
            // double quadform = (normalized_mu.T() * m_cov.I() * normalized_mu)[0];
            double quadform = inner_product(normalized_mu, normalized_mu, m_cov.I());
            double norm = std::pow(std::sqrt(2 * PI), -n) / std::sqrt(m_cov.det());
            return norm * std::exp(-0.5 * quadform);
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

        void fit(const samples &data)
        {
            auto n = data.size();
            auto sum_m = std::accumulate(data.begin(), data.end(), MatrixS<N, 1>());
            m_mean = sum_m / n;
            MatrixS<N, N> sum_s;
            for (size_t i = 0; i < n; i++)
            {
                MatrixS<N, 1> tpx = data.at(i) - m_mean;
                sum_s += tpx * tpx.T();
            }
            m_cov = sum_s / n;
        }

    private:
        MatrixS<N, 1> m_mean;
        MatrixS<N, N> m_cov;
        mutable std::mt19937 m_gen;
    };

    template <size_t N>
    using MVN = MultiGaussianDistribution<N>;

    template <size_t N>
    std::ostream &operator<<(std::ostream &os, const MultiGaussianDistribution<N> &self)
    {
        os << "MultiNormalDistribution<" << N << ">:\n"
           << "mean = \t" << self.mean() << "\n"
           << "cov  = \t" << self.covariance();
        return os;
    }

    template <size_t N, size_t K>
    class MixedGaussianDistribution
    {
    public:
        using dist = MultiGaussianDistribution<N>;
        using samples = std::vector<MatrixS<N, 1>>;

        MixedGaussianDistribution() : m_prior(), m_gen(std::random_device{}()) {}

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
                sum += m_prior[idx];
                ++idx;
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

        double prior(size_t idx) const
        {
            return m_prior.at(idx);
        }

        void setcomp(size_t idx, const dist &d, double p)
        {
            m_guassian[idx] = d;
            m_prior[idx] = p;
        }

        void fit(const samples &data)
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
                                             { return y0 + std::log(pdf(x)); });
                residual = fabs(last_p - tpp);
                last_p = tpp;
                ++its;
            }
        }

        template <size_t L>
        std::enable_if_t<(L < N), std::vector<MatrixS<N - L, 1>>>
        predict(std::vector<MatrixS<L, 1>> &known, const std::array<size_t, L> &idx) const
        {
            std::array<bool, N> xory;
            for (size_t i = 0; i < N; i++)
            {
                xory[i] = std::any_of(idx.cbegin(), idx.cend(), [i](size_t a)
                                      { return a == i; });
            }

            std::vector<MatrixS<N - L, 1>> result;
            std::vector<MatrixS<L, 1>> u_x_k;
            std::vector<MatrixS<L, L>> cov_xx_k;
            std::vector<MatrixS<N - L, 1>> u_y_k;
            std::vector<MatrixS<N - L, L>> reg_k;

            for (size_t k = 0; k < K; k++)
            {
                const auto &u_total = m_guassian[k].mean();
                const auto &cov_total = m_guassian[k].covariance();
                MatrixS<L, 1> u_x;
                MatrixS<N - L, 1> u_y;
                MatrixS<L, L> cov_xx;
                MatrixS<N - L, L> cov_yx;

                size_t idx1 = 0;
                size_t idx2 = 0;
                for (size_t i = 0; i < N; i++)
                {
                    if (xory[i])
                    {
                        u_x[idx1] = u_total[i];
                        size_t jdx1 = 0;
                        for (size_t j = 0; j < N; j++)
                        {
                            if (xory[j])
                            {
                                cov_xx(idx1, jdx1) = cov_total(i, j);
                                jdx1++;
                            }
                        }
                        idx1++;
                    }
                    else
                    {
                        u_y[idx2] = u_total[i];
                        size_t jdx2 = 0;
                        for (size_t j = 0; j < N; j++)
                        {
                            if (xory[j])
                            {
                                cov_yx(idx2, jdx2) = cov_total(i, j);
                                jdx2++;
                            }
                        }
                        idx2++;
                    }
                }
                u_x_k.push_back(std::move(u_x));
                u_y_k.push_back(std::move(u_y));
                reg_k.push_back(cov_yx * cov_xx.I());
                cov_xx_k.push_back(std::move(cov_xx));
            }

            for (const auto &elem : known)
            {
                MatrixS<K, 1> posterior_coff{};
                for (size_t i = 0; i < K; i++)
                {
                    MVN<L> marginal_gmm(u_x_k[i], cov_xx_k[i]);
                    posterior_coff[i] = m_prior[i] * marginal_gmm.pdf(elem);
                }
                posterior_coff /= sum(posterior_coff.data(), K);

                MatrixS<N - L, 1> tmp;
                for (size_t i = 0; i < K; i++)
                {
                    MatrixS<N - L, 1> u = u_y_k[i] + reg_k[i] * (elem - u_x_k[i]);
                    tmp = tmp + posterior_coff[i] * u;
                }
                result.push_back(std::move(tmp));
            }

            return result;
        }

    private:
        std::array<dist, K> m_guassian;
        std::array<double, K> m_prior;
        mutable std::mt19937 m_gen;
        size_t ITMAX = 200;
    };

    template <size_t N, size_t K>
    std::ostream &operator<<(std::ostream &os, const MixedGaussianDistribution<N, K> &self)
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
    using GMM = MixedGaussianDistribution<N, K>;
}
#endif