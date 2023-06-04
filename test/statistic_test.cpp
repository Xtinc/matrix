#include "gtest/gtest.h"
#include "signals.hpp"

using namespace ppx;

bool compare_vector(const std::vector<double> &expected, const std::vector<double> &input)
{
    return std::equal(expected.begin(), expected.end(), input.begin(), [](double a, double b)
                      { return fabs(a - b) < EPS_SP; });
}

class STA_TestCase : public ::testing::Test
{
public:
    STA_TestCase() : rand(0, 0.1)
    {
        for (int i = 0; i < 100; ++i)
        {
            double sample_point = -PI + (double)i / 50 * PI;
            datas1.emplace_back(sin(sample_point) + rand());
        }
        for (int i = 0; i < 8; i++)
        {
            datas2.emplace_back(i);
        }
    };
    MultiNormalDistribution<1> rand;
    std::vector<double> datas1;
    std::vector<double> datas2;

    bool is_vec_equal(const std::vector<double> &a, const std::vector<double> &b)
    {
        return std::equal(a.begin(), a.end(), b.begin(), [](double a, double b)
                          { return ppx::details::is_same(a, b); });
    }
};

TEST_F(STA_TestCase, NormalDistribution)
{
    MultiNormalDistribution<1> dis(3.3, 5);
    for (size_t i = 0; i < 200; i++)
    {
        // std::cout << dis.pdf(-1 + i * 0.01) << std::endl;
    }
}

TEST_F(STA_TestCase, UDFFilter)
{
    Filter flt;
    flt.coff_a() = {1.0, -0.9};
    flt.coff_b() = {0.1};
}

TEST_F(STA_TestCase, FIRFilter)
{
    std::vector<double> expected{0.0679901667389975, 0.864019666522005, 0.0679901667389975};
    FIRFilter<3, FreqProperty::LowPass> flt(0.1);
    EXPECT_TRUE(compare_vector(expected, flt.coff_b()));
}

TEST_F(STA_TestCase, IIRFilter)
{
    std::vector<double> expected_b{8.57655707325943e-06, 5.14593424395566e-05, 0.000128648356098892, 0.000171531141465189, 0.000128648356098892, 5.14593424395566e-05, 8.57655707325943e-06};
    std::vector<double> expected_a{1, -4.78713549885213, 9.64951772872191, -10.4690788925439, 6.44111188100806, -2.12903875003045, 0.295172431349155};
    IIRFilter<6, FreqProperty::LowPass> flt(0.1);
    EXPECT_TRUE(compare_vector(expected_a, flt.coff_a()));
    EXPECT_TRUE(compare_vector(expected_b, flt.coff_b()));
}

TEST_F(STA_TestCase, MovAvgFlt)
{
    MovAvgFilter<2> flt2;
    std::vector<double> results;
    std::vector<double> expects{0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
    for (auto i : datas2)
    {
        results.push_back(flt2(i));
    }
    EXPECT_TRUE(is_vec_equal(expects, results));
    results.clear();

    MovAvgFilter<5> flt5;
    expects = {0, 0.2, 0.6, 1.2, 2, 3, 4, 5};
    for (auto i : datas2)
    {
        results.push_back(flt5(i));
    }
    EXPECT_TRUE(is_vec_equal(expects, results));
    results.clear();
}