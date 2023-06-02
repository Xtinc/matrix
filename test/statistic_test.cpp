#include "gtest/gtest.h"
#include "signals.hpp"

using namespace ppx;

class STA_TestCase : public ::testing::Test
{
public:
    STA_TestCase()
    {
        for (int i = 0; i < 20 * PI; ++i)
        {
            datas1.emplace_back(10 * sin((double)i / 100) + sin((double)i / 10));
        }
        for (int i = 0; i < 8; i++)
        {
            datas2.emplace_back(i);
        }
    };
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
        std::cout << dis.pdf(-1 + i * 0.01) << std::endl;
    }
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

TEST_F(STA_TestCase, FIRFilter_TEST)
{
}