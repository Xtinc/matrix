#include "gtest/gtest.h"
#include "signals.hpp"

using namespace ppx;

bool compare_vector(const std::vector<double> &expected, const std::vector<double> &input, double threshold = 10e-5)
{
    return std::equal(expected.begin(), expected.end(), input.begin(), [threshold](double a, double b)
                      { auto diff=fabs(a - b);
                      if(diff<threshold){
                        return true;
                      }else{
                        std::cout<<"not equal between "<<a<<":"<<b<<std::endl;
                        return false;
                      } });
}

void test_filter(const Filter &flt, const std::vector<double> &input, const std::vector<double> &expected)
{
    std::vector<double> results, velocity, expects_vel;
    for (auto i : input)
    {
        results.push_back(flt(i));
        velocity.push_back(flt.diff());
    }
    expects_vel.push_back(0.0);
    for (size_t i = 1; i < expected.size(); i++)
    {
        expects_vel.push_back(expected[i] - expected[i - 1]);
    }

    EXPECT_TRUE(compare_vector(expected, results));
    EXPECT_TRUE(compare_vector(velocity, expects_vel));
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
    std::vector<double> results{0.0, 0.0680, 1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000};
    FIRFilter<3, FreqProperty::LowPass> flt(0.1, FIRType::Hamming, false);
    EXPECT_TRUE(compare_vector(expected, flt.coff_b()));
    test_filter(flt, datas2, results);
}

TEST_F(STA_TestCase, IIRFilter)
{
    std::vector<double> expected_b{8.57655707325943e-06, 5.14593424395566e-05, 0.000128648356098892, 0.000171531141465189, 0.000128648356098892, 5.14593424395566e-05, 8.57655707325943e-06};
    std::vector<double> expected_a{1, -4.78713549885213, 9.64951772872191, -10.4690788925439, 6.44111188100806, -2.12903875003045, 0.295172431349155};
    std::vector<double> results{0, 8.57655707325942e-06, 0.000109669597409407, 0.000699540295571364, 0.00299783621920436, 0.00980009466680698, 0.0262688145029736, 0.0604915116363023};
    IIRFilter<6, FreqProperty::LowPass> flt(0.1, false);
    EXPECT_TRUE(compare_vector(expected_a, flt.coff_a()));
    EXPECT_TRUE(compare_vector(expected_b, flt.coff_b()));
    test_filter(flt, datas2, results);
}

TEST_F(STA_TestCase, MovAvgFlt)
{
    MovAvgFilter<2> flt1;
    std::vector<double> expects{0, 1.0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
    test_filter(flt1, datas2, expects);

    MovAvgFilter<2> flt2(false);
    expects = {0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
    test_filter(flt2, datas2, expects);

    MovAvgFilter<5> flt3;
    expects = {0, 1, 2, 3, 4, 3, 4, 5};
    test_filter(flt3, datas2, expects);

    MovAvgFilter<5> flt4(false);
    expects = {0, 0.2, 0.6, 1.2, 2, 3, 4, 5};
    test_filter(flt4, datas2, expects);
}