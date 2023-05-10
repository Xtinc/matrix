#include "gtest/gtest.h"
#include "signals.hpp"

using namespace ppx;

class STA_TestCase : public ::testing::Test
{
public:
    STA_TestCase()
    {
        std::normal_distribution<double> dis;
        std::default_random_engine eng;
        for (int i = 0; i < 200 * PI; ++i)
        {
            datas.emplace_back(10 * sin((double)i / 100) + dis(eng));
            //            datas.emplace_back(10 * sin((double) i / 100) + sin(i));
        }
    };
    std::vector<double> datas;
};

TEST_F(STA_TestCase, MovAvgFlt)
{
    FIRFilter<19> flt(0.01, 0.0, FreqProperty::LowPass, WindowType::Hamming);
    ButterWorthFilter<3> but(0.1, 0.0, FreqProperty::LowPass);
    // SGFilter<2, 12, 2> flt;
    for (auto &&elem : datas)
    {
        std::cout << elem << std::endl;
    }
    std::cout << "_____________________________" << std::endl;

    for (auto &&elem : datas)
    {
        std::cout << but(elem) << std::endl;
    }

    // std::cout << flt.coff() << std::endl;
}