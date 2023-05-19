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
    // FIRFilter<19> flt(0.01, 0.0, FreqProperty::LowPass, FIRType::Hamming);
    IIRFilter<5, FreqProperty::BandStop> iir(0.1, 0.2);
    Filter flt(FreqProperty::LowPass);
    flt.coff_a() = {1, -3.98454311961234, 6.43486709027587, -5.25361517035227, 2.16513290972413, -0.359928245063557};
    flt.coff_b() = {5.97957803700031e-05, 0.000298978901850016, 0.000597957803700031, 0.000597957803700031, 0.000298978901850016, 5.97957803700031e-05};

    // SGFilter<2, 12, 2> flt;
    for (auto &&elem : datas)
    {
        std::cout << elem << std::endl;
    }
    std::cout << "_____________________________" << std::endl;

    for (auto &&elem : datas)
    {
        std::cout << flt(elem) << std::endl;
    }

    // std::cout << flt.coff() << std::endl;
}