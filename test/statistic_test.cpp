#include "gtest/gtest.h"
#include "statistics.hpp"

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
    FIRFilter<51> flt(0.42, 0.8, FIR::BandStop, HammingWindow);
    std::cout << flt.coff() << std::endl;
}