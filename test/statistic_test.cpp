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

TEST_F(STA_TestCase, FIRFilter_LP)
{
    std::vector<double> expected{0.0679901667389975, 0.864019666522005, 0.0679901667389975};
    std::vector<double> results{0.0, 0.0680, 1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000};
    FIRFilter<3, FreqProperty::LowPass> flt(0.1, FIRType::Hamming, false);
    EXPECT_TRUE(compare_vector(expected, flt.coff_b()));
    test_filter(flt, datas2, results);
}

TEST_F(STA_TestCase, FIRFilter_HP)
{
    std::vector<double> expected{6.65228083586735e-18, 0.00930428314501813, -0.0475777661344174, 0.122363546361145, -0.202246558429840, 0.237015691859159, -0.202246558429840, 0.122363546361145, -0.0475777661344174, 0.00930428314501813, 6.65228083586735e-18};
    std::vector<double> results{0, 6.65228083586735e-18, 0.00930428314501815, -0.0289691998443812, 0.0551208635273642, -0.0630356315307308, 0.0558235652703332, -0.0275637963584431};
    FIRFilter<11, FreqProperty::HighPass> flt(0.8, FIRType::Hamming, false);
    EXPECT_TRUE(compare_vector(expected, flt.coff_b()));
    test_filter(flt, datas2, results);
}

TEST_F(STA_TestCase, FIRFilter_BS)
{
    std::vector<double> expected{-0.00234328376094847, 0.00699997933273476, 0.0242756533581922, -0.0102228877948111, 0.962581077729665, -0.0102228877948111, 0.0242756533581922, 0.00699997933273476, -0.00234328376094847};
    std::vector<double> results{0, -0.00234328376094847, 0.00231341181083782, 0.0312457607408164, 0.0499552218759838, 1.03124576074082, 2.00231341181084, 2.99765671623905};
    FIRFilter<9, FreqProperty::BandStop> flt(0.4, 0.45, FIRType::Hamming, false);
    EXPECT_TRUE(compare_vector(expected, flt.coff_b()));
    test_filter(flt, datas2, results);
}

TEST_F(STA_TestCase, IIRFilter_LP)
{
    std::vector<double> expected_b{8.57655707325943e-06, 5.14593424395566e-05, 0.000128648356098892, 0.000171531141465189, 0.000128648356098892, 5.14593424395566e-05, 8.57655707325943e-06};
    std::vector<double> expected_a{1, -4.78713549885213, 9.64951772872191, -10.4690788925439, 6.44111188100806, -2.12903875003045, 0.295172431349155};
    std::vector<double> results{0, 8.57655707325942e-06, 0.000109669597409407, 0.000699540295571364, 0.00299783621920436, 0.00980009466680698, 0.0262688145029736, 0.0604915116363023};
    IIRFilter<6, FreqProperty::LowPass> flt(0.1, false);
    EXPECT_TRUE(compare_vector(expected_a, flt.coff_a()));
    EXPECT_TRUE(compare_vector(expected_b, flt.coff_b()));
    test_filter(flt, datas2, results);
}

TEST_F(STA_TestCase, IIRFilter_HP)
{
    std::vector<double> expected_b{0.107998403757355, -0.755988826301485, 2.26796647890446, -3.77994413150743, 3.77994413150743, -2.26796647890446, 0.755988826301485, -0.107998403757355};
    std::vector<double> expected_a{1, -2.78251380035720, 3.96680789197413, -3.40515039448060, 1.87586271598869, -0.650945698531145, 0.130852451560490, -0.0116627280491896};
    std::vector<double> results{0, 0.107998403757355, -0.239484969915386, -0.0147951165394243, 0.196590022304780, 0.127623335139356, -0.0635571388890397, -0.155957236529292};
    IIRFilter<7, FreqProperty::HighPass> flt(0.3, false);
    EXPECT_TRUE(compare_vector(expected_a, flt.coff_a()));
    EXPECT_TRUE(compare_vector(expected_b, flt.coff_b()));
    test_filter(flt, datas2, results);
}

TEST_F(STA_TestCase, IIRFilter_BS)
{
    std::vector<double> expected_b{0.809066862358596, -6.22565762049704, 26.1944285139720, -74.9689223087762, 160.994069793412, -271.116809361367, 367.666609439546, -406.316381735042, 367.666609439546, -271.116809361367, 160.994069793412, -74.9689223087762, 26.1944285139720, -6.22565762049704, 0.809066862358596};
    std::vector<double> expected_a{1, -7.46207919232620, 30.4441268607578, -84.4951944910069, 175.975358068302, -287.432627049382, 378.108892303000, -405.377130759820, 355.900301944903, -254.658587105406, 146.751918266687, -66.3241915216776, 22.4931625871694, -5.18935019671270, 0.654589187766792};
    std::vector<double> results{0, 0.809066862358596, 1.42979710302688, 2.20823887531852, 3.29062017183502, 4.54433190080803, 5.72722677360764, 6.71090557266292};
    IIRFilter<7, FreqProperty::BandStop> flt(0.3, 0.33, false);
    EXPECT_TRUE(compare_vector(expected_a, flt.coff_a()));
    EXPECT_TRUE(compare_vector(expected_b, flt.coff_b()));
    test_filter(flt, datas2, results);
}

TEST_F(STA_TestCase, IIRFilter_BP)
{
    std::vector<double> expected_b{1.22964988731449e-06, 0, -8.60754921120146e-06, 0, 2.58226476336044e-05, 0, -4.30377460560073e-05, 0, 4.30377460560073e-05, 0, -2.58226476336044e-05, 0, 8.60754921120146e-06, 0, -1.22964988731449e-06};
    std::vector<double> expected_a{1, -5.78651995575604, 19.9640294916800, -47.7010481156334, 87.8396722590763, -128.653305217857, 154.475693506059, -153.050400350683, 126.160440824548, -85.8045460115290, 47.8358196887074, -21.2060646062098, 7.24441026706683, -1.71334929459713, 0.241873198978882};
    std::vector<double> results{0, 1.22964988731449e-06, 9.57469338616748e-06, 2.59367881852580e-05, 5.29337085399868e-06, -0.000132313491922598, -0.000292346101087556, -3.40707974209247e-05};
    IIRFilter<7, FreqProperty::BandPass> flt(0.3, 0.4, false);
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