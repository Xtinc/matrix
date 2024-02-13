#include "plog.h"
#include "gtest/gtest.h"
#include <math.h>
#include <list>

#define LOGD(X) LOG_CH(111) << X

class PLOG_TestCase : public ::testing::Test
{
public:
    PLOG_TestCase() = default;

    static void SetUpTestSuite()
    {
        ppx::log_options opts;
        opts.FILE.on = true;
        opts.FILE.max_size_mb = 10;
        opts.FILE.max_size_all = 50;
        ppx::initialize_log(opts);
    }
};

TEST_F(PLOG_TestCase, quick_log)
{
    constexpr int argv = 3;
    constexpr char test_data[56] = "1234567890abcdefghigklmnopqrstuvwxyz_zyw_zsummer_log4z";
    LOGD("char:" << 'c'
                 << ", unsigned char:" << (unsigned char)'c'
                 << ", short:" << (short)-1
                 << ", unsigned short:" << (unsigned short)-1
                 << ", int:" << (int)-1
                 << ", unsigned int:" << (unsigned int)-1
                 << ", long:" << (long)-1
                 << ", unsigned long:" << (unsigned long)-1
                 << ", long long:" << (long long)-1
                 << ", unsigned long long:" << (unsigned long long)-1
                 << ", float:" << (float)-1.234567
                 << ", double:" << (double)-2.34566
                 << ", double:" << pow(2, 52) - 1.0
                 << ", double:" << pow(2, 52) * -1000
                 << ", double:" << pow(2, 52) / 1000
                 << ", double:" << pow(2, 52) / -1000
                 << ", double:" << pow(2, -58)
                 << ", double:" << pow(2, -16) * -1
                 << ", std::string:" << std::string("fffff")
                 << ", int *:" << (int *)argv
                 << ", const int *:" << (const int *)argv
                 << ", constant:" << 1000
                 << ", constant:" << 100.12345678
                 << ", bool:" << true
                 << ", show const char* data:" << test_data);
}

TEST_F(PLOG_TestCase, rand_log)
{
    for (size_t j = 0; j < 1000; j++)
    {
        for (size_t i = 0; i < 100; i++)
        {
            int r = rand() % 14;
            switch (r)
            {
            case 0:
                LOG_CH(111) << 'c';
                break;
            case 1:
                LOG_CH(111) << UINT8_MAX;
                break;
            case 2:
                LOG_CH(111) << INT16_MIN;
                break;
            case 3:
                LOG_CH(111) << INT16_MAX;
                break;
            case 4:
                LOG_CH(111) << UINT16_MAX;
                break;
            case 5:
                LOG_CH(111) << INT32_MIN;
                break;
            case 6:
                LOG_CH(111) << INT32_MAX;
                break;
            case 7:
                LOG_CH(111) << UINT32_MAX;
                break;
            case 8:
                LOG_CH(111) << INT64_MIN;
                break;
            case 9:
                LOG_CH(111) << INT64_MAX;
                break;
            case 10:
                LOG_CH(111) << UINT64_MAX;
                break;
            case 11:
                LOG_CH(111) << (float)pow(rand() % 100, rand() % 20) * ((rand() % 2 == 0 ? -1.0 : 1.0));
                break;
            case 12:
                LOG_CH(111) << (double)pow(rand() % 100, rand() % 200) * ((rand() % 2 == 0 ? -1.0 : 1.0));
                break;
            default:
                LOG_CH(111) << "8";
                break;
            }
        }
    }
}