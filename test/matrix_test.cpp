#include "gtest/gtest.h"
#include "liegroup.hpp"
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUG__) && (defined(__i386__) || defined(__x86_64__))
#include <cpuid.h>
#include <x86intrin.h>
#endif

using namespace ppx;

#ifdef PPX_USE_AVX
// check cpu info
struct cpuidregs
{
    uint32_t EAX;
    uint32_t EBX;
    uint32_t ECX;
    uint32_t EDX;
};

static inline uint32_t getcpuid(uint32_t r_eax, uint32_t r_ecx, cpuidregs *r_out)
{
#ifdef _MSC_VER
    int cpuid_info[4];
    cpuid_info[0] = cpuid_info[1] = cpuid_info[2] = cpuid_info[3] = 0;
    __cpuidex(cpuid_info, r_eax, r_ecx);
#else
    uint32_t cpuid_info[4];
    cpuid_info[0] = cpuid_info[1] = cpuid_info[2] = cpuid_info[3] = 0;
    __get_cpuid_count(r_eax, r_ecx, &cpuid_info[0], &cpuid_info[1],
                      &cpuid_info[2], &cpuid_info[3]);
#endif
    r_out->EAX = cpuid_info[0];
    r_out->EBX = cpuid_info[1];
    r_out->ECX = cpuid_info[2];
    r_out->EDX = cpuid_info[3];

    uint32_t rc = cpuid_info[0] | cpuid_info[1] | cpuid_info[2] | cpuid_info[3];
    return rc;
}

static inline void xgetbv(uint32_t r_ecx, uint32_t *r_eax, uint32_t *r_edx)
{
    uint64_t x = _xgetbv(r_ecx);

    *r_eax = (uint32_t)(x & 0xFFFFFFFF);
    *r_edx = (uint32_t)((x & 0xFFFFFFFF00000000) >> 32);
}

class CpuidInfo
{
public:
    class CacheInfo
    {
    public:
        enum class Type
        {
            Unknown,
            Data,
            Instruction,
            Unified
        };

    private:
        uint32_t m_Level = 0;
        Type m_Type = Type::Unknown;
        uint32_t m_Size = 0;

    public:
        uint32_t GetLevel() const { return m_Level; }
        uint32_t GetSize() const { return m_Size; }
        Type GetType() const { return m_Type; }

        // These are defined in CacheInfo.cpp
        CacheInfo(uint32_t level, uint32_t type, uint32_t size)
        {
            m_Level = level;
            m_Size = size;

            switch (type)
            {
            case 1:
                m_Type = Type::Data;
                break;

            case 2:
                m_Type = Type::Instruction;
                break;

            case 3:
                m_Type = Type::Unified;
                break;

            default:
                m_Type = Type::Unknown;
                break;
            }
        }

        std::string GetTypeString() const
        {
            const char *s = "";
            switch (m_Type)
            {
            case Type::Data:
                s = "Data";
                break;

            case Type::Instruction:
                s = "Instruction";
                break;

            case Type::Unified:
                s = "Unified";
                break;

            default:
                s = "Unknown";
                break;
            }
            return s;
        }
    };

private:
    bool m_Loaded;
    uint32_t m_MaxEax;                             // Max EAX for basic CPUID
    uint32_t m_MaxEaxExt;                          // Max EAX for extended CPUID
    uint64_t m_FeatureFlags;                       // Processor feature flags
    std::vector<CpuidInfo::CacheInfo> m_CacheInfo; // Processor cache information
    char m_VendorId[13]{};                         // Processor vendor ID string
    char m_ProcessorBrand[49]{};                   // Processor brand string
    bool m_OsXsave;                                // XSAVE is enabled for app use
    bool m_OsAvxState;                             // AVX state is enabled by OS
    bool m_OsAvx512State;                          // AVX-512 state is enabled by OS

    void InitProcessorBrand()
    {
        if (m_MaxEaxExt >= 0x80000004)
        {
            cpuidregs r2{}, r3{}, r4{};
            auto *p = (uint32_t *)m_ProcessorBrand;

            getcpuid(0x80000002, 0, &r2);
            getcpuid(0x80000003, 0, &r3);
            getcpuid(0x80000004, 0, &r4);

            p[0] = r2.EAX;
            p[1] = r2.EBX;
            p[2] = r2.ECX;
            p[3] = r2.EDX;
            p[4] = r3.EAX;
            p[5] = r3.EBX;
            p[6] = r3.ECX;
            p[7] = r3.EDX;
            p[8] = r4.EAX;
            p[9] = r4.EBX;
            p[10] = r4.ECX;
            p[11] = r4.EDX;
        }
        else
        {
            strcpy(m_ProcessorBrand, "Unknown");
        }
    }

    void LoadInfo0()
    {
        // Get MaxEax and VendorID
        cpuidregs r{};
        getcpuid(0, 0, &r);
        m_MaxEax = r.EAX;

        auto *p = (uint32_t *)m_VendorId;
        p[0] = r.EBX;
        p[1] = r.EDX;
        p[2] = r.ECX;

        // Get MaxEaxExt
        getcpuid(0x80000000, 0, &r);
        m_MaxEaxExt = r.EAX;

        // Initialize processor brand string
        InitProcessorBrand();
    }

    void LoadInfo1()
    {
        cpuidregs r{};
        if (m_MaxEax < 1)
        {
            return;
        }
        getcpuid(1, 0, &r);

        //
        // Decode r.ECX flags
        //

        // CPUID.(EAX=01H, ECX=00H):ECX.SSE3[bit 0]
        if (r.ECX & (0x1 << 0))
        {
            m_FeatureFlags |= (uint64_t)FF::SSE3;
        }

        // CPUID.(EAX=01H, ECX=00H):ECX.PCLMULQDQ[bit 1]
        if (r.ECX & (0x1 << 1))
        {
            m_FeatureFlags |= (uint64_t)FF::PCLMULQDQ;
        }

        // CPUID.(EAX=01H, ECX=00H):ECX.SSSE3[bit 9]
        if (r.ECX & (0x1 << 9))
        {
            m_FeatureFlags |= (uint64_t)FF::SSSE3;
        }

        // CPUID.(EAX=01H, ECX=00H):ECX.SSE4.1[bit 19]
        if (r.ECX & (0x1 << 19))
        {
            m_FeatureFlags |= (uint64_t)FF::SSE4_1;
        }

        // CPUID.(EAX=01H, ECX=00H):ECX.SSE4.2[bit 20]
        if (r.ECX & (0x1 << 20))
        {
            m_FeatureFlags |= (uint64_t)FF::SSE4_2;
        }

        // CPUID.(EAX=01H, ECX=00H):ECX.MOVBE[bit 22]
        if (r.ECX & (0x1 << 22))
        {
            m_FeatureFlags |= (uint64_t)FF::MOVBE;
        }

        // CPUID.(EAX=01H, ECX=00H):ECX.POPCNT[bit 23]
        if (r.ECX & (0x1 << 23))
        {
            m_FeatureFlags |= (uint64_t)FF::POPCNT;
        }

        // CPUID.(EAX=01H, ECX=00H):ECX.AESNI[bit 25]
        if (r.ECX & (0x1 << 25))
        {
            m_FeatureFlags |= (uint64_t)FF::AESNI;
        }

        // CPUID.(EAX=01H, ECX=00H):ECX.RDRAND[bit 30]
        if (r.ECX & (0x1 << 30))
        {
            m_FeatureFlags |= (uint64_t)FF::RDRAND;
        }

        //
        // Decode r.RDX flags
        //

        // CPUID.(EAX=01H, ECX=00H):EDX.MMX[bit 23]
        if (r.EDX & (0x1 << 23))
        {
            m_FeatureFlags |= (uint64_t)FF::MMX;
        }

        // CPUID.(EAX=01H, ECX=00H):EDX.FXSR[bit 24]
        if (r.EDX & (0x1 << 24))
        {
            m_FeatureFlags |= (uint64_t)FF::FXSR;
        }

        // CPUID.(EAX=01H, ECX=00H):EDX.SSE[bit 25]
        if (r.EDX & (0x1 << 25))
        {
            m_FeatureFlags |= (uint64_t)FF::SSE;
        }

        // CPUID.(EAX=01H, ECX=00H):EDX.SSE2[bit 26]
        if (r.EDX & (0x1 << 26))
        {
            m_FeatureFlags |= (uint64_t)FF::SSE2;
        }
    }

    void LoadInfo2()
    {
        cpuidregs r{};
        if (m_MaxEax < 7)
        {
            return;
        }
        getcpuid(7, 0, &r);

        //
        // Decode EBX flags
        //

        // CPUID.(EAX=07H, ECX=00H):EBX.BMI1[bit 3]
        if (r.EBX & (0x1 << 3))
        {
            m_FeatureFlags |= (uint64_t)FF::BMI1;
        }

        // CPUID.(EAX=07H, ECX=00H):EBX.BMI2[bit 8]
        if (r.EBX & (0x1 << 8))
        {
            m_FeatureFlags |= (uint64_t)FF::BMI2;
        }

        // CPUID.(EAX=07H, ECX=00H):EBX.ERMSB[bit 9]
        // ERMSB = Enhanced REP MOVSB/STOSB
        if (r.EBX & (0x1 << 9))
        {
            m_FeatureFlags |= (uint64_t)FF::ERMSB;
        }

        // CPUID.(EAX=07H, ECX=00H):EBX.RDSEED[bit 18]
        if (r.EBX & (0x1 << 18))
        {
            m_FeatureFlags |= (uint64_t)FF::RDSEED;
        }

        // CPUID.(EAX=07H, ECX=00H):EBX.ADX[bit 19]
        if (r.EBX & (0x1 << 19))
        {
            m_FeatureFlags |= (uint64_t)FF::ADX;
        }

        // CPUID.(EAX=07H, ECX=00H):EBX.CLWB[bit 24]
        if (r.EBX & (0x1 << 24))
        {
            m_FeatureFlags |= (uint64_t)FF::CLWB;
        }

        //
        // Decode ECX flags
        //

        // CPUID.(EAX=07H, ECX=00H):ECX.PREFETCHWT1[bit 0]
        if (r.ECX & (0x1 << 0))
        {
            m_FeatureFlags |= (uint64_t)FF::PREFETCHWT1;
        }

        // CPUID.(EAX=07H, ECX=00H):ECX.GFNI[bit 8]
        if (r.ECX & (0x1 << 8))
        {
            m_FeatureFlags |= (uint64_t)FF::GFNI;
        }
    }

    void LoadInfo3()
    {
        cpuidregs r{};
        if (m_MaxEaxExt < 0x80000001)
        {
            return;
        }
        getcpuid(0x80000001, 0, &r);

        // CPUID.(EAX=80000001H, ECX=00H):ECX.LZCNT[bit 5]
        if (r.ECX & (0x1 << 5))
        {
            m_FeatureFlags |= (uint64_t)FF::LZCNT;
        }

        // CPUID.(EAX=80000001H, ECX=00H):ECX.PREFETCHW[bit 8]
        if (r.ECX & (0x1 << 8))
        {
            m_FeatureFlags |= (uint64_t)FF::PREFETCHW;
        }
    }

    void LoadInfo4()
    {
        cpuidregs r_eax01h{};
        cpuidregs r_eax07h{};
        cpuidregs r_eax07h_ecx01h{};
        if (m_MaxEax < 7)
        {
            return;
        }
        getcpuid(1, 0, &r_eax01h);
        getcpuid(7, 0, &r_eax07h);
        getcpuid(7, 1, &r_eax07h_ecx01h);

        // Test CPUID.(EAX=01H, ECX=00H):ECX.OSXSAVE[bit 27] to verify use of XGETBV
        m_OsXsave = (r_eax01h.ECX & (0x1 << 27)) != 0;
        if (m_OsXsave)
        {
            // Use XGETBV to obtain following information
            // AVX state is enabled by OS if (XCR0[2:1] == '11b') is true
            // AVX-512 state is enabled by OS if (XCR0[7:5] == '111b') is true

            uint32_t xgetbv_eax, xgetbv_edx;

            xgetbv(0, &xgetbv_eax, &xgetbv_edx);
            m_OsAvxState = (((xgetbv_eax >> 1) & 0x03) == 0x03);

            if (m_OsAvxState)
            {
                // CPUID.(EAX=01H, ECX=00H):ECX.AVX[bit 28]
                if (r_eax01h.ECX & (0x1 << 28))
                {
                    m_FeatureFlags |= (uint64_t)FF::AVX;

                    //
                    // Decode ECX flags
                    //

                    // CPUID.(EAX=07H, ECX=00H):EBX.AVX2[bit 5]
                    if (r_eax07h.EBX & (0x1 << 5))
                    {
                        m_FeatureFlags |= (uint64_t)FF::AVX2;
                    }

                    // CPUID.(EAX=07H, ECX=00H):ECX.VAES[bit 9]
                    if (r_eax07h.ECX & (0x1 << 9))
                    {
                        m_FeatureFlags |= (uint64_t)FF::VAES;
                    }

                    // CPUID.(EAX=07H, ECX=00H):ECX.VPCLMULQDQ[bit 10]
                    if (r_eax07h.ECX & (0x1 << 10))
                    {
                        m_FeatureFlags |= (uint64_t)FF::VPCLMULQDQ;
                    }

                    // CPUID.(EAX=01H, ECX=00H):ECX.FMA[bit 12]
                    if (r_eax01h.ECX & (0x1 << 12))
                    {
                        m_FeatureFlags |= (uint64_t)FF::FMA;
                    }

                    // CPUID.(EAX=01H, ECX=00H):ECX.F16C[bit 29]
                    if (r_eax01h.ECX & (0x1 << 29))
                    {
                        m_FeatureFlags |= (uint64_t)FF::F16C;
                    }

                    //
                    // Decode EAX flags (subleaf 1)
                    //

                    // CPUID.(EAX=07H, ECX=01H):EAX.AVX_VNNI[bit 4]
                    if (r_eax07h_ecx01h.EAX & (0x1 << 4))
                    {
                        m_FeatureFlags |= (uint64_t)FF::AVX_VNNI;
                    }

                    m_OsAvx512State = (((xgetbv_eax >> 5) & 0x07) == 0x07);

                    if (m_OsAvx512State)
                    {
                        // CPUID.(EAX=07H, ECX=00H):EBX.AVX512F[bit 16]
                        if (r_eax07h.EBX & (0x1 << 16))
                        {
                            m_FeatureFlags |= (uint64_t)FF::AVX512F;

                            //
                            // Decode EBX flags
                            //

                            // CPUID.(EAX=07H, ECX=00H):EBX.AVX512DQ[bit 17]
                            if (r_eax07h.EBX & (0x1 << 17))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512DQ;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EBX.AVX512_IFMA[bit 21]
                            if (r_eax07h.EBX & (0x1 << 21))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_IFMA;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EBX.AVX512PF[bit 26]
                            if (r_eax07h.EBX & (0x1 << 26))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512PF;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EBX.AVX512ER[bit 27]
                            if (r_eax07h.EBX & (0x1 << 27))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512ER;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EBX.AVX512CD[bit 28]
                            if (r_eax07h.EBX & (0x1 << 28))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512CD;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EBX.AVX512BW[bit 30]
                            if (r_eax07h.EBX & (0x1 << 30))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512BW;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EBX.AVX512VL[bit 31]
                            if (r_eax07h.EBX & (0x1 << 31))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512VL;
                            }

                            //
                            // Decode ECX flags
                            //

                            // CPUID.(EAX=07H, ECX=00H):ECX.AVX512_VBMI[bit 1]
                            if (r_eax07h.ECX & (0x1 << 1))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_VBMI;
                            }

                            // CPUID.(EAX=07H, ECX=00H):ECX.AVX512_VBMI2[bit 6]
                            if (r_eax07h.ECX & (0x1 << 6))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_VBMI2;
                            }

                            // CPUID.(EAX=07H, ECX=00H):ECX.AVX512_VNNI[bit 11]
                            if (r_eax07h.ECX & (0x1 << 11))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_VNNI;
                            }

                            // CPUID.(EAX=07H, ECX=00H):ECX.AVX512_BITALG[bit 12]
                            if (r_eax07h.ECX & (0x1 << 12))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_BITALG;
                            }

                            // CPUID.(EAX=07H, ECX=00H):ECX.AVX512_VPOPCNTDQ[bit 14]
                            if (r_eax07h.ECX & (0x1 << 14))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_VPOPCNTDQ;
                            }

                            //
                            // Decode EDX flags
                            //

                            // CPUID.(EAX=07H, ECX=00H):EDX.AVX512_4FMAPS[bit 2]
                            if (r_eax07h.EDX & (0x1 << 2))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_4FMAPS;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EDX.AVX512_4VNNIW[bit 3]
                            if (r_eax07h.EDX & (0x1 << 3))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_4VNNIW;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EDX.AVX512_VP2INTERSECT[bit 8]
                            if (r_eax07h.EDX & (0x1 << 8))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_VP2INTERSECT;
                            }

                            // CPUID.(EAX=07H, ECX=00H):EDX.AVX512_FP16[bit 23]
                            if (r_eax07h.EDX & (0x1 << 23))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_FP16;
                            }

                            //
                            // Decode EAX flags (subleaf 1)
                            //

                            // CPUID.(EAX=07H, ECX=01H):EAX.AVX512_BF16[bit 5]
                            if (r_eax07h_ecx01h.EAX & (0x1 << 5))
                            {
                                m_FeatureFlags |= (uint64_t)FF::AVX512_BF16;
                            }
                        }
                    }
                }
            }
        }
    }

    void LoadInfo5()
    {
        if (m_MaxEax < 4)
        {
            return;
        }

        bool done = false;
        uint32_t index = 0;
        while (!done)
        {
            cpuidregs r{};
            getcpuid(4, index, &r);
            uint32_t cache_type = r.EAX & 0x1f;
            uint32_t cache_level = ((r.EAX >> 5) & 0x3);

            if (cache_type == 0)
            {
                done = true;
            }
            else
            {
                uint32_t ways = ((r.EBX >> 22) & 0x3ff) + 1;
                uint32_t partitions = ((r.EBX >> 12) & 0x3ff) + 1;
                uint32_t line_size = (r.EBX & 0xfff) + 1;
                uint32_t sets = r.ECX + 1;
                uint32_t cache_size = ways * partitions * line_size * sets;

                CacheInfo ci{cache_level, cache_type, cache_size};
                m_CacheInfo.push_back(ci);
                index++;
            }
        }
    }

public:
    enum class FF : uint64_t
    {
        FXSR = (uint64_t)1 << 0,
        MMX = (uint64_t)1 << 1,
        MOVBE = (uint64_t)1 << 2,
        SSE = (uint64_t)1 << 3,
        SSE2 = (uint64_t)1 << 4,
        SSE3 = (uint64_t)1 << 5,
        SSSE3 = (uint64_t)1 << 6,
        SSE4_1 = (uint64_t)1 << 7,
        SSE4_2 = (uint64_t)1 << 8,
        PCLMULQDQ = (uint64_t)1 << 9,
        POPCNT = (uint64_t)1 << 10,
        PREFETCHW = (uint64_t)1 << 11,
        PREFETCHWT1 = (uint64_t)1 << 12,
        RDRAND = (uint64_t)1 << 13,
        RDSEED = (uint64_t)1 << 14,
        ERMSB = (uint64_t)1 << 15,
        AVX = (uint64_t)1 << 16,
        AVX2 = (uint64_t)1 << 17,
        F16C = (uint64_t)1 << 18,
        FMA = (uint64_t)1 << 19,
        BMI1 = (uint64_t)1 << 20,
        BMI2 = (uint64_t)1 << 21,
        LZCNT = (uint64_t)1 << 22,
        ADX = (uint64_t)1 << 23,
        AVX512F = (uint64_t)1 << 24,
        AVX512ER = (uint64_t)1 << 25,
        AVX512PF = (uint64_t)1 << 26,
        AVX512DQ = (uint64_t)1 << 27,
        AVX512CD = (uint64_t)1 << 28,
        AVX512BW = (uint64_t)1 << 29,
        AVX512VL = (uint64_t)1 << 30,
        AVX512_IFMA = (uint64_t)1 << 31,
        AVX512_VBMI = (uint64_t)1 << 32,
        AVX512_4FMAPS = (uint64_t)1 << 33,
        AVX512_4VNNIW = (uint64_t)1 << 34,
        AVX512_VPOPCNTDQ = (uint64_t)1 << 35,
        AVX512_VNNI = (uint64_t)1 << 36,
        AVX512_VBMI2 = (uint64_t)1 << 37,
        AVX512_BITALG = (uint64_t)1 << 38,
        AVX512_BF16 = (uint64_t)1 << 39,
        AVX512_VP2INTERSECT = (uint64_t)1 << 40,
        CLWB = (uint64_t)1 << 41,
        GFNI = (uint64_t)1 << 42,
        AESNI = (uint64_t)1 << 43,
        VAES = (uint64_t)1 << 44,
        VPCLMULQDQ = (uint64_t)1 << 45,
        AVX_VNNI = (uint64_t)1 << 46,
        AVX512_FP16 = (uint64_t)1 << 47,
    };

    CpuidInfo() : m_Loaded(false), m_MaxEax(0), m_MaxEaxExt(0), m_FeatureFlags(0),
                  m_OsXsave(false), m_OsAvxState(false), m_OsAvx512State(false)
    {
        m_VendorId[0] = '\0';
        m_VendorId[sizeof(m_VendorId) - 1] = '\0';
        m_ProcessorBrand[0] = '\0';
        m_ProcessorBrand[sizeof(m_ProcessorBrand) - 1] = '\0';
    }

    const std::vector<CpuidInfo::CacheInfo> &GetCacheInfo() const
    {
        return m_CacheInfo;
    }

    bool GetFF(FF flag) const
    {
        return (m_FeatureFlags & (uint64_t)flag) != 0;
    }

    std::string Brand() const { return m_ProcessorBrand; }

    std::string Vendor() const { return m_VendorId; }

    void LoadInfo()
    {
        if (m_Loaded)
        {
            return;
        }
        LoadInfo0();
        LoadInfo1();
        LoadInfo2();
        LoadInfo3();
        LoadInfo4();
        LoadInfo5();
        m_Loaded = true;
    }

    friend std::ostream &operator<<(std::ostream &os, const CpuidInfo &self)
    {
        const char nl = '\n';
        std::cout << "\n----- Processor Info  -----" << nl
                  << "Processor vendor: " << self.Vendor() << nl
                  << "Processor brand:  " << self.Brand() << nl;

        std::cout << "\n----- Cache Info  -----" << nl;
        for (const auto &x : self.GetCacheInfo())
        {
            uint32_t cache_size = x.GetSize();
            uint32_t cache_size_kb = cache_size / 1024;

            std::cout << "Cache L" << x.GetLevel() << ": "
                      << cache_size_kb << " KB - "
                      << x.GetTypeString() << nl;
        }
        std::cout << "\n----- Processor CPUID Feature Flags -----" << nl
                  << "FMA:                 " << std::boolalpha << self.GetFF(CpuidInfo::FF::FMA) << nl
                  << "AVX:                 " << self.GetFF(CpuidInfo::FF::AVX) << nl
                  << "AVX2:                " << self.GetFF(CpuidInfo::FF::AVX2) << nl
                  << "AVX512F:             " << self.GetFF(CpuidInfo::FF::AVX512F) << nl
                  << "AVX512CD:            " << self.GetFF(CpuidInfo::FF::AVX512CD) << nl
                  << "AVX512DQ:            " << self.GetFF(CpuidInfo::FF::AVX512DQ) << nl
                  << "AVX512BW:            " << self.GetFF(CpuidInfo::FF::AVX512BW) << nl
                  << "AVX512VL:            " << self.GetFF(CpuidInfo::FF::AVX512VL) << nl
                  << "AVX512_IFMA:         " << self.GetFF(CpuidInfo::FF::AVX512_IFMA) << nl
                  << "AVX512_VBMI:         " << self.GetFF(CpuidInfo::FF::AVX512_VBMI) << nl
                  << "AVX512_VNNI:         " << self.GetFF(CpuidInfo::FF::AVX512_VNNI) << nl
                  << "AVX512_VPOPCNTDQ:    " << self.GetFF(CpuidInfo::FF::AVX512_VPOPCNTDQ) << nl
                  << "AVX512_VBMI2:        " << self.GetFF(CpuidInfo::FF::AVX512_VBMI2) << nl
                  << "AVX512_BITALG:       " << self.GetFF(CpuidInfo::FF::AVX512_BITALG) << nl
                  << "AVX512_BF16:         " << self.GetFF(CpuidInfo::FF::AVX512_BF16) << nl
                  << "AVX512_VP2INTERSECT: " << self.GetFF(CpuidInfo::FF::AVX512_VP2INTERSECT) << nl
                  << "AVX512_FP16:         " << self.GetFF(CpuidInfo::FF::AVX512_FP16) << nl;
        return os;
    }
};

void list_cpu_infos()
{
    static CpuidInfo cpu_infos;
    cpu_infos.LoadInfo();
    std::cout << cpu_infos << "\n";
}

#endif

class MatrixS_TestCase : public ::testing::Test
{
public:
    MatrixS_TestCase() = default;
};

TEST_F(MatrixS_TestCase, ctor)
{
    MatrixS<1, 1> A = {};
    EXPECT_EQ(A[0], 0);
    MatrixS<2, 2> B{1, 2, 3, 4};
    EXPECT_EQ(B(1, 1), 4);
    MatrixS<2, 2> C(std::array<int, 4>{1, 2, 3, 4});
    MatrixS<2, 2> D(std::vector<int>{1, 2, 3, 4});
    EXPECT_EQ(B, C);
    EXPECT_EQ(C, D);
    EXPECT_EQ(sizeof(MatrixS<2, 2>), (sizeof(double)) * 4);
    EXPECT_EQ(sizeof(MatrixS<20, 20>), sizeof(std::vector<double>));
    MatrixS<3, 2> E{1, 2, 3, 4, 5, 6};
    MatrixS<3, 2> F{{1, 2, 3},
                    {4, 5, 6}};
    MatrixS<3, 2> G{{{1, 4}, {2, 5}, {3, 6}}, Ori::Row};
    EXPECT_EQ(E, F);
    EXPECT_EQ(F, G);
    MatrixS<3, 2> H{{1, 2, 3, 4},
                    {4, 5, 6}};
    EXPECT_EQ(G, H);
}

TEST_F(MatrixS_TestCase, expr)
{
    MatrixS<4, 4> a = {1, 2, 3, 4, 5, 6, 7, 8};
    a.sub<4, 2, false>(0, 2) = {1, 2, 3, 4, 5, 6, 7, 8};
    MatrixS<4, 4> result{1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8};
    EXPECT_EQ(a, result);
    a.sub<4, 1, false>(0, 0) = a.sub<4, 1>(0, 0) - a.sub<4, 1>(0, 2);
    a.sub<4, 1, false>(0, 1) = a.sub<4, 1>(0, 1) - a.sub<4, 1>(0, 3);
    a.sub<4, 2>(0, 2) = a.sub<4, 2>(0, 2) - a.sub<4, 2>(0, 2);
    a.sub<4, 2>(0, 2) = a.sub<4, 2>(0, 2) - a.sub<4, 2>(0, 2);
    EXPECT_EQ(a, decltype(a){});
    MatrixS<3, 3> x = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    MatrixS<3, 3> y = x.T();
    MatrixS<3, 3> z;
    z.sub<3, 3>(0, 0) =
        (x * 2 + y * 2 + x * y) / 2 - MatrixS<3, 3>{35.0, 45.0, 55.0, 45.0, 56.5, 68.0, 55.0, 68.0, 81.0};
    EXPECT_EQ(z, decltype(z){});
    MatrixS<3, 3> w = Sqrt(pwmul(x, x));
    EXPECT_EQ(w, x);
    MatrixS<3, 3> v = Abs(pwdiv(pwmul(x, x.sub<3, 3>(0, 0)), 1 + x.sub<3, 3>(0, 0) - 1));
    EXPECT_EQ(w, v);
}

TEST_F(MatrixS_TestCase, func)
{
    MatrixS<2, 2> x = {1, 2, 4, 3};
    MatrixS<2, 2> result{-0.6, 0.4, 0.8, -0.2};
    EXPECT_NEAR(x.det(), -5.0, ppx::EPS_DP);
    EXPECT_EQ(x.I(), result);
    MatrixS<1, 1> c = {3};
    EXPECT_NEAR(c.I()[0], 1.0 / 3.0, ppx::EPS_DP);
}

TEST_F(MatrixS_TestCase, fac)
{
    MatrixS<4, 3> X{3.5, 1.6, 3.7, 4.3,
                    2.7, -5.7, -0.8, -9.8,
                    -3.1, -6.0, 1.9, 6.9};
    MatrixS<4, 1> Y{1, 1, 1, 1};
    MatrixS<3, 1> result{0.271349846985587, -0.030388140654028, -0.062228084118990};
    auto x1 = linsolve<Factorization::QR>(X, Y);
    auto x2 = linsolve<Factorization::SVD>(X, Y);
    EXPECT_EQ(x1.x, result);
    EXPECT_EQ(x2.x, result);

    MatrixS<4, 3> u{1, 4, 7, 11,
                    2, 5, 8, 1,
                    3, 6, 9, 5};
    SVD<4, 3> ru(u);
    EXPECT_EQ(ru.u * ru.w.diag() * ru.v.T(), u);

    MatrixS<4, 4> g{2, -1, 1, 4,
                    -1, 2, -1, 1,
                    1, -1, 2, -2,
                    4, 1, -2, 3};
    EigenValue<4> e1(g, true);
    MatrixS<4, 1> result2{-2.674416988218162, 1.0, 3.945104914279287, 6.729312073938871};
    EXPECT_EQ(e1.d, result2);
    MatrixS<4, 4> result3{-0.657702449551592, 0.127000127000191, 0.452113001308500, 0.588975627376541,
                          -0.196417934159944, 0.762000762001143, -0.611122207034796, 0.0854662618752970,
                          0.367242010902251, 0.635000635000952, 0.642924048180056, -0.220354639725665,
                          0.627678889578617, 0.0, -0.0936597586400323, 0.772817611852141};
    result3 = result3.T() * (-1);
    EXPECT_EQ(e1.z, result3);
}

TEST_F(MatrixS_TestCase, cat)
{
    MatrixS<2, 2> x = {1, 2, 3, 4};
    MatrixS<2, 2> y = {5, 6, 7, 8};
    MatrixS<2, 4> expect1{1, 2, 3, 4, 5, 6, 7, 8};
    MatrixS<4, 2> expect2{1, 2, 5, 6, 3, 4, 7, 8};
    EXPECT_EQ(expect1, concat<Ori::Col>(x, y));
    EXPECT_EQ(expect2, concat<Ori::Row>(x, y));

    MatrixS<2, 3> a{1, 2, 3, 4, 5, 6};
    MatrixS<2, 1> b{9, 9};
    MatrixS<1, 3> c{-1, -1, -1};
    MatrixS<2, 4> expect3{1, 2, 3, 4, 5, 6, 9, 9};
    MatrixS<3, 3> expect4{1, 2, -1, 3, 4, -1, 5, 6, -1};
    EXPECT_EQ(expect3, concat(a, b));
    EXPECT_EQ(expect4, concat(a, c));
    MatrixS<3, 1> expect5{9, 9, 1};
    MatrixS<3, 1> expect6{1, 9, 9};
    MatrixS<1, 4> expect7{-1, -1, -1, 3};
    EXPECT_EQ(expect5, concat(b, 1));
    EXPECT_EQ(expect6, concat(1, b));
    EXPECT_EQ(expect7, concat(c, 3));
}

#ifdef PPX_USE_AVX

template <size_t M, size_t N, size_t L>
MatrixS<M, L> matmul(const MatrixS<M, N> &self, const MatrixS<N, L> &other)
{
    MatrixS<M, L> result;
    for (size_t k = 0; k < N; k++)
    {
        for (size_t j = 0; j < L; j++)
        {
            for (size_t i = 0; i < M; i++)
            {
                result(i, j) += self(i, k) * other(k, j);
            }
        }
    }
    return result;
}

TEST_F(MatrixS_TestCase, cpuinfo)
{
    list_cpu_infos();
}

TEST_F(MatrixS_TestCase, avx_aligned)
{
    // manually aligned
    alignas(32) double b[4]{};
    EXPECT_TRUE(avxt::chkalign(b, 32));
    // on heap. not aligned
    // std::vector<double> c{1, 2, 3, 4};
    // EXPECT_FALSE(avxt::chkalign(c.data(), 32));
}

TEST_F(MatrixS_TestCase, avx_sum)
{
    MatrixS<20, 20> A;
    for (size_t i = 0; i < A.size(); i++)
    {
        random(A, -1e3, 1e3);
        EXPECT_NEAR(sum(A.data(), i),
                    std::accumulate(A.cbegin(), A.cbegin() + i, 0.0),
                    EPS_DP * 1e3 * A.size());
    }
}

TEST_F(MatrixS_TestCase, avx_norm2)
{
    MatrixS<101, 1> A;
    for (size_t i = 0; i < A.size(); i++)
    {
        random(A, -1e3, 1e3);
        EXPECT_NEAR(norm2(A),
                    sqrt(std::inner_product(A.cbegin(), A.cend(), A.cbegin(), 0.0)),
                    EPS_DP * 1e3 * A.size());
    }
}

TEST_F(MatrixS_TestCase, avx_full_mul)
{
    MatrixS<4, 4> A, B;
    random(A, -1e3, 1e3);
    random(B, -1e3, 1e3);
    EXPECT_EQ(A * B, matmul(A, B));

    MatrixS<4, 3> C;
    MatrixS<3, 4> D;
    random(C, -1e3, 1e3);
    random(D, -1e3, 1e3);
    EXPECT_EQ(C * D, matmul(C, D));

    MatrixS<4, 2> E;
    MatrixS<2, 4> F;
    random(E, -1e3, 1e3);
    random(F, -1e3, 1e3);
    EXPECT_EQ(E * F, matmul(E, F));

    MatrixS<4, 1> G;
    MatrixS<1, 4> H;
    random(G, -1e3, 1e3);
    random(H, -1e3, 1e3);
    EXPECT_EQ(G * H, matmul(G, H));

    MatrixS<16, 7> I;
    MatrixS<7, 16> J;
    random(I, -1e3, 1e3);
    random(J, -1e3, 1e3);
    EXPECT_EQ(I * J, matmul(I, J));

    MatrixS<80, 33> K;
    MatrixS<33, 20> L;
    random(K, -1e3, 1e3);
    random(L, -1e3, 1e3);
    EXPECT_EQ(K * L, matmul(K, L));

    MatrixS<100, 104> M;
    MatrixS<104, 100> N;
    random(M, -1e3, 1e3);
    random(N, -1e3, 1e3);
    EXPECT_EQ(M * N, matmul(M, N));
}

TEST_F(MatrixS_TestCase, avx_mask_mul)
{
    MatrixS<3, 3> A, B;
    random(A, -1e3, 1e3);
    random(B, -1e3, 1e3);
    EXPECT_EQ(A * B, matmul(A, B));

    MatrixS<3, 4> C;
    MatrixS<4, 3> D;
    random(C, -1e3, 1e3);
    random(D, -1e3, 1e3);
    EXPECT_EQ(C * D, matmul(C, D));

    MatrixS<3, 2> E;
    MatrixS<2, 3> F;
    random(E, -1e3, 1e3);
    random(F, -1e3, 1e3);
    EXPECT_EQ(E * F, matmul(E, F));

    MatrixS<3, 1> G;
    MatrixS<1, 3> H;
    random(G, -1e3, 1e3);
    random(H, -1e3, 1e3);
    EXPECT_EQ(G * H, matmul(G, H));

    MatrixS<3, 7> I;
    MatrixS<7, 3> J;
    random(I, -1e3, 1e3);
    random(J, -1e3, 1e3);
    EXPECT_EQ(I * J, matmul(I, J));

    MatrixS<3, 33> K;
    MatrixS<33, 3> L;
    random(K, -1e3, 1e3);
    random(L, -1e3, 1e3);
    EXPECT_EQ(K * L, matmul(K, L));

    MatrixS<3, 104> M;
    MatrixS<104, 3> N;
    random(M, -1e3, 1e3);
    random(N, -1e3, 1e3);
    EXPECT_EQ(M * N, matmul(M, N));
}

TEST_F(MatrixS_TestCase, avx_transpose)
{
    MatrixS<4, 4> A({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}}, Ori::Col);
    MatrixS<4, 4> B({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}}, Ori::Row);
    EXPECT_EQ(A.T(), B);
}

#endif

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}