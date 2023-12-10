#ifdef PPX_USE_AVX

#include <immintrin.h>
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUG__)
#include <cpuid.h>
#include <x86intrin.h>
#else
#error Unknown target
#endif
#include <vector>
#include <cstring>
#include <string>
#include <iostream>

namespace ppx
{
    namespace avxt
    {
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

        // avx related
        template <typename T>
        bool chkalign(const T *p, size_t alignment)
        {
            if (p == nullptr)
            {
                return false;
            }
            if (((uintptr_t)p % alignment) != 0)
            {
                return false;
            }
            return true;
        }

        inline double sum(const double *start, size_t n)
        {
            auto q_sum = _mm256_setzero_pd();
            size_t i = 0;
            for (; n - i >= 4; i += 4)
            {
                q_sum = _mm256_add_pd(_mm256_loadu_pd(start + i), q_sum);
            }
            // Peform reduction, final sum in low-order element of temp3
            auto temp = _mm_add_pd(_mm256_extractf128_pd(q_sum, 0),
                                   _mm256_extractf128_pd(q_sum, 1));
            // Process remaining elements
            double sum{};
            _mm_store_sd(&sum, _mm_hadd_pd(temp, temp));

            for (; i < n; i++)
            {
                sum += *(start + i);
            }
            return sum;
        }

        inline double inrpdt(const double *a, const double *b, size_t n)
        {
            auto q_sum = _mm256_setzero_pd();
            size_t i = 0;
            for (; n - i >= 4; i += 4)
            {
                auto q_a = _mm256_loadu_pd(a + i);
                auto q_b = _mm256_loadu_pd(b + i);
                q_sum = _mm256_fmadd_pd(q_a, q_b, q_sum);
            }
            // Peform reduction, final sum in low-order element of temp3
            auto temp = _mm_add_pd(_mm256_extractf128_pd(q_sum, 0),
                                   _mm256_extractf128_pd(q_sum, 1));
            // Process remaining elements
            double sum{};
            _mm_store_sd(&sum, _mm_hadd_pd(temp, temp));

            for (; i < n; i++)
            {
                auto na = *(a + i);
                auto nb = *(b + i);
                sum += na * nb;
            }
            return sum;
        }

        template <size_t M, size_t N, size_t L>
        inline void gemm(const double *a, const double *b, double *c)
        {
            constexpr uint64_t ZR = 0;
            constexpr uint64_t MV = 0x8000000000000000;
            alignas(32) const uint64_t c_Mask0[4]{ZR, ZR, ZR, ZR};
            alignas(32) const uint64_t c_Mask1[4]{MV, ZR, ZR, ZR};
            alignas(32) const uint64_t c_Mask2[4]{MV, MV, ZR, ZR};
            alignas(32) const uint64_t c_Mask3[4]{MV, MV, MV, ZR};

            const uint64_t *c_MaskMovLUT[8]{c_Mask0, c_Mask1, c_Mask2, c_Mask3};
            constexpr size_t left_cols = M % 4;
            auto res_mask = _mm256_load_si256((__m256i *)c_MaskMovLUT[left_cols]);

            for (size_t j = 0; j < L; j++)
            {
                size_t i = 0;
                while (i + 4 <= M)
                {
                    auto quad_c = _mm256_setzero_pd();
                    for (size_t k = 0; k < N; k++)
                    {
                        quad_c = _mm256_fmadd_pd(_mm256_loadu_pd(a + i + k * M),
                                                 _mm256_broadcast_sd(b + k + j * N),
                                                 quad_c);
                    }
                    _mm256_storeu_pd(c + i + j * M, quad_c);
                    i += 4;
                }
                if (left_cols)
                {
                    auto quad_c = _mm256_setzero_pd();
                    for (size_t k = 0; k < N; k++)
                    {
                        quad_c = _mm256_fmadd_pd(_mm256_maskload_pd(a + i + k * M, res_mask),
                                                 _mm256_broadcast_sd(b + k + j * N),
                                                 quad_c);
                    }
                    _mm256_maskstore_pd(c + i + j * M, res_mask, quad_c);
                }
            }
        }
    }

    inline void list_cpu_infos()
    {
        static avxt::CpuidInfo cpu_infos;
        cpu_infos.LoadInfo();
        std::cout << cpu_infos << "\n";
    }
}
#endif