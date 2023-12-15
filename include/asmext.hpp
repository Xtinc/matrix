#ifndef VVERY_SIMPLE_ASMEXT_HEADER
#define VVERY_SIMPLE_ASMEXT_HEADER

#include <immintrin.h>
#include <cinttypes>

namespace ppx
{
#if defined(PPX_USE_AVX)
    namespace avxt
    {
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
            // Perform reduction, final sum in low-order element of temp3
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

        inline void transpose4x4(const double *a, double *b)
        {
            auto q_a = _mm256_loadu_pd(a);
            auto q_b = _mm256_loadu_pd(a + 4);
            auto q_c = _mm256_loadu_pd(a + 8);
            auto q_d = _mm256_loadu_pd(a + 12);

            auto temp0 = _mm256_unpacklo_pd(q_a, q_b);
            auto temp1 = _mm256_unpackhi_pd(q_a, q_b);
            auto temp2 = _mm256_unpacklo_pd(q_c, q_d);
            auto temp3 = _mm256_unpackhi_pd(q_c, q_d);

            _mm256_storeu_pd(b, _mm256_permute2f128_pd(temp0, temp2, 0x20));
            _mm256_storeu_pd(b + 4, _mm256_permute2f128_pd(temp1, temp3, 0x20));
            _mm256_storeu_pd(b + 8, _mm256_permute2f128_pd(temp0, temp2, 0x31));
            _mm256_storeu_pd(b + 12, _mm256_permute2f128_pd(temp1, temp3, 0x31));
        }
    }
#endif
}
#endif