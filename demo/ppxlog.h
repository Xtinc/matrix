#ifndef VVERY_SIMPLE_PPXLOG_HEADER
#define VVERY_SIMPLE_PPXLOG_HEADER

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <sstream>
#include <tuple>
#include <type_traits>

namespace ppx
{
    struct bit_t
    {
        explicit constexpr bit_t(int b) : m_bit_idx(b) {}
        explicit constexpr operator int() const { return m_bit_idx; }

    private:
        int m_bit_idx;
    };

    constexpr bit_t operator"" _bit(unsigned long long int b) { return bit_t{static_cast<int>(b)}; }

    template <typename T, typename Tag, typename Cond = std::enable_if_t<std::is_integral<T>::value>>
    struct bit_flag
    {
        constexpr bit_flag(bit_flag const &rhs) noexcept = default;
        constexpr bit_flag(bit_flag &&rhs) noexcept = default;
        constexpr bit_flag() noexcept : m_val(0) {}
        explicit constexpr bit_flag(T const val) noexcept : m_val(val) {}
        constexpr bit_flag(bit_t const bit) noexcept : m_val(static_cast<T>(T{1} << static_cast<int>(bit))) {}
        explicit constexpr operator T() const noexcept { return m_val; }
        explicit constexpr operator bool() const noexcept { return m_val != 0; }

        static constexpr bit_flag all()
        {
            return bit_flag(~T{0});
        }

        T val() const
        {
            return m_val;
        }

        bool constexpr operator==(bit_flag const f) const noexcept
        {
            return m_val == f.m_val;
        }

        bool constexpr operator!=(bit_flag const f) const noexcept
        {
            return m_val != f.m_val;
        }

        bit_flag &operator|=(bit_flag const f) noexcept
        {
            m_val |= f.m_val;
            return *this;
        }

        bit_flag &operator&=(bit_flag const f) noexcept
        {
            m_val &= f.m_val;
            return *this;
        }

        bit_flag &operator^=(bit_flag const f) noexcept
        {
            m_val ^= f.m_val;
            return *this;
        }

        constexpr friend bit_flag operator|(bit_flag const lhs, bit_flag const rhs) noexcept
        {
            return bit_flag(lhs.m_val | rhs.m_val);
        }

        constexpr friend bit_flag operator&(bit_flag const lhs, bit_flag const rhs) noexcept
        {
            return bit_flag(lhs.m_val & rhs.m_val);
        }

        constexpr friend bit_flag operator^(bit_flag const lhs, bit_flag const rhs) noexcept
        {
            return bit_flag(lhs.m_val ^ rhs.m_val);
        }

        constexpr bit_flag operator~() const noexcept
        {
            return bit_flag(~m_val);
        }

        bit_flag &operator=(bit_flag const &rhs) noexcept = default;
        bit_flag &operator=(bit_flag &&rhs) noexcept = default;

    private:
        T m_val;
    };

    template <typename T, typename Tag>
    std::ostream &operator<<(std::ostream &os, bit_flag<T, Tag> val)
    {
        return os << static_cast<T>(val);
    }
    template <typename Tag>
    std::ostream &operator<<(std::ostream &os, bit_flag<std::uint8_t, Tag> val)
    {
        return os << static_cast<std::uint16_t>(static_cast<std::uint8_t>(val));
    }

    using LogLevel = bit_flag<std::uint8_t, struct LogLevel_tag>;

    constexpr LogLevel CH01 = 0_bit;
    constexpr LogLevel CH02 = 1_bit;
    constexpr LogLevel CH03 = 2_bit;
    constexpr LogLevel CH04 = 3_bit;
    constexpr LogLevel CH05 = 4_bit;
    constexpr LogLevel CH06 = 5_bit;
    constexpr LogLevel CH07 = 6_bit;
    constexpr LogLevel CH08 = 7_bit;

    namespace details
    {
        // log below
        struct string_literal_t
        {
            explicit string_literal_t(const char *s) : m_s(s) {}
            char const *m_s;
        };

        template <typename T, typename Tuple>
        struct TupleIndex;

        template <typename T, typename... Types>
        struct TupleIndex<T, std::tuple<T, Types...>>
        {
            static constexpr const std::size_t value = 0;
        };

        template <typename T, typename U, typename... Types>
        struct TupleIndex<T, std::tuple<U, Types...>>
        {
            static constexpr const std::size_t value = 1 + TupleIndex<T, std::tuple<Types...>>::value;
        };

        using SupportedTypes =
            std::tuple<char, uint32_t, uint64_t, int32_t, int64_t, double, details::string_literal_t, char *>;

        class LogLine
        {
        public:
            LogLine(LogLevel level, char const *file, char const *function, uint32_t line);

            void stringify(std::ostream &os);

            LogLine &operator<<(char arg);
            LogLine &operator<<(int32_t arg);
            LogLine &operator<<(uint32_t arg);
            LogLine &operator<<(int64_t arg);
            LogLine &operator<<(uint64_t arg);
            LogLine &operator<<(double arg);
            LogLine &operator<<(const std::string &arg);
            template <typename T>
            LogLine &operator<<(const std::vector<T> &arg)
            {
                *this << '[';
                for (const auto &i : arg)
                {
                    *this << i << ' ';
                }
                *this << ']';
                return *this;
            }
            template <typename T1, typename T2>
            LogLine &operator<<(const std::map<T1, T2> &arg)
            {
                *this << '[';
                for (const auto &p : arg)
                {
                    *this << p.first << ' ' << p.second;
                }
                *this << ']';
                return *this;
            }
            template <typename T, size_t N>
            LogLine &operator<<(const std::array<T, N> &arg)
            {
                *this << '[';
                for (size_t i = 0; i < N; ++i)
                {
                    *this << arg[i] << ' ';
                }
                *this << ']';
                return *this;
            }
            template <size_t N>
            LogLine &operator<<(const char (&arg)[N])
            {
                encode(string_literal_t(arg));
                return *this;
            }
            template <typename Arg>
            std::enable_if_t<std::is_same<Arg, const char *>::value, LogLine &>
            operator<<(Arg const &arg)
            {
                encode(arg);
                return *this;
            }
            template <typename Arg>
            std::enable_if_t<std::is_same<Arg, char *>::value, LogLine &>
            operator<<(const Arg &arg)
            {
                encode(arg);
                return *this;
            }
            template <typename Arg>
            LogLine &operator<<(const Arg &arg)
            {
                std::stringstream ss;
                ss << arg;
                *this << ss.str();
                return *this;
            }

        private:
            char *buffer();

            template <typename Arg>
            void encode(Arg arg);

            template <typename Arg>
            void encode(Arg arg, uint8_t type_id);

            void encode(char *arg);
            void encode(const char *arg);
            void encode(details::string_literal_t arg);
            void encode_c_string(const char *arg, size_t length);
            void resize_buffer_if_needed(size_t additional_bytes);
            void stringify(std::ostream &os, char *start, char const *end);

        private:
            size_t m_bytes_used;
            size_t m_buffer_size;
            std::unique_ptr<char[]> m_heap_buffer;
            char m_stack_buffer[256 - 2 * sizeof(size_t) - sizeof(decltype(m_heap_buffer)) - 8 /* Reserved */];
        };

        struct PsuedoLog
        {
            bool operator==(LogLine &);
        };
    }

    void set_log_level(LogLevel level);

    bool is_logged(LogLevel level);

    void initialize_log(const std::string &log_directory, const std::string &log_file_name, uint32_t log_file_roll_size_mb);
}

#define PPX_LOG(LEVEL) ppx::details::PsuedoLog() == ppx::details::LogLine(LEVEL, __FILE__, __FUNCTION__, __LINE__)
#define LOG_CH(NUM) ppx::is_logged(ppx::CH##NUM) && PPX_LOG(ppx::CH##NUM)

#endif