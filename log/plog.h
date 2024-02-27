#ifndef VVERY_SIMPLE_PPXLOG_HEADER
#define VVERY_SIMPLE_PPXLOG_HEADER

#include <type_traits>
#include <utility>
#include <sstream>
#include <memory>

// Operating System Evaluation
#if (defined(_WIN32) || defined(_WIN64))
#define PLOG_OS_WINDOWS 1
#else
#define PLOG_OS_WINDOWS 0
#endif
#if (defined(__linux) || defined(__linux__))
#define PLOG_OS_LINUX 1
#else
#define PLOG_OS_LINUX 0
#endif
#if (defined(__QNX__) || defined(__QNXNTO__))
#define PLOG_OS_QNX 1
#else
#define PLOG_OS_QNX 0
#endif
#if (defined(__APPLE__))
#define PLOG_OS_MAC 1
#else
#define PLOG_OS_MAC 0
#endif
// above are *nix
#if ((PLOG_OS_LINUX || PLOG_OS_MAC || PLOG_OS_QNX) && (!PLOG_OS_WINDOWS))
#define PLOG_OS_UNIX 1
#else
#define PLOG_OS_UNIX 0
#endif
// code are not guaranteed except these OSs. let it makes compiler errors.
#if (!(PLOG_OS_UNIX || PLOG_OS_WINDOWS))
#error unknow Operating system for PLOG!
#endif

// Helper macro for declaring functions as having similar signature to printf.
// This allows the compiler to catch format errors at compile-time.
#if defined(__clang__) || defined(__GNUC__)
#define PLOG_PRINTF_CHECK(fmtarg, firstvararg) __attribute__((__format__(__printf__, fmtarg, firstvararg)))
#define PLOG_PRINT_STRING_TYPE const char *
#elif defined(_MSC_VER)
#define PLOG_PRINTF_CHECK(fmtarg, firstvararg)
#define PLOG_PRINT_STRING_TYPE _In_z_ _Printf_format_string_ const char *
#else
#define PLOG_PRINTF_CHECK(fmtarg, firstvararg)
#define PLOG_PRINT_STRING_TYPE const char *
#endif

namespace ppx
{
    namespace details
    {
        template <typename... Ts>
        struct make_noid
        {
            typedef void type;
        };

        template <typename... Ts>
        using noid_t = typename make_noid<Ts...>::type;

        template <typename T, typename U = noid_t<>>
        struct is_overloaded_stream_impl : public std::false_type
        {
        };

        template <typename T>
        struct is_overloaded_stream_impl<T, noid_t<decltype(std::ostringstream() << std::declval<T>())>> : public std::true_type
        {
        };

        template <typename T>
        static constexpr bool is_overloaded_stream()
        {
            return is_overloaded_stream_impl<T>::value;
        }

        void psuedo_log_flush();
    }

    inline const char *filename(const char *path)
    {
        for (auto ptr = path; *ptr; ++ptr)
        {
            if (*ptr == '/' || *ptr == '\\')
            {
                path = ptr + 1;
            }
        }
        return path;
    }

    inline const char *funcname(const char *func)
    {
        for (auto ptr = func + 1; *ptr; ++ptr)
        {
            if (*ptr == ':' && *(ptr - 1) == ':')
            {
                func = ptr + 1;
            }
        }
        return func;
    }

    struct bit_t
    {
        explicit constexpr bit_t(int b) : m_bit_idx(b) {}
        explicit constexpr operator int() const { return m_bit_idx; }

    private:
        int m_bit_idx;
    };

    constexpr bit_t operator"" _bit(unsigned long long int b) { return bit_t{static_cast<int>(b)}; }

    template <typename T, typename Tag, typename Cond = typename std::enable_if<std::is_integral<T>::value>::type>
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

    using LogLevel = bit_flag<std::uint16_t, struct LogLevel_tag>;

    constexpr LogLevel CH000 = 0_bit;
    constexpr LogLevel CH001 = 1_bit;
    constexpr LogLevel CH002 = 2_bit;
    constexpr LogLevel CH003 = 3_bit;
    constexpr LogLevel CH010 = 4_bit;
    constexpr LogLevel CH011 = CH001 | CH010;
    constexpr LogLevel CH012 = CH002 | CH010;
    constexpr LogLevel CH013 = CH003 | CH010;
    constexpr LogLevel CH020 = 5_bit;
    constexpr LogLevel CH021 = CH001 | CH020;
    constexpr LogLevel CH022 = CH002 | CH020;
    constexpr LogLevel CH023 = CH003 | CH020;
    constexpr LogLevel CH030 = 6_bit;
    constexpr LogLevel CH031 = CH001 | CH030;
    constexpr LogLevel CH032 = CH002 | CH030;
    constexpr LogLevel CH033 = CH003 | CH030;
    constexpr LogLevel CH100 = 7_bit;
    constexpr LogLevel CH101 = CH001 | CH100;
    constexpr LogLevel CH102 = CH002 | CH100;
    constexpr LogLevel CH103 = CH003 | CH100;
    constexpr LogLevel CH110 = CH010 | CH100;
    constexpr LogLevel CH111 = CH011 | CH100;
    constexpr LogLevel CH112 = CH012 | CH100;
    constexpr LogLevel CH113 = CH013 | CH100;
    constexpr LogLevel CH120 = CH020 | CH100;
    constexpr LogLevel CH121 = CH021 | CH100;
    constexpr LogLevel CH122 = CH022 | CH100;
    constexpr LogLevel CH123 = CH023 | CH100;
    constexpr LogLevel CH130 = CH030 | CH100;
    constexpr LogLevel CH131 = CH031 | CH100;
    constexpr LogLevel CH132 = CH032 | CH100;
    constexpr LogLevel CH133 = CH033 | CH100;
    constexpr LogLevel CH200 = 8_bit;
    constexpr LogLevel CH201 = CH001 | CH200;
    constexpr LogLevel CH202 = CH002 | CH200;
    constexpr LogLevel CH203 = CH003 | CH200;
    constexpr LogLevel CH210 = CH010 | CH200;
    constexpr LogLevel CH211 = CH011 | CH200;
    constexpr LogLevel CH212 = CH012 | CH200;
    constexpr LogLevel CH213 = CH013 | CH200;
    constexpr LogLevel CH220 = CH020 | CH200;
    constexpr LogLevel CH221 = CH021 | CH200;
    constexpr LogLevel CH222 = CH022 | CH200;
    constexpr LogLevel CH223 = CH023 | CH200;
    constexpr LogLevel CH230 = CH030 | CH200;
    constexpr LogLevel CH231 = CH031 | CH200;
    constexpr LogLevel CH232 = CH032 | CH200;
    constexpr LogLevel CH233 = CH033 | CH200;
    constexpr LogLevel CH300 = 9_bit;
    constexpr LogLevel CH301 = CH001 | CH300;
    constexpr LogLevel CH302 = CH002 | CH300;
    constexpr LogLevel CH303 = CH003 | CH300;
    constexpr LogLevel CH310 = CH010 | CH300;
    constexpr LogLevel CH311 = CH011 | CH300;
    constexpr LogLevel CH312 = CH012 | CH300;
    constexpr LogLevel CH313 = CH013 | CH300;
    constexpr LogLevel CH320 = CH020 | CH300;
    constexpr LogLevel CH321 = CH021 | CH300;
    constexpr LogLevel CH322 = CH022 | CH300;
    constexpr LogLevel CH323 = CH023 | CH300;
    constexpr LogLevel CH330 = CH030 | CH300;
    constexpr LogLevel CH331 = CH031 | CH300;
    constexpr LogLevel CH332 = CH032 | CH300;
    constexpr LogLevel CH333 = CH033 | CH300;

    struct log_options
    {
        struct SCRN_t
        {
            bool on = true;
        } SCRN;

        struct SOCK_t
        {
            bool on = false;
            std::string ip = "127.0.0.1";
            uint16_t port = 9998;
        } SOCK;

        struct FILE_t
        {
            bool on = false;
            bool compressed = false;
            uint32_t max_size_mb = 100;
            uint32_t max_size_all = 100 * 12;
            std::string directory = "";
            std::string rootname = "default";
        } FILE;

        LogLevel LVL = CH111;
    };

    class LogLine
    {
    private:
        LogLine();

        friend void details::psuedo_log_flush();

    public:
        struct string_literal_t
        {
            explicit string_literal_t(const char *s) : m_s(s) {}
            char const *m_s;
        };

        LogLine(LogLevel level, char const *file, char const *function, uint32_t line);

        LogLine(LogLevel level, char const *file, char const *function, uint32_t line, char const *ctx);

        void stringify(std::ostream &os, LogLevel mask = LogLevel::all(), unsigned rsh = 0);

        LogLevel lvl() const;

        LogLine &operator<<(char arg);
        LogLine &operator<<(int32_t arg);
        LogLine &operator<<(uint32_t arg);
        LogLine &operator<<(int64_t arg);
        LogLine &operator<<(uint64_t arg);
        LogLine &operator<<(double arg);
        LogLine &operator<<(const std::string &arg);

        template <size_t N>
        LogLine &operator<<(const char (&arg)[N])
        {
            encode_c_string(arg, (std::max)(N, size_t(1)) - 1);
            return *this;
        }

        template <typename Arg>
        typename std::enable_if<std::is_same<Arg, const char *>::value, LogLine &>::type
        operator<<(Arg const &arg)
        {
            encode(arg);
            return *this;
        }

        template <typename Arg>
        typename std::enable_if<std::is_same<Arg, char *>::value, LogLine &>::type
        operator<<(const Arg &arg)
        {
            encode(arg);
            return *this;
        }

        template <typename Arg>
        typename std::enable_if<details::is_overloaded_stream<Arg>() &&
                                    !std::is_same<Arg, const char *>::value &&
                                    !std::is_same<Arg, char *>::value,
                                LogLine &>::type
        operator<<(const Arg &arg)
        {
            std::ostringstream ss;
            ss << arg;
            *this << ss.str();
            return *this;
        }

    private:
        char *buffer();

        template <typename Arg>
        void encode(Arg arg)
        {
            *reinterpret_cast<Arg *>(buffer()) = arg;
            m_bytes_used += sizeof(Arg);
        }
        template <typename Arg>
        void encode(Arg arg, uint8_t type_id)
        {
            resize_buffer_if_needed(sizeof(Arg) + sizeof(uint8_t));
            encode<uint8_t>(type_id);
            encode<Arg>(arg);
        }

        void encode(char *arg);
        void encode(char const *arg);
        void encode(string_literal_t arg);
        void encode_c_string(char const *arg, size_t length);
        void resize_buffer_if_needed(size_t additional_bytes);
        void stringify(std::ostream &os, char *start, char const *end);

    public:
        bool m_ctrl_bytes;
        size_t m_bytes_used;
        size_t m_buffer_size;
        std::unique_ptr<char[]> m_heap_buffer;
        char m_stack_buffer[256 - 2 * sizeof(size_t) - sizeof(decltype(m_heap_buffer)) - 16];
    };

    namespace details
    {
        struct PsuedoLog
        {
            bool operator==(LogLine &);
        };

        PLOG_PRINTF_CHECK(5, 6)
        void psuedo_log_fmt(LogLevel level, char const *file, char const *function,
                            uint32_t line, PLOG_PRINT_STRING_TYPE format, ...);
    }

    bool is_logged(LogLevel level);

    void initialize_log(const log_options &opts = log_options{});
}

#define PPX_LOG(LEVEL) ppx::details::PsuedoLog() == ppx::LogLine(LEVEL, ppx::filename(__FILE__), ppx::funcname(__FUNCTION__), __LINE__)
#define LOG_CH(NUM) ppx::is_logged(ppx::CH##NUM) && PPX_LOG(ppx::CH##NUM)
#define LOG_FMT(NUM, ...) ppx::is_logged(ppx::CH##NUM)                                                                                                  \
                              ? ppx::details::psuedo_log_fmt(ppx::CH##NUM, ppx::filename(__FILE__), ppx::funcname(__FUNCTION__), __LINE__, __VA_ARGS__) \
                              : (void)0
#define LOG_IF(NUM, COND, ...) ppx::is_logged(ppx::CH##NUM) && (cond)                                                                                        \
                                   ? ppx::details::psuedo_log_fmt(ppx::CH##NUM, ppx::filename(__FILE__), ppx::funcname(__FUNCTION__), __LINE__, __VA_ARGS__) \
                                   : (void)0
#define FLUSH_LOG ppx::details::psuedo_log_flush()
#endif