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
    enum class LogLevel : uint8_t
    {
        INFO,
        WARN,
        CRIT
    };

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
                    *this << i << ' ';
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
            void encode(int *arg);
            void encode(const int *arg);
            void encode(details::string_literal_t arg);
            void encode_c_string(const char *arg, size_t length);
            void resize_buffer_if_needed(size_t additional_bytes);
            void stringify(std::ostream &os, char *start, char const *const end);

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
#define LOG_INFO ppx::is_logged(ppx::LogLevel::INFO) && PPX_LOG(ppx::LogLevel::INFO)
#define LOG_WARN ppx::is_logged(ppx::LogLevel::WARN) && PPX_LOG(ppx::LogLevel::WARN)
#define LOG_CRIT ppx::is_logged(ppx::LogLevel::CRIT) && PPX_LOG(ppx::LogLevel::CRIT)

#endif