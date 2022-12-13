#ifndef VVERY_SIMPLE_PPXLOG_HEADER
#define VVERY_SIMPLE_PPXLOG_HEADER

#include <cstdint>
#include <memory>
#include <string>
#include <type_traits>

namespace ppx
{
    enum class LogLevel : uint8_t
    {
        INFO,
        WARN,
        CRIT
    };

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
        LogLine &operator<<(std::string const &arg);

        template <size_t N>
        LogLine &operator<<(const char (&arg)[N])
        {
            encode(string_literal_t(arg));
            return *this;
        }

        template <typename Arg>
        typename std::enable_if<std::is_same<Arg, char const *>::value, LogLine &>::type
        operator<<(Arg const &arg)
        {
            encode(arg);
            return *this;
        }

        template <typename Arg>
        typename std::enable_if<std::is_same<Arg, char *>::value, LogLine &>::type
        operator<<(Arg const &arg)
        {
            encode(arg);
            return *this;
        }

        struct string_literal_t
        {
            explicit string_literal_t(char const *s) : m_s(s) {}
            char const *m_s;
        };

    private:
        char *buffer();

        template <typename Arg>
        void encode(Arg arg);

        template <typename Arg>
        void encode(Arg arg, uint8_t type_id);

        void encode(char *arg);
        void encode(char const *arg);
        void encode(string_literal_t arg);
        void encode_c_string(char const *arg, size_t length);
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

    void set_log_level(LogLevel level);

    bool is_logged(LogLevel level);

    struct RingBufferLogger
    {
        RingBufferLogger(uint32_t ring_buffer_size_mb_ = 4) : ring_buffer_size_mb(ring_buffer_size_mb_) {}
        uint32_t ring_buffer_size_mb;
    };

    void initialize(RingBufferLogger ngl, std::string const &log_directory, std::string const &log_file_name, uint32_t log_file_roll_size_mb);
}

#define PPX_LOG(LEVEL) ppx::PsuedoLog() == ppx::LogLine(LEVEL, __FILE__, __FUNCTION__, __LINE__)
#define LOG_INFO ppx::is_logged(ppx::LogLevel::INFO) && PPX_LOG(ppx::LogLevel::INFO)
#define LOG_WARN ppx::is_logged(ppx::LogLevel::WARN) && PPX_LOG(ppx::LogLevel::WARN)
#define LOG_CRIT ppx::is_logged(ppx::LogLevel::CRIT) && PPX_LOG(ppx::LogLevel::CRIT)

#endif