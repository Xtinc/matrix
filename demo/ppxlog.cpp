#include "ppxlog.h"
#include <cstring>
#include <chrono>
#include <ctime>
#include <thread>
#include <atomic>
#include <fstream>
#include <cinttypes>

namespace ppx
{
    struct RingBufferLogger
    {
        RingBufferLogger(uint32_t ring_buffer_size_mb_ = 4) : ring_buffer_size_mb(ring_buffer_size_mb_) {}
        uint32_t ring_buffer_size_mb;
    };

    uint64_t timestamp_now()
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    }

    // I want [2022-12-13 00:01:23.528514], not efficient
    void format_timestamp(std::ostream &os, uint64_t timestamp)
    {
        std::time_t time_t = timestamp / 1000000;
        auto gmtime = std::localtime(&time_t);
        char buffer[32];
        strftime(buffer, 32, "%Y-%m-%d %T.", gmtime);
        char microseconds[7];
        sprintf(microseconds, "%06" PRIu64, timestamp % 1000000);
        os << '[' << buffer << microseconds << ']';
    }

    std::thread::id this_thread_id()
    {
        static thread_local const auto id = std::this_thread::get_id();
        return id;
    }

    int getlog2(LogLevel lev)
    {
        unsigned n = lev.val();
        int count = 0;
        while (n)
        {
            n = n >> 1;
            count++;
        }
        return count;
    }

    template <typename Arg>
    char *decode(std::ostream &os, char *b, Arg *)
    {
        Arg arg = *reinterpret_cast<Arg *>(b);
        os << arg;
        return b + sizeof(Arg);
    }

    template <>
    char *decode(std::ostream &os, char *b, details::string_literal_t *)
    {
        details::string_literal_t s = *reinterpret_cast<details::string_literal_t *>(b);
        os << s.m_s;
        return b + sizeof(details::string_literal_t);
    }

    template <>
    char *decode(std::ostream &os, char *b, char **)
    {
        while (*b != '\0')
        {
            os << *b;
            ++b;
        }
        return ++b;
    }

    namespace details
    {
        template <typename Arg>
        void LogLine::encode(Arg arg)
        {
            *reinterpret_cast<Arg *>(buffer()) = arg;
            m_bytes_used += sizeof(Arg);
        }

        template <typename Arg>
        void LogLine::encode(Arg arg, uint8_t type_id)
        {
            resize_buffer_if_needed(sizeof(Arg) + sizeof(uint8_t));
            encode<uint8_t>(type_id);
            encode<Arg>(arg);
        }

        LogLine::LogLine(LogLevel level, char const *file, char const *function, uint32_t line)
            : m_bytes_used(0), m_buffer_size(sizeof(m_stack_buffer))
        {
            encode<uint64_t>(timestamp_now());
            encode<std::thread::id>(this_thread_id());
            encode<details::string_literal_t>(details::string_literal_t(file));
            encode<details::string_literal_t>(details::string_literal_t(function));
            encode<uint32_t>(line);
            encode<LogLevel>(level);
        }

        void LogLine::stringify(std::ostream &os)
        {
            char *b = !m_heap_buffer ? m_stack_buffer : m_heap_buffer.get();
            char const *const end = b + m_bytes_used;
            auto timestamp = *reinterpret_cast<uint64_t *>(b);
            b += sizeof(uint64_t);
            auto threadid = *reinterpret_cast<std::thread::id *>(b);
            b += sizeof(std::thread::id);
            auto file = *reinterpret_cast<details::string_literal_t *>(b);
            b += sizeof(details::string_literal_t);
            auto function = *reinterpret_cast<details::string_literal_t *>(b);
            b += sizeof(details::string_literal_t);
            auto line = *reinterpret_cast<uint32_t *>(b);
            b += sizeof(uint32_t);
            auto loglevel = *reinterpret_cast<LogLevel *>(b);
            b += sizeof(LogLevel);

            format_timestamp(os, timestamp);

            os << '[' << getlog2(loglevel) << ']'
               << '[' << threadid << ']'
               << '[' << file.m_s << ':' << function.m_s << ':' << line << "] ";

            stringify(os, b, end);

            os << "\n";

            if (loglevel & CH01)
            {
                os.flush();
            }
        }

        void LogLine::stringify(std::ostream &os, char *start, char const *const end)
        {
            if (start == end)
            {
                return;
            }

            int type_id = static_cast<int>(*start);
            start++;

            switch (type_id)
            {
            case 0:
                stringify(os, decode(os, start, static_cast<std::tuple_element<0, SupportedTypes>::type *>(nullptr)), end);
                return;
            case 1:
                stringify(os, decode(os, start, static_cast<std::tuple_element<1, SupportedTypes>::type *>(nullptr)), end);
                return;
            case 2:
                stringify(os, decode(os, start, static_cast<std::tuple_element<2, SupportedTypes>::type *>(nullptr)), end);
                return;
            case 3:
                stringify(os, decode(os, start, static_cast<std::tuple_element<3, SupportedTypes>::type *>(nullptr)), end);
                return;
            case 4:
                stringify(os, decode(os, start, static_cast<std::tuple_element<4, SupportedTypes>::type *>(nullptr)), end);
                return;
            case 5:
                stringify(os, decode(os, start, static_cast<std::tuple_element<5, SupportedTypes>::type *>(nullptr)), end);
                return;
            case 6:
                stringify(os, decode(os, start, static_cast<std::tuple_element<6, SupportedTypes>::type *>(nullptr)), end);
                return;
            case 7:
                stringify(os, decode(os, start, static_cast<std::tuple_element<7, SupportedTypes>::type *>(nullptr)), end);
                return;
            }
        }

        char *LogLine::buffer()
        {
            return !m_heap_buffer ? &m_stack_buffer[m_bytes_used] : &(m_heap_buffer.get())[m_bytes_used];
        }

        void LogLine::resize_buffer_if_needed(size_t additional_bytes)
        {
            size_t const required_size = m_bytes_used + additional_bytes;

            if (required_size <= m_buffer_size)
            {
                return;
            }

            if (!m_heap_buffer)
            {
                m_buffer_size = std::max(static_cast<size_t>(512), required_size);
                m_heap_buffer.reset(new char[m_buffer_size]);
                memcpy(m_heap_buffer.get(), m_stack_buffer, m_bytes_used);
                return;
            }
            else
            {
                m_buffer_size = std::max(static_cast<size_t>(2 * m_buffer_size), required_size);
                std::unique_ptr<char[]> new_heap_buffer(new char[m_buffer_size]);
                memcpy(new_heap_buffer.get(), m_heap_buffer.get(), m_bytes_used);
                m_heap_buffer.swap(new_heap_buffer);
            }
        }

        void LogLine::encode(const char *arg)
        {
            if (arg != nullptr)
            {
                encode_c_string(arg, strlen(arg));
            }
        }

        void LogLine::encode(char *arg)
        {
            if (arg != nullptr)
            {
                encode_c_string(arg, strlen(arg));
            }
        }

        void LogLine::encode_c_string(char const *arg, size_t length)
        {
            if (length == 0)
            {
                return;
            }

            resize_buffer_if_needed(1 + length + 1);
            char *b = buffer();
            auto type_id = details::TupleIndex<char *, SupportedTypes>::value;
            *reinterpret_cast<uint8_t *>(b++) = static_cast<uint8_t>(type_id);
            memcpy(b, arg, length + 1);
            m_bytes_used += 1 + length + 1;
        }

        void LogLine::encode(details::string_literal_t arg)
        {
            encode<details::string_literal_t>(arg, details::TupleIndex<details::string_literal_t, SupportedTypes>::value);
        }

        LogLine &LogLine::operator<<(const std::string &arg)
        {
            encode_c_string(arg.c_str(), arg.length());
            return *this;
        }

        LogLine &LogLine::operator<<(int32_t arg)
        {
            encode<int32_t>(arg, details::TupleIndex<int32_t, SupportedTypes>::value);
            return *this;
        }

        LogLine &LogLine::operator<<(uint32_t arg)
        {
            encode<uint32_t>(arg, details::TupleIndex<uint32_t, SupportedTypes>::value);
            return *this;
        }

        LogLine &LogLine::operator<<(int64_t arg)
        {
            encode<int64_t>(arg, details::TupleIndex<int64_t, SupportedTypes>::value);
            return *this;
        }

        LogLine &LogLine::operator<<(uint64_t arg)
        {
            encode<uint64_t>(arg, details::TupleIndex<uint64_t, SupportedTypes>::value);
            return *this;
        }

        LogLine &LogLine::operator<<(double arg)
        {
            encode<double>(arg, details::TupleIndex<double, SupportedTypes>::value);
            return *this;
        }

        LogLine &LogLine::operator<<(char arg)
        {
            encode<char>(arg, details::TupleIndex<char, SupportedTypes>::value);
            return *this;
        }
    }

    struct SpinLock
    {
        SpinLock(std::atomic_flag &flag) : m_flag(flag)
        {
            while (m_flag.test_and_set(std::memory_order_acquire))
            {
            };
        }

        ~SpinLock()
        {
            m_flag.clear(std::memory_order_release);
        }

    private:
        std::atomic_flag &m_flag;
    };

    class RingBuffer
    {
    public:
        struct alignas(64) Item
        {
            Item() : flag{ATOMIC_FLAG_INIT}, written(0), logline(CH01, nullptr, nullptr, 0)
            {
            }

            std::atomic_flag flag;
            char written;
            char padding[256 - sizeof(std::atomic_flag) - sizeof(char) - sizeof(details::LogLine)];
            details::LogLine logline;
        };

        RingBuffer(const size_t size)
            : m_size(size), m_ring(static_cast<Item *>(std::malloc(size * sizeof(Item)))), m_write_index(0), m_read_index(0)
        {
            for (size_t i = 0; i < m_size; ++i)
            {
                new (&m_ring[i]) Item();
            }
            static_assert(sizeof(Item) == 256, "Unexpected size != 256");
        }

        ~RingBuffer()
        {
            for (size_t i = 0; i < m_size; ++i)
            {
                m_ring[i].~Item();
            }
            std::free(m_ring);
        }

        void push(details::LogLine &&logline)
        {
            unsigned int write_index = m_write_index.fetch_add(1, std::memory_order_relaxed) % m_size;
            Item &item = m_ring[write_index];
            SpinLock spinlock(item.flag);
            item.logline = std::move(logline);
            item.written = 1;
        }

        bool try_pop(details::LogLine &logline)
        {
            Item &item = m_ring[m_read_index % m_size];
            SpinLock spinlock(item.flag);
            if (item.written == 1)
            {
                logline = std::move(item.logline);
                item.written = 0;
                ++m_read_index;
                return true;
            }
            return false;
        }

        RingBuffer(const RingBuffer &) = delete;
        RingBuffer &operator=(const RingBuffer &) = delete;

    private:
        const size_t m_size;
        Item *m_ring;
        std::atomic<unsigned int> m_write_index;
        char pad[64];
        unsigned int m_read_index;
    };

    class FileWriter
    {
    public:
        FileWriter(std::string const &log_directory, std::string const &log_file_name, uint32_t log_file_roll_size_mb)
            : m_log_file_roll_size_bytes(log_file_roll_size_mb * 1024 * 1024), m_name(log_directory + log_file_name)
        {
            roll_file();
        }

        void write(details::LogLine &logline, std::chrono::time_point<std::chrono::system_clock> &start)
        {
            auto pos = m_os->tellp();
            logline.stringify(*m_os);
            m_bytes_written += m_os->tellp() - pos;
            auto bytes = m_bytes_written - r_bytes_written;
            if (bytes > 1000 * 1000)
            {
                auto now = std::chrono::system_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
                *m_os << "logs 1MegaBytes in " << duration << "ms, " << bytes / duration << "kb/s.\n";
                r_bytes_written = m_bytes_written;
                start = now;
                m_bytes_written += 36;
            }
            if (m_bytes_written > m_log_file_roll_size_bytes)
            {
                roll_file();
            }
        }

    private:
        void roll_file()
        {
            if (m_os)
            {
                m_os->flush();
                m_os->close();
            }

            if (m_file_number > 10)
            {
                m_file_number = 0;
            }

            m_bytes_written = 0;
            r_bytes_written = 0;
            m_os.reset(new std::ofstream());
            std::string log_file_name = m_name;
            log_file_name.append("_");
            log_file_name.append(std::to_string(++m_file_number));
            log_file_name.append(".txt");
            m_os->open(log_file_name, std::ofstream::out | std::ofstream::trunc);
        }

    private:
        uint32_t m_file_number = 0;
        std::streamoff m_bytes_written = 0;
        std::streamoff r_bytes_written = 0;
        uint32_t const m_log_file_roll_size_bytes;
        std::string const m_name;
        std::unique_ptr<std::ofstream> m_os;
    };

    class Logger
    {
    public:
        Logger(RingBufferLogger ngl, std::string const &log_directory,
               std::string const &log_file_name, uint32_t log_file_roll_size_mb)
            : m_state(State::INIT), m_buffer_base(new RingBuffer(std::max(1u, ngl.ring_buffer_size_mb) * 1024 * 4)),
              m_file_writer(log_directory, log_file_name, std::max(1u, log_file_roll_size_mb)), m_thread(&Logger::pop, this)
        {
            m_state.store(State::READY, std::memory_order_release);
        }

        ~Logger()
        {
            m_state.store(State::SHUTDOWN);
            m_thread.join();
        }

        void add(details::LogLine &&logline)
        {
            m_buffer_base->push(std::move(logline));
        }

        void pop()
        {
            // Wait for constructor to complete and pull all stores done there to this thread / core.
            while (m_state.load(std::memory_order_acquire) == State::INIT)
            {
                std::this_thread::sleep_for(std::chrono::microseconds(50));
            }

            details::LogLine logline(CH01, nullptr, nullptr, 0);
            auto now = std::chrono::system_clock::now();

            while (m_state.load() == State::READY)
            {
                if (m_buffer_base->try_pop(logline))
                {
                    m_file_writer.write(logline, now);
                }
                else
                {
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                }
            }

            // Pop and log all remaining entries
            while (m_buffer_base->try_pop(logline))
            {
                m_file_writer.write(logline, now);
            }
        }

    private:
        enum class State
        {
            INIT,
            READY,
            SHUTDOWN
        };

        std::atomic<State> m_state;
        std::unique_ptr<RingBuffer> m_buffer_base;
        FileWriter m_file_writer;
        std::thread m_thread;
    };

    std::unique_ptr<Logger> logger;
    std::atomic<Logger *> atomic_logger;

    bool details::PsuedoLog::operator==(details::LogLine &logline)
    {
        atomic_logger.load(std::memory_order_acquire)->add(std::move(logline));
        return true;
    }

    void initialize_log(const std::string &log_directory,
                        const std::string &log_file_name, uint32_t log_file_roll_size_mb)
    {
        logger.reset(new Logger(RingBufferLogger(), log_directory, log_file_name, log_file_roll_size_mb));
        atomic_logger.store(logger.get());
    }

    std::atomic<uint8_t> loglevel = {0};

    void set_log_level(LogLevel level)
    {
        loglevel.store(level.val());
    }

    bool is_logged(LogLevel level)
    {
        return (bool)(level & LogLevel(loglevel.load()));
    }
}