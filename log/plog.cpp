#include "plog.h"
#include <sys/stat.h>
#include <inttypes.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <vector>
#include <array>
#include <list>
#include <queue>
#include <tuple>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>
#include <fstream>
#include <iostream>
#if PLOG_OS_WINDOWS
#include <WinSock2.h>
#include <WS2tcpip.h>
#include <Windows.h>
#include <io.h>
#include <strsafe.h>
#else
#include <stdarg.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <dirent.h>
#endif

#if PLOG_OS_WINDOWS
#define PLOG_SAFE_WRITE ::_write
#define PLOG_SAFE_FILENO ::_fileno
#define PLOG_LOCALTIME(a, b) localtime_s(b, a)
#define PLOG_SAFE_SOCKADDR SOCKADDR
#define PLOG_SAFE_TCPPROTOCAL IPPROTO_TCP
#define PLOG_SAFE_INVALID_SOCKET INVALID_SOCKET
#else
#define PLOG_SAFE_WRITE ::write
#define PLOG_SAFE_FILENO ::fileno
#define PLOG_LOCALTIME(a, b) localtime_r(a, b)
#define PLOG_SAFE_SOCKADDR sockaddr
#define PLOG_SAFE_TCPPROTOCAL 0
#define PLOG_SAFE_INVALID_SOCKET -1
#endif

#if defined(_MSC_VER)
#if _MSVC_LANG > 201103L
#define PLOG_SAFE_MAKE_UNIQUE(PTR, TYPE, ...) PTR = std::make_unique<TYPE>(__VA_ARGS__)
#else
#define PLOG_SAFE_MAKE_UNIQUE(PTR, TYPE, ...) PTR.reset(new TYPE(__VA_ARGS__))
#endif
#else
#if __cplusplus > 201103L
#define PLOG_SAFE_MAKE_UNIQUE(PTR, TYPE, ...) PTR = std::make_unique<TYPE>(__VA_ARGS__)
#else
#define PLOG_SAFE_MAKE_UNIQUE(PTR, TYPE, ...) PTR.reset(new TYPE(__VA_ARGS__))
#endif
#endif

#if defined(__GNUC__)
#define PLOG_LIKELY(x) __builtin_expect(x, 1)
#define PLOG_UNLIKELY(x) __builtin_expect(x, 0)
#define PLOG_DONT_WARN_FMT_START _Pragma("GCC diagnostic push") \
    _Pragma("GCC diagnostic ignored \"-Wformat-truncation=\"")
#define PLOG_DONT_WARN_FMT_END _Pragma("GCC diagnostic pop")
#else
#define PLOG_LIKELY(x) x
#define PLOG_UNLIKELY(x) x
#define PLOG_DONT_WARN_FMT_START
#define PLOG_DONT_WARN_FMT_END
#endif

namespace ppx
{
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
        std::tuple<char, uint32_t, uint64_t, int32_t, int64_t, double, LogLine::string_literal_t, char *>;

#if PLOG_OS_QNX
    static std::once_flag start_flag;
#endif
    static constexpr size_t MAX_LINE_LENGTH = 1024 * 50;
    static constexpr size_t kMaxfilenum = 9;
    static constexpr LogLevel kScreenChannel = CH001 | CH002 | CH003;
    static constexpr LogLevel kSocketChannel = CH010 | CH020 | CH030;
    static constexpr LogLevel kDiskFileChannel = CH100 | CH200 | CH300;
    static constexpr std::array<char, 3> kCustomLabels{'I', 'W', 'E'};
    static std::vector<LogLine> MainThreadInnerMsg;

    PLOG_PRINTF_CHECK(1, 0)
    char *vtextprintf(const char *format, va_list vlist)
    {
        static thread_local char log_line_bufs[MAX_LINE_LENGTH + 8]{};
        auto fmt_result = vsnprintf(log_line_bufs, MAX_LINE_LENGTH, format, vlist);
        return fmt_result > -1 && fmt_result <= static_cast<int>(MAX_LINE_LENGTH) ? log_line_bufs : nullptr;
    }

    PLOG_PRINTF_CHECK(5, 6)
    void inner_log_fmt(LogLevel level, char const *file, char const *function,
                       uint32_t line, PLOG_PRINT_STRING_TYPE format, ...)
    {
        va_list vlist;
        va_start(vlist, format);
        auto *ctx = vtextprintf(format, vlist);
        va_end(vlist);
        MainThreadInnerMsg.emplace_back(level, file, function, line, ctx);
    }

#define PLOG_INNER_IMSG(...) inner_log_fmt(CH111, ppx::filename(__FILE__), ppx::funcname(__FUNCTION__), __LINE__, __VA_ARGS__)
#define PLOG_INNER_WMSG(...) inner_log_fmt(CH222, ppx::filename(__FILE__), ppx::funcname(__FUNCTION__), __LINE__, __VA_ARGS__)
#define PLOG_INNER_EMSG(...) inner_log_fmt(CH333, ppx::filename(__FILE__), ppx::funcname(__FUNCTION__), __LINE__, __VA_ARGS__)

    // file function related to file system. should be replace by <filesystem> after C++17
#if PLOG_OS_WINDOWS
    static constexpr const char *kFilePathSeparator = "\\";
#else
    static constexpr const char *kFilePathSeparator = "/";
#endif

    bool PathExist(const std::string &pathstr, bool considerFile = false)
    {
        const char *path = pathstr.c_str();
        if (!path)
        {
            return false;
        }
#if PLOG_OS_UNIX
        (void)considerFile;
        struct stat st;
        return stat(path, &st) == 0;
#else
        DWORD fileType = GetFileAttributesA(path);
        if (fileType == INVALID_FILE_ATTRIBUTES)
        {
            return false;
        }
        return considerFile || (fileType & FILE_ATTRIBUTE_DIRECTORY) != 0;
#endif
    }

    std::vector<std::string> DirectoryFileList(const std::string &directory)
    {
        std::vector<std::string> filelist;
#if PLOG_OS_WINDOWS
        TCHAR szDir[MAX_PATH];
        StringCchCopy(szDir, MAX_PATH, directory.c_str());
        StringCchCat(szDir, MAX_PATH, TEXT("\\*"));
        WIN32_FIND_DATA ffd;
        auto hFind = FindFirstFile(szDir, &ffd);
        if (INVALID_HANDLE_VALUE == hFind)
        {
            return filelist;
        }
        do
        {
            if (!(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
            {
                filelist.emplace_back(ffd.cFileName);
            }
        } while (FindNextFile(hFind, &ffd) != 0);
        auto dwError = GetLastError();
        if (dwError != ERROR_NO_MORE_FILES)
        {
            return filelist;
        }
        FindClose(hFind);
#else
        DIR *dir;
        struct dirent *diread;
        struct stat fileinfo;
        if ((dir = opendir(directory.c_str())) != nullptr)
        {
            while ((diread = readdir(dir)) != nullptr)
            {
#if PLOG_OS_QNX
                if (diread != nullptr &&
                    ::stat((directory + kFilePathSeparator + diread->d_name).c_str(), &fileinfo) == 0)
                {
                    if (S_ISREG(fileinfo.st_mode))
                    {
                        filelist.push_back(diread->d_name);
                    }
                }
#else
                (void)(fileinfo);
                if (diread->d_type == DT_REG)
                {
                    filelist.push_back(diread->d_name);
                }
#endif
            }
            closedir(dir);
        }
#endif
        return filelist;
    }

    void RotateFile(const std::string &directory, size_t max_directory_size)
    {
        if (!PathExist(directory))
        {
            return;
        }
        auto filelist = DirectoryFileList(directory);
        using finfo = std::tuple<std::string, size_t, time_t>;
        std::list<finfo> loglist;
        struct stat fileinfo;
        for (const auto &filename : filelist)
        {
            if (filename.find(".log") != std::string::npos)
            {
                std::string fname(directory);
                fname.append(kFilePathSeparator).append(filename);
                if (!::stat(fname.c_str(), &fileinfo))
                {
                    auto tmp = std::make_tuple(fname, fileinfo.st_size, fileinfo.st_mtime);
                    auto it = std::lower_bound(loglist.cbegin(), loglist.cend(), tmp, [](const finfo &a, const finfo &b)
                                               { return std::get<2>(a) > std::get<2>(b); });
                    loglist.insert(it, tmp);
                }
            }
        }
        size_t total_size = 0;
        for (const auto &log : loglist)
        {
            auto fsize = std::get<1>(log);
            if (total_size + fsize > max_directory_size)
            {
                std::remove(std::get<0>(log).c_str());
                PLOG_INNER_WMSG("remove file: %s", std::get<0>(log).c_str());
            }
            else
            {
                total_size += fsize;
            }
        }
    }

    uint64_t timestamp_now()
    {
#if PLOG_OS_QNX
        static uint64_t high_res_start_time = 0;
        static uint64_t high_res_start_tick = 0;
        std::call_once(start_flag, []()
                       { 
                        auto now_count = std::chrono::system_clock::now().time_since_epoch();
                        high_res_start_time=std::chrono::duration_cast<std::chrono::microseconds>(now_count).count();
                        uint64_t dTicksPerUs=SYSPAGE_ENTRY(qtime)->cycles_per_sec*0.000001;
                        high_res_start_tick=static_cast<uint64_t>(ClockCycles()/dTicksPerUs); });
        dTicksPerUs = SYSPAGE_ENTRY(qtime)->cycles_per_sec * 0.000001;
        auto clock_count = static_cast<uint64_t>(ClockCycles() / dTicksPerUs) - high_res_start_tick;
        return high_res_start_tick * 1000 + clock_count;
#else
        return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
#endif
    }

    uint64_t format_timestamp(char *buf, uint64_t timestamp)
    {
        auto cur_t = static_cast<std::time_t>(timestamp / 1000000);
        static thread_local std::time_t last_t = 0;
        static thread_local tm gmtime{};
        if (cur_t != last_t)
        {
            PLOG_LOCALTIME(&cur_t, &gmtime);
            last_t = cur_t;
        }
        PLOG_DONT_WARN_FMT_START
        snprintf(buf, 22, "[%04d-%02d-%02d %02d:%02d:%02d.", gmtime.tm_year + 1900, gmtime.tm_mon + 1,
                 gmtime.tm_mday, gmtime.tm_hour, gmtime.tm_min, gmtime.tm_sec);
        PLOG_DONT_WARN_FMT_END
        return timestamp % 1000000;
    }

    size_t this_thread_id()
    {
        static thread_local const auto id = std::hash<std::thread::id>{}(std::this_thread::get_id()) % 8209 + 1000;
        return id;
    }

    size_t this_process_id()
    {
#if PLOG_OS_WINDOWS
        static const size_t id = GetCurrentProcessId() % 8209 + 1000;
        return id;
#else
        static const size_t id = static_cast<uint32_t>(::getpid()) % 8209 + 1000;
        return id;
#endif
    }

    // decode
    template <typename Arg>
    char *decode(std::ostream &os, char *b, Arg *)
    {
        Arg arg = *reinterpret_cast<Arg *>(b);
        os << arg;
        return b + sizeof(Arg);
    }

    template <>
    char *decode(std::ostream &os, char *b, LogLine::string_literal_t *)
    {
        LogLine::string_literal_t s = *reinterpret_cast<LogLine::string_literal_t *>(b);
        os << s.m_s;
        return b + sizeof(LogLine::string_literal_t);
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

    // LogLine

    LogLine::LogLine() : m_ctrl_bytes(true), m_bytes_used(0), m_buffer_size(sizeof(m_stack_buffer)) {}

    LogLine::LogLine(LogLevel level, char const *file, char const *function, uint32_t line)
        : m_ctrl_bytes(false), m_bytes_used(0), m_buffer_size(sizeof(m_stack_buffer))
    {
        encode<LogLevel>(level);
        encode<uint64_t>(timestamp_now());
        encode<size_t>(this_thread_id());
        encode<string_literal_t>(string_literal_t(file));
        encode<string_literal_t>(string_literal_t(function));
        encode<uint32_t>(line);
    }

    LogLine::LogLine(LogLevel level, char const *file, char const *function, uint32_t line, char const *ctx)
        : LogLine(level, file, function, line)
    {
        encode(ctx);
    }

    LogLevel LogLine::lvl() const
    {
        const char *b = !m_heap_buffer ? m_stack_buffer : m_heap_buffer.get();
        return *reinterpret_cast<LogLevel *>(const_cast<char *>(b));
    }

    void LogLine::stringify(std::ostream &os, LogLevel mask, unsigned rsh)
    {
        if (PLOG_UNLIKELY(m_ctrl_bytes))
        {
            os.flush();
            return;
        }

        char *b = !m_heap_buffer ? m_stack_buffer : m_heap_buffer.get();
        char const *const end = b + m_bytes_used;
        auto loglvl = *reinterpret_cast<LogLevel *>(b);
        b += sizeof(LogLevel);
        auto timestamp = *reinterpret_cast<uint64_t *>(b);
        b += sizeof(uint64_t);
        auto threadid = *reinterpret_cast<size_t *>(b);
        b += sizeof(size_t);
        auto file = *reinterpret_cast<string_literal_t *>(b);
        b += sizeof(string_literal_t);
        auto function = *reinterpret_cast<string_literal_t *>(b);
        b += sizeof(string_literal_t);
        auto line = *reinterpret_cast<uint32_t *>(b);
        b += sizeof(uint32_t);

        auto buflen = strlen(file.m_s) + strlen(function.m_s);
        auto label = kCustomLabels.at((loglvl & mask).val() >> rsh);

        if (PLOG_UNLIKELY(buflen > 470))
        {
            char buf[30]{};
            auto micro_sec = format_timestamp(buf, timestamp);
            std::snprintf(buf + 21, 9, "%06" PRIu64 "] ", micro_sec);
            os << buf << '[' << label << threadid << ' ' << file.m_s << ':' << line << "][" << function.m_s << "] ";
        }
        else
        {
            char buf[512]{};
            auto micro_sec = format_timestamp(buf, timestamp);
            std::snprintf(buf + 21, 491, "%06" PRIu64 "][%c%zu %s:%d][%s] ", micro_sec, label, threadid, file.m_s, line, function.m_s);
            os << buf;
        }
        stringify(os, b, end);
        os << "\n";
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
        auto required_size = m_bytes_used + additional_bytes;

        if (required_size <= m_buffer_size)
        {
            return;
        }

        if (!m_heap_buffer)
        {
            m_buffer_size = (std::max)(static_cast<size_t>(512), required_size);
            m_heap_buffer.reset(new char[m_buffer_size]);
            memcpy(m_heap_buffer.get(), m_stack_buffer, m_bytes_used);
            return;
        }
        else
        {
            m_buffer_size = (std::max)(static_cast<size_t>(2 * m_buffer_size), required_size);
            std::unique_ptr<char[]> new_heap_buffer(new char[m_buffer_size]);
            memcpy(new_heap_buffer.get(), m_heap_buffer.get(), m_bytes_used);
            m_heap_buffer.swap(new_heap_buffer);
        }
    }

    void LogLine::encode(const char *arg)
    {
        if (PLOG_LIKELY(arg != nullptr))
        {
            encode_c_string(arg, strlen(arg));
        }
    }

    void LogLine::encode(char *arg)
    {
        if (PLOG_LIKELY(arg != nullptr))
        {
            encode_c_string(arg, strlen(arg));
        }
    }

    void LogLine::encode_c_string(char const *arg, size_t length)
    {
        if (PLOG_UNLIKELY(length == 0))
        {
            return;
        }

        resize_buffer_if_needed(length + 2);
        char *b = buffer();
        auto type_id = TupleIndex<char *, SupportedTypes>::value;
        *reinterpret_cast<uint8_t *>(b++) = static_cast<uint8_t>(type_id);
        memcpy(b, arg, length + 1);
        m_bytes_used += length + 2;
    }

    void LogLine::encode(string_literal_t arg)
    {
        encode<string_literal_t>(arg, TupleIndex<string_literal_t, SupportedTypes>::value);
    }

    LogLine &LogLine::operator<<(const std::string &arg)
    {
        encode_c_string(arg.c_str(), (std::min)(strlen(arg.c_str()), arg.length()));
        return *this;
    }

    LogLine &LogLine::operator<<(int32_t arg)
    {
        encode<int32_t>(arg, TupleIndex<int32_t, SupportedTypes>::value);
        return *this;
    }

    LogLine &LogLine::operator<<(uint32_t arg)
    {
        encode<uint32_t>(arg, TupleIndex<uint32_t, SupportedTypes>::value);
        return *this;
    }

    LogLine &LogLine::operator<<(int64_t arg)
    {
        encode<int64_t>(arg, TupleIndex<int64_t, SupportedTypes>::value);
        return *this;
    }

    LogLine &LogLine::operator<<(uint64_t arg)
    {
        encode<uint64_t>(arg, TupleIndex<uint64_t, SupportedTypes>::value);
        return *this;
    }

    LogLine &LogLine::operator<<(double arg)
    {
        encode<double>(arg, TupleIndex<double, SupportedTypes>::value);
        return *this;
    }

    LogLine &LogLine::operator<<(char arg)
    {
        encode<char>(arg, TupleIndex<char, SupportedTypes>::value);
        return *this;
    }

    LogLine CopyLogLine(const LogLine &other)
    {
        LogLine line(CH001, nullptr, nullptr, 0);
        line.m_buffer_size = other.m_buffer_size;
        line.m_bytes_used = other.m_bytes_used;
        if (other.m_heap_buffer)
        {
            line.m_heap_buffer.reset(new char[other.m_buffer_size]);
            memcpy(line.m_heap_buffer.get(), other.m_heap_buffer.get(), other.m_bytes_used);
        }
        memcpy(line.m_stack_buffer, other.m_stack_buffer, sizeof(other.m_stack_buffer));
        return line;
    }

    // basic infrastructure
    using CodeType = std::uint16_t;
    static constexpr CodeType dms = (std::numeric_limits<CodeType>::max)();

    class SearchTree
    {
        struct Node
        {
            CodeType first;
            CodeType left;
            CodeType right;
            char c;
            explicit Node(char _c) : first(dms), left(dms), right(dms), c(_c) {}
        };

        void reset()
        {
            vn.clear();
            constexpr int minc = (std::numeric_limits<char>::min)();
            constexpr int maxc = (std::numeric_limits<char>::max)();
            for (auto c = minc; c <= maxc; ++c)
            {
                vn.emplace_back(c);
            }
        }

    public:
        SearchTree()
        {
            constexpr int minc = (std::numeric_limits<char>::min)();
            constexpr int maxc = (std::numeric_limits<char>::max)();
            CodeType k{0};
            for (auto ic = minc; ic <= maxc; ic++)
            {
                initials[static_cast<unsigned char>(ic)] = k++;
            }
            vn.reserve(dms);
            reset();
        }

        CodeType search_and_insert(CodeType i, char c)
        {
            if (vn.size() == dms)
            {
                reset();
            }

            if (i == dms)
            {
                return search_initials(c);
            }

            const CodeType vn_size = static_cast<CodeType>(vn.size());
            CodeType ci{vn[i].first}; // Current Index

            if (ci != dms)
            {
                for (;;)
                {
                    if (c < vn[ci].c)
                    {
                        if (vn[ci].left == dms)
                        {
                            vn[ci].left = vn_size;
                            break;
                        }
                        else
                        {
                            ci = vn[ci].left;
                        }
                    }
                    else if (c > vn[ci].c)
                    {
                        if (vn[ci].right == dms)
                        {
                            vn[ci].right = vn_size;
                            break;
                        }
                        else
                        {
                            ci = vn[ci].right;
                        }
                    }
                    else
                    {
                        return ci;
                    }
                }
            }
            else
            {
                vn[i].first = vn_size;
            }
            vn.emplace_back(c);
            return dms;
        }

        CodeType search_initials(char c) const
        {
            return initials[static_cast<unsigned char>(c)];
        }

    private:
        std::vector<Node> vn;
        std::array<CodeType, 1u << CHAR_BIT> initials;
    };

    class FdOutBuf : public std::streambuf
    {
    private:
        static constexpr int bufferSize = 1400;

        const int fd;
        const bool sock_type;
        int tol_idx{0};
        char buffer[bufferSize]{};

    public:
        explicit FdOutBuf(int _fd, bool _isSock)
            : fd(_fd), sock_type(_isSock)
        {
            setp(buffer, buffer + (bufferSize - 1));
        }

        void clear()
        {
            sync();
        }

    protected:
        int flushBuffer()
        {
            auto ret = static_cast<int>(pptr() - pbase());
            if (sock_type)
            {
                ret = send(fd, buffer, ret, 0);
                if (PLOG_UNLIKELY(ret < 0))
                {
                    return EOF;
                }
            }
            else
            {
                if (PLOG_UNLIKELY(ret != PLOG_SAFE_WRITE(fd, buffer, ret)))
                {
                    return EOF;
                }
            }
            tol_idx += ret;
            pbump(-ret);
            return ret;
        }

        int_type overflow(int_type c) override
        {
            if (PLOG_LIKELY(c != EOF))
            {
                *pptr() = c;
                pbump(1);
            }
            if (flushBuffer() == EOF)
            {
                return EOF;
            }
            return c;
        }

        int sync() override
        {
            if (flushBuffer() == EOF)
            {
                return -1;
            }
            return 0;
        }

        std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,
                               std::ios_base::openmode mode) override
        {
            return mode == std::ios_base::out && way == std::ios_base::cur ? std::streampos(std::streamoff(tol_idx))
                                                                           : std::streampos(std::streamoff(-1));
        }
    };

    class FdCprdBuf : public std::streambuf
    {
    private:
        static const int bufferSize = 1280;

        const int fd;
        int cur_idx{0};
        int tol_idx{0};
        CodeType i{dms};
        CodeType buffer[bufferSize]{};
        SearchTree dictionary;

        bool write_buffer(CodeType code)
        {
            if (PLOG_UNLIKELY(cur_idx == bufferSize))
            {
                if (PLOG_SAFE_WRITE(fd, reinterpret_cast<const char *>(&buffer),
                                    bufferSize * sizeof(CodeType)) != bufferSize * sizeof(CodeType))
                {
                    return false;
                }
                cur_idx = 0;
                tol_idx += bufferSize * sizeof(CodeType);
            }
            buffer[cur_idx++] = code;
            return true;
        }

    public:
        explicit FdCprdBuf(int _fd) : fd(_fd) {}

        void clear()
        {
            auto num = PLOG_SAFE_WRITE(fd, reinterpret_cast<const char *>(&buffer),
                                       cur_idx * sizeof(CodeType));
            tol_idx += num;
        }

    protected:
        int_type overflow(int_type c) override
        {
            if (PLOG_LIKELY(c != EOF))
            {
                const CodeType temp{i};
                if ((i = dictionary.search_and_insert(temp, c)) == dms)
                {
                    if (!write_buffer(temp))
                    {
                        return EOF;
                    }
                    i = dictionary.search_initials(c);
                }
            }
            else
            {
                if (i != dms && !write_buffer(i))
                {
                    return EOF;
                }
            }
            return c;
        }

        std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,
                               std::ios_base::openmode mode) override
        {
            return mode == std::ios_base::out && way == std::ios_base::cur ? std::streampos(std::streamoff(tol_idx))
                                                                           : std::streampos(std::streamoff(-1));
        }
    };

    class osockstream : public std::ostream
    {
    private:
#if PLOG_OS_WINDOWS
        SOCKET sock_fd;
#else
        int sock_fd;
#endif
        std::unique_ptr<FdOutBuf> buf;

    public:
        osockstream(const std::string &ipaddr, uint16_t port)
            : std::ostream(nullptr), sock_fd(0)
        {
            sockaddr_in sa;
            memset(&sa, 0, sizeof(sa));
            sa.sin_family = AF_INET;
            inet_pton(AF_INET, ipaddr.c_str(), &sa.sin_addr);
            sa.sin_port = htons(port);

            sock_fd = socket(AF_INET, SOCK_STREAM, PLOG_SAFE_TCPPROTOCAL);
            if (sock_fd != PLOG_SAFE_INVALID_SOCKET)
            {
                PLOG_SAFE_MAKE_UNIQUE(buf, FdOutBuf, (int)sock_fd, true);
                timeval tv{4, 0};
                setsockopt(sock_fd, SOL_SOCKET, SO_SNDTIMEO, (char *)&tv, sizeof(tv));
            }

            if (buf)
            {
                rdbuf(buf.get());
            }

            if (connect(sock_fd, (const PLOG_SAFE_SOCKADDR *)&sa, sizeof(sa)) != 0)
            {
                setstate(std::ios::badbit);
            }
        }

        ~osockstream() override
        {
            if (good())
            {
                buf->clear();
            }

#if PLOG_OS_WINDOWS
            shutdown(sock_fd, SD_SEND);
            closesocket(sock_fd);
#else
            shutdown(sock_fd, SHUT_WR);
            close(sock_fd);
#endif
        }
    };

    class ofcprstream : public std::ostream
    {
    private:
        FILE *fptr;
        std::unique_ptr<FdCprdBuf> buf;

    public:
        ofcprstream(const char *filename) : std::ostream(nullptr)
        {
            fptr = ::fopen(filename, "wb");
            if (fptr)
            {
                PLOG_SAFE_MAKE_UNIQUE(buf, FdCprdBuf, PLOG_SAFE_FILENO(fptr));
                rdbuf(buf.get());
            }
        }

        ~ofcprstream() override
        {
            if (fptr)
            {
                buf->clear();
                ::fclose(fptr);
                fptr = nullptr;
            }
        }
    };

    class ofbufstream : public std::ostream
    {
    private:
        FILE *fptr;
        std::unique_ptr<FdOutBuf> buf;

    public:
        ofbufstream(const char *filename) : std::ostream(nullptr)
        {
            fptr = ::fopen(filename, "w");
            if (fptr)
            {
                PLOG_SAFE_MAKE_UNIQUE(buf, FdOutBuf, PLOG_SAFE_FILENO(fptr), false);
                rdbuf(buf.get());
            }
        }

        ~ofbufstream() override
        {
            if (fptr)
            {
                buf->clear();
                ::fclose(fptr);
                fptr = nullptr;
            }
        }
    };

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

    class BufferLeaf
    {
    public:
        struct Item
        {
            Item(LogLine &&logline) : logline(std::move(logline)) {}
            char padding[256 - sizeof(LogLine)];
            LogLine logline;
        };

        static constexpr const size_t SIZE = 32768;

        BufferLeaf() : m_buffer(static_cast<Item *>(std::malloc(SIZE * sizeof(Item))))
        {
            for (size_t i = 0; i <= SIZE; ++i)
            {
                m_write_state[i].store(0, std::memory_order_relaxed);
            }
            static_assert(sizeof(Item) == 256, "Unexpected size != 256");
        }

        ~BufferLeaf()
        {
            unsigned int write_count = m_write_state[SIZE].load(std::memory_order_acquire);
            for (size_t i = 0; i < write_count; ++i)
            {
                m_buffer[i].~Item();
            }
            std::free(m_buffer);
        }

        // Returns true if we need to switch to next buffer
        bool push(LogLine &&logline, unsigned int const write_index)
        {
            new (&m_buffer[write_index]) Item(std::move(logline));
            m_write_state[write_index].store(1, std::memory_order_release);
            return m_write_state[SIZE].fetch_add(1, std::memory_order_acquire) + 1 == SIZE;
        }

        bool try_pop(LogLine &logline, unsigned int const read_index)
        {
            if (m_write_state[read_index].load(std::memory_order_acquire))
            {
                Item &item = m_buffer[read_index];
                logline = std::move(item.logline);
                return true;
            }
            return false;
        }

        BufferLeaf(BufferLeaf const &) = delete;
        BufferLeaf &operator=(BufferLeaf const &) = delete;

    private:
        Item *m_buffer;
        std::atomic<unsigned int> m_write_state[SIZE + 1];
    };

    class BufferTree
    {
    public:
        BufferTree(BufferTree const &) = delete;
        BufferTree &operator=(BufferTree const &) = delete;

        BufferTree() : m_current_read_buffer{nullptr}, m_write_index(0), m_flag{ATOMIC_FLAG_INIT}, m_read_index(0), m_queue_len(0)
        {
            setup_next_write_buffer();
        }

        void push(LogLine &&logline)
        {
            if (PLOG_UNLIKELY(m_queue_len > 11))
            {
                return;
            }
            unsigned int write_index = m_write_index.fetch_add(1, std::memory_order_relaxed);
            if (write_index < BufferLeaf::SIZE)
            {
                if (m_current_write_buffer.load(std::memory_order_acquire)->push(std::move(logline), write_index))
                {
                    setup_next_write_buffer();
                }
            }
            else
            {
                while (m_write_index.load(std::memory_order_acquire) >= BufferLeaf::SIZE)
                    ;
                push(std::move(logline));
            }
        }

        bool try_pop(LogLine &logline)
        {
            if (m_current_read_buffer == nullptr)
            {
                m_current_read_buffer = get_next_read_buffer();
            }

            BufferLeaf *read_buffer = m_current_read_buffer;

            if (PLOG_UNLIKELY(read_buffer == nullptr))
            {
                return false;
            }

            if (read_buffer->try_pop(logline, m_read_index))
            {
                m_read_index++;
                if (m_read_index == BufferLeaf::SIZE)
                {
                    m_read_index = 0;
                    m_current_read_buffer = nullptr;
                    SpinLock spinlock(m_flag);
                    m_buffers.pop();
                    m_queue_len--;
                }
                return true;
            }

            return false;
        }

    private:
        void setup_next_write_buffer()
        {
            std::unique_ptr<BufferLeaf> next_write_buffer(new BufferLeaf());
            m_current_write_buffer.store(next_write_buffer.get(), std::memory_order_release);
            SpinLock spinlock(m_flag);
            m_buffers.push(std::move(next_write_buffer));
            m_write_index.store(0, std::memory_order_relaxed);
            m_queue_len++;
        }

        BufferLeaf *get_next_read_buffer()
        {
            SpinLock spinlock(m_flag);
            return m_buffers.empty() ? nullptr : m_buffers.front().get();
        }

    private:
        std::queue<std::unique_ptr<BufferLeaf>> m_buffers;
        std::atomic<BufferLeaf *> m_current_write_buffer;
        BufferLeaf *m_current_read_buffer;
        std::atomic<unsigned int> m_write_index;
        std::atomic_flag m_flag;
        unsigned int m_read_index;
        unsigned int m_queue_len;
    };

    struct FileWriter
    {
    public:
        FileWriter(std::string const &log_directory, std::string const &log_file_name,
                   uint32_t log_file_roll_size_mb, bool enable_compressed)
            : m_log_file_roll_size_bytes(log_file_roll_size_mb * 1024 * 1024),
              m_name(log_directory + log_file_name), m_enable_cprd(enable_compressed)
        {

            roll_file();
        }

        void write(LogLine &logline)
        {
            logline.stringify(*m_os, kDiskFileChannel, 8);
            if (m_os->tellp() > m_log_file_roll_size_bytes)
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
            }

            if (m_file_number > kMaxfilenum)
            {
                m_file_number = 0;
            }

            if (m_enable_cprd)
            {
                PLOG_SAFE_MAKE_UNIQUE(m_os, ofcprstream, (m_name + "_" + std::to_string(m_file_number++) + ".logz").c_str());
            }
            else
            {
                PLOG_SAFE_MAKE_UNIQUE(m_os, ofbufstream, (m_name + "_" + std::to_string(m_file_number++) + ".log").c_str());
            }
        }

    private:
        uint32_t m_file_number = 0;
        const uint32_t m_log_file_roll_size_bytes;
        const std::string m_name;
        const bool m_enable_cprd;
        std::unique_ptr<std::ostream> m_os;
    };

    class ScreenWriter
    {
    public:
        void write(LogLine &logline)
        {
            logline.stringify(std::cout, kScreenChannel, 2);
        }
    };

    class SocketWriter
    {
        uint16_t port;
        const std::string ipaddr;
        std::unique_ptr<osockstream> m_os;

    public:
        const std::string m_uid;

    public:
        SocketWriter(const std::string &_ipaddr, uint16_t _port)
            : port(_port), ipaddr(_ipaddr)
        {
#if PLOG_OS_WINDOWS
            WSADATA wsaData = {0};
            auto iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
            if (iResult == 0)
            {
                PLOG_SAFE_MAKE_UNIQUE(m_os, osockstream, ipaddr, port);
            }
#else
            PLOG_SAFE_MAKE_UNIQUE(m_os, osockstream, ipaddr, port);
#endif
        }

        ~SocketWriter()
        {
            m_os.reset();
        }

        void write(LogLine &logline)
        {
            static auto reConnTP = std::chrono::system_clock::now();
            if (m_os->good())
            {
                logline.stringify(*m_os, kSocketChannel, 5);
            }
            else
            {
                auto now = std::chrono::system_clock::now();
                if (now - reConnTP > std::chrono::seconds(4))
                {
                    reConnTP = now;
                    m_os.reset(new osockstream(ipaddr, port));
                }
            }
        }
    };

    template <typename T>
    class Logger
    {
        enum class State
        {
            INIT,
            READY,
            SHUTDOWN
        };

    public:
        template <typename... Ts>
        Logger(Ts... args)
            : m_state(State::INIT),
              m_writer(std::forward<Ts>(args)...),
              m_thread(&Logger::pop, this)
        {
            m_state.store(State::READY, std::memory_order_release);
        }

        ~Logger()
        {
            m_state.store(State::SHUTDOWN, std::memory_order_release);
            m_thread.join();
        }

        void add(LogLine &&logline)
        {
            m_buffer.push(std::move(logline));
        }

        void pop()
        {
            while (m_state.load(std::memory_order_acquire) == State::INIT)
            {
                std::this_thread::sleep_for(std::chrono::microseconds(50));
            }

            LogLine logline(CH000, nullptr, nullptr, 0);

            while (m_state.load() == State::READY)
            {
                if (m_buffer.try_pop(logline))
                {
                    m_writer.write(logline);
                }
                else
                {
                    std::this_thread::sleep_for(std::chrono::milliseconds(4));
                }
            }

            while (m_buffer.try_pop(logline))
            {
                m_writer.write(logline);
            }
        }

    private:
        std::atomic<State> m_state;
        T m_writer;
        BufferTree m_buffer;
        std::thread m_thread;
    };

    namespace
    {
        using FLogger = Logger<FileWriter>;
        using SLogger = Logger<SocketWriter>;
        using TLogger = Logger<ScreenWriter>;
        std::unique_ptr<FLogger> flogger;
        std::unique_ptr<SLogger> slogger;
        std::unique_ptr<TLogger> tlogger;
        std::atomic<FLogger *> atomic_flogger;
        std::atomic<SLogger *> atomic_slogger;
        std::atomic<TLogger *> atomic_tlogger;
        std::atomic<uint16_t> loglevel{1};
    }

    bool details::PsuedoLog::operator==(LogLine &logline)
    {
        auto level = logline.lvl();

        if (level & kDiskFileChannel)
        {
            auto lg = atomic_flogger.load(std::memory_order_acquire);
            if (lg)
            {
                lg->add(CopyLogLine(logline));
            }
        }
        if (level & kScreenChannel)
        {
            auto lg = atomic_tlogger.load(std::memory_order_acquire);
            if (lg)
            {
                lg->add(CopyLogLine(logline));
            }
        }
        if (level & kSocketChannel)
        {
            auto lg = atomic_slogger.load(std::memory_order_acquire);
            if (lg)
            {
                lg->add(std::move(logline));
            }
        }

        return true;
    }

    void details::psuedo_log_fmt(LogLevel level, char const *file, char const *function,
                                 uint32_t line, PLOG_PRINT_STRING_TYPE format, ...)
    {
        va_list vlist;
        va_start(vlist, format);
        auto *ctx = vtextprintf(format, vlist);
        va_end(vlist);

        if (level & kDiskFileChannel)
        {
            auto lg = atomic_flogger.load(std::memory_order_acquire);
            if (lg)
            {
                lg->add(LogLine(level, file, function, line, ctx));
            }
        }
        if (level & kSocketChannel)
        {
            auto lg = atomic_slogger.load(std::memory_order_acquire);
            if (lg)
            {
                lg->add(LogLine(level, file, function, line, ctx));
            }
        }
        if (level & kScreenChannel)
        {
            auto lg = atomic_tlogger.load(std::memory_order_acquire);
            if (lg)
            {
                lg->add(LogLine(level, file, function, line, ctx));
            }
        }
    }

    void details::psuedo_log_flush()
    {
        auto lg = atomic_tlogger.load(std::memory_order_acquire);
        if (lg)
        {
            lg->add(LogLine{});
        }
    }

    void set_log_level(LogLevel level)
    {
        loglevel.store(level.val(), std::memory_order_release);
    }

    bool options_legality_check(const log_options &opts)
    {
        if (opts.SOCK.on)
        {
            struct sockaddr_in sa;
            if (inet_pton(AF_INET, opts.SOCK.ip.c_str(), &(sa.sin_addr)) != 1)
            {
                PLOG_INNER_WMSG("%s is an invalid ip, no socket will be created.", opts.SOCK.ip.c_str());
                return false;
            }
        }

        if (opts.FILE.on)
        {
            if (opts.FILE.directory.empty())
            {
                PLOG_INNER_IMSG("storage directory is empty, path to program runtime will be used.");
            }
        }

        return true;
    }

    void initialize_log(const log_options &input_opts)
    {
        log_options opts;
        if (options_legality_check(input_opts))
        {
            opts = input_opts;
        }

        if (opts.SCRN.on)
        {
            opts.LVL = opts.LVL | kScreenChannel;
            PLOG_SAFE_MAKE_UNIQUE(tlogger, TLogger);
            atomic_tlogger.store(tlogger.get(), std::memory_order_release);
            PLOG_INNER_IMSG("teminal output is enabled.");
        }
        else
        {
            opts.LVL = opts.LVL & (~kScreenChannel);
            PLOG_INNER_WMSG("teminal output is disabled.");
        }

        if (opts.SOCK.on)
        {
            opts.LVL = opts.LVL | kSocketChannel;
            PLOG_SAFE_MAKE_UNIQUE(slogger, SLogger, opts.SOCK.ip, opts.SOCK.port);
            atomic_slogger.store(slogger.get(), std::memory_order_release);
            std::stringstream ss;
            PLOG_INNER_IMSG("net server output is enabled, ip address of log server is %s:%hu", opts.SOCK.ip.c_str(), opts.SOCK.port);
        }
        else
        {
            opts.LVL = opts.LVL & (~kSocketChannel);
            PLOG_INNER_IMSG("net server output is disabled.");
        }

        if (opts.FILE.on)
        {
            opts.LVL = opts.LVL | kSocketChannel;
            if (opts.FILE.directory.empty())
            {
                opts.FILE.directory = std::string(".") + kFilePathSeparator;
            }
            RotateFile(opts.FILE.directory, opts.FILE.max_size_all * 1024 * 1024);
            char buffer[17]{};
            auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            tm gmtime{};
            PLOG_LOCALTIME(&now, &gmtime);
            strftime(buffer, 17, "_%Y%m%d%H%M%S", &gmtime);
            opts.FILE.rootname += buffer;
            PLOG_SAFE_MAKE_UNIQUE(flogger, FLogger, opts.FILE.directory, opts.FILE.rootname, (std::max)(1u, opts.FILE.max_size_mb), opts.FILE.compressed);
            atomic_flogger.store(flogger.get(), std::memory_order_release);
            PLOG_INNER_IMSG("disk output is enabled, storage directory is %s, max size of a single file is %u Mib,  max capacity of directory is %u Mib, in %s format.",
                            opts.FILE.directory.c_str(), opts.FILE.max_size_mb, opts.FILE.max_size_all, opts.FILE.compressed ? "compressed" : "normal text");
        }
        else
        {
            opts.LVL = opts.LVL & (~kDiskFileChannel);
            PLOG_INNER_IMSG("disk output is disabled.");
        }

        set_log_level(opts.LVL);
        PLOG_INNER_IMSG("log level set to %hu", opts.LVL.val());

        for (const auto &msg : MainThreadInnerMsg)
        {
            if (tlogger)
            {
                tlogger->add(CopyLogLine(msg));
            }
            if (slogger)
            {
                slogger->add(CopyLogLine(msg));
            }
            if (flogger)
            {
                flogger->add(CopyLogLine(msg));
            }
        }
    }

    bool is_logged(LogLevel level)
    {
        return (bool)(level & LogLevel(loglevel.load(std::memory_order_acquire)));
    }
}
