#pragma once

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h" // For console logging
#include "spdlog/sinks/basic_file_sink.h" // For file logging
#include <sys/stat.h> // 用于跨平台创建目录
#include <direct.h> // Windows 平台使用的目录操作

namespace EnSC {

class Logger {
public:
    // 创建目录辅助函数 (C++14 兼容)
    static bool createDirectory(const std::string& path) {
        #ifdef _WIN32
            return _mkdir(path.c_str()) == 0 || errno == EEXIST;
        #else
            return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
        #endif
    }

    // Initializes the logging system (e.g., creates a console logger)
    static void init() {
        try {
            // Create a color console logger
            // %^ marks the start of color range, %$ marks the end
            // [%Y-%m-%d %H:%M:%S.%e] [%n] [%l] %v%$ - Timestamp, logger name, level, message
            auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            console_sink->set_pattern("[%^%Y-%m-%d %H:%M:%S.%e%$] [%n] [%l] %v");
            
            // 根据不同构建类型设置不同的日志级别和输出方式
            #if defined(_DEBUG) || defined(DEBUG) || defined(CMAKE_BUILD_TYPE_DEBUG)
                // Debug 模式: 显示所有调试信息，同时输出到控制台和文件
                // 创建 logs 目录 (C++14 兼容方式)
                createDirectory("logs");
                
                // 创建文件日志接收器
                auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/enscdynamics_debug.log", true);
                file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%t] [%l] %v");
                
                // 组合接收器
                spdlog::sinks_init_list sink_list = { console_sink, file_sink };
                auto logger = std::make_shared<spdlog::logger>("EnSCDynamics", sink_list.begin(), sink_list.end());
                
                // 设置默认日志记录器
                spdlog::set_default_logger(logger);
                spdlog::set_level(spdlog::level::debug);
                spdlog::flush_on(spdlog::level::debug);
                spdlog::info("Logger 初始化 - Debug 模式 (显示 debug 及以上级别日志，同时输出到控制台和文件: logs/enscdynamics_debug.log)");

            #elif defined(CMAKE_BUILD_TYPE_RELWITHDEBINFO) || defined(RELWITHDEBINFO) || defined(RelWithDebInfo) || defined(_RELWITHDEBINFO) || defined(RelW) || defined(_RelW)
                // RelWithDebInfo 模式: 显示重要的调试信息，同时输出到控制台和文件
                // 创建 logs 目录 (C++14 兼容方式)
                createDirectory("logs");
                
                // 创建文件日志接收器
                auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/enscdynamics_relwithdebinfo.log", true);
                file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%t] [%l] %v");
                
                // 组合接收器
                spdlog::sinks_init_list sink_list = { console_sink, file_sink };
                auto logger = std::make_shared<spdlog::logger>("EnSCDynamics", sink_list.begin(), sink_list.end());
                
                // 设置默认日志记录器
                spdlog::set_default_logger(logger);
                spdlog::set_level(spdlog::level::debug);
                spdlog::flush_on(spdlog::level::debug);
                spdlog::info("Logger 初始化 - RelWithDebInfo 模式 (显示 debug 及以上级别日志，同时输出到控制台和文件: logs/enscdynamics_relwithdebinfo.log)");
                
            #else
                // Release 模式: 不显示任何日志信息，只使用控制台接收器
                auto logger = std::make_shared<spdlog::logger>("EnSCDynamics", console_sink);
                spdlog::set_default_logger(logger);
                spdlog::set_level(spdlog::level::off);
                spdlog::flush_on(spdlog::level::critical); // 只在关键错误时刷新
                // 初始化消息也不会显示，因为日志级别设置为 off
            #endif

        } catch (const spdlog::spdlog_ex& ex) {
            std::cerr << "Log initialization failed: " << ex.what() << std::endl;
        }
    }

    // Shuts down the logging system
    static void shutdown() {
        spdlog::info("Logger shutting down.");
        spdlog::shutdown(); // Release all spdlog resources
    }
};

} // namespace EnSC 