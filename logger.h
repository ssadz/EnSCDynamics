#pragma once

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h" // For console logging

namespace EnSC {

class Logger {
public:
    // Initializes the logging system (e.g., creates a console logger)
    static void init() {
        try {
            // Create a color console logger
            // %^ marks the start of color range, %$ marks the end
            // [%Y-%m-%d %H:%M:%S.%e] [%n] [%l] %v%$ - Timestamp, logger name, level, message
            auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            console_sink->set_pattern("[%^%Y-%m-%d %H:%M:%S.%e%$] [%n] [%l] %v");

            // You can add more sinks here (e.g., file sink)
            // auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/enscdynamics.log", true);
            // file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%t] [%l] %v");

            // Combine sinks if needed
            // spdlog::sinks_init_list sink_list = { console_sink, file_sink };
            // auto logger = std::make_shared<spdlog::logger>("EnSCDynamics", sink_list.begin(), sink_list.end());

            // For simplicity, just use the console sink for now
            auto logger = std::make_shared<spdlog::logger>("EnSCDynamics", console_sink);

            // Set the default logger globally
            spdlog::set_default_logger(logger);
            
            // 根据不同构建类型设置不同的日志级别
#if defined(_DEBUG) || defined(DEBUG)
            // Debug 模式: 显示所有调试信息
            spdlog::set_level(spdlog::level::debug);
            spdlog::flush_on(spdlog::level::debug);
            spdlog::info("Logger 初始化 - Debug 模式 (显示 debug 及以上级别日志)");
#elif defined(CMAKE_BUILD_TYPE_RELWITHDEBINFO) || defined(RELWITHDEBINFO)
            // RelWithDebInfo 模式: 显示重要的调试信息
            spdlog::set_level(spdlog::level::debug);
            spdlog::flush_on(spdlog::level::debug);
            spdlog::info("Logger 初始化 - RelWithDebInfo 模式 (显示 debug 及以上级别日志)");
#else
            // Release 模式: 不显示任何日志信息
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