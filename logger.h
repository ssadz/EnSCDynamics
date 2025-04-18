#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <memory>
#include <string>
#include <vector>

namespace EnSC {

/**
 * @brief 简单的日志工具类，封装spdlog功能
 */
class Logger {
public:
    /**
     * @brief 初始化日志系统
     * @param logName 日志文件名前缀
     */
    static void init(const std::string& logName = "EnSCDynamics") {
        // 创建控制台输出接收器
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::info);
        
        std::vector<spdlog::sink_ptr> sinks{console_sink};
        
        // 在Debug和RelWithDebInfo模式下添加文件日志
#ifdef ENABLE_FILE_LOGGING
        try {
            auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logName + ".log", true);
            file_sink->set_level(spdlog::level::trace);
            sinks.push_back(file_sink);
        } catch (const spdlog::spdlog_ex& ex) {
            std::cerr << "日志文件无法创建: " << ex.what() << std::endl;
        }
#endif

        // 在Release模式下只记录错误日志
#ifdef ENABLE_ERROR_LOGGING_ONLY
        try {
            auto error_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logName + "_error.log", true);
            error_file_sink->set_level(spdlog::level::err);
            sinks.push_back(error_file_sink);
            
            // 在Release模式下降低控制台输出级别
            console_sink->set_level(spdlog::level::warn);
        } catch (const spdlog::spdlog_ex& ex) {
            std::cerr << "错误日志文件无法创建: " << ex.what() << std::endl;
        }
#endif

        // 创建并设置默认记录器
        auto logger = std::make_shared<spdlog::logger>(logName, sinks.begin(), sinks.end());
        spdlog::set_default_logger(logger);
        
        // 设置日志格式
        spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [thread %t] %v");
        
        // 设置刷新级别
        spdlog::flush_on(spdlog::level::warn);
    }
    
    /**
     * @brief 关闭日志系统
     */
    static void shutdown() {
        spdlog::shutdown();
    }
};

} // namespace EnSC 