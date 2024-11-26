function set_loglevel(loglevel)
    if lowercase(loglevel) == "debug"
        logger = ConsoleLogger(stderr, Logging.Debug)
    elseif lowercase(loglevel) == "warning"
        logger = ConsoleLogger(stderr, Logging.Warn)
    else
        logger = ConsoleLogger(stderr, Logging.Info)
    end
    global_logger(logger)
end
