init_log_file <- function(log_file) {
  cat("", file = log_file)
}

init_warning_file <- function(warnings_file) {
  cat("", file = warnings_file)
}

log_info <- function(msg, log_file) {
  line <- sprintf("[%s] INFO  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

log_warn <- function(msg, log_file) {
  line <- sprintf("[%s] WARN  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", file = log_file, append = TRUE)
}

log_error <- function(msg, log_file) {
  line <- sprintf("[%s] ERROR %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", file = log_file, append = TRUE)
}