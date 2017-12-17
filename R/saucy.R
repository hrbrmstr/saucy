#' Saucy
#'
#' @param path path
#' @param mode mode
#' @param timeout timeout
#' @param rpt repeat
#' @export
#' @examples
#' saucy::saucy(system.file("extdata", "graphfile", package="saucy"))
saucy <- function(path, mode = c("", "cnf", "digraph", "shatter"),
                  timeout = 0, rpt = 1) {

  mode <- match.arg(trimws(tolower(mode)), c("", "cnf", "digraph", "shatter"))

  path <- normalizePath(path.expand(path))
  if (!file.exists(path)) stop("File not found", call.=FALSE)

  tf <- tempfile()
  tfc <- file(tf, open="wt")
  sink(tfc)
  sink(tfc, type="message")
  ret <- saucy_int(path, mode, timeout, rpt)
  sink()

  x <- readLines(tf)
  unlink(tf)

  ret$printed_output <- x

  class(ret) <- c("saucy", "list")

  ret

}

#' @param x object
#' @param ... unused
#' @rdname saucy
#' @export
print.saucy <- function(x, ...) {

  orig_x <- x

  cat(x$printed_output, sep="\n");
  cat("\n")

  x$printed_output <- NULL
  x$input_file <- basename(x$input_file)

  cat(sprintf("%20s: %s", names(x), unname(x)), sep="\n")

  invisible(orig_x)

}
