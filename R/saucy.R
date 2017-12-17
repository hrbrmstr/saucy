#' Searching for Automorphisms in Underlying CNF, yes?
#'
#' @md
#' @param path file with graph
#' @param mode using `""` auto-chooses
#' @param timeout timeout (Default: 0)
#' @param rpt repeat (mostly for benchmarking)
#' @references [source](https://github.com/andandandand/saucy-repo/blob/master/saucy.c)
#' @export
#' @examples
#' saucy::saucy(system.file("extdata", "graphfile", package="saucy"))
#' saucy::saucy(system.file("extdata", "graphfile2", package="saucy"))
saucy <- function(path, mode = c("", "cnf", "digraph", "shatter"),
                  timeout = 0, rpt = 1) {

  mode <- match.arg(trimws(tolower(mode)), c("", "cnf", "digraph", "shatter"))

  path <- normalizePath(path.expand(path))
  if (!file.exists(path)) stop("File not found", call.=FALSE)

  tf <- tempfile()
  tfc <- file(tf, open="wt")
  sink(tfc)
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
