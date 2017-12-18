#' Identifies Symmetries in CNF Instances
#'
#' @md
#' @param path file with graph
#' @references [Shatter](http://www.aloul.net/Tools/shatter/)
#' @export
#' @examples
#' saucy::shatter(system.file("extdata", "example.cnf", package="saucy"))
shatter <- function(path) {

  path <- normalizePath(path.expand(path))
  if (!file.exists(path)) stop("File not found", call.=FALSE)

  sbp_file <- tempfile()
  g_file <- tempfile()

  tf <- tempfile()
  tfc <- file(tf, open="wt")
  sink(tfc)
  ret <- shatter_int(path, sbp_file, g_file)
  sink()

  x <- readLines(tf, warn = FALSE)
  s <- readLines(sbp_file, warn = FALSE)
  g <- readLines(g_file, warn = FALSE)

  close(tfc)
  unlink(tf)
  unlink(sbp_file)
  unlink(g_file)

  as.list(unlist(
    lapply(
      lapply(
        lapply(g, strsplit, " = "),
        unlist
      ),
      function(x) setNames(as.numeric(x[2]), gsub("[- ]+", "_", x[1]))
    )
  )) -> gl

  ret <- append(ret, gl)
  ret$symmetry_breaking_predicates <- s
  ret$printed_output <- x

  class(ret) <- c("shatter", "list")

  ret

}

#' @param x object
#' @param ... unused
#' @rdname saucy
#' @export
print.shatter <- function(x, ...) {

  orig_x <- x

  x$input_file <- basename(x$input_file)
  sbp <- x$symmetry_breaking_predicates

  x$symmetry_breaking_predicates <- sprintf("%d symmetry breaking predicates", length(sbp))

  cat(sprintf("%29s: %s", names(x), unname(x)), sep="\n")

  invisible(orig_x)

}
