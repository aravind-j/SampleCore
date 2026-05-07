#' Allocation of Accessions on the Basis of Cluster/Group Size
#'
#' Estimate the number of accessions to be allocated from each cluster/group in
#' the entire collection to construct a core collection on the basis of
#' cluster/group size. The following strategies are implemented. \loadmathjax
#' \itemize{
#' \item{Constant}
#' \item{Proportional}
#' \item{Logarithmic}
#' \item{Square root}}
#'
#' The different methods to determine the number of entries from each group or
#' clusters implemented in \code{allocate.basic} are as follows.
#'
#' These are different methods which estimate the number of entries only on the
#' basis of total number of accessions in each cluster/group.
#'
#' \insertCite{brown_core_1989;textual}{SampleCore} proposed the constant (C),
#' proportional (P) and logarithmic (L) methods and later a similar square root
#' method was proposed by
#' \insertCite{huaman_selecting_1999;textual}{SampleCore}.
#'
#'   \subsection{Constant method}{From an entire collection of size \mjseqn{N},
#'   to construct a core set of sample size \mjseqn{n}, the number of entries to
#'   be selected from the \mjseqn{i}th group among \mjseqn{1 \cdots g} groups
#'   (\mjseqn{n_{i}}) is estimated as below.
#'
#'   \mjsdeqn{n_{i} = \frac{n}{g} \times N}
#'
#'   }
#'
#'   \subsection{Proportional method}{Here the number of entries to be selected
#'   is proportional to the cluster/group size (\mjseqn{N_{i}}) as
#'   below.
#'
#'   \mjsdeqn{n_{i} = n \times \frac{N_{i}}{\sum_{i=1}^{g}N_{i}}}
#'
#'   \mjsdeqn{n_{i} = n \times \frac{N_{i}}{N}}
#'
#'   }
#'
#'   \subsection{Logarithmic method}{Here the number of entries to be selected
#'   is proportional to the logarithm of the cluster/group size
#'   (\mjseqn{N_{i}}) as below.
#'
#'   \mjsdeqn{n_{i} = n \times
#'   \frac{\log{(N_{i})}}{\sum_{i=1}^{g}\log{(N_{i})}}}
#'
#'   }
#'
#'   \subsection{Square root method}{Here the number of entries to
#'   be selected is proportional to the square root of the cluster/group size
#'   (\mjseqn{N_{i}}) as below.
#'
#'   \mjsdeqn{n_{i} = n \times \frac{\sqrt{N_{i}}}{\sum_{i=1}^{g}\sqrt{N_{i}}}}
#'
#'   }
#'
#' @template general-arg
#' @template log-arg
#' @param method The allocation method. Either \code{"const"} for constant or
#'   \code{"prop"} for proportional or \code{"log"} for logarithmic or
#'   \code{"sqrt"} for square root allocation.
#'
#' @returns
#'
#' @import mathjaxr
#' @importFrom Rdpack reprompt
#' @export
#'
#' @references
#'
#' \insertAllCited
#'
#' @examples
allocate.basic <- function(data,
                           names, group,
                           method = c("const", "prop", "log", "sqrt"),
                           log.base = exp(1),
                           size) {

  # Check stats ----

  checks.sample.core(data = data, size = size,
                     names = names, group = group,
                     dist.mat = NULL,
                     quantitative = NULL,
                     qualitative = NULL,
                     log.base = log.base,
                     mode = "alloc")

  method <- match.arg(method)

  # Prepare data ----

  data[, group] <- droplevels(data[, group])

  # Basic group stats ----

  gps <- levels(data[, group])
  gpsize <- summary(data[, group])
  gpcount <- length(gps)

  tcount <- nrow(data)

  # Constant allocation ----

  if (method == "const") {

    # Candidate sample sizes
    smpsize_c <- ceiling((size * tcount) / gpcount)
    smpsize_f <- floor((size * tcount) / gpcount)

    # Choose value closer to target total sample size
    if (abs((smpsize_c * gpcount) - (size * tcount)) >=
        abs((smpsize_f * gpcount) - (size * tcount))) {
      smpsize <- smpsize_f
    } else {
      smpsize <- smpsize_c
    }

    # Assign sample size to groups
    out <- rep(smpsize, length(gps))
    names(out) <- names(gpsize)

    # Check if sample size exceeds the group size
    # ie. group is too small to support the required sample size
    if (any(gpsize < smpsize)) {
      # Problem groups
      sml_gpsize <- gpsize[gpsize < smpsize]

      sml_gpsize <- paste(names(sml_gpsize), " (=", sml_gpsize,")",
                          sep = "", collapse = ", ")

      warning('For the following groups in "', group,
              '" the sampling count estimated with "size" = ', size,
              " (", smpsize, ")", " exceeds the group count.\n", sml_gpsize,
              "\nRestricting the sampling count to group counts",
              " for these groups.")

      out <- ifelse(gpsize < smpsize, gpsize, smpsize)
    }

  }

  if (method != "const") {

  # Freq estimation ----

  ## "prop" ----
  if (method == "prop") {
    freq <- gpsize / sum(gpsize)
  }

  ## "log" ----
  if (method == "log") {
    log_gpsize <- log(gpsize, base = log.base)
    freq <- log_gpsize / sum(log_gpsize)
  }

  ## "sqrt" ----
  if (method == "sqrt") {
    freq <- sqrt(gpsize) / sum(sqrt(gpsize))
  }

  # Get output ---
    out <- setNames(round(freq * size * tcount), names(gpsize))
  }

  return(out)

}
