
#' Diversity Indices
#'
#' Compute the following diversity indices. \itemize{ \item{Simpson's and
#' related indices} \itemize{ \item{Simpson's Index (\mjseqn{d})
#' \insertCite{simpson_measurement_1949,peet_measurement_1974}{SampleCore}}
#' \item{Simpson's Index of Diversity or Gini's Diversity Index or Gini-Simpson
#' Index or Nei's Diversity Index or Nei's Variation Index (\mjseqn{D})
#' \insertCite{gini_variabilita_1912,gini_variabilita_1912-2,greenberg_measurement_1956,berger_diversity_1970,nei_analysis_1973,peet_measurement_1974}{SampleCore}}
#' \item{Maximum Simpson's Index of Diversity or Maximum Nei's
#' Diversity/Variation Index (\mjseqn{D_{max}})
#' \insertCite{hennink_interpretation_1990}{SampleCore}} \item{Simpson's
#' Reciprocal Index or Hill's \mjseqn{N_{2}} (\mjseqn{D_{R}})
#' \insertCite{williams_patterns_1964,hill_diversity_1973}{SampleCore}}
#' \item{Relative Simpson's Index of Diversity or Relative Nei's
#' Diversity/Variation Index (\mjseqn{D'})
#' \insertCite{hennink_interpretation_1990}{SampleCore}} }
#' \item{Shannon-Weaver and related indices} \itemize{ \item{Shannon or
#' Shannon-Weaver or Shannon-Weiner Diversity Index (\mjseqn{H})
#' \insertCite{shannon_mathematical_1949,peet_measurement_1974}{SampleCore}}
#' \item{Maximum Shannon-Weaver Diversity Index (\mjseqn{H_{max}})
#' \insertCite{hennink_interpretation_1990}{SampleCore}} \item{Relative
#' Shannon-Weaver Diversity Index or Shannon Equitability Index (\mjseqn{H'})
#' \insertCite{hennink_interpretation_1990}{SampleCore}} } \item{McIntosh
#' Diversity Index} \itemize{ \item{McIntosh Diversity Index (\mjseqn{D_{Mc}})
#' \insertCite{mcintosh_index_1967,peet_measurement_1974}{SampleCore}} } }
#' \loadmathjax
#'
#' @section Details: The diversity indices and the corresponding statistical
#'   tests implemented in \code{diversity.indices} are as follows.
#'
#'   \subsection{Simpson's and related indices}{Simpson's index (\mjseqn{d})
#'   which estimates the probability that two accessions randomly selected will
#'   belong to the same class of a qualitative trait, is computed as follows
#'   \insertCite{simpson_measurement_1949,peet_measurement_1974}{SampleCore}.
#'
#'   \mjsdeqn{d = \sum_{i = 1}^{k}p_{i}^{2}}
#'
#'   Where, \mjseqn{p_{i}} denotes the proportion/fraction/frequency of
#'   accessions in the \mjseqn{i}th class of a qualitative trait and \mjseqn{k}
#'   is the number of classes for the qualitative trait.
#'
#'   The value of \mjseqn{d} can range from 0 to 1 with 0 representing maximum
#'   diversity and 1, no diversity.
#'
#'   The complement of \mjseqn{d} (\mjseqn{d} is subtracted from 1) is called
#'   the Simpson's index of diversity (\mjseqn{D})
#'   \insertCite{greenberg_measurement_1956,berger_diversity_1970,peet_measurement_1974,hennink_interpretation_1990}{SampleCore}
#'    originally suggested by
#'   \insertCite{gini_variabilita_1912,gini_variabilita_1912-2;textual}{SampleCore}
#'    and described in literature as Gini's diversity index or Gini-Simpson
#'   index. It is the same as Nei's diversity index or Nei's variation index
#'   \insertCite{nei_analysis_1973,hennink_interpretation_1990}{SampleCore}.
#'   Greater the value of \mjseqn{D}, greater the diversity with a range from 0
#'   to 1.
#'
#'   \mjsdeqn{D = 1 - d}
#'
#'   The maximum value of \mjseqn{D}, \mjseqn{D_{max}} occurs when accessions
#'   are uniformly distributed across the classes in the qualitative trait and
#'   is computed as follows
#'   \insertCite{hennink_interpretation_1990}{SampleCore}.
#'
#'   \mjsdeqn{D_{max} = 1 - \frac{1}{k}}
#'
#'   Reciprocal of \mjseqn{d} gives the Simpson's reciprocal index
#'   (\mjseqn{D_{R}})
#'   \insertCite{williams_patterns_1964,hennink_interpretation_1990}{SampleCore}
#'    and can range from 1 to \mjseqn{k}. This was also described in
#'   \insertCite{hill_diversity_1973;textual}{SampleCore} as (\mjseqn{N_{2}}).
#'
#'   \mjsdeqn{D_{R} = \frac{1}{d}}
#'
#'   Relative Simpson's index of diversity or Relative Nei's diversity/variation
#'   index (\mjseqn{H'}) \insertCite{hennink_interpretation_1990}{SampleCore}
#'   is defined as follows \insertCite{peet_measurement_1974}{SampleCore}.
#'
#'   \mjsdeqn{D' = \frac{D}{D_{max}}}
#'
#'   }
#'
#'   \subsection{Shannon-Weaver and related indices}{An index of information
#'   \mjseqn{H}, was described by
#'   \insertCite{shannon_mathematical_1949;textual}{SampleCore} as follows.
#'
#'   \mjsdeqn{H = -\sum_{i=1}^{k}p_{i} \log_{2}(p_{i})}
#'
#'   \mjseqn{H} is described as Shannon or Shannon-Weaver or Shannon-Weiner
#'   diversity index in literature.
#'
#'   Alternatively, \mjseqn{H} is also computed using natural logarithm instead
#'   of logarithm to base 2.
#'
#'   \mjsdeqn{H = -\sum_{i=1}^{k}p_{i} \ln(p_{i})}
#'
#'   The maximum value of \mjseqn{H} (\mjseqn{H_{max}}) is \mjseqn{\ln(k)}. This
#'   value occurs when each class for a qualitative trait has the same
#'   proportion of accessions.
#'
#'   \mjsdeqn{H_{max} = \log_{2}(k)\;\; \textrm{OR} \;\; H_{max} = \ln(k)}
#'
#'   The relative Shannon-Weaver diversity index or Shannon equitability index
#'   (\mjseqn{H'}) is the Shannon diversity index (\mjseqn{I}) divided by the
#'   maximum diversity (\mjseqn{H_{max}}).
#'
#'   \mjsdeqn{H' = \frac{H}{H_{max}}}
#'
#'   }
#'
#'   \subsection{McIntosh Diversity Index}{A similar index of diversity was
#'   described by \insertCite{mcintosh_index_1967;textual}{SampleCore} as
#'   follows (\mjseqn{D_{Mc}}) \insertCite{peet_measurement_1974}{SampleCore}.
#'
#'   \mjsdeqn{D_{Mc} = \frac{N - \sqrt{\sum_{i=1}^{k}n_{i}^2}}{N - \sqrt{N}}}
#'
#'   Where, \mjseqn{n_{i}} denotes the number of accessions in the \mjseqn{i}th
#'   class for a qualitative trait and \mjseqn{N} is the total number of
#'   accessions so that \mjseqn{p_{i} = {n_{i}}/{N}}.}
#'
#' @param x The qualitative trait data as a vector of class factor.
#' @param base The logarithm base to be used for computation of Shannon-Weaver
#'   Diversity Index (\mjseqn{I}). Default is 2.
#'
#' @return The diversity index value.
#'
#' @export
#'
#' @references
#'
#' \insertAllCited{}
#'
diversity.indices <- function(x, index = c("simpson", "simpson.cmp",
                                           "simpson.max", "simpson.inv",
                                           "simpson.rel",
                                           "shannon", "shannon.max",
                                           "shannon.rel",
                                           "mcintosh"),
                              base = 2) {

  index <- match.arg(index)

  x <- droplevels(x)

  count <- as.vector(table(x))
  total.count <- sum(count, na.rm = TRUE)

  k <- length(count)

  prob <- count / total.count

  # Simpson's Index (d) ----
  if (index == "simpson") {
    d <- sum(prob ^ 2)
    out <- d
  }

  # Simpson's Index of Diversity ----
  if (index == "simpson.cmp") {
    D <- 1 - d
    out <- D
  }

  # Max Simpson's index of diversity ----
  if (index == "simpson.max") {
    D.max <- 1 - (1 / k)
    out <- D.max
  }

  # Reciprocal Simpson's Index ----
  if (index == "simpson.inv") {
    D.inv <- 1 / d
    out <- D.inv
  }

  # Relative Simpson's Index ----
  if (index == "simpson.rel") {
    D.rel <- D / D.max
    out <- D.rel
  }

  #-----------------------------------------------------------------------------

  # Shannon-Weaver Diversity Index (H) ----
  if (index == "shannon") {
    I <- - sum(prob * log(prob, base = base))
    out <- I
  }

  # Maximum Shannon-Weaver Diversity Index ----
  if (index == "shannon.max") {
    I.max <- log(k)
    out <- I.max
  }

  # Relative Shannon-Weaver Diversity Index  ----
  if (index == "shannon.rel") {
    I.rel <- I / I.max
    out <- I.rel
  }

  #-----------------------------------------------------------------------------

  # McIntosh Index ----
  if (index == "mcintosh") {
    D.Mc <- (total.count - sqrt(sum(count ^ 2))) /
      (total.count - sqrt(total.count))
    out <- D.Mc
  }

  #-----------------------------------------------------------------------------

  name(out) <- index
  return(out)

}
