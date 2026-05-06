#' @section Methods:
#'
#'   \subsection{Diversity method}{From an entire collection of size \mjseqn{N},
#'   to construct a core set of sample size \mjseqn{n}, the number of entries to
#'   be selected from the \mjseqn{i}th group among \mjseqn{1 \cdots g} groups
#'   (\mjseqn{n_{i}}) is estimated as below.
#'
#'    \mjsdeqn{n_{i} = n \times \frac{D_{i}}{\sum_{i=1}^{g}D_{i}}}
#'
#'    Where, \mjseqn{D_{i}} is a measure of the extent of diversity present in
#'    the \mjseqn{i}th cluster.
#'
#'    }
#'
#'    \subsection{Diversity and proportional method}{Here the number of entries
#'    to be selected is proportional to the diversity of the cluster/group
#'    (\mjseqn{D_{i}}) weighted by the the cluster/group size (\mjseqn{N_{i}}).
#'
#'     \mjsdeqn{n_{i} = n \times \frac{N_{i}D_{i}}{\sum_{i=1}^{g}N_{i}D_{i}}}
#'
#'    }
#'
#'    \subsection{Diversity and logarithmic method}{Here the number of entries
#'    to be selected is proportional to the diversity of the cluster/group
#'    (\mjseqn{D_{i}}) weighted by the logarithm of the cluster/group size
#'    (\mjseqn{N_{i}}).
#'
#'    \mjsdeqn{n_{i} = n \times
#'    \frac{\log(N_{i})D_{i}}{\sum_{i=1}^{g}\log(N_{i})D_{i}}}
#'
#'    }
#'
#'    \subsection{Diversity and square root method}{Here the number of entries
#'    to be selected is proportional to the diversity of the cluster/group
#'    (\mjseqn{D_{i}}) weighted by the square root of the cluster/group size
#'    (\mjseqn{N_{i}}).
#'
#'    \mjsdeqn{n_{i} = n \times
#'    \frac{\sqrt{N_{i}}D_{i}}{\sum_{i=1}^{g}\sqrt{N_{i}}D_{i}}}
#'
#'    }
#'
