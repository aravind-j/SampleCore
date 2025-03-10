% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.size.R
\name{sample.size}
\alias{sample.size}
\title{Sample Size Estimation}
\usage{
sample.size(
  data,
  names,
  group,
  dist,
  quantitative,
  qualitative,
  size = 0.2,
  sample.method = c("const", "prop", "log", "sqrt", "diversity", "diverstiy.prop",
    "diversity.sqrt", "diversity.log")
)
}
\arguments{
\item{data}{The data as a data frame object. The data frame should possess
one row per individual and columns with the individual names and multiple
trait/character data.}

\item{names}{Name of column with the accession names as a character string.}

\item{group}{Name of column with the accession group/cluster names as a
character string.}

\item{dist}{A precomputed distance matrix of distance measures between the
accessions in \code{data}.}

\item{size}{The desired core set size proportion.}

\item{sample.method}{}
}
\description{
Estimate the number of accessions to be sampled from each group/cluster in
the entire collection to construct a core collection. The following
strategies are implemented. \loadmathjax
\itemize{ \item{Basic methods \insertCite{brown_Core_1989,huaman_Selecting_1999}{SampleCore}}{
\itemize{
\item{Constant}
\item{Proportional}
\item{Logarithmic}
\item{Square root}}}
\item{Diversity dependent methods \insertCite{yonezawa_Sampling_1995,schoen_Conservation_1993,bisht_Assessment_1999,mahajan_Sampling_1999,franco_sampling_2005}{SampleCore}}{ \itemize{
\item{Diversity}
\item{Diversity & Proportional}
\item{Diversity & Logarithmic}
\item{Diversity & Square root}}}
}
}
\section{Details}{
 The different methods to determine the number of entries
  from each group or clusters implemented in \code{sample.size} are as
  follows.

  \subsection{Basic methods}{These are different methods which estimate the
  number of entries only on the basis of total number of accessions in
  each group/cluster. \insertCite{brown_Core_1989;textual}{SampleCore}
  proposed the constant (C), proportional (P) and logarithmic (L) methods
  and later a similar square root method was proposed by
  \insertCite{huaman_Selecting_1999;textual}{SampleCore}.

  \subsection{Constant method}{From an entire collection of size \mjseqn{N},
  to construct a core set of sample size \mjseqn{n}, the number of entries to
  be selected from the \mjseqn{i}th group among \mjseqn{1 \cdots g} groups
  (\mjseqn{n_{i}}) is estimated as below.

  \mjsdeqn{n_{i} = \frac{n}{g} \times N}

  }

  \subsection{Proportional method}{Here the number of entries to be selected
  is proportional to the group/cluster size (\mjseqn{N_{i}}) as
  below.

  \mjsdeqn{n_{i} = n \times \frac{N_{i}}{\sum_{i=1}^{g}N_{i}}}

  \mjsdeqn{n_{i} = n \times \frac{N_{i}}{N}}

  }

  \subsection{Logarithmic method}{Here the number of entries to be selected
  is proportional to the logarithm of the group/cluster size
  (\mjseqn{N_{i}}) as below.

  \mjsdeqn{n_{i} = n \times
  \frac{\log{(N_{i})}}{\sum_{i=1}^{g}\log{(N_{i})}}}

  }

  \subsection{Square root method}{Here the number of entries to
  be selected is proportional to the square root of the group/cluster size
  (\mjseqn{N_{i}}) as below.

  \mjsdeqn{n_{i} = n \times \frac{\sqrt{N_{i}}}{\sum_{i=1}^{g}\sqrt{N_{i}}}}

  }

  }

  \subsection{Diversity dependent methods}{These are different methods which
  estimate the number of entries on the basis of how diverse are the
  accessions in each group/cluster. There are several methods proposed on the
  basis of diversity indices such as genetic multiplicity (G) dependent
  method based on the range of genetic diversity
  \insertCite{yonezawa_Sampling_1995}{SampleCore}, H strategy based on Nei's
  gene diversity \insertCite{nei_analysis_1973}{SampleCore} and a method
  based on the pooled Shannon diversity index
  \insertCite{bisht_Assessment_1999,mahajan_Sampling_1999}{SampleCore}.
  Similarly, measures such as expected proportion of heterozygous loci per
  individual and effective number of alleles have also been employed as a
  diversity measure for determining sample size
  \insertCite{franco_sampling_2006}{SampleCore}.
  \insertCite{franco_sampling_2005;textual}{SampleCore} proposed a method
  based on mean Gower's distance
  \insertCite{gower_general_1971}{SampleCore} which was also extended to
  other distance measure averages named D Allocation strategy
  \insertCite{franco_sampling_2006}{SampleCore}. These methods were also
  combined with the proportional and logarithmic methods. For example, the GP
  and GL strategy of \insertCite{bisht_Assessment_1999;textual}{SampleCore}
  and \insertCite{mahajan_Sampling_1999;textual}{SampleCore} as well as the
  NY and LD allocation methods of
  \insertCite{franco_sampling_2005;textual}{SampleCore}.

  \subsection{Diversity method}{From an entire collection of size \mjseqn{N},
  to construct a core set of sample size \mjseqn{n}, the number of entries to
  be selected from the \mjseqn{i}th group among \mjseqn{1 \cdots g} groups
  (\mjseqn{n_{i}}) is estimated as below.

   \mjsdeqn{n_{i} = n \times \frac{D_{i}}{\sum_{i=1}^{g}D_{i}}}

   Where, \mjseqn{D_{i}} is a measure of the extent of diversity present in
   the \mjseqn{i}th cluster.

   \mjseqn{D} can be either
   1) Range of a diversity index 2) Pooled value of a diversity index or 3)
   Mean genetic distance.

   }

    \subsection{Diversity and proportional method}{Here the number of entries
   to be selected is proportional to the diversity of the group/cluster
   (\mjseqn{D_{i}}) weighted by the the group/cluster size (\mjseqn{N_{i}}).

    \mjsdeqn{n_{i} = n \times \frac{N_{i}D_{i}}{\sum_{i=1}^{g}N_{i}D_{i}}}

   }

   \subsection{Diversity and logarithmic method}{Here the number of entries
   to be selected is proportional to the diversity of the group/cluster
   (\mjseqn{D_{i}}) weighted by the logarithm of the group/cluster size
   (\mjseqn{N_{i}}).

   \mjsdeqn{n_{i} = n \times
   \frac{\log(N_{i})D_{i}}{\sum_{i=1}^{g}\log(N_{i})D_{i}}}

   }

   \subsection{Diversity and square root method}{Here the number of entries
   to be selected is proportional to the diversity of the group/cluster
   (\mjseqn{D_{i}}) weighted by the square root of the group/cluster size
   (\mjseqn{N_{i}}).

   \mjsdeqn{n_{i} = n \times
   \frac{\sqrt{N_{i}}D_{i}}{\sum_{i=1}^{g}\sqrt{N_{i}}D_{i}}}

   }

  }
}

\references{
\insertAllCited
}
