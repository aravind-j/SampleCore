# check if all levels are present in factors


#' Common checks for all functions
#'
#' Not exported. Strictly internal
#'
#' @keywords internal
#'
#' @template general-arg
#' @template qualquant-arg
#
checks.evaluate.core <- function(data, names, quantitative = NULL,
                                 qualitative = NULL,
                                 method, size, group
                                 # selected
                                 ) {

  # Declare nulls ----

  if (missing(quantitative)) {
    quantitative <- NULL
  }

  if (missing(qualitative)) {
    qualitative <- NULL
  }

  if (length(c(quantitative, qualitative)) == 1) {
    stop("Only one trait specified.")
  }

  # data ----

  # check if 'data' is a data frame object
  if (!is.data.frame(data)) {
    stop('"data" should be a data frame object.')
  }

  if (any(c("tbl_dataf", "tbl") %in% class(data))) {
    warning('"data" is of type tibble\nCoercing to data frame.')
    data <- as.data.frame(data)
  }

  # Argument lengths and type ----

  # check if 'names' argument is character vector of unit length
  if (!(is.character(names) && length(names) == 1)) {
    stop('"names" should be a character vector of unit length.')
  }

  # check if 'group' argument is character vector of unit length
  if (!(is.character(group) && length(group) == 1)) {
    stop('"group" should be a character vector of unit length.')
  }

  # check if 'quantitative' argument is a character vector
  if (!is.null(quantitative)) {
    if (!is.character(quantitative)) {
      stop('"quantitative" should be a character vector.')
    }
  }

  # check if 'qualitative' argument is a character vector
  if (!is.null(qualitative)) {
    if (!is.character(qualitative)) {
      stop('"qualitative" should be a character vector.')
    }
  }

  # check if 'group' column is of type factor
  if (!is.factor(data[, group])) {
    stop('"group" column in "data" should be of type factor.')
  }

  # check if 'size' argument is numeric vector of unit length
  if (!is.null(size)) {
    if (!(is.numeric(size) && length(size) == 1)) {
      stop('"size" should be a numeric vector of unit length.')
    }
  }

  # Columns exist in data ----

  # check if 'names' column is present in 'data'
  if (!(names %in% colnames(data))) {
    stop(paste('Column ', names,
               ' specified as the "names" column is not present in "data".',
               sep = ""))
  }

  # check if 'quantitative' columns are present in 'data'
  if (!is.null(quantitative)) {
    if (FALSE %in% (quantitative %in% colnames(data)))  {
      stop(paste('The following column(s) specified in "quantitative" ',
                 'not present in "data":\n',
                 paste(quantitative[!(quantitative %in% colnames(data))],
                       collapse = ", "),
                 sep = ""))
    }
  }

  # check if 'qualitative' columns are present in 'data'
  if (!is.null(qualitative)) {
    if (FALSE %in% (qualitative %in% colnames(data)))  {
      stop(paste('The following column(s) specified in "qualitative" ',
                 'not present in "data":\n',
                 paste(qualitative[!(qualitative %in% colnames(data))],
                       collapse = ", "),
                 sep = ""))
    }
  }

  # check if overlap exists between 'quantitative' and 'qualitative'
  if ((!is.null(quantitative)) && (!is.null(qualitative))) {
    if (length(intersect(quantitative, qualitative)) != 0) {
      stop(paste('The following column(s) is/are specified in both ',
                 '"quantitative" and "qualitative":\n',
                 paste(intersect(quantitative, qualitative),
                       collapse = ", "),
                 sep = ""))
    }
  }

  # check if 'group' column is present in 'data'
  if (!(group %in% colgroup(data))) {
    stop(paste('Column ', group,
               ' specified as the "group" column is not present in "data".',
               sep = ""))
  }

  # Check column types ----

  # check if 'names' column is of type character
  if (!is.character(data[, names])) {
    stop('"names" column in "data" should be of type character.')
  }

  # check if 'quantitative' columns are of type numeric/integer
  if (!is.null(quantitative)) {
    intquantcols <-
      unlist(lapply(data[, quantitative],
                    function(x) FALSE %in% (is.vector(x, mode = "integer") |
                                              is.vector(x, mode = "numeric"))))
    if (TRUE %in% intquantcols) {
      stop(paste('The following "quantitative" column(s) in "data" are not ',
                 'of type numeric:\n',
                 paste(names(intquantcols[intquantcols]), collapse = ", ")))
    }
  }

  # check if 'qualitative' columns are of type factor
  if (!is.null(qualitative)) {
    intqualcols <- unlist(lapply(data[, qualitative],
                                 function(x) is.factor(x)))
    if (FALSE %in% intqualcols) {
      stop(paste('The following "qualitative" column(s) in "data" are not ',
                 'of type factor:\n',
                 paste(names(intqualcols[!intqualcols]), collapse = ", ")))
    }
  }

  # Missing values ----

  # check for missing values
  missvcols <- unlist(lapply(data[, quantitative],
                             function(x) TRUE %in% is.na(x)))
  if (TRUE %in% missvcols) {
    warning(paste('The following column(s) in "data" have missing values:\n',
                  paste(names(missvcols[missvcols]), collapse = ", ")))
  }

  # Duplication ----

  # check for duplication in names
  if (any(duplicated(data[, names]))) {
    stop('Duplicated entries exist in "names" column.')
  }

  # Size ----

  # check if 'size' is a proportion between 0 and 1
  if (size <= 0 || size >= 1) {
    stop('"size" should be a proportion between 0 and 1.')
  }
}
