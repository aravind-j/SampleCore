# check if all levels are present in factors


#' Common checks for all functions
#'
#' Not exported. Strictly internal
#'
#' @keywords internal
#'
#' @template general-arg
#' @template qualquant-arg
#' @template dist-arg
#' @template log-arg
#' @template sel-arg
#' @template size-arg
#
checks.sample.core <- function(data, names,
                               size, group,
                               quantitative = NULL,
                               qualitative = NULL,
                               dist.mat = NULL,
                               log.base = NULL,
                               alloc,
                               always.selected = NULL,
                               mode = c("alloc", "sel")) {

  # Declare nulls ----

  if (missing(quantitative)) {
    quantitative <- NULL
  }

  if (missing(qualitative)) {
    qualitative <- NULL
  }

  if (missing(dist.mat)) {
    dist.mat <- NULL
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
  if (mode == "alloc") {
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
  if (!(group %in% colnames(data))) {
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

  # Distance matrix ----

  if (!is.null(dist.mat)) {

    # check if 'dist.mat' is a distance matrix
    if (!inherits(dist.mat, "dist")) {
      stop('"dist.mat" should be a distance matrix of class "dist".')
    }
    dist_labels <- labels(dist.mat) # No dups as it is a distance matrix

    if (is.null(dist_labels)) {
      stop('Labels are missing in "dist.mat".')
    }

    # Check for mismatch between elements of 'data' and 'dist.mat'
    if (!setequal(data[, names], labels(dist.mat))) {
      only_in_data <- setdiff(data[, names], labels(dist.mat))
      only_in_dist <- setdiff(labels(dist.mat), data[, names])

      stop(
        paste0(
          'Mismatch between "dist.mat" labels and entries in "names" column.\n',
          'Only in "data": ', paste(only_in_data, collapse = ', '), '\n',
          'Only in "dist.mat": ', paste(only_in_dist, collapse = ', ')
        )
      )
    }
  }

  # Size ----

  # check if 'size' is a proportion between 0 and 1
  if (mode == "alloc") {
    if (size <= 0 || size >= 1) {
      stop('"size" should be a proportion between 0 and 1.')
    }
  }

  # Log base ----

  if (!is.null(log.base)) {

    if (!is.numeric(log.base) || length(log.base) != 1 || is.na(log.base)) {
      stop('"log.base" must be a single numeric value.')
    }

    if (log.base <= 0) {
      stop('"base" must be positive.')
    }

    if (log.base == 1) {
      stop('"log.base" of 1 is undefined.')
    }

  }

  # Always selected ----

  if (mode == "sel" & !is.null(always.selected)) {
    # check if 'always.selected' is a character vector
    if (!is.character(always.selected)) {
      stop('"always.selected" should be a character vector.')
    }

    # check if always.selected is present in the entire set
    if (any(!(always.selected %in% data[, names]))) {
      alsel_miss <- always.selected[!(always.selected %in% data[, names])]
      stop(paste('The following entry/entries specified in "always.selected" ',
                 'are not present in "data":\n',
                 paste(alsel_miss, collapse = ", "),
                 sep = ""))
    }
  }

  # Allocation vector ----

  if (mode == "sel") {

    nm <- names(alloc)
    lv <- levels(data[, group])

    if (!setequal(nm, lv)) {

      only_in_alloc <- setdiff(nm, lv)
      only_in_group <- setdiff(lv, nm)

      msg <- paste0(
        'Mismatch between "alloc" names and levels of "group" column in data.\n\n',
        if (length(only_in_alloc) > 0)
          paste0(
            'Present in "alloc" but not in "group" column levels: ',
            paste(shQuote(only_in_alloc), collapse = ', '),
            '\n'
          ) else '',
        if (length(only_in_group) > 0)
          paste0(
            'Present in "group" column levels but not in "alloc": ',
            paste(shQuote(only_in_group), collapse = ', '),
            '\n'
          ) else ''
      )

      stop(msg, call. = FALSE)

    }
  }

}




