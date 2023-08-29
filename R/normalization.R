#' @import data.table

# processing --------------------------------------------------------------

# hybridization control elutions
.hce <- c('seq.2171.12', 'seq.2178.55', 'seq.2194.91', 'seq.2229.54',
  'seq.2249.25', 'seq.2273.34', 'seq.2288.7', 'seq.2305.52', 'seq.2312.13',
  'seq.2359.65', 'seq.2430.52', 'seq.2513.7')

# functions, helper -------------------------------------------------------

# @export
derive_intraplate_ref <- \(ds, v, ref_name) {
  ds[SampleType == 'Calibrator', (ref_name) := median(get(v)), .(PlateId, variable)]
}

fill_column <- function(x){
  res <- x %>% unique() %>% na.omit()
  if (length(res) != 1) {
    stop('fill column issue: ', length(res), paste(res, collapse = ' '))
  }
  return(res)
}

calc_ratios <- function(ds, value_name, ref_name) {
  # processing
  ds[, .(ratio = unique(get(ref_name)) / median(get(value_name))), variable] %>%
    .$ratio
}

# functions, pipeline ---------------------------------------------------------------

#' Normalize data using pipeline
#'
#' @param ds Dataset as data.table
#' @param value_name String of variable name containing value to be normalized
#' @param steps List of functions that will be used to calculate scale factors
#' @param step_names Names of each step provided by steps parameter
#' @param refs List containing a string or function that provides references to be used by each step
#' @param extra_params Optional parameter to provide extra parameters to functions given by steps parameter
#' @examples
#'
#' normalization_pipeline(
#'   ana.raw,
#'   value,
#'   steps = list(hyb_norm, intraplate_norm, plate_norm, calibration_norm),
#'   step_names = c('H', 'HI', 'HIP', 'HIPC'),
#'   refs = list('ref.all', derive_intraplate_ref, 'ref.calibrator', 'ref.calibrator')
#' )
#'
#' @export

normalization_pipeline <- function(
    ds,
    value_name,
    steps,
    step_names,
    refs,
    extra_params = list()
    ) {
  # derive
  n_steps <- step_names %>% length()
  exts <- paste0('_', step_names)
  prev_step <- exts %>%
    head(n_steps - 1) %>%
    c('', .)

  # main loop
  for(i in seq_along(exts)) {
    # get entries
    ext <- exts[[i]]
    prv <- prev_step[[i]]
    f <- steps[[i]]
    ref <- refs[[i]]
    nm <- step_names[[i]]

    v <- paste0(value_name, prv)

    # handle reference
    if (class(ref) == 'function') {
      ref_name <- paste0('ref', ext)
      ref(ds, v, ref_name)
    } else if (class(ref) == 'character') {
      ref_name <- ref
    } else {
      stop('invalid class of reference')
    }

    # get extra params
    parms <- list(ds, v, ref_name, ext) %>%
      c(extra_params[[nm]])

    # calcuate SF
    do.call('f', parms)

    # calculate value
    ds[, paste0(value_name, ext) := get(v) * get(paste0('sf', ext))]

    # return
    NULL
  }
  return(NULL)
}


# function, normalization steps -------------------------------------------

#' Hybridization normalization
#'
#' @param ds data.table
#' @param v String of name of variable for values
#' @param ref String of name of variable for references
#' @param ext String giving extension given to scale factor variable name
#' @importFrom magrittr %>%
#' @export

hyb_norm <- function(ds, v, ref, ext) {
  sf_name <- paste0('sf', ext)
  ds[variable %in% .hce, (sf_name) := median(get(ref) / get(v)), .(SlideId, SampleId, PlateId)]
  ds[, (sf_name) := get(sf_name) %>% fill_column(), .(SlideId, SampleId, PlateId)]
}

#' Intra-plate normalization
#'
#' @param ds data.table
#' @param v String of name of variable for values
#' @param ref String of name of variable for references
#' @param ext String giving extension given to scale factor variable name
#' @importFrom magrittr %>%
#' @export

intraplate_norm <- function(ds, v, ref, ext) {
  sf_name <- paste0('sf', ext)
  ds[SampleType == 'Calibrator',
    (sf_name) := median(get(ref) / get(v)),
          .(PlateId, SlideId, Dilution)]
  ds[SampleType != 'Calibrator', (sf_name) := 1]
}

#' Plate Scaling
#'
#' @param ds data.table
#' @param v String of name of variable for values
#' @param ref String of name of variable for references
#' @param ext String giving extension given to scale factor variable name
#' @importFrom magrittr %>%
#' @export

plate_norm <- function(ds, v, ref, ext) {
  sf_name <- paste0('sf', ext)
  ds[SampleType == 'Calibrator', (sf_name) := .SD %>% calc_ratios(v, ref) %>% median(), .(PlateId)]
  ds[, (sf_name) := get(sf_name) %>% fill_column(), .(PlateId)]
}

#' Calibration
#'
#' @param ds data.table
#' @param v String of name of variable for values
#' @param ref String of name of variable for references
#' @param ext String giving extension given to scale factor variable name
#' @importFrom magrittr %>%
#' @export

calibration_norm <- function(ds, v, ref, ext) {
  sf_name <- paste0('sf', ext)
  ds[SampleType == 'Calibrator', (sf_name) := get(ref) / median(get(v)), .(PlateId, variable)]
  ds[, (sf_name) := get(sf_name) %>% fill_column(), .(PlateId, variable)]
}

#' ANML
#'
#' @param ds data.table
#' @param v String of name of variable for values
#' @param ref String of name of variable for references
#' @param ext String giving extension given to scale factor variable name
#' @param mad String giving variable name of median absolute deviation of references
#' @param iteration Integer giving number of times to iterate
#' @importFrom magrittr %>%
#' @export

anml_norm <- function(ds, v, ref, ext, mad = 'mad', iterations = 100, verbose = FALSE) {
  sf_name <- paste0('sf', ext)

  ds[!(SampleType %in% c('Sample', 'QC')), (sf_name) := 1]
  ds[SampleType %in% c('Sample', 'QC'), (sf_name) := {
    if (verbose) cat('\r ANML: progress', (.GRP / .NGRP) %>% multiply_by(100) %>% round(1))
    anml_helper(value, get(ref), get(mad.ref), iterations = iterations, tolerance = 1e-4)
  }, .(Dilution, SampleId, SlideId)]
}

#' Sample level median fold change normalization (alternative to ANML)
#'
#' @param ds data.table
#' @param v String of name of variable for values
#' @param ref String of name of variable for references
#' @param ext String giving extension given to scale factor variable name
#' @param mad String giving variable name of median absolute deviation of references
#' @param iteration Integer giving number of times to iterate
#' @importFrom magrittr %>%
#' @export
sample_norm <- function(ds, v, ref, ext) {
  sf_name <- paste0('sf', ext)
  ds[, (sf_name) := {
    cat('\rprogress', (.GRP / .NGRP) %>% multiply_by(100) %>% round(1))
    fifelse(SampleType %in% c('Sample', 'QC'),
                            median(get(ref) / get(v)), 1)
    }, .(Dilution, SampleId, SlideId)]
}

