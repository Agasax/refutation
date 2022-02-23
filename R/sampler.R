#' Refutation estimation function
#'
#' @param effect_size Effect size type, "RR" or "OR"
#' @param current_RCT_mean Mean effect for the current study
#' @param current_RCT_lb Lower bound for the confidence interval for the effect from the current study
#' @param current_RCT_ub Upper bound for the confidence interval for the effect from the current study
#' @param event_percentage Base rate/risk for the neutral hypothetical study
#' @param cutoff Lower limit for ROPE
#' @param total_future_patients Size of hypothetica neutral future trial
#' @importFrom dplyr tibble
#' @importFrom purrr map
#' @return
#' @export
#'
#' @examples
#' refutation_sampler(effect_size = "RR" ,  # effect_size
#' current_RCT_mean = 0.7 ,                    # current_RCT_mean
#' current_RCT_lb = 0.51 ,                     # current_RCT_lb
#' current_RCT_ub = 0.96 ,                     # current_RCT_ub
#' event_percentage = 45 ,               # event_percentage
#' cutoff = 0.9 )                     # lower limit for ROPE
refutation_sampler= function(effect_size = c("RR", "OR"),
         current_RCT_mean,
         current_RCT_lb,
         current_RCT_ub,
         event_percentage,
         cutoff,
         total_future_patients)
{
  DT = dplyr::tibble(patients = seq(from = 2, to = 5000, by = 2))

  # Apply function
  output =
    DT |>
    dplyr::mutate(probability = purrr::map(
      patients,
      # use this columns as "total_future_patients"

      ~ posterior_samples(
        effect_size = effect_size,
        current_RCT_mean = current_RCT_mean,
        current_RCT_lb = current_RCT_lb,
        current_RCT_ub = current_RCT_ub,
        event_percentage = event_percentage,
        cutoff = cutoff,
        total_future_patients = .
      )
    ))
  return(output)
}
