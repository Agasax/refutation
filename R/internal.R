


#' Posterior sample function for
#'
#' @param effect_size Effect size type, "RR" or "OR"
#' @param current_RCT_mean Mean effect for the current study
#' @param current_RCT_lb Lower bound for the confidence interval for the effect from the current study
#' @param current_RCT_ub Upper bound for the confidence interval for the effect from the current study
#' @param event_percentage Base rate/risk for the neutral hypothetical study
#' @param cutoff Lower limit for ROPE
#' @param total_future_patients
#'
#' @return A list containing posterior sample parameters
#'
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @keywords internal
#'
#' @examples
posterior_samples =
  function(effect_size = c("RR", "OR"),
                             current_RCT_mean,
                             current_RCT_lb,
                             current_RCT_ub,
                             event_percentage,
                             cutoff,
                             total_future_patients){

  ## Current RCT data
  prior.mean = log(current_RCT_mean)
  prior.sd = (log(current_RCT_ub) - log(current_RCT_lb))/(2*qnorm(0.975))
  prior.var = prior.sd^2

  ## Future RCT data
  risk = event_percentage/100
  # 1:1 allocation (thus, divide total by 2)
  patients_per_arm = round(total_future_patients)/2

  # Calculate sampling variance (sigma^2)

  if (effect_size == "RR"){

     events_control = risk*patients_per_arm
     # Because there is no effect (mean RR = 1), the number of events are the same
     events_exp = events_control

     data.var = 1/events_control - 1/patients_per_arm + 1/events_exp - 1/patients_per_arm
  }

  if (effect_size == "OR"){

  # Experimental arm
  a = risk*patients_per_arm       # Events
  b = patients_per_arm - a       # No events
  # Control arm
  c = risk*patients_per_arm       # Events
  d = patients_per_arm - c       # No events

  data.var = 1/(a + 1/2) + 1/(b + 1/2) + 1/(c + 1/2) + 1/(d + 1/2)
  }

  # No effect (effect size = 1)
  data.mean = log(1)

  ### Normal conjugate analyses
  post.mean.numerator = prior.mean/prior.var + data.mean/data.var
  post.mean.denominator = 1/prior.var + 1/data.var
  post.mean =  post.mean.numerator/post.mean.denominator
  post.sd = sqrt(1/(1/prior.var + 1/data.var))

  # Calculate posterior probability of interest (Pr(ROPE))
  # ROPE, region of practical equivalence

  upper_cutoff = 1/cutoff

   prob_upper_cutoff = # eg, Pr(< log(1.11))
    pnorm(log(upper_cutoff),
          mean = post.mean, sd = post.sd)

  prob_lower_cutoff = # eg, Pr(< log(0.9))
    pnorm(log(cutoff),
          mean = post.mean, sd = post.sd)


  rope_prob = # eg, Pr(> log(0.9) & < log(1.11))

# https://stackoverflow.com/questions/46031346/probability-of-a-column-between-a-range-for-a-normal-distribution/46033017#46033017
    prob_upper_cutoff - prob_lower_cutoff

  rope_prob = 100*rope_prob

  d = list(prior.mean = prior.mean,
           prior.sd = prior.sd,
           data.mean = data.mean,
           data.sd = sqrt(data.var),
           post.mean = post.mean,
           post.sd = post.sd,
           probability = rope_prob)

  return(d)

}


