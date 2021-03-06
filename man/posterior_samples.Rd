% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{posterior_samples}
\alias{posterior_samples}
\title{Posterior sample function for}
\usage{
posterior_samples(
  effect_size = c("RR", "OR"),
  current_RCT_mean,
  current_RCT_lb,
  current_RCT_ub,
  event_percentage,
  cutoff,
  total_future_patients
)
}
\arguments{
\item{effect_size}{Effect size type, "RR" or "OR"}

\item{current_RCT_mean}{Mean effect for the current study}

\item{current_RCT_lb}{Lower bound for the confidence interval for the effect from the current study}

\item{current_RCT_ub}{Upper bound for the confidence interval for the effect from the current study}

\item{event_percentage}{Base rate/risk for the neutral hypothetical study}

\item{cutoff}{Lower limit for ROPE}

\item{total_future_patients}{}
}
\value{
A list containing posterior sample parameters
}
\description{
Posterior sample function for
}
\keyword{internal}
