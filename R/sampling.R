## ----setup, include=FALSE-------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----packages-------------------------------------------------------------------------------------------------
# Ensures the package "pacman" is installed
if (!require("pacman")) install.packages("pacman")

pacman::p_load(dplyr,
               flextable,
               tidytable,
               ggplot2,
               purrr,
               metafor,
               ggdist,
               cowplot, 
               gridGraphics)


## ----function-------------------------------------------------------------------------------------------------
posterior_samples = function(effect_size = c("RR", "OR"),
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

# Vectorize the function to be able to apply it below
posterior_samples = base::Vectorize(posterior_samples)



## ----user function inputs-------------------------------------------------------------------------------------
# Define inputs
ef = "RR"         # effect_size
mean = 0.7        # current_RCT_mean
lb = 0.51         # current_RCT_lb
ub = 0.96         # current_RCT_ub
percentage = 45   # event_percentage
cut = 0.9         # cutoff

prob_cutoff = 80  # probability cutoff defined as "null"


## ----run function---------------------------------------------------------------------------------------------

# Create vector from 2 to 20k
# "by 2" because we are creating future RCTs with equal allocation (1:1) 
# thus the increase in number of total patients must be by 2
DT = dplyr::tibble(patients = seq(from = 2, to = 5000, by = 2))

# Apply function
output = 
  DT |> 
  dplyr::mutate(
  probability = purrr::map(patients, # use this columns as "total_future_patients"
                           
                          ~posterior_samples(
                            effect_size = ef,
                            current_RCT_mean = mean,
                            current_RCT_lb = lb,
                            current_RCT_ub = ub,
                            event_percentage = percentage,
                            cutoff = cut,
                            total_future_patients = .)
                          )
  )



## ----calculate posterior probabilities------------------------------------------------------------------------
probability = data.frame(patients = output$patients)

# Extract calculated posterior probabilities
for (i in 1:nrow(output)) {
  
  # function to extract elements from stored lists in the "output" tibble
  p = purrr::pluck(output,
                   2, #2nd column
                   i, # row
                   7) # probability (7th list element)
  
  probability[i, "p"] = p
  
}
 
# Find how many patients are required to achieve the probability cutoff
new_df = subset(probability, p <= prob_cutoff)

cutoff = tail(new_df, 1)



## ----probability plot-----------------------------------------------------------------------------------------
p_prob = 
  probability |> 
  ggplot() +
  aes(x = patients, y = p) +
  geom_line() +
  geom_point(data = cutoff, color = "firebrick", size = 3) +
  scale_y_continuous(breaks = seq(0, 100, 20),
                     limits = c(0, 100)) +
  annotate("text", x = cutoff[[1,1]] + 1500, y = prob_cutoff - 10,
           label = paste0("At least ", cutoff[[1,1]], " patients"),
           color = "firebrick") +
  labs(x = "Required Future Patients",
       y = "Posterior Probability (%)",
       title = "Posterior Probabilities with Different Amounts of Patients\n") +
  theme_bw()


## ----forest plot----------------------------------------------------------------------------------------------

cutoff_output = 
  output |> 
  dplyr::filter(patients == cutoff[["patients"]])

extract_fun = function(element){purrr::pluck(cutoff_output,
                                             2, # 2nd column
                                             1, # access single list
                                             element)}

forest_data = data.frame(label = c("Current RCT",
                                   paste0("Future RCT (", cutoff$patients, " patients)")
                                   ),
                         yi = c(extract_fun(1), # "prior.mean"
                                extract_fun(3)), # "data.mean"
                         sei = c(extract_fun(2), # "prior.sd"
                                 extract_fun(4))) # "data.sd"

forest_fun = function(){
  with(forest_data,
       metafor::forest(x = yi, sei = sei,
                       order = desc(yi),
                       slab = label,
                       ylim = c(-1, 5),
                       rows = 1:2,
                       at = log(seq(0.5, 1.5, 0.25)),
                       pch = 19, # circles
                       atransf=exp,
                       refline = c(log(cut), log(1/cut)),
                       xlab = paste0(ef, " (log scale)")) # label used before
       )
  
  abline(h=0) # horizontal line
  
  metafor::addpoly(x = extract_fun(5), # "post.mean"
                   sei = extract_fun(6), # "post.sd"
                   
                   mlab = "Combined Results", # label
                   row=-0.5, # location
                   atransf=exp,
                   efac = 1.5 # polygon size
                   )
  
  text(x = -0.1, y = 4.5,
       paste0("Assumed Baseline Risk in Future RCT: ", percentage, "%")
       ) 
  
  text(x = -0.1, y = 3.5,
       paste0("ROPE: between ", cut, " and ", round(1/cut, 2), " ", ef)
       ) 
}

forest_fun()

p_forest = grDevices::recordPlot()


## ----posterior distribution plot------------------------------------------------------------------------------
p_dist =
  data.frame(mean = extract_fun(5), # "post.mean"
           sd = extract_fun(6)) |>  # "post.sd"
  ggplot() +
  aes(y = 1,
      xdist = distributional::dist_normal(mean, sd),
      fill = after_stat(log(cut) < x & x < log(1/cut))) +
  ggdist::stat_slab() +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  geom_vline(xintercept = log(c(cut, 1/cut)), linetype = 2 , color = "gray30") +
  annotate("text", y = 1.8, x = log(1.05),
           label = paste0(prob_cutoff, "%"), size = 6.5) +
  
  scale_x_continuous(breaks = log(seq(0.6, 1.2, 0.1)),
                     labels = seq(0.6, 1.2, 0.1)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(x = log(c(0.8, 1.2))) +
  labs(x = paste0(ef, " (log scale)"),
       y = "Density",
       title = "Combined Results\n(Posterior Distribution)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(), # To remove axis + ticks
        axis.text.y = element_blank()) 


## ----combined plot, fig.height= 5, fig.width=6----------------------------------------------------------------
# https://wilkelab.org/cowplot/articles/plot_grid.html

# first align the top-row plot (p_forest) with the left-most plot of the
# bottom row (p_dist)
plots = cowplot::align_plots(p_forest, p_dist, align = 'v', axis = 'l')

# then build the bottom row
bottom_row = cowplot::plot_grid(plots[[2]], p_prob, 
                                 label_size = 12,
                                rel_widths = c(1, 2))

# then combine with the top row for final plot
plot_grid(plots[[1]], bottom_row, label_size = 12, ncol = 1)



## ----2by2 table-----------------------------------------------------------------------------------------------
dplyr::tribble(
  ~" ", ~"event", ~"no event", ~"total",
  "Control group", "a", "b", "nc",
  "Experimental group", "c", "d", "ne"
) |> 
  flextable::flextable() |> 
  flextable::autofit()


## -------------------------------------------------------------------------------------------------------------
sessionInfo()

