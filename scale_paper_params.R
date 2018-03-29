initialise_user_global_params <- function(){
  
  global_params = list()
  
  global_params$scenario_subset = 'all'
  global_params$user_simulated_ecology_params_file = 'scale_paper_params.R'  # path to file
  global_params$number_of_cores = 'all'
  
  # Where simulation outputs will be written
  global_params$simulation_folder = 'default'
  
  # The number of realizations to run
  global_params$realisation_num = 1
  
  # Makes a single pdf at the end of the simulation showing the locatons of all offsets
  global_params$write_offset_layer = FALSE
  
  # Create an animation of the outputs
  global_params$write_movie = FALSE
  
  return(global_params)
}



initialise_user_simulation_params <- function(){ 
  
  simulation_params = list()
  
  # what subset of features to use in the simulation
  simulation_params$features_to_use_in_simulation = 1
  
  # The total number of layers to use in the offset calcuation (iterating from the start)
  simulation_params$features_to_use_in_offset_calc = 1
  
  simulation_params$features_to_use_in_offset_intervention = 1
  
  # How long to run the simulaton in years
  simulation_params$time_steps = 50
  
  # generate the intervention vector stochastically with parameters in format (time_steps, start, end, total_number, standard_deviation)
  simulation_params$intervention_vec = generate_stochastic_intervention_vec(time_steps = simulation_params$time_steps, 
                                                                            intervention_start = 1, 
                                                                            intervention_end = simulation_params$time_steps, 
                                                                            intervention_num = 50, 
                                                                            sd = 1)
  
  # The maxoimum number of parcels can be selected to offset a single development
  
  simulation_params$max_offset_parcel_num = 10
  
  # Stops the offset from delivering any further gains once it has acheived the gains required
  simulation_params$limit_offset_restoration = TRUE
  
  
  # The probability per parcel of it being illegally cleared, every parcel gets set to this number - set to zero to turn off
  simulation_params$unregulated_loss_prob = 0.002
  
  # Exclude parcels with less than this number of pixels.
  simulation_params$site_screen_size = 50
  
  # The mean and the standard deviation of a normal distribution from which to sample the restoration parameters from
  simulation_params$restoration_rate = 0.02
  
  simulation_params$restoration_rate_std = 0.005
  #   c('net_gains', 'restoration_gains', 'avoided_condition_decline', 'avoided_loss',
  #     'protected_condition', 'current_condition', 'restored_condition')
  
  simulation_params$offset_action_params = list(c('net_gains', 'restore'),
                                                c('restoration_gains', 'restore'),
                                                c('avoided_loss', 'maintain'))
  
  # This is the equivalent of offset_calc_type for the dev parcel. Options
  # are: 'current_condition' - losses are calcuated relative to the value of
  # the site at the time of the intervention 
  # 'future_condition' - is the do nothing trjectory of the development site.
  simulation_params$dev_calc_type = list('future_condition', 'current_condition')
  
  # Track accumulated credit from previous exchanges (eithger in current or
  # previous time step) and use them to allow developments to proceed if the
  # credit is large enough. FALSE means ignore any exces credit from offset exchanges
  simulation_params$allow_developments_from_credit = TRUE
  
  # How the development parcels are selected options are 'random' or
  # 'weighted'. Note tha weighted requires an additonal weighting layer. If
  # you are running on your own data you need to specify the weights file in
  # intialise_routines.R  (or put the files in simulation_inputs)
  simulation_params$development_selection_type = 'random'  
  
  # The time horizon in which the offset gains need to equal the devlopment impact
  simulation_params$offset_time_horizon = list(15, 30)
  
  # Include future legal developments in calculating contribution of avoided
  # losses to the impact of the offset. This increases the impact of the
  # offset (due to future losses that are avoided)
  simulation_params$include_potential_developments_in_offset_calc = list(TRUE, FALSE)
  
  # Include future stochastic developments in calculating contribution of avoided losses
  # to the impact of the offset. This increases the impact of the
  # offset (due to future losses that are avoided)
  simulation_params$include_unregulated_loss_in_offset_calc = list(TRUE, FALSE)
  
  simulation_params$dev_counterfactual_adjustment = 'as_offset'
  # The development impacts is multiplied by this factor (irrespective of how
  # they were caluclated) and the offset impact then needs to match this
  # multiplied development impact
  simulation_params$offset_multiplier = 1
  
  
  return(simulation_params)
  
}


initialise_user_simulated_ecology_params <- function(){
  
  # Construct the static initial landscape 
  
  simulated_ecology_params = list()
  
  #how many feature layers to generate
  simulated_ecology_params$feature_num = 1
  
  # Number of pixels in (y, x) for the feature layes 
  simulated_ecology_params$ecology_size = c(500, 500)
  
  simulated_ecology_params$mean_decline_rates = rep(list(-1e-2), simulated_ecology_params$feature_num)
  
  #set this parameter to zero to yield no noise
  simulated_ecology_params$decline_rate_std = rep(list(1e-3), simulated_ecology_params$feature_num)
  
  # Numnber of parcels in x (but total size varies)
  simulated_ecology_params$parcel_num_x = 50
  
  # Numnber of parcels in y (but total size varies)
  simulated_ecology_params$parcel_num_y = 50
  
  #how much the site dimensions should vary
  
  simulated_ecology_params$site_width_variation_param = 1
  # Minimum allowable initial ecological value of smallest ecological element
  # (pixel) ie min value to sample from
  simulated_ecology_params$min_initial_eco_val = 20
  
  # Max allowable initial ecological value of largest element (pixel) ie max
  # value to sample from
  simulated_ecology_params$max_initial_eco_val = 90
  
  simulated_ecology_params$occupation_ratio = list(1) 
  # Mow much initial variation in pixels per land parcel (this is the width of
  # uniform dist) used to add noise to each pixel. Eg if the pixel has a vlaue
  # of 35, a new value will be sampled from between 35-45
  simulated_ecology_params$initial_eco_noise = 10
  
  return(simulated_ecology_params)
}


initialise_user_output_params <- function(){
  output_params = list()
  output_params$output_plot_folder = vector()
  output_params$plot_type = 'impacts' # can be 'outcomes'  or 'impacts',
  output_params$output_plot = TRUE # switch to choose wheteher plots are output
  output_params$output_csv_file = TRUE #switch to choose whether impacts are exported as csv
  output_params$realisation_num = 'all' # 'all'  or number to plot
  output_params$write_pdf = FALSE
  output_params$sets_to_plot = 5 # example site to plot
  output_params$scenario_vec = 19 #c(1,4,7,10, 8, 2,3,5,6,9,11,12 ) #1:12
  output_params$site_impact_col_vec = c('darkgreen', 'red', 'black')
  output_params$program_col_vec = c('darkgreen', 'red', 'black') 
  output_params$cfac_col = 'blue' 
  output_params$landscape_col = 'black'
  output_params$lwd_vec = c(3, 0.5)
  
  output_params$plot_subset_type = 'all' #c('offset_action_type') # 'offset_calc_type', 'offset_action_type', offset_time_horizon'
  output_params$plot_subset_param = 'all' #c('maintain') # 'net_gains', 'restore', 15
  
  output_params$site_impact_lwd = 0.5
  output_params$site_outcome_lwd_vec = c(0.5)
  output_params$program_lwd_vec = c(3, 0.5)
  output_params$program_outcome_lwd_vec = c(3, 0.5)
  output_params$landscape_lwd_vec  = c(3)
  output_params$landscape_outcome_lwd_vec = c(3)
  
  output_params$string_width = 3 # how many digits are used to store scenario index and realisation index
  output_params$nx = 3 
  output_params$ny = 4
  
  output_params$site_outcome_plot_lims_set = rep(list(c(0, 3e4)), length(output_params$scenario_vec))
  output_params$program_outcome_plot_lims_set = rep(list(c(0e6, 1e7)), length(output_params$scenario_vec))
  output_params$landscape_outcome_plot_lims_set = rep(list(c(0, 2e7)), length(output_params$scenario_vec))
  
  output_params$site_impact_plot_lims_set = rep(list(c(-5e3, 5e3)), length(output_params$scenario_vec))
  output_params$program_impact_plot_lims_set = rep(list(c(-1e6, 1e6)), length(output_params$scenario_vec)) 
  output_params$landscape_impact_plot_lims_set = rep(list(c(-1e6, 0)), length(output_params$scenario_vec))
  
  
  
  return(output_params)
  
}