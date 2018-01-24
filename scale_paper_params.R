initialise_user_global_params <- function(){
  
  global_params = list()
  
  global_params$overwrite_default_params = TRUE
  
  global_params$user_simulated_ecology_params_file = 'user_params/user_simulated_ecology_params.R'  # path to file
  
  global_params$number_of_cores = 1
  # Where simulation outputs will be written
  global_params$simulation_folder = paste0(path.expand('~'), '/offset_data/simulated/')
  
  # The number of realizations to run
  global_params$realisation_num = 1
  
  # Makes a single pdf at the end of the simulation showing the locatons of all offsets
  global_params$write_offset_layer = TRUE
  
  # Create an animation of the outputs
  global_params$write_movie = TRUE
 
  return(global_params)
}



initialise_user_simulation_params <- function(){ 
  
  simulation_params = list()
  
  # what subset of features to use in the simulation
  simulation_params$features_to_use_in_simulation = 1
  
  # The total number of layers to use in the offset calcuation (iterating from the start)
  simulation_params$features_to_use_in_offset_calc = 1
  
  simulation_params$features_to_use_in_offset_intervention = 1
  
  # The total number of parcels that will be developed
  simulation_params$total_dev_num = 100
  
  # The time step at which development starts
  simulation_params$dev_start = 1
  
  # The time at which development ends
  simulation_params$dev_end = 50
  
  # How long to run the simulaton in years
  simulation_params$time_steps = 50
  
  # The maxoimum number of parcels can be selected to offset a single development
  
  simulation_params$max_offset_parcel_num = 10
  
  # Stops the offset from delivering any further gains once it has acheived the gains required
  simulation_params$limit_offset_restoration = TRUE
  
  # The probability per parcel of it being illegally cleared, every parcel gets set to this number - set to zero to turn off
  simulation_params$illegal_clearing_prob = 0
  
  # Exclude parcels with less than this number of pixels.
  simulation_params$site_screen_size = 50
  
  # The mean and the standard deviation of a normal distribution from which to sample the restoration parameters from
  simulation_params$restoration_rate = 0.02
  
  simulation_params$restoration_rate_std = 0.005
#   c('net_gains', 'restoration_gains', 'avoided_condition_decline', 'avoided_loss',
#     'protected_condition', 'current_condition', 'restored_condition')
  
  simulation_params$offset_action_params = list(c('net_gains', 'restore'), 
                                                c('avoided_condition_decline', 'maintain'))

  # This is the equivalent of offset_calc_type for the dev parcel. Options
  # are: 'current_condition' - losses are calcuated relative to the value of
  # the site at the time of the intervention 
  # 'future_condition' - is the do nothing trjectory of the development site.
  simulation_params$dev_calc_type = 'future_condition'    #'future_condition', 'current_condition' 
  
  # Track accumulated credit from previous exchanges (eithger in current or
  # previous time step) and use them to allow developments to proceed if the
  # credit is large enough. FALSE means ignore any exces credit from offset exchanges
  simulation_params$allow_developments_from_credit = TRUE
  
  # How the development parcels are selected options are 'random' or
  # 'weighted'. Note tha weighted requires an additonal weighting layer. If
  # you are running on your own data you need to specify the weights file in
  # intialise_routines.R  (or put the files in simulation_inputs)
  simulation_params$development_selection_type = 'random'  
  
  # Whether to use banking. FALSE - means perform offsets simultaneously with development, TRUE -
  # means perform offset banking prior to development according to offset bank
  # parameters
  simulation_params$use_offset_bank = c(FALSE)
  
  # The time at which the offset in the bank offsets are first are implemented and start acurring grains, 
  simulation_params$offset_bank_start = 1 
  
  # The time at which no more offsets are added to the bank. The number of
  # offsets per time step is determined as follows: First the mean number
  # number per time step is determined, then sampling is done around this
  # mean number using a normal distribution such that the total number of
  # developments will always equal the total number (Note sd for this
  # distribution is set in the code the currently isn't user settable)
  simulation_params$offset_bank_end = 1 
  
  # THe number parcels to include in banking scheme. These are randomly selected.
  simulation_params$offset_bank_num = 200 
  
  # Options are 'credit' or 'parcel_set'. 'credit' means there is accumulated
  # gain that is subtracted as parcels are developed. 'parcel_set' one or more
  # parcels in the bank are traded for one development site. If there is left
  # over credit (and allow_developments_from_credit is set to TRUE) then this excess credit is used on subsequent developments
  simulation_params$offset_bank_type = c('credit') #c('parcel_set', 'credit')     
  
  # The time horizon in which the offset gains need to equal the devlopment impact
  simulation_params$offset_time_horizon = list(15)
  
  # Include illegal clearing in the calculating the contribution of avoided
  # losses to the impact of the development. 
  # simulation_params$include_illegal_clearing_in_dev_calc = simulation_params$include_illegal_clearing_in_offset_calc
  
  # Include future legal developments in calculating contribution of avoided
  # losses to the impact of the offset. This increases the impact of the
  # offset (due to future losses that are avoided)
  simulation_params$include_potential_developments_in_offset_calc = list(TRUE)
  
  # Include future illegal developments in calculating contribution of avoided losses
  # to the impact of the offset. This increases the impact of the
  # offset (due to future losses that are avoided)
  simulation_params$include_illegal_clearing_in_offset_calc = list(TRUE)
  
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
  simulated_ecology_params$feature_num = 3
  
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
  
  simulated_ecology_params$occupation_ratio = list(0.80, 0.5, 0.2) 
  # Mow much initial variation in pixels per land parcel (this is the width of
  # uniform dist) used to add noise to each pixel. Eg if the pixel has a vlaue
  # of 35, a new value will be sampled from between 35-45
  simulated_ecology_params$initial_eco_noise = 10
  
  # Defining multiple regions eg different states where different polcies can apply 
  simulated_ecology_params$region_num_x = 1
  
  # Defining multiple regions eg different states where different rules can apply 
  simulated_ecology_params$region_num_y = 1
  
  return(simulated_ecology_params)
}


initialise_user_plot_params <- function(){
  plot_params = list()
  plot_params$output_plot_folder = vector()
  plot_params$plot_type = 'impacts' # can be 'outcomes'  or 'impacts',
  plot_params$output_type = 'scenarios' # set to plot through 'features', 'scenarios' or 'site_sets'
  plot_params$realisation_num = 1 # 'all' or number to plot
  plot_params$write_pdf = FALSE
  plot_params$sets_to_plot = 5 # example site to plot
  plot_params$scenario_vec = 'all' #c(1,4,7,10, 8, 2,3,5,6,9,11,12 ) #1:12
  plot_params$site_impact_col_vec = c('darkgreen', 'red', 'black')
  plot_params$program_col_vec = c('darkgreen', 'red', 'black') 
  plot_params$cfac_col = 'blue' 
  plot_params$landscape_col = 'black'
  plot_params$lwd_vec = c(3, 0.5)
  
  plot_params$plot_subset_type = c('dev_calc_type') # 'offset_calc', 'time_horizon'
  plot_params$plot_subset_param = c('future_condition')
  
  plot_params$site_impact_lwd = 0.5
  plot_params$site_outcome_lwd_vec = c(0.5)
  plot_params$program_lwd_vec = c(3, 0.5)
  plot_params$program_outcome_lwd_vec = c(3, 0.5)
  plot_params$landscape_lwd_vec  = c(3)
  plot_params$landscape_outcome_lwd_vec = c(3)
  
  plot_params$string_width = 3 # how many digits are used to store scenario index and realisation index
  plot_params$nx = 3 
  plot_params$ny = 4
  
  plot_params$site_outcome_plot_lims_set = rep(list(c(0, 3e4)), length(plot_params$scenario_vec))
  plot_params$program_outcome_plot_lims_set = rep(list(c(0e6, 1e7)), length(plot_params$scenario_vec))
  plot_params$landscape_outcome_plot_lims_set = rep(list(c(0, 2e7)), length(plot_params$scenario_vec))
  
  plot_params$site_impact_plot_lims_set = rep(list(c(-5e3, 5e3)), length(plot_params$scenario_vec))
  plot_params$program_impact_plot_lims_set = rep(list(c(-1e6, 1e6)), length(plot_params$scenario_vec)) 
  plot_params$landscape_impact_plot_lims_set = rep(list(c(-1e6, 0)), length(plot_params$scenario_vec))
  
  
  
  return(plot_params)
}

