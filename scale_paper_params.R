initialise_user_global_params <- function(){
  
  global_params = list()
  
#  global_params$user_simulated_ecology_params_file = 'scale_paper_params.R'  # path to file
  global_params$number_of_cores = 1

  global_params$simulation_folder = paste0(path.expand('~'), '/offset_data/scale_paper/')
  
  global_params$feature_raster_files = 'default' #paste0(global_params$simulation_folder, 'simulation_inputs/feature_001.tif')
  global_params$planning_units_raster = 'default' #paste0(global_params$simulation_folder, 'simulation_inputs/planning_units.tif')
  
  # a vector composed of what features to use in the simulation e.g. first 3 features is 1:3 NOT 3
  global_params$features_to_use_in_simulation = 1

  # Where simulation outputs will be written

  global_params$time_steps = 50
  
  # The number of realizations to run
  global_params$realisation_num = 1
  

  
  # if a file is supplied set this to false to use values in provided list of probabilities, otherwise set to true for equal probability of site development 
  global_params$overwrite_dev_probability_list = FALSE
  
  # if a file is supplied set this to false to use values in provided list of probabilities, otherwise set to true for equal probability of site offset 
  global_params$overwrite_offset_probability_list = FALSE
  
  # NOTE - set this to TRUE for initial simulation (when running with simulated data) - 
  # setting to false on subsequent runs will skip data builds and hence speed up simulations
  global_params$build_simulated_feature_layers = TRUE
  # if a file is supplied set this to false to use values in provided list of dynamics, otherwise set to true for on the fly dynamics calculations
  global_params$overwrite_management_dynamics = TRUE
  # if a file is supplied set this to false to use values in provided list of dynamics, otherwise set to true for on the fly dynamics calculations
  global_params$overwrite_feature_dynamics = TRUE
  # if a file is supplied set this to false to use values in provided raster layer of condition classes, otherwise set to true for on the fly condition class calculations
  global_params$overwrite_condition_classes = FALSE
  global_params$overwrite_features = FALSE
  global_params$overwrite_site_characteristics = FALSE
  
  return(global_params)
}



initialise_user_simulation_params <- function(global_params){ 
  
  simulation_params = list()
  
  simulation_params$offset_calc_type = list(c('net_gains'),
                                                c('restoration_gains'),
                                                c('avoided_condition_decline'))

  # This is the equivalent of offset_calc_type for the dev site. Options
  # are: 'current_condition' - losses are calcuated relative to the value of
  # the site at the time of the intervention 
  # 'future_condition' - is the do nothing trjectory of the development site.
  simulation_params$dev_calc_type = list('future_condition')    #'future_condition', 'current_condition' 
  
  
  # The total number of layers to use in the offset calcuation (iterating from the start)
  simulation_params$features_to_use_in_offset_calc = list(global_params$features_to_use_in_simulation)
  
  simulation_params$features_to_offset = list(global_params$features_to_use_in_simulation) 
  
  simulation_params$features_to_clear = list(global_params$features_to_use_in_simulation)  
  
  # simulation_params$use_transform_metric = list(FALSE)
  
  # The maximum number of sites can be selected to offset a single development
  
  simulation_params$max_offset_site_num = list(20)

  simulation_params$use_uncoupled_offsets = list(FALSE)
  
  simulation_params$uncoupled_offset_type = list('credit')
  
  simulation_params$uncoupled_offset_selection_type = list('stochastic')
  simulation_params$development_selection_type = list('stochastic')
  simulation_params$offset_selection_type = list('greedy')
  
  #define when the interventions are set to take place, in this case force to occur once per year
  simulation_params$development_control = list(build_stochastic_intervention(global_params$time_steps,
                                                                             intervention_start = 1,
                                                                             intervention_end = global_params$time_steps,
                                                                             intervention_num = 500,
                                                                             sd = 1))

  # Stops the offset from delivering any further gains once it has acheived the gains required
  simulation_params$limit_offset_restoration = list(TRUE)
  
  # The probability per site of it being unregulatedly cleared, every site gets set to this number - set to zero to turn off
  simulation_params$unregulated_loss_prob = list(1e-5)
  
  # Exclude sites with less than this number of elements
  simulation_params$min_site_screen_size = list(5)
  
  # setting (0-1] ignore sites with size above this number of elements 
  simulation_params$max_site_screen_size_quantile = list(0.99)
# 
 
  # Track accumulated credit from previous exchanges (eithger in current or
  # previous time step) and use them to allow developments to proceed if the
  # credit is large enough. FALSE means ignore any exces credit from offset exchanges
  simulation_params$develop_from_credit = list(TRUE)
  
  simulation_params$use_credit_in_offset = list(TRUE)
  
  # How the development sites are selected options are 'random' or
  # 'weighted'. Note tha weighted requires an additonal weighting layer. If
  # you are running on your own data you need to specify the weights file in
  # intialise_routines.R  (or put the files in simulation_inputs)
  
  # The time horizon in which the offset gains need to equal the devlopment impact
  simulation_params$offset_time_horizon = list(20)
  
  # Include future legal developments in calculating contribution of avoided
  # losses to the impact of the offset. This increases the impact of the
  # offset (due to future losses that are avoided)
  simulation_params$include_potential_developments_in_offset_calc = list(FALSE)
  
  # Include future unregulated developments in calculating contribution of avoided losses
  # to the impact of the offset. This increases the impact of the
  # offset (due to future losses that are avoided)
  simulation_params$include_unregulated_loss_in_offset_calc = list(FALSE)
  
  # Include unregulated clearing in the calculating the contribution of avoided
  # losses to the impact of the development. 
  simulation_params$include_unregulated_loss_in_dev_calc = list(simulation_params$include_unregulated_loss_in_offset_calc)
  
  simulation_params$dev_counterfactual_adjustment = list('as_offset')
  # The development impacts is multiplied by this factor (irrespective of how
  # they were caluclated) and the offset impact then needs to match this
  # multiplied development impact
  simulation_params$offset_multiplier = list(1)
  
  
  return(simulation_params)
  
}



build_logistic_dynamics_set <- function(logistic_params_set, condition_class_bounds, time_vec){
  
  
  # dynamics_set = lapply(seq_along(logistic_params_set), 
  #                       function(i) lapply(1:nrow(logistic_params_set[[i]]),
  #                                          function(j) lapply(1:ncol(logistic_params_set[[i]]),
  #                                                             function(k) logistic_projection(site_vals = logistic_params_set[[i]][i, j], 
  #                                                                                             x_min = condition_class_bounds[[i]][[j]][1], 
  #                                                                                             x_max = condition_class_bounds[[i]][[j]][3], 
  #                                                                                             current_dec_rate = logistic_params_set[[i]][[j]][[k]][2], 
  #                                                                                             time_vec = time_vec))))
  # 
  
  dynamics_set = lapply(seq_along(logistic_params_set), 
                        function(i) lapply(seq_along(logistic_params_set[[i]]),
                                           function(j) lapply(1:ncol(logistic_params_set[[i]][[j]]),
                                                              function(k) logistic_projection(site_vals = logistic_params_set[[i]][[j]][[k]][1], 
                                                                                              x_min = condition_class_bounds[[i]][j, names(condition_class_bounds[[i]]) == 'lower_bound'], 
                                                                                              x_max = condition_class_bounds[[i]][j, names(condition_class_bounds[[i]]) == 'upper_bound'], 
                                                                                              current_dec_rate = logistic_params_set[[i]][[j]][[k]][2], 
                                                                                              time_vec = time_vec))))
  dynamics_set = lapply(seq_along(logistic_params_set), 
                        function(i) lapply(seq_along(logistic_params_set[[i]]),
                                           function(j) as.data.frame(setNames(dynamics_set[[i]][[j]], c('lower_bound', 'best_estimate', 'upper_bound')))))
  
  return(dynamics_set)
}


logistic_projection <- function(site_vals, x_min, x_max, current_dec_rate, time_vec){
  
  t_sh = -1/current_dec_rate * log( ((site_vals - x_min)/(x_max - site_vals)))
  
  # define logistic curve given logistic parameter set.
  eco_projected = x_min + (x_max - x_min)/(1 + exp(-current_dec_rate*(time_vec - t_sh)))
  
  return(eco_projected)
}


initialise_user_feature_params <- function(global_params, simulation_params){
   
  feature_params = list()
 
  # How many feature layers to generate
  feature_params$simulated_feature_num = length(global_params$features_to_use_in_simulation)
  
  # The number of pixels in (y, x) for the feature layers 
  feature_params$feature_layer_size = c(1000, 1000)
  
  # Numnber of rectangiular sites for a given feature in x, y (but total size varies), third param is control the standard deviation for the sampling. 
  # If zero all sites the same size, greater values mean greater variation in sizes of the sites.
  feature_params$site_num_characteristics = c(100, 100, 5)
  
  # Numnber of features in x, y (but total size varies), third param is controld the sd deviation for the sampling. 
  # if zero all features the same size
  feature_params$feature_num_characteristics = c(30, 30, 10)
  

  # Specifies the proprotion of rectangular sites that are retained in the landscape for that feature. 
  # As sites are randomly removed, it them means there is more variation in terms of how much each given 
  # feature overlaps a land parcel. The more you crank up this value the more similar the sites become in value.
  feature_params$occupation_ratio = rep(list(0.3), feature_params$simulated_feature_num) 
  
  # Creates a list of data frams per feature that describe the condition class charateristics per feature.
  # For this analysis we just have one feature and one condition class. 
  feature_params$condition_class_bounds = vector('list', feature_params$simulated_feature_num)
  
  # Sets the bounds for the condition classes. For this analyis swe have single condition class with the bounds zero to one 
  for (feature_ind in seq(feature_params$simulated_feature_num)){
    feature_params$condition_class_bounds[[feature_ind]] = data.frame(feature = feature_ind, condition_class = 1, lower_bound = c(0), upper_bound = c(1))
  }
  feature_params$management_condition_class_bounds = feature_params$condition_class_bounds
  
  # Decline reate to use in the logistic condition curve
  decline_rate = -0.02
  
  # Per feature per condition class decline rates for unmanaged vegetation. 
  # This analysis is only using one feature and one conditoin class.
  # The lower bound means constant zero. Upper bound is constant 1. 
  # If don't have lower/upper 0/1 the for example, if you set lower bound to 0.3, ever cell with condition value 
  # less than that 0.3 will have the same decline values. 
  # Don't change lower and upper zero and 1.
  background_logistic_params_set = rep(list(list(data.frame(lower_bound = c(0, decline_rate), 
                                                            mean = c(0.5, decline_rate), 
                                                            upper_bound = c(1, decline_rate)))), 
                                       feature_params$simulated_feature_num)
  
  # Build 3 logistic curves and using the two params to improvement due to management. 
  # The first param is starting value, the 2nd param is the shape variable which determines the improvement rate 
  # Just change the rates here, don't change starting value, if make too small value takes ages to reach higher values
  # Incease values of the second param to increase restoration rate. Make 2nd param of the lower and upper bound bigger to make more variation.
  management_logistic_params_set = rep(list(list(data.frame(lower_bound = c(0.01, 0.04),
                                                            mean = c(0.01, 0.05), 
                                                            upper_bound = c(0.01, 0.06)))), 
                                       feature_params$simulated_feature_num)  
  
  # Length of dynamics curves use to model logistic change.
  feature_params$dynamics_time = 0:200
  
  # This functions builds the logistic curves for background the and management dynamics. 
  # Can put a browser() in here and see the curves.
  # This returns a dataframe, if you want you can replace with with your own csv file to have bespoke curves.
  feature_params$background_dynamics_bounds <- build_logistic_dynamics_set(background_logistic_params_set, 
                                                                           feature_params$condition_class_bounds,
                                                                           feature_params$dynamics_time)
  
  feature_params$management_dynamics_bounds <- build_logistic_dynamics_set(management_logistic_params_set, 
                                                                           feature_params$condition_class_bounds,
                                                                           feature_params$dynamics_time)
  
  
  return(feature_params)
}


dynamics_characteristics <- function(feature_num, condition_class_bounds){

  mean_logistic_rate = -0.02
  mean_restoration_rate = 0.04
  background_logistic_params_set = rep(list(list(list(c(0, mean_logistic_rate), c(0.5, mean_logistic_rate), c(1, mean_logistic_rate)))), feature_num)
  
  management_logistic_params_set = rep(list(list(list(c(0.01, 0.04), c(0.01, 0.05), c(0.01, 0.06)))), feature_num)
  
  dynamics_characteristics$dynamics_time = 0:200
  
  dynamics_characteristics$background_dynamics_bounds <- build_logistic_dynamics_set(background_logistic_params_set, 
                                                                   condition_class_bounds,
                                                                   dynamics_characteristics$dynamics_time)
  
  dynamics_characteristics$management_dynamics_bounds <- build_logistic_dynamics_set(management_logistic_params_set, 
                                                                   condition_class_bounds,
                                                                   dynamics_characteristics$dynamics_time)
  
  dynamics_characteristics$initial_condition_class_bounds = lapply(seq_along(dynamics_characteristics$background_dynamics_bounds), 
                                                         function(i) lapply(seq_along(dynamics_characteristics$background_dynamics_bounds[[i]]), 
                                                                            function(j) c(dynamics_characteristics$background_dynamics_bounds[[i]][[j]]$lower_bound[1], 
                                                                                          dynamics_characteristics$background_dynamics_bounds[[i]][[j]]$best_estimate[1], 
                                                                                          dynamics_characteristics$background_dynamics_bounds[[i]][[j]]$upper_bound[1])))
  
}
# 
# initialise_user_output_params <- function(){
#   output_params = list()
#   output_params$plot_type = 'impacts' # can be 'outcomes'  or 'impacts',
#   output_params$output_type = 'plot' #switch to choose whether impacts are exported as csv
#   output_params$realisation_num = 'all'; 10 #'all' # 'all'  or number to plot
#   output_params$write_pdf = TRUE
#   output_params$sets_to_plot = 5 # example site to plot
#   output_params$scenario_vec = 'all' #19 #c(1,4,7,10, 8, 2,3,5,6,9,11,12 ) #1:12
#   output_params$site_impact_col_vec = c('darkgreen', 'red', 'black')
#   output_params$program_col_vec = c('darkgreen', 'red', 'black') 
#   output_params$cfac_col = 'blue' 
#   output_params$landscape_col = 'black'
#   output_params$lwd_vec = c(3, 0.5)
#   output_params$print_dev_offset_sites = TRUE
# 
#   # Plot subset of the data. For example:
#   # output_params$plot_subset_type = c('offset_time_horizon', 'dev_calc_type' )
#   # output_params$plot_subset_param = c('30', 'future_condition')  
# 
#   # Note the two elemets of the vector for
#   # simulation_params$offset_action_params can be accessed by the variable
#   # offset_calc_type and offset_action_type. For example
#   #     output_params$plot_subset_type = c('offset_calc_type', 'offset_action_type' )
#   #     output_params$plot_subset_param = c('net_gains', 'restore')
# 
# 
#   # The conditions to apply to the to plots
#   my_subset_conditions <- matrix( ncol=2, byrow=TRUE, c(
# 
#     'offset_time_horizon',                             '15',
#     'dev_calc_type',                                   'future_condition' # current_condition, future_condition
#     ,'include_unregulated_loss_in_offset_calc',        'FALSE'
#     ,'include_potential_developments_in_offset_calc',  'TRUE'             # this should change for the dev too to be the same
#     #,'offset_calc_type',                               'avoided_condition_decline' # net_gains  avoided_condition_decline restoration_gains
# 
#     ) )
# 
#  
#   #output_params$plot_subset_type = 'all' #c('offset_action_type') # 'offset_calc_type', 'offset_action_type', offset_time_horizon'
#   #output_params$plot_subset_param = 'all' #c('maintain') # 'net_gains', 'restore', 15
#   
#   #output_params$plot_subset_type = my_subset_conditions[,1]
#   #output_params$plot_subset_param = my_subset_conditions[,2]
# 
#   output_params$site_impact_lwd = 0.5
#   output_params$site_outcome_lwd_vec = c(0.5)
#   output_params$program_lwd_vec = c(3, 0.5)
#   output_params$program_outcome_lwd_vec = c(3, 0.5)
#   output_params$landscape_lwd_vec  = c(3)
#   output_params$landscape_outcome_lwd_vec = c(3)
# 
#   output_params$nx = 3 
#   output_params$ny = 4
# 
#   
#   output_params$site_outcome_plot_lims_set = rep(list(list(c(0, 3e4))), length(output_params$scenario_vec))
#   output_params$program_outcome_plot_lims_set = rep(list(list(c(0e6, 1e7))), length(output_params$scenario_vec))
#   output_params$landscape_outcome_plot_lims_set = rep(list(list(c(0, 2e7))), length(output_params$scenario_vec))
#   
# 
#   output_params$site_impact_plot_lims_set = rep(list(list(c(-7e3, 1.1e4))), length(output_params$scenario_vec))
#   output_params$program_impact_plot_lims_set = rep(list(list(c(-3e5, 3e5))), length(output_params$scenario_vec)) 
#   output_params$landscape_impact_plot_lims_set = rep(list(list(c(-25e5, 0))), length(output_params$scenario_vec))
# 
#   
#   
#   return(output_params)
#   
# }


initialise_user_output_params <- function(global_params){
  output_params = list()
  output_params$plot_type = 'outcomes' # can be 'outcomes'  or 'impacts',
  output_params$output_type = 'plot' #switch to choose whether impacts are exported as csv
  output_params$realisation_num = 'all';  #'all' # 'all'  or number to plot
  output_params$write_pdf = FALSE
  output_params$sets_to_plot = 6 # example site to plot
  output_params$scenario_vec = 3 #19 #c(1,4,7,10, 8, 2,3,5,6,9,11,12 ) #1:12
  output_params$site_impact_col_vec = c('darkgreen', 'red', 'black')
  output_params$program_col_vec = c('darkgreen', 'red', 'black') 
  output_params$cfac_col = 'blue' 
  output_params$landscape_col = 'black'
  output_params$lwd_vec = c(3, 0.5)
  output_params$print_dev_offset_sites = TRUE

  # Plot subset of the data. For example:
  # output_params$plot_subset_type = c('offset_time_horizon', 'dev_calc_type' )
  # output_params$plot_subset_param = c('30', 'future_condition')  
  
  # Note the two elemets of the vector for
  # simulation_params$offset_action_params can be accessed by the variable
  # offset_calc_type and offset_action_type. For example
  #     output_params$plot_subset_type = c('offset_calc_type', 'offset_action_type' )
  #     output_params$plot_subset_param = c('net_gains', 'restore')
  
  
  # The conditions to apply to the to plots
  my_subset_conditions <- matrix( ncol=2, byrow=TRUE, c(
    
    'offset_time_horizon',                             '15',
    'dev_calc_type',                                   'future_condition' # current_condition, future_condition
    ,'include_unregulated_loss_in_offset_calc',        'FALSE'
    ,'include_potential_developments_in_offset_calc',  'TRUE'             # this should change for the dev too to be the same
    #,'offset_calc_type',                               'avoided_condition_decline' # net_gains  avoided_condition_decline restoration_gains
    
  ) )
  
  
  #output_params$plot_subset_type = 'all' #c('offset_action_type') # 'offset_calc_type', 'offset_action_type', offset_time_horizon'
  #output_params$plot_subset_param = 'all' #c('maintain') # 'net_gains', 'restore', 15
  
  #output_params$plot_subset_type = my_subset_conditions[,1]
  #output_params$plot_subset_param = my_subset_conditions[,2]
  
  output_params$site_impact_lwd = 0.5
  output_params$site_outcome_lwd_vec = c(0.5)
  output_params$program_lwd_vec = c(3, 0.5)
  output_params$program_outcome_lwd_vec = c(3, 0.5)
  output_params$landscape_lwd_vec  = c(3)
  output_params$landscape_outcome_lwd_vec = c(3)
  
  output_params$nx = 3 
  output_params$ny = 4
  
  
  # output_params$site_outcome_plot_lims_set = rep(list(list(c(0, 3e3))), max(output_params$scenario_vec))
  # output_params$program_outcome_plot_lims_set = rep(list(list(c(0, 25e3))), max(output_params$scenario_vec))
  # output_params$landscape_outcome_plot_lims_set = rep(list(list(c(0, 50e3))), max(output_params$scenario_vec))

  output_params$site_outcome_plot_lims_set = rep(list(rep(list(c(0, 10e1)), max(global_params$features_to_use_in_simulation))), max(output_params$scenario_vec))
  output_params$program_outcome_plot_lims_set = rep(list(rep(list(c(0, 0.3e5)), max(global_params$features_to_use_in_simulation))), max(output_params$scenario_vec)) 
  output_params$landscape_outcome_plot_lims_set = rep(list(rep(list(c(0.1e5, 0.3e5)), max(global_params$features_to_use_in_simulation))), max(output_params$scenario_vec))

  output_params$site_impact_plot_lims_set = rep(list(rep(list(c(-1e1, 1e1)), max(global_params$features_to_use_in_simulation))), max(output_params$scenario_vec))
  output_params$program_impact_plot_lims_set = rep(list(rep(list(c(-1e4, 1e4)), max(global_params$features_to_use_in_simulation))), max(output_params$scenario_vec)) 
  output_params$landscape_impact_plot_lims_set = rep(list(rep(list(c(-1e4, 1e4)), max(global_params$features_to_use_in_simulation))), max(output_params$scenario_vec))
  
  
  
  return(output_params)
  
}
