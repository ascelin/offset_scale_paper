
# to run:
#    source('run_offset_simulation_with_simulated_data.R')


library(offsetsim)


source('scale_paper_params.R')

user_simulation_params = initialise_user_simulation_params()
user_global_params = initialise_user_global_params()
user_simulated_ecology_params = initialise_user_simulated_ecology_params()
user_output_params <- initialise_user_output_params()

osim.run(user_global_params, user_simulation_params, user_simulated_ecology_params, loglevel = 'INFO')
#include run_number for specified run folder - leave to automatically select latest

current_simulation_folder = find_current_run_folder(user_global_params$simulation_folder)
#current_simulation_folder = '/Users/ascelin/analysis/offset_simulator/new_package_runs/simulation_runs/00008'

# current_simulation_folder = '/Users/ascelin/analysis/offset_simulator/new_package_runs/March/00002' # old code
# current_simulation_folder = '/Users/ascelin/analysis/offset_simulator/new_package_runs/March/00003' # mmew code

osim.output(user_output_params, current_simulation_folder, loglevel = 'INFO')
