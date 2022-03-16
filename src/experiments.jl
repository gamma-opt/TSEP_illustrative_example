src_link  = "/Users/nikitabelyak/Dropbox (Aalto)/TSEP/illustrative example"
cd(src_link)
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using CSV, DataFrames, JuMP, Gurobi, PyCall, Plots, Ipopt, Statistics, XLSX
#np = pyimport("numpy")

# set unique envinronment for Gurobi
const GRB_ENV = Gurobi.Env()

# scaling factor for the monetary parameters 
scaling_factor = 1.0

# scaling factor for the objective coefficients 

obg_scaling_factor = 1.0

# maximum number of decimal digits when roudning 
#( minimum value should correspond to the biggest numebr of digits used in scaling factors)

max_digits = max(1, Int(max(log10(scaling_factor), log10(obg_scaling_factor))))
max_digits = 20


data_src_link  =  src_link * "/2022 data"

include(src_link*"/src/utils/data_preprocessing_functions.jl")
include(src_link*"/src/types/parameters.jl")
include(src_link*"/src/utils/parameters_initialisation.jl")
include(src_link*"/src/utils/models_generation.jl")
include(src_link*"/src/utils/models_generation_new.jl")
include(src_link*"/src/utils/data_postprocessing_functions.jl")

single_level_problem = single_level_problem_generation(input_parameters)
#io = open(data_src_link*"/single_level_problem_generation.txt" ,"w")
#println(io,single_level_problem)
#close(io)
optimize!(single_level_problem)
objective_value(single_level_problem)
print_output(data_src_link*"/optimisation_results" , input_parameters, objective_value(single_level_problem), value.(single_level_problem[:l_plus]), value.(single_level_problem[:g_VRES_plus]) , value.(single_level_problem[:g_conv_plus]), value.(single_level_problem[:g_VRES]), value.(single_level_problem[:g_conv]), value.(single_level_problem[:f]), value.(single_level_problem[:q]), incentives, input_parameters.gen_budget, "single_level", scaling_factor)
g_e_plus = sum(value.(single_level_problem[:g_VRES_plus]))
g_conv_plus = sum(value.(single_level_problem[:g_conv_plus]))
g_e = sum(value.(single_level_problem[:g_VRES]))*scaling_factor
g_conv = sum(value.(single_level_problem[:g_conv]))*scaling_factor
consumption = sum(value.(single_level_problem[:q]))*scaling_factor
g_e_plus*100/(g_e_plus+g_conv_plus)
value.(single_level_problem[:q])


#write_xlsx_output(data_src_link*"/optimisation_results", input_parameters, value.(single_level_problem[:g_VRES_plus]), value.(single_level_problem[:g_conv_plus]), value.(single_level_problem[:l_plus]), "centralised")
bi_level_problem_cournot = bi_level_problem_generation_new(input_parameters, "cournot")
#JuMP.fix.(bi_level_problem_cournot[:l_plus][:,:], 0; force = true)
optimize!(bi_level_problem_cournot)
objective_value(bi_level_problem_cournot)
print_output(data_src_link*"/optimisation_results", input_parameters, objective_value(bi_level_problem_cournot), value.(bi_level_problem_cournot[:l_plus]), value.(bi_level_problem_cournot[:g_VRES_plus]) , value.(bi_level_problem_cournot[:g_conv_plus]), value.(bi_level_problem_cournot[:g_VRES]), value.(bi_level_problem_cournot[:g_conv]), value.(bi_level_problem_cournot[:f]), value.(bi_level_problem_cournot[:q]), incentives, input_parameters.gen_budget, "bi_level_cournot", scaling_factor)
g_e_plus = sum(value.(bi_level_problem_cournot[:g_VRES_plus]))
g_conv_plus = sum(value.(bi_level_problem_cournot[:g_conv_plus]))
g_e = sum(value.(bi_level_problem_cournot[:g_VRES]))*scaling_factor
g_conv = sum(value.(bi_level_problem_cournot[:g_conv]))*scaling_factor
consumption = sum(value.(bi_level_problem_cournot[:q]))*scaling_factor
g_e_plus*100/(g_e_plus+g_conv_plus)



#check_duality_gap(bi_level_problem_cournot, input_parameters, "cournot")

bi_level_problem_perfect = bi_level_problem_generation_new(input_parameters, "perfect")
optimize!(bi_level_problem_perfect)
objective_value(bi_level_problem_perfect)
print_output(data_src_link*"/optimisation_results", input_parameters, objective_value(bi_level_problem_perfect), value.(bi_level_problem_perfect[:l_plus]), value.(bi_level_problem_perfect[:g_VRES_plus]) , value.(bi_level_problem_perfect[:g_conv_plus]), value.(bi_level_problem_perfect[:g_VRES]), value.(bi_level_problem_perfect[:g_conv]), value.(bi_level_problem_perfect[:f]), value.(bi_level_problem_perfect[:q]), incentives, input_parameters.gen_budget, "bi_level_perfect", scaling_factor)
g_e_plus = sum(value.(bi_level_problem_perfect[:g_VRES_plus]))
g_conv_plus = sum(value.(bi_level_problem_perfect[:g_conv_plus]))
g_e = sum(value.(bi_level_problem_perfect[:g_VRES]))*scaling_factor
g_conv = sum(value.(bi_level_problem_perfect[:g_conv]))*scaling_factor
consumption = sum(value.(bi_level_problem_perfect[:q]))*scaling_factor
g_e_plus*100/(g_e_plus+g_conv_plus)

#check_duality_gap(bi_level_problem_perfect, input_parameters, "perfect")





using Dualization
primal_lower_level_problem = primal_problem_generation(input_parameters)
dual_lower_level_problem_1 = dualize(primal_lower_level_problem)
set_optimizer(dual_lower_level_problem_1, () -> Gurobi.Optimizer(GRB_ENV))
optimize!(dual_lower_level_problem_1)
dual_1 = objective_value(dual_lower_level_problem_1)

dual_lower_level_problem_2 = dual_problem_generation(input_parameters)
optimize!(dual_lower_level_problem_2)
dual_2 = objective_value(dual_lower_level_problem_2)

dual_1 - dual_2