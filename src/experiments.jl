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
print_output(data_src_link*"/optimisation_results" , input_parameters, objective_value(single_level_problem), value.(single_level_problem[:l_plus]), value.(single_level_problem[:g_VRES_plus]) , value.(single_level_problem[:g_conv_plus]), value.(single_level_problem[:g_VRES]), value.(single_level_problem[:g_conv]), value.(single_level_problem[:f]), "single_level", scaling_factor)
g_e_plus = sum(value.(single_level_problem[:g_VRES_plus]))
g_conv_plus = sum(value.(single_level_problem[:g_conv_plus]))
g_e = sum(value.(single_level_problem[:g_VRES]))*scaling_factor
g_conv = sum(value.(single_level_problem[:g_conv]))*scaling_factor
consumption = sum(value.(single_level_problem[:q]))*scaling_factor
g_e_plus*100/(g_e_plus+g_conv_plus)
value.(single_level_problem[:q])


#write_xlsx_output(data_src_link*"/optimisation_results", input_parameters, value.(single_level_problem[:g_VRES_plus]), value.(single_level_problem[:g_conv_plus]), value.(single_level_problem[:l_plus]), "centralised")

bi_level_problem_cournot = bi_level_problem_generation(input_parameters, "cournot")
optimize!(bi_level_problem_cournot)
objective_value(bi_level_problem_cournot)
print_output(data_src_link*"/optimisation_results", input_parameters, objective_value(bi_level_problem_cournot), value.(bi_level_problem_cournot[:l_plus]), value.(bi_level_problem_cournot[:g_VRES_plus]) , value.(bi_level_problem_cournot[:g_conv_plus]), value.(bi_level_problem_cournot[:g_VRES]), value.(bi_level_problem_cournot[:g_conv]), value.(bi_level_problem_cournot[:f]), "bi_level_cournot", scaling_factor)
g_e_plus = sum(value.(bi_level_problem_cournot[:g_VRES_plus]))
g_conv_plus = sum(value.(bi_level_problem_cournot[:g_conv_plus]))
g_e = sum(value.(bi_level_problem_cournot[:g_VRES]))*scaling_factor
g_conv = sum(value.(bi_level_problem_cournot[:g_conv]))*scaling_factor
consumption = sum(value.(bi_level_problem_cournot[:q]))*scaling_factor
g_e_plus*100/(g_e_plus+g_conv_plus)
value.(bi_level_problem_cournot[:θ]).*scaling_factor
q = value.(bi_level_problem_cournot[:q])
g_conv = value.(bi_level_problem_cournot[:g_conv])
g_conv_plus = value.(bi_level_problem_cournot[:g_conv_plus])
g_VRES = value.(bi_level_problem_cournot[:g_VRES])

g_VRES_plus = value.(bi_level_problem_cournot[:g_VRES_plus])

value.(bi_level_problem_cournot[:g_VRES_plus])
market = "cournot"

cournot_primal_obj =  sum( sum( inpa.scen_prob[s] *
    (inpa.id_intercept[s,t,n]*q[s,t,n] 
        - 0.5* inpa.id_slope[s,t,n]* q[s,t,n]^2

        - (market == "perfect" ? 0 : 
            (0.5* inpa.id_slope[s,t,n] * 
                sum( 
                    (sum( g_conv[s,t,n,i,e]  for e in 1:inpa.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:inpa.num_VRES))^2
                for i in 1:inpa.num_prod)
            )
        )
    
        - sum( 
            (inpa.conv.operational_costs[n,i,e] 
            + inpa.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
        for i in 1:inpa.num_prod, e in 1:inpa.num_conv)

        #- sum(inpa.transm.transmissio_costs[n,m]*f[s,t,n,m] 
        #for n in 1:inpa.num_nodes, m in 1:inpa.num_nodes)
        #- 1E-1 * sum( f[s,t,n,m] for m in 1:inpa.num_nodes)
        )

    for t in 1:inpa.num_time_periods, s in 1:inpa.num_scen)

    -sum( 

        sum( 
            inpa.vres.maintenance_costs[n,i,r]*
            (inpa.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
            + 
            inpa.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
        for r in 1:inpa.num_VRES)
    
        +
        
        sum( 
            inpa.conv.maintenance_costs[n,i,e]*
            (inpa.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
            + 
            inpa.conv.investment_costs[n,i,e]*g_conv_plus[n,i,e]
        for e in 1:inpa.num_conv) 
         
    for i in 1:inpa.num_prod)

    #- sum(

        #inpa.transm.maintenance_costs[n,m]*
        #(inpa.transm.installed_capacities[n,m] + l_plus[n,m])
        #+
        #inpa.transm.investment_costs[n,m]*l_plus[n,m]

    #for m in inpa.num_nodes)

for n in inpa.num_nodes )/scaling_factor 

β_f_1 = value.(bi_level_problem_cournot[:β_f_1])
β_f_2 = value.(bi_level_problem_cournot[:β_f_2])
β_conv = value.(bi_level_problem_cournot[:β_conv])
β_VRES = value.(bi_level_problem_cournot[:β_VRES])
β_up_conv = value.(bi_level_problem_cournot[:β_up_conv])
β_down_conv =value.(bi_level_problem_cournot[:β_down_conv])

cournot_dual_obj =  (sum( inpa.scen_prob[s] * 
        ( 0.5* inpa.id_slope[s,t,n]*
            ( q[s,t,n]^2 +  (market == "perfect" ? 0 : 
                                    (sum( 
                                        (sum( g_conv[s,t,n,i,e]  for e in 1:inpa.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:inpa.num_VRES))^2 
                                    for i in 1:inpa.num_prod)
                                    )
                            )
            )              
        )
    for n in 1:inpa.num_nodes, t in 1:inpa.num_time_periods, s in 1:inpa.num_scen)
    
    + sum(inpa.time_periods[t] * inpa.transm.installed_capacities[n,m] * (β_f_1[s,t,n,m] + β_f_2[s,t,n,m]) for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  m in 1:inpa.num_nodes)
    
    + sum(inpa.time_periods[t] * inpa.conv.installed_capacities[n,i,e] * β_conv[s,t,n,i,e] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, e in 1:inpa.num_conv)

    + sum(inpa.time_periods[t] * inpa.vres.availability_factor[s,t,n,r] * inpa.vres.installed_capacities[n,i,r] * β_VRES[s,t,n,i,r] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, r in 1:inpa.num_VRES)

    + sum(inpa.time_periods[t] * inpa.conv.ramp_up[n,i,e] * inpa.conv.installed_capacities[n,i,e] * β_up_conv[s,t,n,i,e] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, e in 1:inpa.num_conv)

    + sum(inpa.time_periods[t] * inpa.conv.ramp_down[n,i,e] * inpa.conv.installed_capacities[n,i,e] * β_down_conv[s,t,n,i,e] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, e in 1:inpa.num_conv)
)/scaling_factor 

cournot_dual_obj - cournot_primal_obj


#write_xlsx_output(data_src_link*"/optimisation_results", input_parameters, value.(bi_level_problem_cournot[:g_VRES_plus]), value.(bi_level_problem_cournot[:g_conv_plus]), value.(bi_level_problem_cournot[:l_plus]), "cournot")

bi_level_problem_perfect = bi_level_problem_generation(input_parameters, "perfect")
optimize!(bi_level_problem_perfect)
value.(bi_level_problem_perfect[:θ])
value.(bi_level_problem_perfect[:β_down_conv])
value.(bi_level_problem_perfect[:β_down_conv])

β_up_conv

input_parameters.conv.operational_costs

objective_value(bi_level_problem_perfect)
print_output(data_src_link*"/optimisation_results", input_parameters, objective_value(bi_level_problem_perfect), value.(bi_level_problem_perfect[:l_plus]), value.(bi_level_problem_perfect[:g_VRES_plus]) , value.(bi_level_problem_perfect[:g_conv_plus]), value.(bi_level_problem_perfect[:g_VRES]), value.(bi_level_problem_perfect[:g_conv]), value.(bi_level_problem_perfect[:f]), "bi_level_perfect", scaling_factor)
g_e_plus = sum(value.(bi_level_problem_perfect[:g_VRES_plus]))
g_conv_plus = sum(value.(bi_level_problem_perfect[:g_conv_plus]))
g_e = sum(value.(bi_level_problem_perfect[:g_VRES]))*scaling_factor
g_conv = sum(value.(bi_level_problem_perfect[:g_conv]))*scaling_factor
consumption = sum(value.(bi_level_problem_perfect[:q]))*scaling_factor
g_e_plus*100/(g_e_plus+g_conv_plus)

write_xlsx_output(data_src_link*"/optimisation_results", input_parameters, value.(bi_level_problem_perfect[:g_VRES_plus]), value.(bi_level_problem_perfect[:g_conv_plus]), value.(bi_level_problem_perfect[:l_plus]), "perfect")

q = value.(bi_level_problem_perfect[:q])
g_conv = value.(bi_level_problem_perfect[:g_conv])
g_conv_plus = value.(bi_level_problem_perfect[:g_conv_plus])
g_VRES = value.(bi_level_problem_perfect[:g_VRES])
g_VRES_plus = value.(bi_level_problem_perfect[:g_VRES_plus])

market = "perfect"
inpa = input_parameters

perfect_primal_obj =  sum( sum( inpa.scen_prob[s] *
    (inpa.id_intercept[s,t,n]*q[s,t,n] 
        - 0.5* inpa.id_slope[s,t,n]* q[s,t,n]^2

        - (market == "perfect" ? 0 : 
            (0.5* inpa.id_slope[s,t,n] * 
                sum( 
                    (sum( g_conv[s,t,n,i,e]  for e in 1:inpa.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:inpa.num_VRES))^2
                for i in 1:inpa.num_prod)
            )
        )
    
        - sum( 
            (inpa.conv.operational_costs[n,i,e] 
            + inpa.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
        for i in 1:inpa.num_prod, e in 1:inpa.num_conv)

        #- sum(inpa.transm.transmissio_costs[n,m]*f[s,t,n,m] 
        #for n in 1:inpa.num_nodes, m in 1:inpa.num_nodes)
        #- 1E-1 * sum( f[s,t,n,m] for m in 1:inpa.num_nodes)
        )

    for t in 1:inpa.num_time_periods, s in 1:inpa.num_scen)

    -sum( 

        sum( 
            inpa.vres.maintenance_costs[n,i,r]*
            (inpa.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
            + 
            inpa.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
        for r in 1:inpa.num_VRES)
    
        +
        
        sum( 
            inpa.conv.maintenance_costs[n,i,e]*
            (inpa.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
            + 
            inpa.conv.investment_costs[n,i,e]*g_conv_plus[n,i,e]
        for e in 1:inpa.num_conv) 
         
    for i in 1:inpa.num_prod)

    #- sum(

        #inpa.transm.maintenance_costs[n,m]*
        #(inpa.transm.installed_capacities[n,m] + l_plus[n,m])
        #+
        #inpa.transm.investment_costs[n,m]*l_plus[n,m]

    #for m in inpa.num_nodes)

for n in inpa.num_nodes )/scaling_factor 

β_f_1 = value.(bi_level_problem_perfect[:β_f_1])
β_f_2 = value.(bi_level_problem_perfect[:β_f_2])
β_conv = value.(bi_level_problem_perfect[:β_conv])
β_VRES = value.(bi_level_problem_perfect[:β_VRES])
β_up_conv = value.(bi_level_problem_perfect[:β_up_conv])
β_down_conv =value.(bi_level_problem_perfect[:β_down_conv])

perfect_dual_obj =  (sum( inpa.scen_prob[s] * 
        ( 0.5* inpa.id_slope[s,t,n]*
            ( q[s,t,n]^2 +  (market == "perfect" ? 0 : 
                                    (sum( 
                                        (sum( g_conv[s,t,n,i,e]  for e in 1:inpa.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:inpa.num_VRES))^2 
                                    for i in 1:inpa.num_prod)
                                    )
                            )
            )              
        )
    for n in 1:inpa.num_nodes, t in 1:inpa.num_time_periods, s in 1:inpa.num_scen)
    
    + sum(inpa.time_periods[t] * inpa.transm.installed_capacities[n,m] * (β_f_1[s,t,n,m] + β_f_2[s,t,n,m]) for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  m in 1:inpa.num_nodes)
    
    + sum(inpa.time_periods[t] * inpa.conv.installed_capacities[n,i,e] * β_conv[s,t,n,i,e] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, e in 1:inpa.num_conv)

    + sum(inpa.time_periods[t] * inpa.vres.availability_factor[s,t,n,r] * inpa.vres.installed_capacities[n,i,r] * β_VRES[s,t,n,i,r] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, r in 1:inpa.num_VRES)

    + sum(inpa.time_periods[t] * inpa.conv.ramp_up[n,i,e] * inpa.conv.installed_capacities[n,i,e] * β_up_conv[s,t,n,i,e] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, e in 1:inpa.num_conv)

    + sum(inpa.time_periods[t] * inpa.conv.ramp_down[n,i,e] * inpa.conv.installed_capacities[n,i,e] * β_down_conv[s,t,n,i,e] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, e in 1:inpa.num_conv)
)/scaling_factor 

perfect_dual_obj - perfect_primal_obj


input_parameters.conv.maintenance_costs
input_parameters.conv.investment_costs
input_parameters.conv.operational_costs

input_parameters.vres.maintenance_costs
input_parameters.vres.investment_costs

print(bi_level_problem_perfect)