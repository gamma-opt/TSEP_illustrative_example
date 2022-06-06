""" 
function print_output(...) prints the output resutling from the model optimisation 
    in the user friendly format that is further saved in the single_level_problem_output.txt file 
    in the folder specified in the arguments of the function 
Arguments of the function:
* src_link::String:                 The link at which the output .txt file to be stored 
* ip::initial_parameters            Initial parameters of the model 
* ob_val::Float64                   The objective function value 
* l_plus::Array{Float64}            The optimal values the transmission lines capacity expansion 
* g_VRES_plus::Array{Float64}       The optimal values the VRES capacity expansion
* g_conv_plus::Array{Float64}       The optimal values conventional sources capacity expansion
* g_VRES::Array{Float64}            The optimal VRES generation values
* g_conv::Array{Float64}            The optimal convetional sources generation values
* f::Array{Float64}                 The optimal energy flow values
* model::String                     The indicator on which model are we considering (single_level, bi_level_cournot or bi_level_perfect)
* sf::Float64                       The scaling factor for the model parameters
"""


function print_output(src_link::String, ip::initial_parameters, ob_val::Float64, l_plus::Array{Float64}, g_VRES_plus::Array{Float64}, g_conv_plus::Array{Float64}, g_VRES::Array{Float64}, g_conv::Array{Float64}, f::Array{Float64}, q::Array{Float64}, incentives::Array{Float64}, gen_budget::Array{Float64}, model::String, sf::Float64)
    
    #rounding all the values 
    round_digits = 3
    ob_val = round(ob_val, digits = round_digits)
    l_plus = round.(l_plus, digits = round_digits)
    g_VRES_plus = round.(g_VRES_plus, digits = round_digits)
    g_conv_plus = round.(g_conv_plus, digits = round_digits)
    g_VRES = round.(g_VRES, digits = round_digits)
    g_conv = round.(g_conv, digits = round_digits)
    f = round.(f, digits = round_digits)
    
    if model == "single_level"
        io = open(src_link*"/single_level_problem_output_"* string(Int(input_parameters.transm.budget_limit)) * "_budget_per_line_" * string(Int(input_parameters.conv.CO2_tax[1])) * "_tax_"* string(incentives[1]) * "_incentives_" * string(gen_budget[1]) * "_gen_exp_budget_com_1_" * string(gen_budget[2])* "_gen_exp_budget_com_2.txt","w")
    elseif model == "bi_level_cournot"
        io = open(src_link*"/bi_level_cournot_problem_output_" * string(Int(input_parameters.transm.budget_limit)) * "_budget_per_line_" * string(Int(input_parameters.conv.CO2_tax[1])) * "_tax_"* string(incentives[1]) * "_incentives_" * string(gen_budget[1]) * "_gen_exp_budget_com_1_" * string(gen_budget[2])* "_gen_exp_budget_com_2.txt" ,"w")
    else 
        io = open(src_link*"/bi_level_perfect_problem_output_" * string(Int(input_parameters.transm.budget_limit)) * "_budget_per_line_" * string(Int(input_parameters.conv.CO2_tax[1])) * "_tax_"* string(incentives[1]) * "_incentives_" * string(gen_budget[1]) * "_gen_exp_budget_com_1_" * string(gen_budget[2])* "_gen_exp_budget_com_2.txt" ,"w")
    end

    println(io, "OBJECTIVE FUNCTION VALUE  $(ob_val*sf)")
    println(io, "\n")


    println(io, "TRANSMISSION CAPACITY EXPANSION")
    println(io, "\n")
    nodes = string.(1:ip.num_nodes)
    df_transm_exp = DataFrame()
    df_transm_exp.node = nodes
    for n = 1:ip.num_nodes
        insertcols!(df_transm_exp, n+1, "node_$n" => l_plus[n, :].*sf )
    end
    println(io, df_transm_exp)
    println(io, "\n")

    println(io, "VRES CAPACITY EXPANSION")
    println(io, "\n")
    for i = 1:ip.num_prod
        println(io, "PRODUCER $(i)")
        df_VRES_exp = DataFrame()
        df_VRES_exp.node = nodes
        for r in 1:ip.num_VRES
            insertcols!(df_VRES_exp, r+1, "VRES_$r" => g_VRES_plus[:,i,r].* sf )
        end
        println(io, df_VRES_exp)
        println(io, "\n")
    end
    println(io, "\n")

    println(io, "CONV CAPACITY EXPANSION")
    println(io, "\n")
    for i = 1:ip.num_prod
        println(io, "PRODUCER $(i)")
        df_conv_exp = DataFrame()
        df_conv_exp.node = nodes
        for e in 1:ip.num_conv
            insertcols!(df_conv_exp, e+1, "conv_$e" => g_conv_plus[:,i,e] .*sf)
        end
        println(io, df_conv_exp)
        println(io, "\n")
    end
    println(io, "\n")

    println(io, "VRES GENERATION")
    println(io, "\n")
    for i = 1:ip.num_prod
        println(io, "PRODUCER $(i)")
        for s = 1:ip.num_scen
            println(io, "scenario $(s)")
            df_VRES = DataFrame()
            df_VRES.node = nodes
            for r in 1:ip.num_VRES
                insertcols!(df_VRES, r+1, "VRES_$r" => vec(sum(g_VRES[s,:,:,i,r], dims = 1).* sf) )
            end
            println(io, df_VRES)
            println(io, "\n")
        end
    end
    println(io, "\n")

    println(io, "CONV GENERATION")
    println(io, "\n")
    for i = 1:ip.num_prod
        println(io, "PRODUCER $(i)")
        for s = 1:ip.num_scen
            println(io, "scenario $(s)")
            df_conv = DataFrame()
            df_conv.node = nodes
            for e in 1:ip.num_conv
                insertcols!(df_conv, e+1, "conv_$e" => vec(sum(g_conv[s,:,:,i,e], dims = 1) .* sf) )
            end
            println(io, df_conv)
            println(io, "\n")
        end
    end
    println(io, "\n")

    println(io, "ENERGY FLOW")
    println(io, "\n")
    for n = 1:ip.num_nodes
        for m = n+1:ip.num_nodes
            println(io, "flow: $n -> $m ")
            df_flow = DataFrame()
            df_flow.scenario = 1:ip.num_scen
            for t in 1:ip.num_time_periods
                insertcols!(df_flow, t+1, "time_$t" => f[:,t,n,m] .* sf )
            end
            println(io, df_flow)
            println(io, "\n")

            println(io, "flow: $m -> $n ")
            df_flow = DataFrame()
            df_flow.scenario = 1:ip.num_scen
            for t in 1:ip.num_time_periods
                insertcols!(df_flow, t+1, "time_$t" => f[:,t,m,n] .* sf)
            end
            println(io, df_flow)
            println(io, "\n")
        end
    end
    println(io, "\n")

    println(io, "TOTAL PRODUCTION/CONSUMPTION")
    println(io, "\n")
    df_total = DataFrame()
    df_total.scenario = 1:ip.num_scen
    
    #summing up the total of vres generation
    total_vres = Array{Float64}(undef, ip.num_scen)
    for s = 1:ip.num_scen
        total_vres[s] = sum(g_VRES[s,:,:,:,:])
    end
    insertcols!(df_total, 2, "VRES_generation" => total_vres[:] .* sf)

    #summing up the total consumption
    total_consump = Array{Float64}(undef, ip.num_scen)
    for s = 1:ip.num_scen
        total_consump[s] = sum(q[s,:,:])
    end


    #calculating the percentage of vres 
    vres_share = Array{Float64}(undef, ip.num_scen)
    for s = 1:ip.num_scen
        vres_share[s] = round(total_vres[s]*sf*100/total_consump[s]*sf, digits = 0)
    end
    insertcols!(df_total, 3, "VRES_share" =>  vres_share[:])

    #summing up the total of conv generation
    total_conv = Array{Float64}(undef, ip.num_scen)
    for s = 1:ip.num_scen
        total_conv[s] = sum(g_conv[s,:,:,:,:])
    end
    insertcols!(df_total, 4, "conv_generation" => total_conv[:] .* sf)

    # inserting total consumption column
    insertcols!(df_total, 5, "consumption" => total_consump[:] .* sf)
    
    println(io, df_total)
    println(io, "\n")

    println(io, "NODAL SLOPE")
    println(io, "\n")
    for n = 1:ip.num_nodes 
        println(io, "NODE $n:")
        df_node = DataFrame()
        df_node.scenario = 1:ip.num_scen
            
            for t = 1:ip.num_time_periods
                insertcols!(df_node, t+1, "time_$t" => ip.id_slope[:,t,n] )
            end

        println(io, df_node)
        println(io, "\n")
    end

    println(io, "NODAL INTERCEPT")
    println(io, "\n")
    for n = 1:ip.num_nodes 
        println(io, "NODE $n:")
        df_node = DataFrame()
        df_node.scenario = 1:ip.num_scen
            
            for t = 1:ip.num_time_periods
                insertcols!(df_node, t+1, "time_$t" => ip.id_intercept[:,t,n] )
            end

        println(io, df_node)
        println(io, "\n")
    end


    println(io, "NODES' PRICES")
    println(io, "\n")

    for n = 1:ip.num_nodes 
        println(io, "NODE $n")
        df_node = DataFrame()
        df_node.scenario = 1:ip.num_scen
            
            for t = 1:ip.num_time_periods
                insertcols!(df_node, t+1, "time_$t" => ip.id_intercept[:,t,n] .- 0.5 .* ip.id_slope[:,t,n] .* q[:,t,n]./ip.time_periods[t] )
            end

        println(io, df_node)
        println(io, "\n")
        
    end

    println(io, "NODAL REVENUES")
    println(io, "\n")

    for n = 1:ip.num_nodes 
        println(io, "NODE $n")
        df_node = DataFrame()
        df_node.scenario = 1:ip.num_scen
        
        n_rev_s = zeros(ip.num_scen)
        for t = 1:ip.num_time_periods
            n_rev_s = n_rev_s.+(ip.id_intercept[:,t,n] .- 0.5 .* ip.id_slope[:,t,n] .* q[:,t,n]./ip.time_periods[t]).*q[:,t,n]
        end
        insertcols!(df_node, 2, "revenue" => n_rev_s)

        println(io, df_node)
        println(io, "\n")
        
    end


    println(io, "TSO'S EXPENCES")
    println(io, "\n")
    df_TSO = DataFrame()


    # TSO INV EXPENSES
    TSO_inv_exp = sum( - sum( ip.transm.investment_costs[n,m]*0.5*l_plus[n,m] for m in 1:ip.num_nodes) for n in 1:ip.num_nodes)/scaling_factor/obg_scaling_factor

    # TSO OP COSTS

    TSO_op_exp = sum(- sum(
        ip.transm.maintenance_costs[n,m]*
        (ip.transm.installed_capacities[n,m] + 0.5 * l_plus[n,m])
        for m in 1:ip.num_nodes)
    for n in 1:ip.num_nodes)

    insertcols!(df_TSO, 1, "INV COSTS" => TSO_inv_exp )
    insertcols!(df_TSO, 2, "OP COSTS" => TSO_op_exp )
    println(io, df_TSO)
    println(io, "\n")


    println(io, "COMPANIES VRES FIXED COSTS")
    println(io, "\n")

    df_vres_costs = DataFrame()
    df_vres_costs.producer = 1:ip.num_prod
    for n in 1:ip.num_nodes
        node_exp = Array{Float64}(undef,ip.num_prod)
        for i = 1:ip.num_prod
            node_exp[i] = sum( ip.vres.maintenance_costs[n,i,r]*(ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) + ip.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r] for r in 1:ip.num_VRES)
        end
        insertcols!(df_vres_costs, n+1, "node_$n" => node_exp )
    end

    println(io, df_vres_costs)
    println(io, "\n")

    println(io, "COMPANIES CONV FIXED COSTS")
    println(io, "\n")

    df_conv_costs = DataFrame()
    df_conv_costs.producer = 1:ip.num_prod
    for n in 1:ip.num_nodes
        node_exp = Array{Float64}(undef,ip.num_prod)
        for i = 1:ip.num_prod
            node_exp[i] = sum( ip.conv.maintenance_costs[n,i,r]*(ip.conv.installed_capacities[n,i,r] + g_conv_plus[n,i,r]) + ip.conv.investment_costs[n,i,r]*g_conv_plus[n,i,r] for r in 1:ip.num_conv) 
        end
        insertcols!(df_conv_costs, n+1, "node_$n" => node_exp )
    end

    println(io, df_conv_costs)
    println(io, "\n")

    println(io, "COMPANIES CONV VARIABLE COSTS")
    println(io, "\n")

    for i = 1:ip.num_prod
        println(io, "Producer $i")
        println(io, "\n")
        df_conv_costs_var = DataFrame()
        df_conv_costs_var.scenario = 1:ip.num_scen
        for n in 1:ip.num_nodes
            node_exp = Array{Float64}(undef,ip.num_scen)
            for s = 1:ip.num_scen
                node_exp[s] = sum( (ip.conv.operational_costs[n,i,e] + ip.conv.CO2_tax[e])*g_conv[s,t,n,i,e] for t in 1:ip.num_time_periods, e in 1:ip.num_conv)
            end
            insertcols!(df_conv_costs_var, n+1, "node_$n" => node_exp )
        end
        println(io, df_conv_costs_var)
        println(io, "\n")
    end



    close(io)
end

# function to calculate the duality gap for the lower-level problem
function check_duality_gap(bi_level_problem::JuMP.Model, inpa::initial_parameters, market::String)
    
    # saving the optimal values 
    q = value.(bi_level_problem[:q])
    g_conv = value.(bi_level_problem[:g_conv])
    g_conv_plus = value.(bi_level_problem[:g_conv_plus])
    g_VRES = value.(bi_level_problem[:g_VRES])
    g_VRES_plus = value.(bi_level_problem[:g_VRES_plus])

    β_f_1 = value.(bi_level_problem[:β_f_1])
    β_f_2 = value.(bi_level_problem[:β_f_2])
    β_conv = value.(bi_level_problem[:β_conv])
    β_VRES = value.(bi_level_problem[:β_VRES])
    β_up_conv = value.(bi_level_problem[:β_up_conv])
    β_down_conv =value.(bi_level_problem[:β_down_conv])
    β_plus = value.(bi_level_problem[:β_plus])

    #calculating the primal bound (lower bound)

    primal_objective = sum(
        sum( inpa.scen_prob[s] * 

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
            )

        for t in 1:inpa.num_time_periods, s in 1:inpa.num_scen)
    
        -sum( 

            sum( 
                inpa.vres.maintenance_costs[n,i,r]*
                (g_VRES_plus[n,i,r]) 
                + 
                inpa.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
            for r in 1:inpa.num_VRES)
        
            +
            
            sum( 
                inpa.conv.maintenance_costs[n,i,e]*
                ( g_conv_plus[n,i,e]) 
                + 
                inpa.conv.investment_costs[n,i,e]*g_conv_plus[n,i,e]
            for e in 1:inpa.num_conv) 
             
        for i in 1:inpa.num_prod)
    
    for n in 1:inpa.num_nodes)

    dual_objective = (
        sum( inpa.scen_prob[s] * 
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

        + sum(round(inpa.time_periods[t] * inpa.vres.availability_factor[s,t,n,r] * inpa.vres.installed_capacities[n,i,r], digits = 4) * β_VRES[s,t,n,i,r] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, r in 1:inpa.num_VRES)

        + sum(inpa.time_periods[t] * inpa.conv.ramp_up[n,i,e] * inpa.conv.installed_capacities[n,i,e] * β_up_conv[s,t,n,i,e] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, e in 1:inpa.num_conv)

        + sum(inpa.time_periods[t] * inpa.conv.ramp_down[n,i,e] * inpa.conv.installed_capacities[n,i,e] * β_down_conv[s,t,n,i,e] for s in 1:inpa.num_scen, t in 1:inpa.num_time_periods, n in 1:inpa.num_nodes,  i in 1:inpa.num_prod, e in 1:inpa.num_conv)
        
        + sum(inpa.gen_budget[i] * β_plus[i] for i in 1:inpa.num_prod)

    )

    println("UB - LB = $(dual_objective - primal_objective)")
    return dual_objective - primal_objective
    
end