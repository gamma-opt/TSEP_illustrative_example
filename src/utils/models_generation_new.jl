function single_level_problem_generation(ip::initial_parameters)    
    # Defining single-level model
    single_level_problem = Model(() -> Gurobi.Optimizer(GRB_ENV))
    #single_level_problem = Model(() -> Ipopt.Optimizer())
    #set_optimizer_attribute(single_level_problem, "OutputFlag", 0)
    set_optimizer_attribute(single_level_problem, "NonConvex", 2)
    #set_optimizer_attribute(single_level_problem, "InfUnbdInfo", 1)
    set_optimizer_attribute(single_level_problem, "Presolve", 0)
    set_optimizer_attribute(single_level_problem, "IntFeasTol", 1E-9)
    set_optimizer_attribute(single_level_problem, "FeasibilityTol", 1E-9)
    set_optimizer_attribute(single_level_problem, "FeasibilityTol", 1E-9)
    set_optimizer_attribute(single_level_problem, "NumericFocus", 2)

    ## VARIABLES

    # Conventional generation related variable
    @variable(single_level_problem, g_conv[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv] >= 0)

    # VRES generation related variable
    @variable(single_level_problem, g_VRES[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES] >= 0)

    # Consumption realted variable
    @variable(single_level_problem, q[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes] >= 0)

    # Energy transmission realted variable
    @variable(single_level_problem, f[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes]) # >= 0)

    # Transmission capacity expansion realted variable
    @variable(single_level_problem, l_plus[ 1:ip.num_nodes, 1:ip.num_nodes]>=0)

    # Conventional energy capacity expansion realted variable
    @variable(single_level_problem, g_conv_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv]>=0)

    # VRES capacity expansion realted variable
    @variable(single_level_problem, g_VRES_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES]>=0)

    ## OBJECTIVE
    @objective(single_level_problem, Max, 
            sum(
                sum( ip.scen_prob[s] * 

                    (ip.id_intercept[s,t,n]*q[s,t,n] 
                    - 0.5* ip.id_slope[s,t,n]* q[s,t,n]^2
                    
                    - sum( 
                        (ip.conv.operational_costs[n,i,e] 
                        + ip.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
                    for i in 1:ip.num_prod, e in 1:ip.num_conv)

                    #- sum(ip.transm.transmissio_costs[n,m]*f[s,t,n,m] 
                    #for n in 1:ip.num_nodes, m in 1:ip.num_nodes)
                    #- 1E-4 * sum( f[s,t,n,m] for m in 1:ip.num_nodes)
                    )

                    

                for t in 1:ip.num_time_periods, s in 1:ip.num_scen)
                
                -sum( 

                    sum( 
                        ip.vres.maintenance_costs[n,i,r]*
                        (ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
                        + 
                        ip.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
                    for r in 1:ip.num_VRES)
                    
                    +
                        
                    sum( 
                        ip.conv.maintenance_costs[n,i,e]*
                        (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
                        + 
                        ip.conv.investment_costs[n,i,e]* g_conv_plus[n,i,e]
                    for e in 1:ip.num_conv) 
                         
                for i in 1:ip.num_prod)
                
                - sum(
                    ip.transm.maintenance_costs[n,m]*
                    (ip.transm.installed_capacities[n,m] + 0.5 * l_plus[n,m])
                    +
                    ip.transm.investment_costs[n,m] * 0.5 * l_plus[n,m]
                for m in 1:ip.num_nodes)
            
            for n in 1:ip.num_nodes)/scaling_factor/obg_scaling_factor
    )

    ## CONSTRAINTS

    @constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes],
        q[s,t,n]/scaling_factor  - sum( 
                    sum(g_conv[s,t,n,i,e] for e in 1:ip.num_conv)    
                    + 
                    sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES) 
                    for i = 1:ip.num_prod) / scaling_factor 
                + sum( f[s,t,n,m] for m in n+1:ip.num_nodes)/scaling_factor - sum(f[s,t,m,n] for m in 1:n-1)/scaling_factor #- sum( 0.98*f[s,t,m,n] for m in 1:n-1)
        == 0
    )
    
    # Transmission bounds 
    @constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        f[s,t,n,m]/scaling_factor - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor <= 0
    )

    @constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        f[s,t,n,m]/scaling_factor + ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor >= 0
    )
    
    # Primal feasibility for the transmission 
    @constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:n],
        f[s,t,n,m]/scaling_factor == 0 
    )

    #@constraint(single_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        #f[s,t,n,m] - sum( 
            #sum(g_conv[s,t,n,i,e] for e in 1:ip.num_conv)    
           #  + 
           # sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES) 
           # for i = 1:ip.num_prod) - sum(f[s,t,m1,n] for m1 in 1:ip.num_nodes)<= 0
    #)

    # Conventional generation bounds
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        g_conv[s,t,n,i,e]/scaling_factor - ip.time_periods[t]*(ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor <= 0
    )

    # Conventional generation budget limits
    @constraint(single_level_problem, [i = 1:ip.num_prod],
        sum(ip.conv.investment_costs[n,i,e] * g_conv_plus[n,i,e] for  e in 1:ip.num_conv, n in 1:ip.num_nodes)/scaling_factor <= ip.conv.budget_limits[i]/scaling_factor
    )

    # VRES generation bounds 
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, r in 1:ip.num_VRES],
        g_VRES[s,t,n,i,r]/scaling_factor - ip.time_periods[t]* ip.vres.availability_factor[s,t,n,r]*(ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r])/scaling_factor <= 0
    )

    # VRES generation budget limits
    @constraint(single_level_problem, [i = 1:ip.num_prod],
        sum(ip.vres.investment_costs[n,i,r] * g_VRES_plus[n,i,r] for  r in 1:ip.num_VRES, n in 1:ip.num_nodes)/scaling_factor <= ip.vres.budget_limits[i]/scaling_factor
    )

    # primal feasibility for the transmission expansion investments
    #@constraint(single_level_problem, [n in 1:ip.num_nodes],
        #l_plus[n,n] == 0
    #)

    @constraint(single_level_problem, [n in 1:ip.num_nodes, m in 1:ip.num_nodes],
        l_plus[n,m]/scaling_factor - l_plus[m,n]/scaling_factor == 0
    )

    # Primal feasibility for the transmission expansion 2 (budget related)
    @constraint(single_level_problem,
        sum(ip.transm.investment_costs[n,m] * l_plus[n,m] for n in 1:ip.num_nodes, m in 1:ip.num_nodes)/scaling_factor - ip.transm.budget_limit/scaling_factor <= 0 
    )



    # Maximum ramp-up rate for conventional units
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
       g_conv[s,t,n,i,e]/scaling_factor - (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor - ip.time_periods[t] * ip.conv.ramp_up[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor <= 0
    )

    # Maximum ramp-down rate for conventional units
    @constraint(single_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor - g_conv[s,t,n,i,e]/scaling_factor - ip.time_periods[t] * ip.conv.ramp_down[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor <= 0
    )
    
    return single_level_problem
end

function bi_level_problem_generation(ip::initial_parameters, market::String)    
    # Defining single-level model
    bi_level_problem = Model(() -> Gurobi.Optimizer(GRB_ENV))
    #bi_level_problem = Model(() -> Ipopt.Optimizer())
    #set_optimizer_attribute(bi_level_problem, "OutputFlag", 0)
    set_optimizer_attribute(bi_level_problem, "NonConvex", 2)
    #set_optimizer_attribute(single_level_problem, "InfUnbdInfo", 1)
    set_optimizer_attribute(bi_level_problem, "Presolve", 0)
    #set_optimizer_attribute(bi_level_problem, "IntFeasTol", 1E-9)
    
    #set_optimizer_attribute(bi_level_problem, "FeasibilityTol", 1E-9)
    set_optimizer_attribute(bi_level_problem, "NumericFocus", 3)
    #set_optimizer_attribute(bi_level_problem, "BarCorrectors", 2000000000)
    #set_optimizer_attribute(bi_level_problem, "BarQCPConvTol", 0.0)
    #set_optimizer_attribute(bi_level_problem, "BarConvTol", 0.0)
    #set_optimizer_attribute(bi_level_problem, "BarHomogeneous", 1)
    #set_optimizer_attribute(bi_level_problem, "Method", 4)

    ## PRIMAL VARIABLES

    # Conventional generation related variable
    @variable(bi_level_problem, g_conv[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv] >= 0)

    # VRES generation related variable
    @variable(bi_level_problem, g_VRES[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES] >= 0)

    # Consumption realted variable
    @variable(bi_level_problem, q[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes]>= 0)

    # Energy transmission realted variable
    @variable(bi_level_problem, f[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes] ) #>= 0)

    # Transmission capacity expansion realted variable
    @variable(bi_level_problem, l_plus[ 1:ip.num_nodes, 1:ip.num_nodes]>=0)

    # Conventional energy capacity expansion realted variable
    @variable(bi_level_problem, g_conv_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv]>=0)

    # VRES capacity expansion realted variable
    @variable(bi_level_problem, g_VRES_plus[ 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES]>=0)

    ## DUAL VARIABLES

    # Shadow price on the power balance
    @variable(bi_level_problem, θ[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes]) #>=0)

    # Shadow price on the power flow primal feasibility constraint
    @variable(bi_level_problem, λ_f[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes])

    # Shadow price on the transmission capacity for the power flow
    #@variable(bi_level_problem, β_f[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes])
    @variable(bi_level_problem, β_f_1[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes] >= 0)
    @variable(bi_level_problem, β_f_2[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_nodes] >= 0)

    # Shadow price on conventional energy capacity
    @variable(bi_level_problem, β_conv[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv] >= 0)

    # Shadow price on vres energy capacity
    @variable(bi_level_problem, β_VRES[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_VRES] >= 0)

    # Shadow price maximum ramp-up rate for the conventional generation
    @variable(bi_level_problem, β_up_conv[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv] >= 0)

    # Shadow price maximum ramp-down rate for the conventional generation
    @variable(bi_level_problem, β_down_conv[1:ip.num_scen, 1:ip.num_time_periods, 1:ip.num_nodes, 1:ip.num_prod, 1:ip.num_conv] >= 0)

    # Shadow price maximum budget limits for convetional generation 
    @variable(bi_level_problem, β_conv_plus[1:ip.num_prod] >= 0)

    # Shadow price maximum budget limits for renewable generation 
    @variable(bi_level_problem, β_vres_plus[1:ip.num_prod] >= 0)

    ## OBJECTIVE

    @objective(bi_level_problem, Max, 
        sum(
            sum( ip.scen_prob[s] * 

                (ip.id_intercept[s,t,n]*q[s,t,n] 
                - 0.5* ip.id_slope[s,t,n]* q[s,t,n]^2

                #- (market == "perfect" ? 0 : 
                    #(0.5* ip.id_slope[s,t,n] * sum( 
                    #    (sum( g_conv[s,t,n,i,e]  for e in 1:ip.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES))^2
                    #for i in 1:ip.num_prod))
                    #)
            
                - sum( 
                    (ip.conv.operational_costs[n,i,e] 
                    + ip.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
                for i in 1:ip.num_prod, e in 1:ip.num_conv)

                #- sum(ip.transm.transmissio_costs[n,m]*f[s,t,n,m] 
                #for n in 1:ip.num_nodes, m in 1:ip.num_nodes)
                #- 1E-1 * sum( f[s,t,n,m] for m in 1:ip.num_nodes)
                )

            for t in 1:ip.num_time_periods, s in 1:ip.num_scen)
        
            -sum( 

                sum( 
                    ip.vres.maintenance_costs[n,i,r]*
                    (ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
                    + 
                    ip.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
                for r in 1:ip.num_VRES)
            
                +
                
                sum( 
                    ip.conv.maintenance_costs[n,i,e]*
                    (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
                    + 
                    ip.conv.investment_costs[n,i,e]*g_conv_plus[n,i,e]
                for e in 1:ip.num_conv) 
                 
                for i in 1:ip.num_prod)
        
            - sum(

                ip.transm.maintenance_costs[n,m]*
                (ip.transm.installed_capacities[n,m] + 0.5 * l_plus[n,m])
                +
                ip.transm.investment_costs[n,m]*0.5*l_plus[n,m]

            for m in 1:ip.num_nodes)
    
        for n in 1:ip.num_nodes)/scaling_factor/obg_scaling_factor
    )

    #@constraint(bi_level_problem,
        #sum(l_plus[n,m] for  n in 1:ip.num_nodes, m in 1:ip.num_nodes) <= 100000
    #)

    #@constraint(bi_level_problem, [n in 1:ip.num_nodes],
        #l_plus[n,n] == 0
    #)

    ## PRIMAL CONSTRAINTS

    # Power balance
    @constraint(bi_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes],
        q[s,t,n]/scaling_factor  - sum( 
                        sum(g_conv[s,t,n,i,e] for e in 1:ip.num_conv)    
                         + 
                        sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES) 
                        for i = 1:ip.num_prod)/scaling_factor 
                        + sum( f[s,t,n,m] for m in n+1:ip.num_nodes)/scaling_factor  - sum(f[s,t,m,n] for m in 1:n-1)/scaling_factor  #- sum( 0.98*f[s,t,m,n] for m in 1:n-1)
        == 0
    )

    # Transmission bounds 
    #@constraint(bi_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
       # f[s,t,n,m] - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m]) <= 0
    #)

    @constraint(bi_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        f[s,t,n,m]/scaling_factor  - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor  <= 0
    )

    @constraint(bi_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes ],
        -f[s,t,n,m]/scaling_factor  - ip.time_periods[t]*(ip.transm.installed_capacities[n,m] + l_plus[n,m])/scaling_factor  <= 0
    )

    # Primal feasibility for the transmission 
    @constraint(bi_level_problem, [s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:n],
        f[s,t,n,m]/scaling_factor  == 0 
    )

    # Primal feasibility for the transmission expansion 1
    @constraint(bi_level_problem, [n in 1:ip.num_nodes, m in 1:ip.num_nodes],
        l_plus[n,m]/scaling_factor  - l_plus[m,n]/scaling_factor  == 0 
    )


    # Primal feasibility for the transmission expansion 2 (buget related)
    @constraint(bi_level_problem,
        sum(ip.transm.investment_costs[n,m] * l_plus[n,m] for n in 1:ip.num_nodes, m in 1:ip.num_nodes)/scaling_factor  - ip.transm.budget_limit/scaling_factor  <= 0 
    )

    # Conventional generation bounds
    @constraint(bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        g_conv[s,t,n,i,e]/scaling_factor  - ip.time_periods[t]*(ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )
    
    # Conventional generation budget limits
    @constraint(bi_level_problem, [ i = 1:ip.num_prod],
        sum(ip.conv.investment_costs[n,i,e] * g_conv_plus[n,i,e] for  e in 1:ip.num_conv, n in 1:ip.num_nodes)/scaling_factor  <= ip.conv.budget_limits[i]/scaling_factor 
    )

    # VRES generation bounds 
    @constraint(bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, r in 1:ip.num_VRES],
        g_VRES[s,t,n,i,r]/scaling_factor  - round(ip.time_periods[t]*ip.vres.availability_factor[s,t,n,r], digits = 4)*(ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r])/scaling_factor  <= 0
    )
    
    # VRES generation budget limits
    @constraint(bi_level_problem, [ i = 1:ip.num_prod],
        sum(ip.vres.investment_costs[n,i,e] * g_VRES_plus[n,i,e] for  e in 1:ip.num_VRES, n in 1:ip.num_nodes)/scaling_factor  <= ip.vres.budget_limits[i]/scaling_factor 
    )


    # Maximum ramp-up rate for conventional units
    @constraint(bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        g_conv[s,t,n,i,e]/scaling_factor  - (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor  - ip.time_periods[t] * ip.conv.ramp_up[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )

    # Maximum ramp-down rate for conventional units
    @constraint(bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        (t == 1 ? 0 : g_conv[s,t-1,n,i,e])/scaling_factor  - g_conv[s,t,n,i,e]/scaling_factor  - ip.time_periods[t] * ip.conv.ramp_down[n,i,e] * (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e])/scaling_factor  <= 0
    )

    ## DUAL CONSTRAINTS
    @constraint(bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes],
        -ip.scen_prob[s] * (ip.id_intercept[s,t,n] - 0.5* ip.id_slope[s,t,n]* q[s,t,n] )/scaling_factor  + θ[s,t,n]/scaling_factor  >= 0  
    )

    #@constraint( bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:ip.num_nodes],
        #-0.02*θ[s,t,n] + β_f[s,t,n,m] + λ_f[s,t,n,m] >= 0  
    #)

    @constraint( bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in n+1:ip.num_nodes],
        θ[s,t,n]/scaling_factor  - θ[s,t,m]/scaling_factor  + β_f_1[s,t,n,m]/scaling_factor  - β_f_2[s,t,n,m]/scaling_factor  == 0  
    )
    
    #@constraint( bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:n-1],
       #-0.98*θ[s,t,n] + β_f_1[s,t,n,m] - β_f_2[s,t,n,m] + λ_f[s,t,n,m] >= 0  
    #)

        
    @constraint( bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m in 1:n],
        β_f_1[s,t,n,m]/scaling_factor  - β_f_2[s,t,n,m]/scaling_factor  + λ_f[s,t,n,m]/scaling_factor  == 0  
    )

    #@constraint( bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, m = n],
       # β_f_1[s,t,n,m] - β_f_2[s,t,n,m] >= 0  
    #)

    @constraint(bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        ip.scen_prob[s] * ( (market == "perfect" ? 0 : (ip.id_slope[s,t,n] * ( sum( g_conv[s,t,n,i,e1]  for e1 in 1:ip.num_conv) + sum(  g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES)))) +  ip.conv.operational_costs[n,i,e] + ip.conv.CO2_tax[e])/scaling_factor 
        - θ[s,t,n]/scaling_factor + β_conv[s,t,n,i,e]/scaling_factor  + β_up_conv[s,t,n,i,e]/scaling_factor  - (t == ip.num_time_periods ? 0 : β_up_conv[s,t+1,n,i,e])/scaling_factor  + (t == ip.num_time_periods ? 0 : β_down_conv[s,t+1,n,i,e])/scaling_factor  - β_down_conv[s,t,n,i,e]/scaling_factor  >= 0

    )
    
    @constraint(bi_level_problem, [ s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes, i = 1:ip.num_prod, r in 1:ip.num_VRES],
        ip.scen_prob[s] * ( (market == "perfect" ? 0 : (ip.id_slope[s,t,n] * ( sum( g_conv[s,t,n,i,e]  for e in 1:ip.num_conv) + sum(  g_VRES[s,t,n,i,r1] for r1 in 1:ip.num_VRES)))))/scaling_factor 
        - θ[s,t,n]/scaling_factor  + β_VRES[s,t,n,i,r]/scaling_factor   >= 0
    )

    @constraint(bi_level_problem, [ n in 1:ip.num_nodes, i = 1:ip.num_prod, e in 1:ip.num_conv],
        ip.conv.maintenance_costs[n,i,e]/scaling_factor  + ip.conv.investment_costs[n,i,e]/scaling_factor + ip.conv.investment_costs[n,i,e]*β_conv_plus[i]/scaling_factor 
        - sum(ip.time_periods[t] * β_conv[s,t,n,i,e] for t in 1:ip.num_time_periods, s in 1:ip.num_scen ) /scaling_factor 
        - sum(ip.time_periods[t] * ip.conv.ramp_up[n,i,e] * β_up_conv[s,t,n,i,e] for t in 1:ip.num_time_periods, s in 1:ip.num_scen )/scaling_factor  
        - sum(ip.time_periods[t] * ip.conv.ramp_down[n,i,e] * β_down_conv[s,t,n,i,e] for t in 1:ip.num_time_periods, s in 1:ip.num_scen )/scaling_factor  
        >= 0
    )

    @constraint(bi_level_problem, [ n in 1:ip.num_nodes, i = 1:ip.num_prod, r in 1:ip.num_VRES],
        ip.vres.maintenance_costs[n,i,r]/scaling_factor + ip.vres.investment_costs[n,i,r]/scaling_factor + ip.vres.investment_costs[n,i,r]*β_vres_plus[i]/scaling_factor 
        - sum(round(ip.time_periods[t] * ip.vres.availability_factor[s,t,n,r], digits = 4) * β_VRES[s,t,n,i,r] for t in 1:ip.num_time_periods, s in 1:ip.num_scen )/scaling_factor 
        >= 0
    )

    ## WEAK DUALITY CONSTRAINT
    @constraint(bi_level_problem, 
        sum(
            sum( ip.scen_prob[s] * 

                (ip.id_intercept[s,t,n]*q[s,t,n] 
                - 0.5* ip.id_slope[s,t,n]* q[s,t,n]^2

                - (market == "perfect" ? 0 : 
                    (0.5* ip.id_slope[s,t,n] * 
                        sum( 
                            (sum( g_conv[s,t,n,i,e]  for e in 1:ip.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES))^2
                        for i in 1:ip.num_prod)
                    )
                )
            
                - sum( 
                    (ip.conv.operational_costs[n,i,e] 
                    + ip.conv.CO2_tax[e])*g_conv[s,t,n,i,e] 
                for i in 1:ip.num_prod, e in 1:ip.num_conv)

                #- sum(ip.transm.transmissio_costs[n,m]*f[s,t,n,m] 
                #for n in 1:ip.num_nodes, m in 1:ip.num_nodes)
                #- 1E-1 * sum( f[s,t,n,m] for m in 1:ip.num_nodes)
                )

            for t in 1:ip.num_time_periods, s in 1:ip.num_scen)
        
            -sum( 

                sum( 
                    ip.vres.maintenance_costs[n,i,r]*
                    (ip.vres.installed_capacities[n,i,r] + g_VRES_plus[n,i,r]) 
                    + 
                    ip.vres.investment_costs[n,i,r]*g_VRES_plus[n,i,r]
                for r in 1:ip.num_VRES)
            
                +
                
                sum( 
                    ip.conv.maintenance_costs[n,i,e]*
                    (ip.conv.installed_capacities[n,i,e] + g_conv_plus[n,i,e]) 
                    + 
                    ip.conv.investment_costs[n,i,e]*g_conv_plus[n,i,e]
                for e in 1:ip.num_conv) 
                 
            for i in 1:ip.num_prod)
        
            #- sum(

                #ip.transm.maintenance_costs[n,m]*
                #(ip.transm.installed_capacities[n,m] + l_plus[n,m])
                #+
                #ip.transm.investment_costs[n,m]*l_plus[n,m]

            #for m in ip.num_nodes)
    
        for n in 1:ip.num_nodes )./scaling_factor 
        -
        (
            sum( ip.scen_prob[s] * 
                ( 0.5* ip.id_slope[s,t,n]*
                    ( q[s,t,n]^2 +  (market == "perfect" ? 0 : 
                                            (sum( 
                                                (sum( g_conv[s,t,n,i,e]  for e in 1:ip.num_conv) + sum(g_VRES[s,t,n,i,r] for r in 1:ip.num_VRES))^2 
                                            for i in 1:ip.num_prod)
                                            )
                                    )
                    )              
                )
            for n in 1:ip.num_nodes, t in 1:ip.num_time_periods, s in 1:ip.num_scen)
            
            + sum(ip.time_periods[t] * ip.transm.installed_capacities[n,m] * (β_f_1[s,t,n,m] + β_f_2[s,t,n,m]) for s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes,  m in 1:ip.num_nodes)
            
            + sum(ip.time_periods[t] * ip.conv.installed_capacities[n,i,e] * β_conv[s,t,n,i,e] for s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes,  i in 1:ip.num_prod, e in 1:ip.num_conv)

            + sum(round(ip.time_periods[t] * ip.vres.availability_factor[s,t,n,r] * ip.vres.installed_capacities[n,i,r], digits = 4) * β_VRES[s,t,n,i,r] for s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes,  i in 1:ip.num_prod, r in 1:ip.num_VRES)

            + sum(ip.time_periods[t] * ip.conv.ramp_up[n,i,e] * ip.conv.installed_capacities[n,i,e] * β_up_conv[s,t,n,i,e] for s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes,  i in 1:ip.num_prod, e in 1:ip.num_conv)

            + sum(ip.time_periods[t] * ip.conv.ramp_down[n,i,e] * ip.conv.installed_capacities[n,i,e] * β_down_conv[s,t,n,i,e] for s in 1:ip.num_scen, t in 1:ip.num_time_periods, n in 1:ip.num_nodes,  i in 1:ip.num_prod, e in 1:ip.num_conv)
            
            + sum(ip.vres.budget_limits[i] * β_vres_plus[i] for i in 1:ip.num_prod)

            + sum(ip.conv.budget_limits[i] * β_conv_plus[i] for i in 1:ip.num_prod)

        )./scaling_factor 
        >= 0
    )
    
    return bi_level_problem
end