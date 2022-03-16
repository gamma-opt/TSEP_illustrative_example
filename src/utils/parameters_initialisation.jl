## defining parameters 

## Constants
const interest_rate = 0.05
 
## Sets and indices 

# Nodes 
Nodes = collect(skipmissing(DataFrame(CSV.File(data_src_link  * "/sets.csv"))[!,1]))
N_nodes = length(Nodes)

# Days allocted to each scenariio 
S_days = collect(skipmissing(DataFrame(CSV.File(data_src_link  * "/sets.csv"))[!,2]))

# Scenarios probabilities
S_prob = collect(skipmissing(DataFrame(CSV.File(data_src_link  * "/sets.csv"))[!,3]))
N_scen = length(S_prob)

# Time periods 
T = collect(skipmissing(DataFrame(CSV.File(data_src_link  * "/sets.csv"))[!,4]))
N_T = length(T)

# VRES sources 
R = collect(skipmissing(DataFrame(CSV.File(data_src_link  * "/sets.csv"))[!,5]))
N_R = length(R)

# Conventional sources 
E = collect(skipmissing(DataFrame(CSV.File(data_src_link  * "/sets.csv"))[!,6]))
N_E = length(E)

# Producers 
I = collect(skipmissing(DataFrame(CSV.File(data_src_link  * "/sets.csv"))[!,7]))
N_I = length(I)


# variable correspondent to the representative week number if each scenario
# that is used to define the characterisitcs of the scenario 
representative_week = 2


## Demand 


# Slope and the intercept of the inverse demand functoion
#id_slope, id_intercept = lif_slope_and_intercept_generation(data_src_link, N_nodes, N_scen, T, hours_for_each_scenario, scaling_factor)
epsilon = -0.3 # demand elasticity 
id_slope = Array{Float64}(undef, N_scen, N_T, N_nodes) 
id_intercept = Array{Float64}(undef, N_scen, N_T, N_nodes) 
for n = 1:N_nodes 
    node_data = (Matrix(DataFrame(CSV.File(data_src_link * "/inverse_demand/" * "demand_and_prices_" * string(n) * ".csv"))[!,2:3]))
    for s = 1:N_scen
        @show sum(S_days[1:s])-S_days[s]+1
        scenario_data = node_data[sum(S_days[1:s])-S_days[s]+1 : sum(S_days[1:s]), :]
        @show scenario_data[:, 1]
        for t = 1:N_T
            id_slope[s,t,n] = - round(scenario_data[t, 2] * 1 / (epsilon * scenario_data[t, 1]), digits = 4)
            id_intercept[s,t,n] =  round(scenario_data[t, 2] * 1 + id_slope[s,t,n] * scenario_data[t, 1], digits = 4)
        end
    end
end
#id_slope =  id_slope.*1000
#id_intercept = id_intercept.*1000

## Transmission lines parameters 

# Installed capacities 
L_max = Matrix(DataFrame(CSV.File(data_src_link * "/transmission_lines/" * "lines_capacity.csv"))[!,:])

#scalling 
#L_max =  L_max./scaling_factor

# Costs
Costs_lines = Matrix(DataFrame(CSV.File(data_src_link * "/transmission_lines/" * "lines_costs.csv"))[!,:])

# Converter costs 
Converter_costs_lines = Matrix(DataFrame(CSV.File(data_src_link * "/transmission_lines/" * "lines_converter_costs.csv"))[!,:])

# Distance
Distance_lines = Matrix(DataFrame(CSV.File(data_src_link * "/transmission_lines/" * "lines_distances.csv"))[!,:])

# Lifetime
Lifetime_lines = Matrix(DataFrame(CSV.File(data_src_link * "/transmission_lines/" * "lines_lifetime.csv"))[!,:])

# Investment costs (considering the expansion)
I_lines = equivalent_annual_cost.(float.(Costs_lines .* Distance_lines .+ Converter_costs_lines),
float.(Lifetime_lines), interest_rate) 
I_lines = replace(I_lines, NaN=>0)

I_lines = round.(I_lines,digits = 4)
#I_lines = I_lines./100

#scalling
#I_lines = I_lines./scaling_factor
#I_lines = zeros(size(I_lines))

# Budget limits (considering the expansion for each line) 
B_lines = Matrix(DataFrame(CSV.File(data_src_link * "/transmission_lines/" * "budget_limits.csv"))[!,:])[1]

#scalling
#B_lines = B_lines./scaling_factor

#@show I_lines
# Maintenance costs 
M_lines = Matrix(DataFrame(CSV.File(data_src_link * "/transmission_lines/" * "lines_maintenance_costs.csv"))[!,:]) .* I_lines
#scalling
M_lines = M_lines./scaling_factor
#TRANSM = transmission_parameters(round.(L_max, digits = max_digits), round.(M_lines, digits = max_digits), round.(I_lines, digits = max_digits), round.(B_lines, digits = max_digits))
TRANSM = transmission_parameters(L_max, M_lines, I_lines, B_lines)


## Budget limits for generation expansion (for each producer)
GEN_BUDGET = Array(DataFrame(CSV.File(data_src_link * "/generation_budget_limits/budget_limits.csv"))[!,1:N_I])

## Conventional generation

# Maximum generation capacity at each node for conventional energy
G_max_E = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    G_max_E[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/generation_capacities/" * "generation_capacities_producer_"* string(i)* ".csv"))[!,2:N_E+1])
end
#scaling 
#G_max_E = G_max_E./scaling_factor

# Maintenace costs at each node for conventional units 
M_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    M_conv[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/fixed_maintenance_costs/" * "fixed_maintenance_costs_producer_"* string(i)* ".csv"))[!,2:N_E+1]) .* 1000
end
#scalling
#M_conv = M_conv./scaling_factor

# Fuel costs at each node for conventional units 
Fuel_costs_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    Fuel_costs_conv[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/fuel_costs/" * "fuel_costs_producer_"* string(i)* ".csv"))[!,2:N_E+1])
end

# Technology efficiency at each node for conventional units
Efficiency_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    Efficiency_conv[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/efficiency/" * "efficiency_producer_"* string(i)* ".csv"))[!,2:N_E+1])
end  

# Variable maintenance costs at each node for conventional units
Var_m_costs_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    Var_m_costs_conv[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/variable_maintenance_costs/" * "variable_maintenance_costs_producer_"* string(i)* ".csv"))[!,2:N_E+1])
end  

# Operational costs at each node for conventional units 
C_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
C_conv = Fuel_costs_conv ./ Efficiency_conv .+ Var_m_costs_conv 
#scalling
#C_conv = C_conv./scaling_factor

# Investment costs at each node for conventional units 
Investment_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    Investment_conv[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/investment_costs/" * "investment_costs_producer_"* string(i)* ".csv"))[!,2:N_E+1])
end



# lifetime at each node for conventional units 
Lifetime_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    Lifetime_conv[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/lifetime/" * "lifetime_producer_"* string(i)* ".csv"))[!,2:N_E+1])
end

# Annualised investment costs at each node for conventional units 
I_conv= Array{Float64}(undef, N_nodes, N_I, N_E )
I_conv = equivalent_annual_cost.(Investment_conv .* 1000, Lifetime_conv, interest_rate) 

I_conv = round.(I_conv, digits = 4)
#scalling 
#I_conv = I_conv./scaling_factor

# Maximum ramp-up rate for conventional units 
R_up_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    R_up_conv[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/ramp_up/" * "ramp_up_producer_"* string(i)* ".csv"))[!,2:N_E+1])
end
#scalling 
#R_up_conv = R_up_conv./scaling_factor

# Maximum ramp-down rate for conventional units 
R_down_conv = Array{Float64}(undef, N_nodes, N_I, N_E )
for i in 1:N_I
    R_down_conv[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/ramp_down/" * "ramp_down_producer_"* string(i)* ".csv"))[!,2:N_E+1])
end
#scalling 
#R_down_conv = R_down_conv./scaling_factor

# Carbon tax for conventional units 
CO2_tax = Matrix(DataFrame(CSV.File(data_src_link * "/conventional_generation_units/" * "carbon_tax.csv"))[!,2:N_E+1])
#CO2_tax = CO2_tax .* 10
#scalling 
#CO2_tax = CO2_tax./scaling_factor
#CONV = conventional_generation_parameters(round.(G_max_E, digits = max_digits), round.(M_conv, digits = max_digits), round.(I_conv, digits = max_digits), round.(B_conv, digits = max_digits), round.(C_conv, digits = max_digits), round.(R_up_conv, digits = max_digits), round.(R_down_conv, digits = max_digits), round.(CO2_tax, digits = max_digits))
CONV = conventional_generation_parameters(G_max_E, M_conv, I_conv, C_conv, R_up_conv, R_down_conv, CO2_tax)

## VRES generation

# Maximum generation capacity at each node for VRES
G_max_VRES = Array{Float64}(undef, N_nodes, N_I, N_R )
for i in 1:N_I
    G_max_VRES[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/VRES_generation_units/installed_generation_capacities/" * "installed_generation_capacities_producer_"* string(i)* ".csv"))[!,2:N_R+1])
end

# Maintenace costs at each node for VRES
M_VRES = Array{Float64}(undef, N_nodes, N_I, N_R )
for i in 1:N_I
    M_VRES[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/VRES_generation_units/maintenance_costs/" * "maintenance_costs_producer_"* string(i)* ".csv"))[!,2:N_R+1]) .* 1000
end
#M_VRES = M_VRES./100
#scalling
#M_VRES = M_VRES./scaling_factor

# Investment costs at each node for VRES units 
Investment_VRES = Array{Float64}(undef, N_nodes, N_I, N_R )
for i in 1:N_I
    Investment_VRES[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/VRES_generation_units/investment_costs/" * "investment_costs_producer_"* string(i)* ".csv"))[!,2:N_R+1])
end
#Investment_VRES = Investment_VRES./100

# lifetime at each node for VRES units 
Lifetime_VRES = Array{Float64}(undef, N_nodes, N_I, N_R)
for i in 1:N_I
    Lifetime_VRES[:,i,:] = Matrix(DataFrame(CSV.File(data_src_link * "/VRES_generation_units/lifetime/" * "lifetime_producer_"* string(i)* ".csv"))[!,2:N_R+1])
end

# Annualised investment costs at each node for VRES units 
I_VRES= Array{Float64}(undef, N_nodes, N_I, N_R )
I_VRES = equivalent_annual_cost.(Investment_VRES .* 1000, Lifetime_VRES, interest_rate) 
#scalling 
#I_VRES = I_VRES./scaling_factor

I_VRES = round.(I_VRES, digits = 4)

#reading the data about vres investemnt costs incentives
incentives = Array(DataFrame(CSV.File(data_src_link * "/VRES_generation_units/incentives.csv"))[!,2])
for n = 1:N_nodes
    I_VRES[n,:,:] = I_VRES[n,:,:].*(1-incentives[n]/100)
end

# Availability factor 
# create a structure to keep the availability factor values
A_VRES = Array{Float64}(undef, N_scen, N_T, N_nodes, N_R)
#A_VRES = availability_factor_generation(data_src_link, N_nodes, N_scen, T, N_R, hours_for_each_scenario)
for n = 1:N_nodes 
    node_data = Matrix(DataFrame(CSV.File(data_src_link * "/VRES_generation_units/availability_factor/" * "hourly_data_node_$n.csv"))[!,2:end])
    for s = 1:N_scen
        @show sum(S_days[1:s])-S_days[s]+1
        scenario_data = node_data[sum(S_days[1:s])-S_days[s]+1 : sum(S_days[1:s]), :]
        for t = 1:N_T
            for r in 1:N_R
                A_VRES[s,t,n,r] = scenario_data[t,r]
            end
        end
    end
end
#scaling
#A_VRES = A_VRES./scaling_factor


#VRES =  VRES_parameters(round.(G_max_VRES, digits = max_digits), round.(M_VRES,digits = max_digits), round.(I_VRES, digits = max_digits), round.(B_VRES, digits = max_digits), round.(A_VRES, digits = max_digits))
VRES =  VRES_parameters(G_max_VRES, M_VRES, I_VRES, A_VRES)

input_parameters = initial_parameters(N_scen, N_nodes, N_T, N_R, N_E, N_I, S_prob, T, id_slope, id_intercept, GEN_BUDGET, VRES, CONV, TRANSM)