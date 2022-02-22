"""Equivalent annual cost (EAC)
# Arguments
* `cost::Real`: Net present value of a project.
* `n::Integer`: Number of payments.
* `r::Real`: Interest rate.
"""
function equivalent_annual_cost(cost::Float64, n::Float64, r::Float64)
    factor = if iszero(r) n else (1 - (1 + r)^(-n)) / r end
    @show factor
    return cost ./ factor
end


"""
Function hours_for_each_scenario_generation generates the hour-strating and ending value for each of the scenarios
based on the week number send as an input and the set containg the info on the number of days considered in each scenario 

Parameters:
* days                          the array that contains the number of days considered in each scenario
* week::Int                     the week number in the month to be considered for generating the scenario 
"""

function hours_for_each_scenario_generation(days::Array{Int64}, week::Int64)
    output = Array{UnitRange{Int64}}(undef,length(days))
    for s = 1:length(days)
        start_hour = sum(days[1:s-1])*24 + (week-1)*7*24
        end_hour =  start_hour+24*7
        output[s] = start_hour:end_hour
    end
    return output
end


"""
Function lif_slope_and_intercept_generation generates the slope and intercept for the inverse demand function based on the input data. 
It extracts the sets of data (demand and the electricity price) for each of the scenarios (with predefined first and last hours for each of the scenarios). 
Then In each of the extracted sets based on the number of hours (N) dedicated to each time period approximates 
each N points of the set with the linear function and returns the slope and the intercept of it.
Parameters:

* data_src_link::String         A link to the folder containing the hourly-based data on the demand and electiricty prices 
* N_nodes::Int64                Number of nodes
* N_scen::Int64                 Number of scenarios
* T::Array{Int64}               The array that contains the number of hours allocated to each time period 
* hours_to_be_extracted         An array with the range of hours defined for each scenario 
* sf::Float64                   Scalling factor for the model parameters (price in case of this function)

"""

function lif_slope_and_intercept_generation(data_src_link::String, N_nodes::Int64, N_scen::Int64, T::Array{Int64}, hours_to_be_extracted::Array{UnitRange{Int64}}, sf::Float64)

    # accumulate the time periods
    T_accumulated = accumulate(+, T)

    # keep the length of the time periods array for simplifying the notation 
    N_t = length(T)

    # define the structures to keep the slope and the intercept for the inverse demand function
    id_slope = Array{Float64}(undef, N_scen, N_T, N_nodes) 
    id_intercept = Array{Float64}(undef, N_scen, N_T, N_nodes) 

    for n = 1:N_nodes

        # read the correspondent to the node demand and price data
        node_data = (Matrix(DataFrame(CSV.File(data_src_link * "/inverse_demand/" * "demand_and_prices_" * string(n) * ".csv"))[!,2:3]))
        # making sure all the elements are of the type Float (were not read like String)
        #node_data[:,1] = -0.5 .* node_data[:,1]

        for s in 1:N_scen
    
        # extract the hourly data that corresponds to the scenario under consideration
        scenario_data = (node_data[hours_to_be_extracted[s], :])
        @show typeof(scenario_data)
        if typeof(scenario_data) == Array{Any,2}
            #convert(Array{Float64}, scenario_data)
            scenario_data = map(x->parse(Float64,x),string.(scenario_data))
        end
        #scalling
        scenario_data[:,2] = scenario_data[:,2] ./ sf
            for t = 1:N_T

                # caluclate the hours within the scenario that correspond to the time period t
                time_period_indexes = (t==1 ? 0 : T_accumulated[t-1])+1 : T_accumulated[t]
                
                # approximate the data set corresponding to the time period with the polynomial of the 1st degree
                #LinearRegression = np.polyfit(scenario_data[time_period_indexes, 1], scenario_data[time_period_indexes, 2], 1)

                # save the correspodent values
                #id_slope[s,t,n] = LinearRegression[1]
                #id_intercept[s,t,n] = LinearRegression[2]

                epsilon = -0.3 # demand elasticity 
                @show scenario_data[time_period_indexes, 2]
                id_slope[s,t,n] = -mean(scenario_data[time_period_indexes, 2]) * 1 / (epsilon * mean(scenario_data[time_period_indexes, 1]))
                id_intercept[s,t,n] = mean(scenario_data[time_period_indexes, 2]) * 1+ id_slope[s,t,n] * mean(scenario_data[time_period_indexes, 1])
            end
        end
    end

    return id_slope, id_intercept
    
end

"""
Function availability_factor_generation generates availability factors for VRES based on the input data. 
It extracts the sets of data (hurly availability factor at each node for each type of VRES) for each of the scenarios 
(with predefined first and last hours for each of the scenarios). 
Then In each of the extracted sets based on the number of hours (N) dedicated to each time period approximates 
each N points of the set with the expected value of the correspodent availibilty factor data in it.
P
arameters:

* data_src_link::String         A link to the folder containing the hourly-based data on the demand and electiricty prices 
* N_nodes::Int64                Number of nodes
* N_scen::Int64                 Number of scenarios
* T::Array{Int64}               The array that contains the number of hours allocated to each time period 
* N_R::Int64                    Number of the VRES types
* hours_to_be_extracted         An array with the range of hours defined for each scenario 

"""
function availability_factor_generation(data_src_link::String, N_nodes::Int64, N_scen::Int64, T::Array{Int64}, N_R::Int64, hours_to_be_extracted::Array{UnitRange{Int64}})
    
    # accumulate the time periods
    T_accumulated = accumulate(+, T)

    # keep the length of the time periods array for simplifying the notation 
    N_t = length(T)

    # create a structure to keep the availability factor values
    A_VRES = Array{Float64}(undef, N_scen, N_T, N_nodes, N_R)
    
    for n = 1:N_nodes

        # read the correspondent to the node hourly VRES availabiliy data
        node_data = Matrix(DataFrame(CSV.File(data_src_link * "/VRES_generation_units/availability_factor/" * "hourly_data_node_$n.csv"))[!,:])
    
        for r in 1:N_R
        
            for s in 1:N_scen
    
            # extract the hourly data that corresponds to the scenario under consideration
            scenario_data = node_data[hours_to_be_extracted[s], r]
        
                for t = 1:N_T

                    # caluclate the hours within the scenario that correspond to the time period t
                    time_period_indexes = (t==1 ? 0 : T_accumulated[t-1])+1 : T_accumulated[t]

                    # save the correspodent expected value calculated for hourly data on the availability of the 
                    # correspondent VRES
                    A_VRES[s,t,n,r] = sum(scenario_data[time_period_indexes])/T[t]
                
                end
            end
        end
    end

    return A_VRES

end

