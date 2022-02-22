
"""
VRES_parameters
Stores VRES realted parameters. Has the following fields:

* installed_capacities::Array{Float64}              Installed capacties at each node from each producer for each type of the VRES sources (MW)
* maintenance_costs::Array{Float64}                 Maintenace costs at each node from each producer for each type of the VRES sources (€/MW) 
* investment_costs::Array{Float64}                  Capacity expansion investment costs at each node from each producer for each type of the VRES sources (€/MW) 
* budget_limits                                     Capacity expnasion budget limits for each of the producers 
                                                    at each of the nodes (€)
* availability_factor::Array{Float64}               Availability factor for each scenario and each time period at each node for each type of the VRES sources 

"""
mutable struct VRES_parameters
    installed_capacities::Array{Float64}
    maintenance_costs::Array{Float64}
    investment_costs::Array{Float64}
    budget_limits::Array{Float64}
    availability_factor::Array{Float64}
end

"""
conventional_generation_parameters
Stores conventional energy sources related parameters. Has the following fields:

* installed_capacities::Array{Float64}              Installed capacties at each node from each producer (MW)
* maintenance_costs::Array{Float64}                 Maintenace costs at each node from each producer (€/MW) 
* investment_costs::Array{Float64}                  Capacity expansion investment costs at each node from each producer (€/MW) 
* budget_limits                                     Capacity expnasion budget limits for each of the producers 
                                                    at each of the nodes (€)
* operational_costs::Array{Float64}                 Operational costs at each node from each producer (€/MWh)
* ramp_up::Array{Float64}                           Maximum ramp-up rate
* ramp_down::Array{Float64}                         Maximum ramp-down rate
* CO2_tax::Array{Float64}                           Carbon tax for conventional generation (€/MWh)
"""
mutable struct conventional_generation_parameters
    installed_capacities::Array{Float64}
    maintenance_costs::Array{Float64}
    investment_costs::Array{Float64}
    budget_limits::Array{Float64}
    operational_costs::Array{Float64}
    ramp_up::Array{Float64}
    ramp_down::Array{Float64}
    CO2_tax::Array{Float64}
end

"""
transmission_parameters
Stores energy transmission realted parameters. Has the following fields:

* installed_capacities::Array{Float64}              Installed capacities for each of the nodes (MW)
* maintenance_costs::Array{Float64}                 Maintenace costs for each of the nodes (€/MW)
* investment_costs::Array{Float64}                  Capacity expansion costs for each node (€/MW)
* budget_limit::Float64                             Capacity expnasion budget limit (€)

"""
mutable struct transmission_parameters
    installed_capacities::Array{Float64}
    maintenance_costs::Array{Float64}
    investment_costs::Array{Float64}
    budget_limit::Float64
end

"""
initial_parameters
Stores attributes for generating JuMP model. Has the following fields:

* num_scen::Int                                     Number of scenarios used in the model
* num_nodes::Int                                    Number of nodes in the model
* num_time_periods::Int                             Number of time periods used in the model 
* num_VRES::Int                                     Number of types of VRES sources
* num_conv::Int                                     Number of types of conventional sources
* num_prod::Int                                     Number of energy producer 
* scen_prob::Vector{Float64}                        Vector of probabilities associated with each scenario
* time_periods::Vector{Float64}                     Vector of the time periods values (h)
* id_slope::Array{Float64}                          Slope of the inverse demand function
* id_intercept::Array{Float64}                      Intercept of the inverse demand function
* vres::VRES_parameters                             VRES related parameters
* conv::conventional_generation_parameters          Conventional generation related parameters 
* transm::transmission_parameters                   Transmission lines related paramters

"""

mutable struct initial_parameters
    # Indices and sets 
    num_scen::Int
    num_nodes::Int
    num_time_periods::Int
    num_VRES::Int
    num_conv::Int
    num_prod::Int

    scen_prob::Vector{Float64}
    time_periods::Vector{Float64}

    # slope of the inverse demand function
    id_slope::Array{Float64}

    # intercept of the inverse demand function
    id_intercept::Array{Float64}

    # VRES parameters 
    vres::VRES_parameters

    # Conventional generation parameters
    conv::conventional_generation_parameters 

    # Transmission related parameters
    transm::transmission_parameters
end
