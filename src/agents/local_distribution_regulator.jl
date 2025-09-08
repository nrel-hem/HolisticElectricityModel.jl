struct LocalDistributionRegulatorOptions{T <: RateDesign, U <: NetMeteringPolicy} <: AbstractRegulatorOptions
    rate_design::T
    net_metering_policy::U
    tou_suffix::AbstractString

    planning_reserve_margin::AbstractFloat
    allowed_return_on_investment::AbstractFloat

    # other options as needed
end

function LocalDistributionRegulatorOptions(
    rate_design::RateDesign, 
    net_metering_policy::NetMeteringPolicy;
    tou_suffix::AbstractString = "NE2025",
    planning_reserve_margin::AbstractFloat = 0.12, 
    allowed_return_on_investment::AbstractFloat = 0.112
)
    return LocalDistributionRegulatorOptions(rate_design, net_metering_policy, tou_suffix, planning_reserve_margin, allowed_return_on_investment)
end


mutable struct LocalDistributionRegulator <: AbstractRegulator
    id::String
    current_year::Symbol

    index_rate_tou::Dimension

    # Cost parameters
    distribution_cost::ParamArray
    administration_cost::ParamArray
    transmission_cost::ParamArray
    interconnection_cost::ParamArray
    system_cost::ParamArray
    "other cost not related to the optimization problem"
    othercost::ParamArray
    "revenue requirement by year for each component"
    revenue_req_my::ParamArray

    # other fields as needed
end

function LocalDistributionRegulator(input_dir::AbstractString, model_data::HEMData; id = DEFAULT_ID)
    # change this to read from input files
    index_rate_tou = Dimension("index_rate_tou", [:peak, :non_peak]; prose_name="index_rate_tou", description="index for time-of-use rates")

    distribution_cost = ParamArray("distribution_cost", (model_data.index_z_local, model_data.index_y), zeros((length(model_data.index_z_local), length(model_data.index_y))) )
    administration_cost = ParamArray("administration_cost", (model_data.index_z_local, model_data.index_y), zeros((length(model_data.index_z_local), length(model_data.index_y))) )
    transmission_cost = ParamArray("transmission_cost", (model_data.index_z_local, model_data.index_y), zeros((length(model_data.index_z_local), length(model_data.index_y))) )
    interconnection_cost = ParamArray("interconnection_cost", (model_data.index_z_local, model_data.index_y), zeros((length(model_data.index_z_local), length(model_data.index_y))) )
    system_cost = ParamArray("system_cost", (model_data.index_z_local, model_data.index_y), zeros((length(model_data.index_z_local), length(model_data.index_y))) )
    othercost = ParamArray("othercost", (model_data.index_z_local, model_data.index_y), zeros((length(model_data.index_z_local), length(model_data.index_y))) )
    # from stage 1 results
    revenue_req_my = ParamArray("revenue_req_my", (model_data.index_z_local, model_data.index_y), zeros((length(model_data.index_z_local), length(model_data.index_y))) )
    return LocalDistributionRegulator(id, first(model_data.index_y), index_rate_tou, distribution_cost, administration_cost, transmission_cost, interconnection_cost, system_cost, othercost, revenue_req_my)
end

get_id(x::LocalDistributionRegulator) = x.id


function solve_agent_problem!(
    regulator::LocalDistributionRegulator,
    regulator_opts::LocalDistributionRegulatorOptions,
    model_data::HEMData,
    hem_opts::HEMOptions{LocalDistributionAndDER},
    agent_store::AgentStore,
    w_iter,
    window_length,
    jump_model,
    export_file_path,
    update_results::Bool,
    output_intermediate_results::Bool
)
    # Implement the logic here
    @info("Solving problem for Local Distribution Regulator: $(regulator.id)")
    return 0.0
end

function save_results(
    regulator::LocalDistributionRegulator,
    regulator_opts::LocalDistributionRegulatorOptions,
    hem_opts::HEMOptions,
    export_file_path::AbstractString)

    @info("Saving results for Local Distribution Regulator: $(regulator.id)")
    # Implement the logic to save results here
end