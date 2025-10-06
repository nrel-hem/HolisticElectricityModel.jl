struct LocalDistributionRateMakerOptions{T <: RateDesign, U <: NetMeteringPolicy} <: AbstractRegulatorOptions
    rate_design::T
    net_metering_policy::U
    tou_suffix::AbstractString

    planning_reserve_margin::AbstractFloat
    allowed_return_on_investment::AbstractFloat

    # other options as needed
end

function LocalDistributionRateMakerOptions(
    rate_design::RateDesign, 
    net_metering_policy::NetMeteringPolicy;
    tou_suffix::AbstractString = "NE2025",
    planning_reserve_margin::AbstractFloat = 0.12, 
    allowed_return_on_investment::AbstractFloat = 0.112
)
    return LocalDistributionRateMakerOptions(rate_design, net_metering_policy, tou_suffix, planning_reserve_margin, allowed_return_on_investment)
end


mutable struct LocalDistributionRateMaker <: AbstractRegulator
    id::String
    current_year::Symbol

    index_rate_tou::Dimension

    # index from stage 1 results
    index_cost_type::Dimension
    index_demand_type::Dimension
    index_green_tech_option::Dimension

    # parameters from stage 1 results
    "cost allocation by year, zone, customer type, and cost type"
    cost_allocation_my::ParamArray
    "net demand (peak and actual) by year, zone, customer type, and green tech option"
    net_demand_peak_my::ParamArray
    "revenue requirement by year and zone"
    revenue_req_my::ParamArray

    # other fields as needed
end

function LocalDistributionRateMaker(input_dir::AbstractString, model_data::HEMData; id = DEFAULT_ID)
    # change this to read from input files
    index_rate_tou = Dimension("index_rate_tou", [:peak, :non_peak]; prose_name="index_rate_tou", description="index for time-of-use rates")

    # read from stage 1 results
    index_cost_type = Dimension("index_cost_type", [:Administration, :Distribution, :Interconnection, :Transmission, :System, :DERA]; prose_name="index_cost_type", description="types of costs")
    index_demand_type = Dimension("index_demand_type", [:Peak, :Actual]; prose_name="index_demand_type", description="types of demand")
    index_green_tech_option = Dimension("index_green_tech_option", [:WithGreenTech, :WithoutGreenTech]; prose_name="index_green_tech_option", description="with or without green technology")

    cost_allocation_my = ParamArray(
        "cost_allocation_my", 
        (model_data.index_y, model_data.index_z, model_data.index_h, index_cost_type), 
        zeros((length(model_data.index_y), length(model_data.index_z), length(model_data.index_h), length(index_cost_type))) 
    )

    net_demand_peak_my = ParamArray(
        "net_demand_peak_my", 
        (model_data.index_y, model_data.index_z, model_data.index_h, index_demand_type, index_green_tech_option), 
        zeros((length(model_data.index_y), length(model_data.index_z), length(model_data.index_h), length(index_demand_type), length(index_green_tech_option))) 
    )

    revenue_req_my = ParamArray("revenue_req_my", (model_data.index_z, model_data.index_y), zeros((length(model_data.index_z), length(model_data.index_y))) )
    return LocalDistributionRateMaker(
        id, 
        first(model_data.index_y), 
        index_rate_tou, 
        index_cost_type, 
        index_demand_type, 
        index_green_tech_option, 
        cost_allocation_my,
        net_demand_peak_my,
        revenue_req_my
    )
end

get_id(x::LocalDistributionRateMaker) = x.id


function solve_agent_problem!(
    rate_maker::LocalDistributionRateMaker,
    rate_maker_opts::LocalDistributionRateMakerOptions,
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
    @info("Solving problem for Local Distribution Rate Maker: $(rate_maker.id)")

    # get the distribution upgrade costs from DistributionUtility

    # calculate retail rates based on cost allocation, net demand and revenue requirement

    return 0.0
end

function save_results(
    rate_maker::LocalDistributionRateMaker,
    rate_maker_opts::LocalDistributionRateMakerOptions,
    hem_opts::HEMOptions,
    export_file_path::AbstractString)

    @info("Saving results for Local Distribution Rate Maker: $(rate_maker.id)")
    # Implement the logic to save results here
end