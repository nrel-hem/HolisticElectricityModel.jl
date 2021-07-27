# This file defines data and functions associated with the customer.

mutable struct PVAdoptionModel
    Shape::ParamArray
    MeanPayback::ParamArray
    Bass_p::ParamVector
    Bass_q::ParamVector

    Rate::ParamArray
end

"""
Constructs a PVAdoptionModel by computing Rate from the other parameters.
"""
function PVAdoptionModel(Shape, 
                         MeanPayback, 
                         Bass_p, 
                         Bass_q)
    return PVAdoptionModel(
        Shape,
        MeanPayback,
        Bass_p,
        Bass_q,
        ParamArray("Rate", Shape.dims, Dict(key => Shape[key]/MeanPayback[key] for key in keys(Shape))) # Rate
    )    
end

mutable struct Customers <: Agent
    # Sets
    index_m::Set1D # behind-the-meter technologies

    # Parameters
    "number of customers of type h"
    gamma::ParamVector
    "demand (MWh per representative agent per hour)"
    d::ParamArray
    x_DG_E::ParamArray
    Opti_DG::ParamArray
    "DER generation by a representative customer h and DER technology m"
    DERGen::ParamArray
    CapEx_DG::ParamArray
    FOM_DG::ParamArray
    rho_DG::ParamArray
    "Annualization factor for net consumer surplus of PV installation"
    delta::ParamScalar

    # Primal Variables
    x_DG_new::ParamArray

    # Auxiliary Variables
    Payback::ParamArray
    MarketShare::ParamArray
    MaxDG::ParamArray
    F::ParamArray
    year::ParamArray
    A::ParamArray
    ConPVNetSurplus::ParamArray

    pv_adoption_model::PVAdoptionModel
end

function Customers(input_filename::AbstractString, model_data::HEMData)
    index_m = read_set(input_filename, "index_m", "index_m",
                       prose_name = "behind-the-meter technologies m")

    gamma = read_param("gamma", input_filename, "Gamma", model_data.index_h, 
        description = "number of customers of type h")
    x_DG_E = read_param("x_DG_E", input_filename, "ExistingDER", index_m, 
        row_indices = [model_data.index_h], description = "existing DG capacity")
    Opti_DG = read_param("Opti_DG", input_filename, "OptimalDER", index_m, 
        row_indices = [model_data.index_h])
    rho_DG = read_param("rho_DG", input_filename, "AvailabilityDER", 
        model_data.index_t, row_indices = [model_data.index_h, index_m])
    # Define total DER generation per individual customer per hour
    DERGen = initialize_param("DERGen", model_data.index_h, model_data.index_t, value=1.0)
    for h in model_data.index_h, t in model_data.index_t
        if sum(rho_DG[h,m,t]*Opti_DG[h,m] for m in index_m) != 0.0
            DERGen[h,t] = sum(rho_DG[h,m,t]*Opti_DG[h,m] for m in index_m)
        else
            DERGen[h,t] = 1.0
        end
    end

    pv_adoption_model = PVAdoptionModel(
        initialize_param("Shape", model_data.index_h, index_m, value=1.7), # Shape
        initialize_param("MeanPayback", model_data.index_h, index_m, value=8.8), # MeanPayback
        ParamVector("Bass_p", model_data.index_h, 
                    Dict(:Residential => 7.7E-07, :Commercial => 6.0E-04, :Industrial => 6.0E-04)), # Bass_p
        ParamVector("Bass_q", model_data.index_h, 
                    Dict(:Residential => 0.663, :Commercial => 0.133, :Industrial => 0.133)), # Bass_q
    )
    
    return Customers(
        index_m,
        gamma,
        read_param("d", input_filename, "Demand", model_data.index_t, 
            row_indices = [model_data.index_h]),
        x_DG_E,
        Opti_DG,
        DERGen,
        read_param("CapEx_DG", input_filename, "CapExDER", index_m, 
            row_indices = [model_data.index_h]),
        read_param("FOM_DG", input_filename, "FOMDER", index_m, 
            row_indices = [model_data.index_h]),
        rho_DG,
        ParamScalar("delta", 0.05),
        initialize_param("x_DG_new", model_data.index_h, index_m),
        initialize_param("Payback", model_data.index_h, index_m),
        initialize_param("MarketShare", model_data.index_h, index_m),
        initialize_param("MaxDG", model_data.index_h, index_m),
        initialize_param("F", model_data.index_h, index_m),
        initialize_param("year", model_data.index_h, index_m),
        initialize_param("A", model_data.index_h, index_m),
        initialize_param("ConPVNetSurplus", model_data.index_h, index_m),
        pv_adoption_model
    )
end

function solve_agent_problem(
        customers::Customers, 
        customers_opts::AgentOptions,
        model_data::HEMData, 
        hem_opts::HEMOptions,
        other_agents::Vector{Agent})

    regulator = get_agent(other_agents, Regulator)

    x_DG_before = copy(customers.x_DG_new)

    adopt_model = customers.pv_adoption_model

    # Calculate payback period of DER
    # The NetProfit represents the energy saving/credit per representative agent per DER technology, assuming the optimal DER technology size
    NetProfit = Dict((h,m) =>
        # value of distributed generation (offset load)
        sum(model_data.omega[t] * regulator.p[h,t] * min(customers.d[h,t], sum(customers.rho_DG[h,m,t] * customers.Opti_DG[h,m] for m in customers.index_m)) *
        (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] for t in model_data.index_t) +
        # value of distributed generation (excess generation)
        sum(model_data.omega[t] * regulator.p_ex[h,t] * max(0, sum(customers.rho_DG[h,m,t] * customers.Opti_DG[h,m] for m in customers.index_m) - customers.d[h,t]) *
        (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] for t in model_data.index_t) -
        # cost of distributed generation 
        customers.FOM_DG[h,m] * customers.Opti_DG[h,m] for h in model_data.index_h, m in customers.index_m)
    
    for h in model_data.index_h, m in customers.index_m
        if NetProfit[h,m] >= 0.0
            customers.Payback[h,m] =
                customers.CapEx_DG[h,m] * customers.Opti_DG[h,m] / NetProfit[h,m]
            # Calculate maximum market share and maximum DG potential (based on WTP curve)
            customers.MarketShare[h,m] = 
                1.0 - Distributions.cdf(Distributions.Gamma(adopt_model.Shape[h,m], 1/adopt_model.Rate[h,m]), customers.Payback[h,m])
            customers.MaxDG[h,m] = 
                customers.MarketShare[h,m]*customers.gamma[h]*customers.Opti_DG[h,m]
            # Calculate the percentage of existing DER (per agent type per DER technology) as a fraction of maximum DG potential
            customers.F[h,m] =
                min(customers.x_DG_E[h,m]/customers.MaxDG[h,m], 1.0)
            # Back out the reference year of DER based on the percentage of existing DER
            customers.year[h,m] =
                -log((1-customers.F[h,m])/(customers.F[h,m]*adopt_model.Bass_q[h]/adopt_model.Bass_p[h]+1))/(adopt_model.Bass_p[h]+adopt_model.Bass_q[h])
            # Calculate incremental DG build
            customers.A[h,m] =
                (1.0 - exp(-(adopt_model.Bass_p[h] + adopt_model.Bass_q[h]) * (customers.year[h,m]+1))) / 
                (1.0 + (adopt_model.Bass_q[h]/adopt_model.Bass_p[h]) * 
                    exp(-(adopt_model.Bass_p[h] + adopt_model.Bass_q[h])*(customers.year[h,m]+1)))
            customers.x_DG_new[h,m] =
                max(0.0, customers.A[h,m]*customers.MaxDG[h,m]-customers.x_DG_E[h,m])
        else
            customers.x_DG_new[h,m] = 0.0
        end
    end

    @info "Original new DG" x_DG_before
    @info "New new DG" customers.x_DG_new

    return compute_difference_one_norm([
        (x_DG_before, customers.x_DG_new)
    ])
end

function save_results(
        customers::Customers, 
        customers_opts::AgentOptions,
        hem_opts::HEMOptions,
        exportfilepath::AbstractString, 
        fileprefix::AbstractString)

    # Primal Variables
    save_param(customers.x_DG_new.values, [:CustomerType, :DERTech], :Capacity_MW, 
               joinpath(exportfilepath, "$(fileprefix)_x_DG.csv"))
end

function welfare_calculation(
    customers::Customers,
    model_data::HEMData, 
    regulator::Agent)

    adopt_model = customers.pv_adoption_model

    """
    Net Consumer Surplus Calculation:

    + Annualized Net Consumer Surplus of New DER installation (including Energy Savings and DER Excess Credits associated with New DERs)

    + Annualized Net Consumer Surplus of Existing PV installation (including Energy Savings and DER Excess Credits associated with Existing DERs)
    (Note: this term is assumed to be a constant carried over from previous years and not quantified)

    + Gross Surplus and energy consumption
    (Note: this term is assumed to be a constant and not quantified (demand is inelastic))

    - Cost of Energy Purchase
    (Note: this is the out-of-pocket payment for purchasing energy from the utility company, therefore, this term double-counted the energy savings already accounted for in the Annualized Net Consumer Surplus)

    - Energy Savings associaed with both new and existing DERs
    (Note: this term is to remove the double-counted energy savings from the terms above)

    Also note that DER Excess Credits are not listed here because they're implictly accounted for in the Annualized Net Consumer Surplus.

    """

    # The NetProfit represents the energy saving/credit per representative agent per DER technology, assuming the optimal DER technology size
    NetProfit = Dict((h,m) =>
        # value of distributed generation (offset load)
        sum(model_data.omega[t] * regulator.p[h,t] * min(customers.d[h,t], sum(customers.rho_DG[h,m,t] * customers.Opti_DG[h,m] for m in customers.index_m)) * (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] 
            for t in model_data.index_t) +
        # value of distributed generation (excess generation)
        sum(model_data.omega[t] * regulator.p_ex[h,t] * max(0, sum(customers.rho_DG[h,m,t] * customers.Opti_DG[h,m] for m in customers.index_m) - customers.d[h,t]) * (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t]
            for t in model_data.index_t) -
        # cost of distributed generation 
        customers.FOM_DG[h,m] * customers.Opti_DG[h,m] for h in model_data.index_h, m in customers.index_m)

    for h in model_data.index_h, m in customers.index_m
        if NetProfit[h,m] >= 0.0
            # Calculate total Net Consumer Surplus of PV installation
            Integral = Dict((h,m) =>
                quadgk(x -> customers.gamma[h]*customers.Opti_DG[h,m] * (1-Distributions.cdf(Distributions.Gamma(adopt_model.Shape[h,m], 1/adopt_model.Rate[h,m]*NetProfit[h,m]/customers.Opti_DG[h,m]),x)), 
                customers.CapEx_DG[h,m], 100*customers.CapEx_DG[h,m], rtol=1e-8))
            # Calculate annualized Net Consumer Surplus of PV installation
            if customers.MaxDG[h,m] == 0.0
                customers.ConPVNetSurplus[h,m] = 0.0
            else
                customers.ConPVNetSurplus[h,m] =
                    customers.delta * customers.x_DG_new[h,m]/customers.MaxDG[h,m]*Integral[h,m][1]
            end
        else
            customers.ConPVNetSurplus[h,m] = 0.0
        end
    end

    # Calculate energy savings associated with both new and existing DER
    EnergySaving = Dict((h,m) =>
        sum(model_data.omega[t] * regulator.p[h,t] * min(customers.d[h,t], sum(customers.rho_DG[h,m,t] * customers.Opti_DG[h,m] for m in customers.index_m)) * (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] 
            for t in model_data.index_t) *  (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]
            for h in model_data.index_h, m in customers.index_m)
    # Calculate out-of-pocket energy costs           
    EnergyCost = Dict((h,m) =>
        sum(model_data.omega[t]*regulator.p[h,t]*(customers.gamma[h]*customers.d[h,t] - 
        min(sum(customers.rho_DG[h,m,t]*customers.Opti_DG[h,m] for m in customers.index_m), customers.d[h,t]) * (customers.rho_DG[h,m,t]*customers.Opti_DG[h,m])/customers.DERGen[h,t] * 
        (customers.x_DG_E[h,m]+customers.x_DG_new[h,m])/customers.Opti_DG[h,m]) for t in model_data.index_t)
        for h in model_data.index_h, m in customers.index_m)
    # Finally, calculate Net Consumer Surplus
    ConNetSurplus = Dict((h,m) =>
        customers.ConPVNetSurplus[h,m] - EnergySaving[h,m] - EnergyCost[h,m] for h in model_data.index_h, m in customers.index_m)
    # Sum of Net Consumer Surplus across customer tpye and DER technology
    TotalConNetSurplus =  sum(ConNetSurplus[h,m] for h in model_data.index_h, m in customers.index_m)

    return customers.ConPVNetSurplus, EnergySaving, EnergyCost, ConNetSurplus, TotalConNetSurplus

end
