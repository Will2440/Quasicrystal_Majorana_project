module SliceSweep

using ..DynamicBandwidth
using ..SliceSearch

export sweep_parameter, sweep_grid, sweep_progress, fixed_mu_critical, physical_parameters

sweep_parameter(cfg) = get(cfg, :sweep_parameter, :rho)

function sweep_grid(cfg)
    parameter = sweep_parameter(cfg)
    if parameter == :rho
        return SliceSearch.build_rho_grid(cfg[:rho_start], cfg[:rho_max],
                                          cfg[:n_rho_slices], cfg[:rho_spacing_mode])
    elseif parameter == :Delta
        return SliceSearch.build_delta_grid(cfg[:Delta_start], cfg[:Delta_max],
            cfg[:n_Delta_slices], cfg[:Delta_spacing_mode];
            concentration=get(cfg, :Delta_exponential_concentration, 4.0))
    end
    error("Unknown sweep_parameter: $parameter (expected :rho or :Delta)")
end

function sweep_progress(value::Real, cfg)
    parameter = sweep_parameter(cfg)
    start = Float64(cfg[parameter == :rho ? :rho_start : :Delta_start])
    stop = Float64(cfg[parameter == :rho ? :rho_max : :Delta_max])
    start == stop && return 0.0
    return clamp((Float64(value) - start) / (stop - start), 0.0, 1.0)
end

function fixed_mu_critical(cfg, gamma::Real)
    sweep_parameter(cfg) == :Delta || return nothing
    t1, t2 = Float64(cfg[:t1_fixed]), Float64(cfg[:t2_fixed])
    isfinite(t1) && isfinite(t2) && t1 > 0 && t2 > 0 ||
        error("t1_fixed and t2_fixed must be finite and positive for a Delta sweep")
    return DynamicBandwidth.mu_critical_estimate(Float64(gamma), t1, t2)
end

function physical_parameters(value::Real, cfg, gamma::Real, fixed_critical=nothing)
    if sweep_parameter(cfg) == :rho
        rho = Float64(value)
        t1, t2 = DynamicBandwidth.tpair_from_rho(rho;
            sweep_mode=cfg[:sweep_mode], t1_fixed=Float64(cfg[:t1_fixed]),
            t2_fixed=Float64(cfg[:t2_fixed]))
        Delta = Float64(cfg[:Delta])
        mu_critical = DynamicBandwidth.mu_critical_estimate(Float64(gamma), t1, t2)
    elseif sweep_parameter(cfg) == :Delta
        t1, t2 = Float64(cfg[:t1_fixed]), Float64(cfg[:t2_fixed])
        rho, Delta = t2 / t1, Float64(value)
        mu_critical = isnothing(fixed_critical) ? fixed_mu_critical(cfg, gamma) : fixed_critical
    else
        error("Unknown sweep_parameter: $(sweep_parameter(cfg))")
    end
    return (rho=rho, Delta=Delta, t1=t1, t2=t2, mu_critical=Float64(mu_critical))
end

end