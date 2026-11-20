module SliceSweep

using ..DynamicBandwidth
using ..SliceSearch

export sweep_parameter, sweep_grid, sweep_progress, fixed_mu_critical, physical_parameters, kappas_for_delta

sweep_parameter(cfg) = get(cfg, :sweep_parameter, :rho)

function kappas_for_delta(cfg, delta::Real)
    mode = get(cfg, :kappa_mode, :user_defined)
    if mode == :user_defined
        values = Float64.(collect(cfg[:kappas]))
        !isempty(values) && all(value -> isfinite(value) && value >= 0, values) ||
            error("kappas must be a nonempty collection of finite nonnegative values")
        return values
    elseif mode == :delta_over_L
        length_sites = cfg[:N]
        length_sites isa Integer && length_sites > 0 || error("N must be a positive integer for kappa_mode=:delta_over_L")
        isfinite(delta) && delta >= 0 || error("Delta must be finite and nonnegative for kappa_mode=:delta_over_L")
        return [Float64(delta) / length_sites]
    end
    error("Unknown kappa_mode: $mode (expected :user_defined or :delta_over_L)")
end

function sweep_grid(cfg)
    parameter = sweep_parameter(cfg)
    if parameter == :rho
        return SliceSearch.build_rho_grid(cfg[:rho_start], cfg[:rho_max],
            cfg[:n_rho_slices], cfg[:rho_spacing_mode];
            concentration=get(cfg, :rho_log_concentration, 4.0))
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