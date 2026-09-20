module DynamicBandwidth

export tpair_from_rho, bandwidth_estimate, mu_critical_estimate

"""
Map rho = t2/t1 to (t1,t2) according to sweep mode.
- :vary_t2 => keep t1 fixed, vary t2 = rho * t1
- :vary_t1 => keep t2 fixed, vary t1 = t2 / rho
"""
function tpair_from_rho(
    rho::Float64;
    sweep_mode::Symbol=:vary_t2,
    t1_fixed::Float64=1.0,
    t2_fixed::Float64=1.0
)
    if sweep_mode == :vary_t2
        t1 = t1_fixed
        t2 = rho * t1
    elseif sweep_mode == :vary_t1
        t2 = t2_fixed
        t1 = t2 / rho
    else
        error("Unknown sweep_mode: $sweep_mode. Use :vary_t2 or :vary_t1")
    end

    return t1, t2
end

"""
Bandwidth estimate used in the paper's model section:
W_gamma ≈ 4 * tbar_gamma,
tbar_gamma = (1 - gamma) * t1 + gamma * t2.
"""
function bandwidth_estimate(gamma::Float64, t1::Float64, t2::Float64)::Float64
    tbar = (1.0 - gamma) * t1 + gamma * t2
    return 4.0 * tbar
end

"""
Critical chemical potential estimate from bandwidth:
|mu_c| ≈ W_gamma / 2.
"""
function mu_critical_estimate(gamma::Float64, t1::Float64, t2::Float64)::Float64
    return 0.5 * bandwidth_estimate(gamma, t1, t2)
end

end # module
