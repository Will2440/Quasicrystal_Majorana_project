using DataFrames
using ProgressMeter

# function prep_df_for_IDOS(df::DataFrame; fixed_values...)
#     """
#     Collect mu-range vectors and corresponding discretised indicators per unique (phi, N, t_n, Delta).
#     Assumes :disc_mp exists (scalar 0/1 per row) and :disc_eigenvalues exists (Vector{Int} per row).
#     Returns DataFrame with columns:
#       :phi, :N, :t_n, :Delta, :mu_values, :disc_mp, :disc_evals_first
#     """
#     filtered_df = filter(row -> all(getproperty(row, k) == v for (k, v) in fixed_values), df)
#     if isempty(filtered_df)
#         error("No results match the specified fixed values.")
#     end

#     rows = Vector{Dict{Symbol,Any}}()

#     for g in DataFrames.groupby(filtered_df, [:phi, :N, :t_n, :Delta])
#         g_sorted = sort(g, :mu)

#         mu_values = [r.mu for r in eachrow(g_sorted)]
#         # assume disc_mp is scalar 0/1 (may be Float64 but integer-like)
#         disc_mp_vals = [Int(r.disc_mp) for r in eachrow(g_sorted)]
#         # assume disc_eigenvalues is a vector and we want its first element
#         disc_eval_first_vals = [Int(r.disc_eigenvalues[1]) for r in eachrow(g_sorted)]

#         eigenvalues_per_mu = [r.eigenvalues for r in eachrow(g_sorted)]

#         first_row = first(g_sorted)

#         push!(rows, Dict(
#             :phi => first_row.phi,
#             :N => first_row.N,
#             :t_n => first_row.t_n,
#             :Delta => first_row.Delta,
#             :mu_values => mu_values,
#             :disc_mp => disc_mp_vals,
#             :disc_evals => disc_eval_first_vals,
#             :eigenvalues => eigenvalues_per_mu
#         ))
#     end

#     return DataFrame(rows)
# end

function compute_idos_df!(df::DataFrame)
    idos_col = Vector{Vector{Float64}}(undef, nrow(df))
    for (i, row) in enumerate(eachrow(df))
        eigenvalues = sort(real(row.eigenvalues))
        N = length(eigenvalues)
        idos = collect(1:N) ./ N
        idos_col[i] = idos
    end
    df.idos = idos_col
    return df
end

function find_idos_plateaus(df::DataFrame; threshold=0.05, fixed_values...)
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    plateau_idos = Vector{Vector{Float64}}(undef, nrow(filtered_df))
    for (i, row) in enumerate(eachrow(filtered_df))
        energies = sort(real(row.eigenvalues))
        idos = row.idos
        gaps = Int[]
        for j in 1:length(energies)-1
            if energies[j+1] - energies[j] > threshold
                push!(gaps, j)
            end
        end
        # IDOS value at the gap is between j and j+1
        plateau_idos[i] = [(j + 0.5) / length(idos) for j in gaps]
    end
    return plateau_idos
end

function compute_plateaus_df!(df::DataFrame; threshold=0.01)
    plateaus_col = Vector{Vector{Float64}}(undef, nrow(df))
    @showprogress for (i, row) in enumerate(eachrow(df))
        energies = sort(real(row.eigenvalues))
        idos = row.idos
        gaps = Int[]
        for j in 1:length(energies)-1
            if energies[j+1] - energies[j] > threshold
                push!(gaps, j)
            end
        end
        plateaus_col[i] = [(j + 0.5) / length(idos) for j in gaps]
    end
    df.plateaus = plateaus_col
    return df
end