module DynamicAdaptive

using Statistics

export initial_cells,
       cell_sample_points,
       collect_candidate_points,
       classify_and_refine,
       estimate_trivial_gap_count_at_rho_max

function initial_cells(cfg::Dict{Symbol,Any})
    rho_lo, rho_hi = cfg[:rho_range]
    u_lo, u_hi = cfg[:u_range]
    n_rho = cfg[:init_rho_bins]
    n_u = cfg[:init_u_bins]

    rho_edges = collect(range(rho_lo, rho_hi, length=n_rho + 1))
    u_edges = collect(range(u_lo, u_hi, length=n_u + 1))

    cells = NamedTuple[]
    cell_id = 1
    for i in 1:n_rho
        for j in 1:n_u
            push!(cells, (
                id = cell_id,
                depth = 0,
                rho_lo = rho_edges[i],
                rho_hi = rho_edges[i+1],
                u_lo = u_edges[j],
                u_hi = u_edges[j+1],
            ))
            cell_id += 1
        end
    end
    return cells
end

function cell_sample_points(cell)
    rho_mid = 0.5 * (cell.rho_lo + cell.rho_hi)
    u_mid = 0.5 * (cell.u_lo + cell.u_hi)

    return [
        (rho=cell.rho_lo, u=cell.u_lo, role=:corner),
        (rho=cell.rho_lo, u=cell.u_hi, role=:corner),
        (rho=cell.rho_hi, u=cell.u_lo, role=:corner),
        (rho=cell.rho_hi, u=cell.u_hi, role=:corner),
        (rho=rho_mid, u=u_mid, role=:center),
    ]
end

function point_key(rho::Float64, u::Float64, digits::Int)
    return string(round(rho, digits=digits), "|", round(u, digits=digits))
end

function collect_candidate_points(
    cells,
    known_results::Dict{String,NamedTuple},
    cfg::Dict{Symbol,Any}
)
    digits = cfg[:point_round_digits]
    todo = Dict{String,NamedTuple}()

    for cell in cells
        for p in cell_sample_points(cell)
            key = point_key(p.rho, p.u, digits)
            if !haskey(known_results, key) && !haskey(todo, key)
                todo[key] = (
                    point_key = key,
                    rho = p.rho,
                    u = p.u,
                )
            end
        end
    end

    return collect(values(todo))
end

function split_cell(cell, next_id::Int)
    rho_mid = 0.5 * (cell.rho_lo + cell.rho_hi)
    u_mid = 0.5 * (cell.u_lo + cell.u_hi)

    children = [
        (id=next_id, depth=cell.depth + 1, rho_lo=cell.rho_lo, rho_hi=rho_mid, u_lo=cell.u_lo, u_hi=u_mid),
        (id=next_id + 1, depth=cell.depth + 1, rho_lo=cell.rho_lo, rho_hi=rho_mid, u_lo=u_mid, u_hi=cell.u_hi),
        (id=next_id + 2, depth=cell.depth + 1, rho_lo=rho_mid, rho_hi=cell.rho_hi, u_lo=cell.u_lo, u_hi=u_mid),
        (id=next_id + 3, depth=cell.depth + 1, rho_lo=rho_mid, rho_hi=cell.rho_hi, u_lo=u_mid, u_hi=cell.u_hi),
    ]

    return children
end

function classify_and_refine(
    cells,
    known_results::Dict{String,NamedTuple},
    cfg::Dict{Symbol,Any}
)
    digits = cfg[:point_round_digits]
    max_depth = cfg[:max_depth]
    rho_tol = cfg[:rho_cell_tol]
    u_tol = cfg[:u_cell_tol]
    gap_refine_tol = cfg[:gap_refine_tol]

    refined = NamedTuple[]
    resolved = NamedTuple[]
    unresolved = NamedTuple[]
    next_id = maximum(c.id for c in cells) + 1

    for cell in cells
        samples = cell_sample_points(cell)

        present = true
        labels = Int[]
        gaps = Float64[]

        for p in samples
            key = point_key(p.rho, p.u, digits)
            if !haskey(known_results, key)
                present = false
                break
            end
            r = known_results[key]
            push!(labels, r.phase_label)
            push!(gaps, r.gap_min)
        end

        if !present
            push!(unresolved, (cell..., reason=:missing_points))
            continue
        end

        unique_corner_labels = unique(labels[1:4])
        mixed_phase = length(unique_corner_labels) > 1
        near_critical = minimum(gaps) <= gap_refine_tol

        rho_width = cell.rho_hi - cell.rho_lo
        u_width = cell.u_hi - cell.u_lo
        splittable = cell.depth < max_depth && rho_width > rho_tol && u_width > u_tol

        should_refine = (mixed_phase || near_critical) && splittable

        if should_refine
            children = split_cell(cell, next_id)
            append!(refined, children)
            next_id += 4
        else
            label = round(Int, mean(labels) >= 0.5)
            conf = minimum(gaps)
            push!(resolved, (cell..., phase_label=label, confidence_gap=conf, mixed_phase=mixed_phase))
        end
    end

    unresolved_base = [
        (
            id=c.id,
            depth=c.depth,
            rho_lo=c.rho_lo,
            rho_hi=c.rho_hi,
            u_lo=c.u_lo,
            u_hi=c.u_hi,
        ) for c in unresolved
    ]

    return (
        next_cells = vcat(refined, unresolved_base),
        resolved_cells = resolved,
        unresolved_cells = unresolved,
        n_refined = length(refined),
    )
end

function find_cell_label(cells_with_labels, rho::Float64, u::Float64)
    best = nothing
    best_depth = -1

    for c in cells_with_labels
        inside = (c.rho_lo <= rho <= c.rho_hi) && (c.u_lo <= u <= c.u_hi)
        if inside && c.depth > best_depth
            best = c
            best_depth = c.depth
        end
    end

    return best === nothing ? missing : best.phase_label
end

function estimate_trivial_gap_count_at_rho_max(
    resolved_cells,
    cfg::Dict{Symbol,Any}
)
    if isempty(resolved_cells)
        return 0
    end

    rho_max = cfg[:rho_range][2]
    u_lo, u_hi = cfg[:u_range]
    n_probe = cfg[:rho_max_probe_points]
    min_run = cfg[:min_gap_run_points]

    us = collect(range(u_lo, u_hi, length=n_probe))
    labels = Vector{Union{Int,Missing}}(undef, length(us))

    for i in eachindex(us)
        labels[i] = find_cell_label(resolved_cells, rho_max, us[i])
    end

    for i in eachindex(labels)
        if labels[i] === missing
            left = findlast(x -> x !== missing, labels[1:i])
            right_local = findfirst(x -> x !== missing, labels[i:end])
            right = right_local === nothing ? nothing : i + right_local - 1

            if left !== nothing && right !== nothing
                labels[i] = (i - left) <= (right - i) ? labels[left] : labels[right]
            elseif left !== nothing
                labels[i] = labels[left]
            elseif right !== nothing
                labels[i] = labels[right]
            else
                labels[i] = 1
            end
        end
    end

    gap_count = 0
    run = 0
    for lbl in labels
        if lbl == 0
            run += 1
        else
            if run >= min_run
                gap_count += 1
            end
            run = 0
        end
    end
    if run >= min_run
        gap_count += 1
    end

    return gap_count
end

end # module
