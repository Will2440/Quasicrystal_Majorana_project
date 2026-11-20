module ResolutionAnalysis

export grid_indices, phase_runs, grid_metrics, coarsest_passing_stride

function coarsest_passing_stride(rows)
    reference = only(filter(row -> row.stride == 1 && row.offset == 0, rows))
    for stride in sort!(unique([row.stride for row in rows if row.stride > 1]); rev=true)
        pair = filter(row -> row.stride == stride, rows)
        length(pair) == 2 || continue
        Set(row.offset for row in pair) == Set((0, stride ÷ 2)) || continue
        all(row -> row.missed_fine_runs == 0 && row.phase0_runs == reference.phase0_runs &&
            row.phase1_runs == reference.phase1_runs && row.transitions == reference.transitions, pair) && return stride
    end
    return nothing
end

function grid_indices(length_fine::Int, stride::Int, offset::Int=0)
    length_fine >= 2 && stride >= 1 && (length_fine - 1) % stride == 0 || error("incompatible grid")
    0 <= offset < stride || error("offset outside stride")
    return sort!(unique!(vcat(1, collect((1 + offset):stride:length_fine), length_fine)))
end

function phase_runs(labels)
    isempty(labels) && return NamedTuple[]
    runs = NamedTuple[]
    start = firstindex(labels)
    for index in (start + 1):lastindex(labels)
        if labels[index] != labels[index - 1]
            push!(runs, (phase=Int(labels[start]), first=start, last=index - 1))
            start = index
        end
    end
    push!(runs, (phase=Int(labels[start]), first=start, last=lastindex(labels)))
    return runs
end

function grid_metrics(labels, stride::Int, offset::Int, mu_step::Real)
    indices = grid_indices(length(labels), stride, offset)
    observed = phase_runs(labels[indices])
    reference = phase_runs(labels)
    missed = count(run -> !any(index -> run.first <= index <= run.last, indices), reference)
    interior = filter(run -> run.phase == 1 && run.first > 1 && run.last < length(labels), reference)
    min_observed_width = isempty(interior) ? missing : minimum((run.last - run.first + 1) * mu_step for run in interior)
    return (n_points=length(indices), nominal_mu_spacing=stride * mu_step,
        phase0_runs=count(run -> run.phase == 0, observed),
        phase1_runs=count(run -> run.phase == 1, observed),
        transitions=length(observed) - 1, missed_fine_runs=missed,
        min_fine_topological_sample_width=min_observed_width)
end

end