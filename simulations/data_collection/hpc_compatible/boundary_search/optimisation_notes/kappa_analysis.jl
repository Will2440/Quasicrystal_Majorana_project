module KappaAnalysis

export positive_windows, agreement_counts, aligned_changes, validate_comparable_configs

function validate_comparable_configs(reference, observed)
    keys_to_compare = (:N, :sequence_bson_path, :target_slope, :slope_tolerance,
        :target_phason, :phason_tolerance, :sequence_selection_mode, :x0_count, :x0_alpha,
        :topological_threshold, :topological_fraction)
    for key in keys_to_compare
        isequal(get(reference, key, nothing), get(observed, key, nothing)) ||
            throw(ArgumentError("configurations differ in $key"))
    end
    return nothing
end

function positive_windows(labels)
    windows = NamedTuple[]
    start = nothing
    for index in eachindex(labels)
        if labels[index] == 1 && start === nothing
            start = index
        end
        if start !== nothing && (labels[index] != 1 || index == lastindex(labels))
            last_positive = labels[index] == 1 ? index : index - 1
            push!(windows, (first=start, last=last_positive,
                lower_censored=start == firstindex(labels), upper_censored=last_positive == lastindex(labels)))
            start = nothing
        end
    end
    return windows
end

function agreement_counts(reference, observed)
    length(reference) == length(observed) || throw(ArgumentError("labels must align"))
    return (reference_topological=count(==(1), reference), reference_trivial=count(==(0), reference),
        topological_lost=count(index -> reference[index] == 1 && observed[index] == 0, eachindex(reference)),
        trivial_changed=count(index -> reference[index] == 0 && observed[index] == 1, eachindex(reference)))
end

function aligned_changes(reference_u, reference, observed_u, observed)
    length(reference_u) == length(reference) == length(observed_u) == length(observed) || throw(ArgumentError("grids and labels must align"))
    all(isapprox.(reference_u, observed_u; rtol=0, atol=1e-14)) || throw(ArgumentError("coordinates do not align"))
    return findall(index -> reference[index] != observed[index], eachindex(reference))
end

end