# Read the file
content = read("processing.jl", String)

# Fix the problematic line - remove [axes(group, 1)] indexing
old_line = "@. norm_values[group_indices] = 2 * (group.mu_specloc_invariant[axes(group, 1)] - vmin) / (vmax - vmin) -"
new_line = "@. norm_values[group_indices] = 2 * (group.mu_specloc_invariant - vmin) / (vmax - vmin) -"

content = replace(content, old_line => new_line)

# Write back
write("processing.jl", content)
println("Fixed processing.jl")
