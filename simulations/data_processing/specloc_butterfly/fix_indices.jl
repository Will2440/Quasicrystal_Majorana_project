# Read the file
content = read("processing.jl", String)

# The issue is that axes(group, 1) gives indices relative to the subset,
# but we need indices relative to the original DataFrame
# Instead, we need to track the row numbers from the GroupedDataFrame

# Let's look for the problematic section and replace it
old_section = """    # Use groupby to process groups while preserving row indices
    for group in groupby(df, by_cols)
        key = Tuple(group[1, col] for col in by_cols)
        group_indices = axes(group, 1)  # Get row indices for this group in original df"""

new_section = """    # Use groupby to process groups while preserving row indices
    for group in groupby(df, by_cols)
        key = Tuple(group[1, col] for col in by_cols)
        group_indices = parentindices(group)[1]  # Get row indices for this group in original df"""

if occursin(old_section, content)
    content = replace(content, old_section => new_section)
    write("processing.jl", content)
    println("Fixed indices issue using parentindices")
else
    println("Could not find section to fix")
end
