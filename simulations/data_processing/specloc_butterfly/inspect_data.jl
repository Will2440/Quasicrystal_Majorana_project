"""
    file name:   inspect_data.jl
    created:     2026-09-22

    overview:
        Quick inspection script to examine the structure of JLD2 files
        and verify that the processing.jl module will work correctly.
    
    usage:
        julia inspect_data.jl
        
    This will:
        - Load the first chunk file
        - Display column names and types
        - Show sample rows
        - Verify presence of key columns (specloc_gap, mu, phi, etc.)
"""

using JLD2
using DataFrames
using Statistics

# Path to the raw data folder
data_folder = "hof_style_slopes_N1000_phason_0.0-1-0.0_nbins1000_npb1_N(500-500-1)_t1(1p0-1p0-1)_t2(2p0-2p0-1)_mu(0p0-4p0-1201)_Delta(0p05-0p05-1)_x0(218p75-281p25-5)_kappa(0p0001-0p0001-1)"
folder_path = normpath(joinpath(@__DIR__, "..", "..", "raw_data", "specloc_butterfly", data_folder))

println("=" ^ 80)
println("SPECLOC BUTTERFLY DATA INSPECTION")
println("=" ^ 80)
println("\nData folder: $folder_path\n")

# Check if folder exists
if !isdir(folder_path)
    println("ERROR: Folder does not exist!")
    println("Expected: $folder_path")
    exit(1)
end

# Get list of .jld2 files
files = filter(f -> endswith(f, ".jld2"), readdir(folder_path))
sort!(files)

println("Found $(length(files)) .jld2 files")
if isempty(files)
    println("ERROR: No .jld2 files found!")
    exit(1)
end

# Load the first file
first_file = files[1]
first_path = joinpath(folder_path, first_file)

println("\nLoading first file: $first_file\n")

data = JLD2.load(first_path)

println("Keys in file: $(keys(data))")

# Try to extract DataFrame
df = nothing
if haskey(data, "result")
    df = data["result"]
    println("Found 'result' key")
elseif length(data) == 1
    df = first(values(data))
    println("Extracted data from single key: $(first(keys(data)))")
else
    # Try to use the whole thing as a DataFrame
    if isa(data, DataFrame)
        df = data
        println("File data is already a DataFrame")
    else
        df = first(values(data))
        println("Extracted data from first value")
    end
end

if isnothing(df)
    println("ERROR: Could not extract DataFrame from file!")
    exit(1)
end

if !isa(df, DataFrame)
    println("ERROR: Extracted data is not a DataFrame!")
    println("Type: $(typeof(df))")
    exit(1)
end

println("\n" * "=" ^ 80)
println("DATAFRAME STRUCTURE")
println("=" ^ 80)

println("\nShape: $(size(df))")
println("Rows: $(nrow(df))")
println("Columns: $(ncol(df))\n")

println("Column names and types:")
for (i, col) in enumerate(names(df))
    T = eltype(df[!, col])
    print(lpad(i, 3), ") $col")
    print(" " ^ max(1, 25 - length(col)))
    println(": $T")
end

println("\n" * "=" ^ 80)
println("DATA SUMMARY")
println("=" ^ 80)

# Show first few rows
println("\nFirst 5 rows:\n")
show(stdout, "text/plain", first(df, 5))

# Check for key columns
println("\n" * "=" ^ 80)
println("KEY COLUMN CHECK")
println("=" ^ 80 * "\n")

key_columns = [
    :N, :t_n, :mu, :phi, :Delta,
    :specloc_gap, :localiser_gap, :specloc_localiser_gap,
    :specloc_winding, :winding, :projector_winding
]

found_cols = intersect(Symbol.(names(df)), key_columns)
missing_cols = setdiff(key_columns, found_cols)

println("✓ Found columns:")
for col in found_cols
    println("  - $col ($(eltype(df[!, col])))")
end

if !isempty(missing_cols)
    println("\n✗ Missing standard column names (but may exist with different names):")
    for col in missing_cols
        println("  - $col")
    end
end

# Show statistics for numeric columns
println("\n" * "=" ^ 80)
println("NUMERIC COLUMN STATISTICS")
println("=" ^ 80 * "\n")

for col in names(df)
    if eltype(df[!, col]) <: Union{Float64, Int64}
        values = df[!, col]
        if !isa(values, Vector)
            continue  # Skip if not a simple vector
        end
        
        if all(x -> !ismissing(x), values)
            min_val = minimum(values)
            max_val = maximum(values)
            mean_val = mean(skipmissing(values))
            
            println("$col:")
            println("  Min:  $min_val")
            println("  Max:  $max_val")
            println("  Mean: $mean_val\n")
        end
    end
end

println("=" ^ 80)
println("INSPECTION COMPLETE")
println("=" ^ 80)
println("\nNotes:")
println("- Verify that key columns are present for plotting")
println("- If column names differ, update the column mapping in processing.jl")
println("- Check that specloc_gap and mu/phi values are reasonable")
