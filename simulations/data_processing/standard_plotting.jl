include("bson_unpacker.jl")

using Plots
using LaTeXStrings
using Measures
using Base.Threads
using ProgressMeter


###########################################################
########### Sec 1: Generic Plotting Functions #############
###########################################################
# # Following functions make plots of the raw data, they are designed to be used for any variable name which are the key names of the dataframe

function plt_eigenvalues_generalised(
    df::DataFrame, x_variable::Symbol;
    filename::Union{Nothing,String}=nothing,
    xlims=nothing,
    ylims=nothing,
    label_fontsize::Int=24,
    tick_fontsize::Int=18,
    plot_mp::Bool=false,
    mp_yrange::Tuple{Float64,Float64}=(-1.0, 0.0),
    plot_rho_delta_tracker::Bool=false,
    tracker_ratio::Union{Nothing,Float64}=nothing,
    tracker_position::Union{Nothing,Float64}=nothing,
    tracker_limits::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    shape_ratio::Float64=0.95,
    fixed_values...
)
    """
    Plots eigenvalues as a function of a specified variable from the DataFrame.

    Parameters:
    df           -- DataFrame containing the data.
    x_variable   -- Symbol representing the variable to plot on the x-axis.
    filename     -- Optional output path. If provided, the plot is saved there.
    ylims        -- Optional tuple for y-axis limits (ymin, ymax). If not provided, auto.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
    """

    width = 800
    height = Int(round(width * shape_ratio))
    main_plot = Plots.scatter(
        legend=false,
        xlabel=L"\mu/t_A",
        ylabel=L"E",
        grid=false,
        size=(width, height),
        framestyle=:box,
        xlabelfontsize=label_fontsize,
        ylabelfontsize=label_fontsize,
        xtickfontsize=tick_fontsize,
        ytickfontsize=tick_fontsize
    )
    if xlims !== nothing
        Plots.xlims!(main_plot, xlims)
    end
    if ylims !== nothing
        Plots.ylims!(main_plot, ylims)
    end

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    x_values = filtered_df[:, x_variable]
    all_eigenvalues = filtered_df[:, :eigenvalues]

    sort_indices = sortperm(x_values)
    x_values = x_values[sort_indices]
    all_eigenvalues = all_eigenvalues[sort_indices]

    num_eigenvalues = length(all_eigenvalues[1])

    for i in 1:num_eigenvalues
        eigenvalue_line = [eigenvalues[i] for eigenvalues in all_eigenvalues]
        plot!(main_plot, x_values, eigenvalue_line, color=:steelblue, alpha=0.8)
    end

    # Plot MP on a secondary y-axis if the raw MP data is available.
    if plot_mp && (:mp in propertynames(filtered_df))
        mp_values = [real(mp[floor(Int, length(mp) / 2) + 1]) for mp in filtered_df[:, :mp]]
        mp_values = mp_values[sort_indices]
        mp_axis = Plots.twinx(main_plot)
        plot!(
            mp_axis,
            x_values,
            mp_values,
            color=:orange,
            alpha=1.0,
            linewidth=2,
            legend=false,
            ylabel=L"\mathcal{M}",
            ylims=mp_yrange,
            ylabelfontsize=label_fontsize,
            ytickfontsize=tick_fontsize,
            grid=false,
            framestyle=:box,
        )
        if xlims !== nothing
            Plots.xlims!(mp_axis, xlims)
        end
    end

    full_plot = main_plot
    if plot_rho_delta_tracker
        if isnothing(tracker_ratio) || isnothing(tracker_position) || isnothing(tracker_limits)
            error("tracker_ratio, tracker_position, and tracker_limits must be provided when plot_rho_delta_tracker=true")
        end

        ratio_min, ratio_max = tracker_limits
        min_label = string(round(ratio_min; digits=3))
        max_label = string(round(ratio_max; digits=3))

        tracker_plot = Plots.scatter(
            [tracker_position],
            [0.0];
            legend=false,
            xlabel=L"\rho/\Delta",
            ylabel="",
            color=:orange,
            markerstrokecolor=:black,
            markersize=10,
            grid=false,
            framestyle=:none,
            xlims=(0.0, 1.0),
            ylims=(-0.5, 0.5),
            xticks=false,
            yticks=false,
            xlabelfontsize=label_fontsize,
            xtickfontsize=tick_fontsize,
            left_margin=5mm,
            right_margin=18mm,
            bottom_margin=5mm,
            top_margin=0mm,
        )
        hline!(tracker_plot, [0.0], color=:black, linewidth=1.5)
        annotate!(tracker_plot, 0.0, -0.22, text(min_label, tick_fontsize, :black, :left))
        annotate!(tracker_plot, 1.0, -0.22, text(max_label, tick_fontsize, :black, :right))

        full_plot = plot(
            main_plot,
            tracker_plot;
            layout=@layout([a{0.84h}; b{0.16h}]),
            size=(width, Int(round(width * (shape_ratio + 0.15)))),
        )
    end

    if isnothing(filename)
        Plots.display(full_plot)
    else
        Plots.savefig(full_plot, filename)
    end

    return full_plot
end

function plt_mp_generalised(df::DataFrame, x_variable::Symbol, filename::String; fixed_values...)
    """
    Plots Majorana Polarization (MP) as a function of a specified variable from the DataFrame.

    Parameters:
    df           -- DataFrame containing the data.
    x_variable   -- Symbol representing the variable to plot on the x-axis.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
    """

    plt = Plots.plot(
        legend=false, 
        xlabel=string(x_variable), 
        ylabel="MP", 
        grid=true,
        size=(800,600)
    )

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    x_values = filtered_df[:, x_variable]
    mp_values = [real(mp[floor(Int, length(mp) / 2) + 1]) for mp in filtered_df[:, :mp]]

    sort_indices = sortperm(x_values)
    x_values = x_values[sort_indices]
    mp_values = mp_values[sort_indices]

    Plots.plot!(x_values, mp_values, line=:solid) #marker=:circle,

    Plots.savefig(filename)
end

function plt_mp_heatmap(df::DataFrame, colour_var::Symbol, x_variable::Symbol, y_variable::Symbol, filename::String, title::String; fixed_values...)
    """
    Plots Majorana Polarization (MP) as a heatmap with specified x and y variables.

    Parameters:
    df           -- DataFrame containing the data.
    x_variable   -- Symbol representing the variable to plot on the x-axis.
    y_variable   -- Symbol representing the variable to plot on the y-axis.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
    """

    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    x_values = unique(filtered_df[!, x_variable]) |> sort
    y_values = unique(filtered_df[!, y_variable]) |> sort

    value_matrix = Matrix{Float64}(undef, length(y_values), length(x_values))

    for row in eachrow(filtered_df)
        x_index = findfirst(x -> x == row[x_variable], x_values)
        y_index = findfirst(y -> y == row[y_variable], y_values)
    
        if !isnothing(x_index) && !isnothing(y_index)
            value = row[colour_var]  # Already Vector{ComplexF64}
    
            # # Ensure it's non-empty and valid
            # if isa(mp_data, Vector{ComplexF64}) && !isempty(mp_data)
            #     mid_index = div(length(mp_data), 2) + 1
            #     mp_matrix[y_index, x_index] = real(mp_data[mid_index])
            # else
            #     mp_matrix[y_index, x_index] = NaN
            # end
            value_matrix[y_index, x_index] = real(value)
        end
    end

    sorted_indices = sortperm(y_values)
    y_values = y_values[sorted_indices]
    value_matrix = value_matrix[sorted_indices, :]

    plt = Plots.heatmap(
        x_values, 
        y_values, 
        value_matrix, 
        xlabel=string(x_variable), 
        ylabel=string(y_variable),
        # xlabel=L"\mu/t_A",
        # ylabel=L"\rho",
        title=title,
        color=:viridis,
        clims=(minimum(value_matrix), maximum(value_matrix)), # (0.0, 1.0),
        colorbar=true,
        grid=false,
        size=(800,600)
    )
    
    # display(plt)
    Plots.savefig(filename)
end

function plt_mp_heatmap_with_contour(
    df::DataFrame, 
    colour_var::Symbol, 
    x_variable::Symbol, 
    y_variable::Symbol, 
    filename::String,
    tol::Float64;
    fixed_values...
)
    """
    Plots Majorana Polarization (MP) as a heatmap with specified x and y variables,
    and adds a single contour line where the value crosses from < -1 + tol to > -1 + tol.

    Parameters:
    df           -- DataFrame containing the data.
    colour_var   -- Symbol for the variable used for heatmap coloring.
    x_variable   -- Symbol for the variable to plot on the x-axis.
    y_variable   -- Symbol for the variable to plot on the y-axis.
    filename     -- String: Path to save the plot.
    tol          -- Float64: Tolerance used for detecting crossings (default 1e-2).
    fixed_values -- Keyword arguments specifying fixed values for filtering the DataFrame.
    """
    # # Filter DataFrame for fixed values
    # filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)
    
    filtered_df = filter(
    row -> all(getproperty(row, key) == value for (key, value) in fixed_values) && getproperty(row, x_variable) != 0.0, df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    # Extract and sort unique x and y values
    x_values = unique(filtered_df[!, x_variable]) |> sort
    y_values = unique(filtered_df[!, y_variable]) |> sort

    # Initialize value matrix
    value_matrix = Matrix{Float64}(undef, length(y_values), length(x_values))

    # Populate the value matrix
    for row in eachrow(filtered_df)
        x_index = findfirst(x -> x == row[x_variable], x_values)
        y_index = findfirst(y -> y == row[y_variable], y_values)
        if !isnothing(x_index) && !isnothing(y_index)
            value_matrix[y_index, x_index] = real(row[colour_var])
        end
    end

    # Ensure data is sorted properly
    sorted_indices = sortperm(y_values)
    y_values = y_values[sorted_indices]
    value_matrix = value_matrix[sorted_indices, :]

    # Define the level for the contour line
    contour_level = -1 + tol

    # Create the heatmap
    heatmap_plot = Plots.heatmap(
        x_values, 
        y_values, 
        value_matrix, 
        xlabel=L"\mu/t_A", #string(x_variable), 
        ylabel=L"\rho", #string(y_variable), 
        color=:viridis,
        labelfontsize = 35,
        tickfontsize = 24,
        left_margin=10mm,
        bottom_margin=5mm,
        right_margin=10mm,
        top_margin=5mm,
        clims=(minimum(value_matrix), maximum(value_matrix)),
        colorbar=true,
        grid=false,
        size=(1600, 1200)
    )

    # Add the contour line for the level
    Plots.contour!(
        x_values, 
        y_values, 
        value_matrix, 
        levels=[contour_level], 
        color=:red, 
        lw=2.0, 
        legend=false
    )

    # Save the plot to the specified file
    Plots.savefig(heatmap_plot, filename)
end

function plt_mp_heatmap_mbs_disc_final(
    df::DataFrame, 
    colour_var::Symbol, 
    x_variable::Symbol, 
    y_variable::Symbol, 
    filename::String; 
    fixed_values...
)
    """
    Plots Majorana Polarization (MP) as a heatmap with specified x and y variables.
    Values equal to 0.0 are shown as grey, while non-zero values are plotted in color.

    Parameters:
    df           -- DataFrame containing the data.
    colour_var   -- Symbol representing the variable used for coloring the heatmap.
    x_variable   -- Symbol for the variable to plot on the x-axis.
    y_variable   -- Symbol for the variable to plot on the y-axis.
    filename     -- String: Path to save the heatmap plot.
    fixed_values -- Keyword arguments specifying fixed values for all other DataFrame columns.
    """

    # Filter the DataFrame for the fixed values
    filtered_df = filter(row -> all(getproperty(row, key) == value for (key, value) in fixed_values), df)

    if isempty(filtered_df)
        error("No results match the specified fixed values.")
    end

    # Extract unique x and y values, sorted
    x_values = unique(filtered_df[!, x_variable]) |> sort
    y_values = unique(filtered_df[!, y_variable]) |> sort

    # Initialize the value matrix
    value_matrix = Matrix{Float64}(undef, length(y_values), length(x_values))

    for row in eachrow(filtered_df)
        x_index = findfirst(x -> x == row[x_variable], x_values)
        y_index = findfirst(y -> y == row[y_variable], y_values)
        if !isnothing(x_index) && !isnothing(y_index)
            value = real(row[colour_var])
            # Set 0.0 values to NaN to exclude them from the color map
            value_matrix[y_index, x_index] = value == 0.0 ? NaN : value
        end
    end

    # Ensure y_values and value_matrix are sorted consistently
    sorted_indices = sortperm(y_values)
    y_values = y_values[sorted_indices]
    value_matrix = value_matrix[sorted_indices, :]

    # Define the colormap and the background color for NaN (grey)
    heatmap_plot = Plots.heatmap(
        x_values,
        y_values,
        value_matrix,
        xlabel=L"\mu/t_A",
        ylabel=L"\rho",
        color=:plasma,
        nan_color=:grey,           # Set NaN (0.0 values) to be displayed as grey
        bg_inside=:grey,
        labelfontsize=35,
        tickfontsize=24,
        left_margin=10mm,
        bottom_margin=5mm,
        right_margin=10mm,
        top_margin=5mm,
        clims=(0.0, 2.0), #(minimum(skipmissing(value_matrix)), maximum(skipmissing(value_matrix))),
        colorbar=true,
        grid=false,
        size=(1600, 1200)
    )

    # Save the plot
    Plots.savefig(heatmap_plot, filename)
end


###########################################################
########### Sec 2: Standard Plotting Processes ############
###########################################################


function standard_plotting(
    df::DataFrame,
    seq_name::String,
    mp_tol::Float64,
    filepath::String;
    N::Int64=50,
    rho::Float64=1.5,
    Delta::Float64=0.1,
)
    rho_safe = replace(string(rho), "." => "p")
    Delta_safe = replace(string(Delta), "." => "p")

    # filename = "$(filepath)_$(seq_name)_N$(N)_rho$(rho_safe)_Delta$(Delta_safe)_mp_line"
    # plt_mp_generalised(df, :mu_t, filename, N=N, rho=rho, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mp"
    plt_mp_heatmap(df, :mp, :mu_t, :rho, filename, "MP for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mp_disc"
    plt_mp_heatmap(df, :mp_disc, :mu_t, :rho, filename, "MP discretised (tol=$mp_tol) for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_ipr"
    # plt_mp_heatmap(df, :ipr, :mu_t, :rho, filename, "IPR for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_ipr_disc"
    # plt_mp_heatmap(df, :ipr_masked, :mu_t, :rho, filename, "IPR phase masked (MP tol=$mp_tol) for $(seq_name) N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mbs_gap"
    # plt_mp_heatmap(df, :maj_gap, :mu_t, :rho, filename, "MBS gap for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mbs_gap_disc"
    # plt_mp_heatmap_mbs_disc_final(df, :maj_gap_masked, :mu_t, :rho, filename, "MBS gap masked (MP tol=$mp_tol) for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    return
end

function standard_plotting_PQC(
    df::DataFrame,
    seq_name::String,
    mp_tol::Float64,
    filepath::String;
    N::Int64=50,
    rho::Float64=1.5,
    sigma::Float64=2.0,
    Delta::Float64=0.1,
)
    rho_safe = replace(string(rho), "." => "p")
    sig_safe = replace(string(sigma), "." => "p")
    Delta_safe = replace(string(Delta), "." => "p")

    # filename = "$(filepath)_$(seq_name)_N$(N)_rho$(rho_safe)_Delta$(Delta_safe)_mp_line"
    # plt_mp_generalised(df, :mu_t, filename, N=N, rho=rho, Delta_t=Delta)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_mp"
    plt_mp_heatmap(df, :mp, :mu_t, :rho, filename, "MP for $(seq_name), sigma=$sigma, N=$N, Delta=$Delta"; N=N, Delta_t=Delta, sigma=sigma)#, seq_name=seq_name)

    filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_mp_disc"
    plt_mp_heatmap(df, :mp_disc, :mu_t, :rho, filename, "MP discretised (tol=$mp_tol) for $(seq_name), sigma=$sigma, N=$N, Delta=$Delta"; N=N, Delta_t=Delta, sigma=sigma)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_ipr"
    # plt_mp_heatmap(df, :ipr, :mu_t, :rho, filename, "IPR for $(seq_name), sigma=$sigma, N=$N, Delta=$Delta"; N=N, Delta_t=Delta, sigma=sigma)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_ipr_disc"
    # plt_mp_heatmap(df, :ipr_masked, :mu_t, :rho, filename, "IPR phase masked (MP tol=$mp_tol) for $(seq_name), sigma=$sigma, N=$N, Delta=$Delta"; N=N, Delta_t=Delta, sigma=sigma)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_sig$(sig_safe)_N$(N)_Delta$(Delta_safe)_mbs_gap"
    # plt_mp_heatmap(df, :maj_gap, :mu_t, :rho, filename, "MBS gap for $(seq_name), N=$N, sigma=$sigma, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    # filename = "$(filepath)_$(seq_name)_N$(N)_Delta$(Delta_safe)_mbs_gap_disc"
    # plt_mp_heatmap_mbs_disc_final(df, :maj_gap_masked, :mu_t, :rho, filename, "MBS gap masked (MP tol=$mp_tol) for $(seq_name), N=$N, Delta=$Delta"; N=N, Delta_t=Delta)#, seq_name=seq_name)

    return
end


###########################################################
###################### Sec 3: Run #########################
###########################################################

if abspath(PROGRAM_FILE) == @__FILE__
    filepath = "change/" # have the / at the end!
    mp_tol = 0.1
    df = unpack_bason_standard(filepath; mp_tol=mp_tol)

    N = 50
    seq_name = "GQC"


    # # One instance
    # # Delta = 1.2
    # plots_savepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance/PQC_N(50-50-1)_t1(1.0-1.0-1)_t2(0.0-15.0-151)_t3(3.0-3.0-1)_mu(0.0-15.0-151)_Delta(0.1-2.0-20)/"
    # # standard_plotting(df, seq_name, mp_tol, filepath; N=N, Delta=Delta)
    # for Delta in collect(range(0.1, 2.0, 20))
    #     standard_plotting_PQC(df, seq_name, mp_tol, filepath; Delta=Delta, sigma=3.0)
    # end


    # # Instances of all Delta
    # plots_savepath = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/np/abundance_data/GQC_N(50-50-1)_t1(1.0-1.0-1_t2(0.0-1000.0-1001)_mu(0.0-3.0-31)_Delta(0.5-0.5-1)/"
    # isdir(plots_savepath) || mkpath(plots_savepath)

    # Delta_range = collect(range(0.1, 2.0, 20))
    # for Delta in Delta_range
    #     standard_plotting(df, seq_name, mp_tol, plots_savepath; N=N, Delta=Delta)
    # end


    # Instances of Delta and sigma for PQC
    Delta_range = collect(range(0.0, 2.0, 21))
    sigma_range = [3.0] #[2.0, 3.0, 4.0]
    seq_name = "PQC"
    for sigma in sigma_range
        for Delta in Delta_range
            standard_plotting_PQC(df, seq_name, mp_tol, filepath; Delta=Delta, sigma=sigma)
        end
    end
end