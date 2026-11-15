include("standard_plotting.jl")
include("jld2_unpacker.jl")

const DATA_PATH = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/raw_data/rho-Delta_vary/hp/GQC_N(500-500-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-501)_Delta(0.01-1.0-101)/"
const OUTPUT_PATH = "/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/results/rho-Delta_vary/hp/GQC_N(500-500-1)_t1(1.0-1.0-1)_t2(1.5-1.5-1)_mu(0.0-3.0-501)_Delta(0.01-1.0-101)/"
const MP_TOL = 1e-5

function float_label(value::Real)
	return replace(string(value), "." => "p")
end

function has_column(df::DataFrame, column::Symbol)
	return column in Symbol.(names(df))
end

function unique_fixed_values(df::DataFrame)
	fixed_values = Dict{Symbol,Any}()

	for column in (:N, :mu_t, :sequence_name)
		if has_column(df, column)
			unique_values = unique(df[!, column])
			if length(unique_values) == 1
				fixed_values[column] = only(unique_values)
			end
		end
	end

	return (; pairs(fixed_values)...)
end

function unpack_supported_data(folder_path::String; mp_tol::Float64=MP_TOL)
	bson_files = glob("*.bson", folder_path)
	if !isempty(bson_files)
		return unpack_bason_standard(folder_path; mp_tol=mp_tol)
	end

	jld2_files = glob("*.jld2", folder_path)
	if !isempty(jld2_files)
		df = process_jld2_files_fluid(folder_path)
		df = calc_norms_df!(df)
		df = calc_mp_disc!(df, mp_tol)
		println("DataFrame keynames: $(names(df))")
		return df
	end

	return DataFrame()
end


# New: Plot eigenvalues vs mu for each (rho, Delta)
function save_eigenvalue_plots_by_mu(data_path::String, output_path::String; mp_tol::Float64=MP_TOL, xlims=nothing, ylims=nothing, label_fontsize::Int=24, tick_fontsize::Int=18, plot_mp::Bool=false, mp_yrange::Tuple{Float64,Float64}=(-1.0, 0.0), plot_rho_delta_tracker::Bool=false, tracker_rho::Float64=1.5, shape_ratio::Float64=0.95)
	isdir(output_path) || mkpath(output_path)

	df = unpack_supported_data(data_path; mp_tol=mp_tol)
	if isempty(df)
		error("No data was unpacked from $(data_path)")
	end

	required_columns = [:mu, :rho, :Delta_t, :eigenvalues]
	missing_columns = [column for column in required_columns if !has_column(df, column)]
	if !isempty(missing_columns)
		error("Missing required DataFrame columns: $(missing_columns)")
	end

	base_fixed_values = unique_fixed_values(df)
	delta_values = sort(unique(df[!, :Delta_t]))
	rho_values = sort(unique(df[!, :rho]))
	delta_min = minimum(delta_values)
	delta_max = maximum(delta_values)
	tracker_limits = plot_rho_delta_tracker ? (tracker_rho / maximum(delta_values), tracker_rho / minimum(delta_values)) : nothing

	println("Saving eigenvalue plots vs mu for $(length(delta_values)) Delta and $(length(rho_values)) rho values to: $(output_path)")

	for rho in rho_values
		rho_folder = joinpath(output_path, "rho_$(float_label(rho))")
		isdir(rho_folder) || mkpath(rho_folder)
		for delta_t in delta_values
			filename = joinpath(rho_folder, "eigenvalues_vs_mu_rho_$(float_label(rho))_Delta_$(float_label(delta_t)).pdf")
			tracker_ratio = plot_rho_delta_tracker ? tracker_rho / delta_t : nothing
			tracker_position = plot_rho_delta_tracker ? (delta_max - delta_t) / (delta_max - delta_min) : nothing
			plt_eigenvalues_generalised(
				df,
				:mu;
				filename=filename,
				xlims=xlims,
				ylims=ylims,
				label_fontsize=label_fontsize,
				tick_fontsize=tick_fontsize,
				plot_mp=plot_mp,
				mp_yrange=mp_yrange,
				plot_rho_delta_tracker=plot_rho_delta_tracker,
				tracker_ratio=tracker_ratio,
				tracker_position=tracker_position,
				tracker_limits=tracker_limits,
				shape_ratio=shape_ratio,
				base_fixed_values...,
				Delta_t=delta_t,
				rho=rho,
			)
			# println("Saved $(filename)")
		end
	end
	return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
	save_eigenvalue_plots_by_mu(
		DATA_PATH, 
		OUTPUT_PATH, 
		xlims=(0.0, 3.0), 
		ylims=(-0.01, 0.01),
		plot_mp=false,
		mp_yrange=(-1.01, 1.01),
		plot_rho_delta_tracker=false,
		tracker_rho=1.5,
		shape_ratio=0.95,
		)
end
