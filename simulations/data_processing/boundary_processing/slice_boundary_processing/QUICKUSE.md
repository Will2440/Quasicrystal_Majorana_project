cd /user/work/hb21877/Quasicrystal_Majorana_project/simulations/data_processing/boundary_processing/slice_boundary_processing

export RUN=/user/work/hb21877/Quasicrystal_Majorana_project/simulations/data_collection/hpc_compatible/dynamic_slice_search/runs/dynamic_mb_slice_search_20261105_163222

julia export_slice_boundaries.jl "$RUN"
julia plot_slice_boundary.jl "$RUN" --save "$RUN/boundary_overview.png" --points
julia plot_slice_progression.jl "$RUN" --save "$RUN/progression_frames" --interval 5