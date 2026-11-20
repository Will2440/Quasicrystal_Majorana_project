# Merged dynamic slice-search run

- First source: `/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/hpc_compatible/dynamic_slice_search/runs/dynamic_mb_slice_search_20261106_093019`
- Priority source: `/Users/Will/Documents/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/hpc_compatible/dynamic_slice_search/runs/dynamic_mb_slice_search_20261106_120131`
- Priority rho interval: [2.25, 4.0] 
- Retained from first source: 48 slices
- Dropped from first source inside priority interval: 32 slices
- Retained from priority source: 79 slices

The slice summary, transition table, and `slices/slice_NNNN` folders are merged and reindexed together in increasing rho order. Nonzero transition track IDs from the priority source are offset to prevent unrelated run-local IDs from being connected by plotting tools. Source plots and restart checkpoints are not copied because they do not represent the merged data.
