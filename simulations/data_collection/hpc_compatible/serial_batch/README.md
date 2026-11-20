# Serial Batch Jobs

## Sweep progress logging

The HPC `UserOptions` includes `logging` and `logging_interval` settings for
the generic normal- and high-precision solvers. Set `logging=false` to silence
progress messages. With logging enabled, a positive `logging_interval` prints
after that many completed parameter points; `0` selects an automatic interval
of approximately 1% of the row's sweep. The generic solvers report the row,
completion count, elapsed time, average time per point, ETA and final runtime.

`main_Delta.jl` currently sets `logging=true` and `logging_interval=1`, so that
run reports every completed point. The other serial-batch entry points use
`logging_interval=0` for approximately 1% updates. Every-point logging makes
the Slurm output file substantially larger for long sweeps; increase the
interval for less frequent output.

## Julia cache compatibility on Slurm

Some array tasks failed during Julia package loading with `Unable to find
compatible target in cached code image` and `Target 0 (znver2): Rejecting this
target due to use of runtime-disabled features`. The cached native code was
incompatible with the runtime on those nodes. The array tasks share a Julia
depot, and different compute nodes or runtime CPU features can make a cache
built for one target unsuitable elsewhere. This happens before the solver runs;
it is not evidence that the affected parameter rows are invalid.

`runscripts/serial_runscript.sh` now sets `JULIA_CPU_TARGET=generic` before
starting Julia. This is an attempted portability fix: it requests a generic
target for Julia's package images, which may sacrifice some CPU-specific
optimizations. It does not repair existing cache files or guarantee that a
custom Julia system image or external library is portable. The change only
affects subsequently submitted jobs.

Before a large array, warm the required package caches once on an allocated
compute node using the same Julia module, depot, CPU target, working directory,
and Julia environment as the array. For the current launcher, run the following
from `runscripts/`, adjusting the partition, account, time, and memory to the
cluster's requirements:

```bash
sbatch --partition=compute --account=phys030424 --ntasks=1 --cpus-per-task=1 --time=00:30:00 --mem=5G \
	--wrap="module add languages/julia && export JULIA_DEPOT_PATH=/user/work/hb21877/.julia JULIA_CPU_TARGET=generic JULIA_NUM_THREADS=1 && julia -e 'using GenericLinearAlgebra, DataFrames, JLD2, KrylovKit, BSON, IterTools; println(\"Package warm-up complete\")'"
```

Wait for this job to complete successfully and inspect its Slurm output before
submitting the array. This imports packages without running the sweep; it can
still take time and may compile packages on the allocated node. The current
serial-batch `Project.toml` has no dependencies, so `Pkg.precompile()` in that
project alone is not a substitute for importing the packages used by the job.
Keep the module version and environment consistent between the warm-up and the
array; repeat the warm-up after changing them or updating dependencies.

If the incompatible-image error persists, capture the complete error, Julia
version, module, node, CPU target and depot, and consult cluster support. Do
not delete the shared depot or its caches while other jobs are using them. A
separate clean, populated depot or a site-supported Julia build may be needed.
First test a small array across representative compute nodes before restarting
the full sweep.

A Slurm output file containing only a CPU-binding message does not establish a
hang. Startup or shared-filesystem contention is possible, but requires
separate investigation using Slurm job status and node diagnostics.
