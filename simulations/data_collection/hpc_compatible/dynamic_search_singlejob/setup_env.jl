#!/usr/bin/env julia

using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.precompile()

println("DynamicSearch single-job environment ready at: ", @__DIR__)
