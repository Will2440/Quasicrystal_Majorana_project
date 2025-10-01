"""
    file name:   param_restriction_cuts.jl
    created:     30/09/2025
    last edited: 30/09/2025

    overview:
        Contains all of the cuts::Vector{Dict{Any}} used to restrict the 2D parameter space for the :restricted solver type in local_machine/main.jl.
        See .png's in auxilliary/param_restriction_examples for the plots produced by these cuts.
"""

include("param_comb_gen.jl")
using .ParamCombGen


# GQC optimised for Delta=2.0 maximum
GQC_D20_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => -2.0,
        :x_range => (2.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(35.0),
        :y_intercept => 5.0,
        :x_range => (2.5, maximum(xs)),
        :y_range => (7.5, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -10.0,
        :y_intercept => 30.0,
        :x_range => (minimum(xs), 2.5),
        :y_range => (7.5, maximum(ys)),
        :cut_which_side => "above"
    )
]

# GQC optimised for Delta=1.0
GQC_D10_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => -2.0,
        :x_range => (1.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => 2.0,
        :x_range => (2.5, maximum(xs)),
        :y_range => (5.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -7.5,
        :y_intercept => 15.0,
        :x_range => (minimum(xs), 2.5),
        :y_range => (5.0, maximum(ys)),
        :cut_which_side => "above"
    )
]

# GQC optimised for Delta=0.1 minimum
GQC_D01_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => -2.0,
        :x_range => (1.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => 2.0,
        :x_range => (1.0, maximum(xs)),
        :y_range => (2.5, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -20.0,
        :y_intercept => 20.0,
        :x_range => (minimum(xs), 2.5),
        :y_range => (4.0, maximum(ys)),
        :cut_which_side => "above"
    )
]

# SQC optimised for Delta=0.1
SQC_D01_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => -2.0,
        :x_range => (1.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient =>  ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => 1.5,
        :x_range => (2.0, maximum(xs)),
        :y_range => (3.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -20.0,
        :y_intercept => 40.0,
        :x_range => (minimum(xs), 2.5),
        :y_range => (4.0, maximum(ys)),
        :cut_which_side => "above"
    )
]

# SQC optimised for Delta=1.0
SQC_D10_cute = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(47.0),
        :y_intercept => -2.5,
        :x_range => (1.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => 2.0,
        :x_range => (2.5, maximum(xs)),
        :y_range => (3.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -10.0,
        :y_intercept => 25.0,
        :x_range => (minimum(xs), 2.5),
        :y_range => (5.0, maximum(ys)),
        :cut_which_side => "above"
    )
]

# SQC optimised for Delta=2.0
SQC_D20_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(47.0),
        :y_intercept => -2.5,
        :x_range => (1.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => 2.0,
        :x_range => (2.5, maximum(xs)),
        :y_range => (8.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -5.0,
        :y_intercept => 20.0,
        :x_range => (minimum(xs), 2.5),
        :y_range => (8.0, maximum(ys)),
        :cut_which_side => "above"
    )
]

# TMQC optimised for Delta=2.0
TMQC_D20_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(35.0),
        :y_intercept => -1.5,
        :x_range => (1.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(40.0),
        :y_intercept => 3.0,
        :x_range => (3.5, maximum(xs)),
        :y_range => (7.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -5.0,
        :y_intercept => 22.0,
        :x_range => (minimum(xs), 3.5),
        :y_range => (7.0, maximum(ys)),
        :cut_which_side => "above"
    )
]

# TMQC optimised for Delta=1.0
TMQC_D10_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(35.0),
        :y_intercept => -1.5,
        :x_range => (1.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(40.0),
        :y_intercept => 2.5,
        :x_range => (2.5, maximum(xs)),
        :y_range => (5.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -5.0,
        :y_intercept => 17.0,
        :x_range => (minimum(xs), 2.5),
        :y_range => (5.0, maximum(ys)),
        :cut_which_side => "above"
    )
]

# TMQC optimised for Delta=0.1
TMQC_D01_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(37.0),
        :y_intercept => -1.5,
        :x_range => (1.5, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => 1.5,
        :x_range => (1.5, maximum(xs)),
        :y_range => (3.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    # Dict(
    #     :gradient => -15.0,
    #     :y_intercept => 17.0,
    #     :x_range => (minimum(xs), 2.5),
    #     :y_range => (5.0, maximum(ys)),
    #     :cut_which_side => "above"
    # )
]

# PQC optimised for Delta=0.5
PQC_D05_uts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(38.0),
        :y_intercept => -1.5,
        :x_range => (6.0, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(40.0),
        :y_intercept => 1.5,
        :x_range => (6.0, maximum(xs)),
        :y_range => (3.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    # Dict(
    #     :gradient => -15.0,
    #     :y_intercept => 17.0,
    #     :x_range => (minimum(xs), 2.5),
    #     :y_range => (5.0, maximum(ys)),
    #     :cut_which_side => "above"
    # )
]

# PQC optimised for Delta=1.5
PQC_D15_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(35.0),
        :y_intercept => -1.5,
        :x_range => (5.0, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => 2.0,
        :x_range => (5.0, maximum(xs)),
        :y_range => (3.0, maximum(ys)),
        :cut_which_side => "above"
    ),
    Dict(
        :gradient => -15.0,
        :y_intercept => 70.0,
        :x_range => (minimum(xs), 5.0),
        :y_range => (7.0, maximum(ys)),
        :cut_which_side => "above"
    )
]

# PQC optimised for Delta=2.0
PQC_D20_cuts = [
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(35.0),
        :y_intercept => -1.5,
        :x_range => (5.0, maximum(xs)),
        :y_range => (0.0, maximum(ys)),
        :cut_which_side => "below"
    ),
    Dict(
        :gradient => ParamCombGen.angle_to_gradient(45.0),
        :y_intercept => 2.0,
        :x_range => (5.0, maximum(xs)),
        :y_range => (7.5, maximum(ys)),
        :cut_which_side => "above"
    ),
    # Dict(
    #     :gradient => -15.0,
    #     :y_intercept => 70.0,
    #     :x_range => (minimum(xs), 5.0),
    #     :y_range => (7.0, maximum(ys)),
    #     :cut_which_side => "above"
    # )
]
