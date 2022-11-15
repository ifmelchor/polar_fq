#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

__precompile__

module PolarF

    using LinearAlgebra
    
    export polar_analysis, polar_run

    include("types.jl")

    include("angles.jl")

    include("polar.jl")

end
