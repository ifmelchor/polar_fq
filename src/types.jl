#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor


struct BaseP
    data :: Array{Float64}
    fs :: Integer
    fq_band :: Tuple{Float64, Float64}
    nwin :: Integer
    lwin :: Integer
    nadv :: Float64
    nini :: Integer
    NW :: Float64
    pad :: Float64
end

struct PolarT
    freq :: Array{Float64}
    polar_degree :: Array{Float64}
    rect :: Array{Float64}
    azimuth :: Array{Float64}
    elev  :: Array{Float64}
    phyHH :: Array{Float64}
    phyVH :: Array{Float64}
end
