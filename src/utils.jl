#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using LinearAlgebra


function rotate_vector(z::Array{ComplexF64})
    
    zrot_(a) = [zk * (cos(a) + sin(a)im) for zk in z]
    alla = range(0, 2*pi, 100)
    allz = [zrot_(a) for a in alla]
    comp = [abs(dot([x.re for x in z_rotk],[x.im for x in z_rotk])) for z_rotk in allz]
    phi_min = alla[findmin(comp)[2]]
    return zrot_(phi_min)

end
