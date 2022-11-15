#!/usr/local/bin julia
# coding=utf-8

# Utility functions for .jl

# Ivan Melchor

using Multitaper


function _crosscorr(data_i::Array{Float64,1}, data_j::Array{Float64,1}, fs::Integer; NW=3.5, pad=1.0)
    
    K = convert(Int64, 2*NW - 1)
    sij = multispec(data_i, data_j, outp=:spec, dt=1/fs, NW=NW, K=K, ctr=true, pad=pad, guts=true)
    Si = sum(sij.coef[1].coef, dims=2)
    Sj = sum(sij.coef[2].coef, dims=2)

    return Si .* conj(Sj)
end


function _csm(data::Array{Float64,2}, fs::Integer, fq_band::Tuple{Float64,Float64}; NW=3.5, pad=1.0)

    # build the cross spectral matrix
    lengt, fftleng, halffreq = Multitaper.output_len(data[1,:], pad)
    freq = fs*range(0,1,length=fftleng+1)[1:halffreq]
    frmin = findmin(abs.(freq.-fq_band[1]))[2]
    frmax = findmin(abs.(freq.-fq_band[2]))[2]
    freq = freq[frmin:frmax]
    nfs = size(freq)[1]

    csm = Array{ComplexF64}(undef, 3, 3, nfs)
    for i in 1:3
        for j in i:3
            csm[i,j,:] = _crosscorr(data[i,:], data[j,:], fs, NW=NW, pad=pad)[frmin:frmax, 1]
            if i != j
                csm[j,i,:] = conj(csm[i,j,:])
            end
        end
    end

    # do singular value decomposition
    csm_svd = map(svd,[csm[:,:,i] for i in 1:size(csm)[3]])

    return freq, csm_svd
end


function polar_analysis(data::Array{Float64,2}, fs::Integer, fq_band::Tuple{Float64,Float64}; NW=3.5, pad=1.0)

    freq, csm_svd = _csm(data, fs, fq_band, NW=NW, pad=pad)

    # compute polar degree
    polar_degree = [(3*sum(s.S.^2) - sum(s.S)^2)/(2*sum(s.S)^2) for s in csm_svd]
    
    # compute rectilinearity of Vt
    z_rot = map(_rotate_vector, [s.Vt[1,:] for s in csm_svd])
    rect = map(_rectiliniarity, z_rot)
    angles = map(_angles, z_rot)
    azimuth = [a[1] for a in angles]
    elevation = [a[2] for a in angles]
    pHH = [a[3] for a in angles]
    pVH = [a[4] for a in angles]

    return PolarT(freq, polar_degree, rect, azimuth, elevation, pHH, pVH)
end


function polar_run(data::Array{Float64,2}, fs::Integer, nwin::Integer, lwin::Integer, nini::Integer, fq_band::Tuple{Float64,Float64}; NW=3.5, pad=1.0)

    #This code processes a higher matrix in overlaps of nadv

    # define base
    base = BaseP(data, fs, fq_band, nwin, lwin, nadv, nini, NW, pad)

    ret = []
    for n in 1:base.nwin
        ninikk = base.nini + base.nadv*(n-1)
        data_n = @view base.data[:, ninikk:ninikk+base.lwin]
        push!(ret, polar_analysis(data_n, base.fs, base.fq_band, NW=base.NW, pd=base.pad))
    end
    
    return ret
end