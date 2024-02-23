using LinearAlgebra


#Given modular Hamiltonian, compute expectation value of exp(iθQ)
function disorder_operator(H::Matrix, θ::Float64)
    v, _ = eigen(Hermitian(H))
    l = length(v)
    ee = 1
    for i in 1:l
        v1 = v[i]
        ee = ee*(1+exp(-v1+1im*θ))/(1+exp(-v1))
    end
    return ee
end

function disorder_operators(H::Matrix, N::Int)
    es = zeros(ComplexF64,N)
    v, _ = eigen(Hermitian(H))
    l = length(v)
    for i in 1:N
        ee = 1
        θ = 2*pi*(i-1)/N
        for j in 1:l
            v1 = v[j]
            ee = ee*(1+exp(-v1+1im*θ))/(1+exp(-v1))
        end
        es[i] = ee
    end
    return es
end

#real part of FCS generating function
function real_fcs(H::Matrix, N::Int)
    es = zeros(N)
    v, _ = eigen(Hermitian(H))
    l = length(v)
    for i in 1:N
        ee = 0
        θ = 2*pi*(i-1)/N
        for j in 1:l
            v1 = v[j]
            ev = exp(-v1)
            ee = ee + log(sqrt(1+2*ev*cos(θ)+ev^2)/(1+ev))
        end
        es[i] = ee
    end
    return es
end

#value of real_fcs in the range of d⊆(-π,π)
function real_fcs_sub(H::Matrix, d::Vector)
    es = zeros(length(d))
    v, _ = eigen(Hermitian(H))
    l = length(v)
    for i in 1:length(d)
        ee = 0
        θ = d[i]
        for j in 1:l
            v1 = v[j]
            ev = exp(-v1)
            ee = ee + log(sqrt(1+2*ev*cos(θ)+ev^2)/(1+ev))
        end
        es[i] = ee
    end
    return es
end

#first derivative of real_fcs
function real_fcs_d1(H::Matrix, N::Int)
    es = zeros(N)
    v, _ = eigen(Hermitian(H))
    l = length(v)
    for i in 1:N
        ee = 0
        θ = 2*pi*(i-1)/N
        for j in 1:l
            v1 = v[j]
            ev = exp(v1)
            ee = ee - (ev*sin(θ))/((ev+cos(θ))^2+sin(θ)^2)
        end
        es[i] = ee
    end
    return es
end

#second derivative
function real_fcs_d2(H::Matrix, d::Vector)
    es = zeros(length(d))
    v, _ = eigen(Hermitian(H))
    l = length(v)
    for i in 1:length(d)
        ee = 0
        θ = d[i]
        for j in 1:l
            v1 = v[j]
            ev = exp(-v1)
            ee = ee - ev*(cos(θ)*(1+ev^2)+2*ev)/(1+ev^2+2ev*cos(θ))^2
        end
        es[i] = ee
    end
    return es
end

#imaginary part of fcs generating function
function imag_fcs(H::Matrix, N::Int)
    es = zeros(N)
    v, _ = eigen(Hermitian(H))
    l = length(v)
    for i in 1:N
        ee = 0
        θ = 2*pi*(i-1)/N
        for j in 1:l
            v1 = v[j]
            ev = exp(v1)
            ee = ee + atan(sin(θ)/(ev+cos(θ)))
        end
        es[i] = ee
    end
    return es
end

function imag_fcs_sub(H::Matrix, d::Vector)
    es = zeros(length(d))
    v, _ = eigen(Hermitian(H))
    l = length(v)
    for i in 1:length(d)
        ee = 0
        θ = d[i]
        for j in 1:l
            v1 = v[j]
            ev = exp(v1)
            ee = ee + atan(sin(θ)/(ev+cos(θ)))
        end
        es[i] = ee
    end
    return es
end
