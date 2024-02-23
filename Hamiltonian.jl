#construct Hamiltonian with M(horizental)xN(vertical) sites, PBC. t1 set to be 1
function Haldane_p(M::Int, N::Int, t2::Float64, μ::Float64, ϕ::Float64)
    return AB(M,N) .+ t2*AA(M,N,ϕ) .+ μ*chem(M,N)
end

#generate connecting AB sublattices in a M x N graph
function AB(M::Int, N::Int)
    H = zeros(ComplexF64, 2*M*N, 2*M*N)
    for i in 1:N
        l = (i-1)*2*M
        m = mod(l+2*M,N*2*M)
        for j in 1:M
            
            
            H[l+2*j-1,l+2*j]=1
            H[l+2*j,l+2*j-1]=1
            
            H[l+2*j,l+mod(2*j+1,2*M)]=1
            H[l+mod(2*j+1,2*M),l+2*j]=1
            
            H[l+2*j,m+2*j-1]=1
            H[m+2*j-1,l+2*j]=1
        end
    end
    return H
end

#generate connecting AA and BB sublattices in a M x N graph, with direction labeled by ϕ
function AA(M::Int, N::Int, ϕ::Float64)
    H = zeros(ComplexF64, 2*M*N, 2*M*N)
    A = exp(1im*ϕ)
    for i in 1:N
        l = (i-1)*2*M
        m = mod(l+2*M,N*2*M)
        for j in 1:M
            H[l+2*j-1,l+mod(2*j+1,2*M)]=A
            H[l+mod(2*j+1,2*M),l+2*j-1]=1/A
            H[l+2*j,l+mod(2*j+1,2*M)+1]=1/A
            H[l+mod(2*j+1,2*M)+1,l+2*j]=A
            
            H[l+2*j-1,m+2*j-1]=1/A
            H[m+2*j-1,l+2*j-1]=A
            H[l+2*j,m+2*j]=A
            H[m+2*j,l+2*j]=1/A
            
            H[l+2*j-1,m+mod(2*j-3,2*M)]=A
            H[m+mod(2*j-3,2*M),l+2*j-1]=1/A
            H[l+2*j,m+mod(2*j-3,2*M)+1]=1/A
            H[m+mod(2*j-3,2*M)+1,l+2*j]=A
        end
    end
    return H
end

function chem(M::Int, N::Int)
    H = zeros(2*M*N, 2*M*N)
    for i in 1:N
        for j in 1:M
            l = (i-1)*2*M
            H[l+2*j-1,l+2*j-1]=1
            H[l+2*j,l+2*j]=-1
        end
    end
    return H
end

#open boundary condition
function Haldane(M::Int, N::Int, t2::Float64, μ::Float64, ϕ::Float64)
    return ABo(M,N) .+ t2*AAo(M,N,ϕ) .+ μ*chem(M,N)
end

function ABo(M::Int, N::Int)
    H = zeros(ComplexF64, 2*M*N, 2*M*N)
    for i in 1:N-1
        l = (i-1)*2*M
        m = l+2M
        for j in 1:M-1
            H[l+2*j-1,l+2*j]=1
            H[l+2*j,l+2*j-1]=1
            
            H[l+2j,l+2j+1]=1
            H[l+2j+1,l+2j]=1
            
            H[l+2*j,m+2*j-1]=1
            H[m+2*j-1,l+2*j]=1
        end
        H[l+2M-1,l+2M] = 1
        H[l+2M,l+2M-1] = 1

        H[l+2M,m+2M-1] = 1
        H[m+2M-1,l+2M] = 1
    end
    l = 2(N-1)M
    for j in 1:M-1
        H[l+2j-1,l+2j] = 1
        H[l+2j,l+2j-1] = 1

        H[l+2j,l+2j+1] = 1
        H[l+2j+1,l+2j] = 1
    end
    H[l+2M-1,l+2M] = 1
    H[l+2M,l+2M-1] = 1
    return H
end

function AAo(M::Int, N::Int, ϕ::Float64)
    H = zeros(ComplexF64, 2*M*N, 2*M*N)
    A = exp(1im*ϕ)
    for i in 1:N-1
        l = (i-1)*2*M
        m = l+2M
        for j in 1:M-1
            H[l+2*j-1,l+2j+1]=A
            H[l+2j+1,l+2*j-1]=1/A
            H[l+2*j,l+2j+2]=1/A
            H[l+2j+2,l+2*j]=A
            
            H[l+2*j-1,m+2*j-1]=1/A
            H[m+2*j-1,l+2*j-1]=A
            H[l+2*j,m+2*j]=A
            H[m+2*j,l+2*j]=1/A
            
            H[l+2*j+1,m+2j-1]=A
            H[m+2j-1,l+2*j+1]=1/A
            H[l+2*j+2,m+2j]=1/A
            H[m+2j,l+2*j+2]=A
        end
        H[l+2M-1,m+2M-1] = 1/A
        H[m+2M-1,l+2M-1] = A
        H[l+2M,m+2M] = A
        H[m+2M,l+2M] = 1/A
    end
    l = 2(N-1)M
    for j in 1:M-1
        H[l+2j-1,l+2j+1] = A
        H[l+2j+1,l+2j-1] = 1/A
        H[l+2*j,l+2j+2]=1/A
        H[l+2j+2,l+2*j]=A
    end
    return H
end

#SSH model Hamiltonian. #unit cells: N. #hopping between unit cell:t.
#PBC
function SSHp(t::Float64, v::Float64, N::Int)
    H = zeros(2*N,2*N)
    for i in 1:N-1
        H[2*i-1,2*i] = v
        H[2*i,2*i-1] = v
        H[2*i,2*i+1] = t
        H[2*i+1,2*i] = t
    end
    H[2*N-1,2*N] = v
    H[2*N,2*N-1] = v
    H[1,2*N] = t
    H[2*N,1] = t
    return H
end

#obc
function SSHo(t::Float64, v::Float64, N::Int)
    H = zeros(2*N,2*N)
    for i in 1:N-1
        H[2*i-1,2*i] = v
        H[2*i,2*i-1] = v
        H[2*i,2*i+1] = t
        H[2*i+1,2*i] = t
    end
    H[2*N-1,2*N] = v
    H[2*N,2*N-1] = v
    return H
end
