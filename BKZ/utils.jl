using LinearAlgebra
function GSO(B, N::Int64)::Tuple{Matrix{Float64},Matrix{Float64}}
    B_p = zeros(Float64,size(B))
    mu = zeros(Float64,(N,N))
    for i in 1:N
        B_p[i,:] = B[i,:]
        mu[i,i] = 1
        for j in 1:i-1
            mu[i,j] = B[i,:]⋅B_p[j,:]/norm(B_p[j,:])^2
            B_p[i,:] -= mu[i,j]*B_p[j,:]
        end
    end
    return B_p,mu
end

function part_Size_reduce(B::Matrix{Int64},mu::Matrix{Float64},i::Int64,j::Int64)::Tuple{Matrix{Int64},Matrix{Float64}}
    if abs(mu[i,j]) > 0.50
        q = round(mu[i,j])
        B[i,:] =B[i,:] -  q*B[j,:]
        for l in 1:j
            mu[i,l] = mu[i,l] -  q*mu[j,l]
        end
    end
    return B,mu
end

function Size_reduce(B::Matrix{Int64})::Tuple{Matrix{Int64},Matrix{Float64}}
    N,_ = size(B)
    B_p,mu = GSO(B,N)
    for i in 2:N
        for j in i-1:-1:1
            part_Size_reduce(B,mu,i,j)
        end
    end
    return B,mu
end
