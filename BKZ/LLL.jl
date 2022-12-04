using LinearAlgebra
include("utils.jl")

function GSOUp_LLL(mu::Matrix{Float64},norm_Bp::Array{Float64},k::Int64,N::Int64)::Tuple{Matrix{Float64},Array{Float64}}
    nu = mu[k,k-1]
    beta = norm_Bp[k] + (nu^2)*norm_Bp[k-1]
    mu[k,k-1] = nu*norm_Bp[k-1]/beta
    norm_Bp[k] = norm_Bp[k]*norm_Bp[k-1]/beta
    norm_Bp[k-1] = beta
    for j in 1:(k-2)
        mu[k-1,j],mu[k,j] = mu[k,j],mu[k-1,j]
    end
    for i in (k+1):N
        t = mu[i,k]
        mu[i,k] = mu[i,k-1] - nu*t
        mu[i,k-1] = t + mu[k,k-1]*mu[i,k]
    end
    return mu,norm_Bp
end
    
function LLL(B::Matrix{Int64},delta::Float64)
    N,_ = size(B)
    B_p,mu = GSO(B,N)
    norm_Bp = zeros(Float64,N)
    for i in 1:N
        norm_Bp[i] = norm(B_p[i,:])^2
    end
    k = 2
    
    while k <= N
        for j in (k-1):-1:1
            B,mu = part_Size_reduce(B,mu,k,j)
        end
        if norm_Bp[k] >= (delta - mu[k,k-1]^2)*norm_Bp[k-1]
            k = k+1
        else
            B[k-1,:],B[k,:] = B[k,:],B[k-1,:]
            mu,norm_Bp = GSOUp_LLL(mu,norm_Bp,k,N)
            k = max(k-1,2)
        end
    end
    return B,mu,norm_Bp
end
