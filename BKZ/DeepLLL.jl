using LinearAlgebra
include("utils.jl")
function GSOUp_DeepLLL(mu::Matrix{Float64},norm_B::Array{Float64},i::Int64,k::Int64,N::Int64)::Tuple{Matrix{Float64},Array{Float64}}
    P = zeros(Float64,N)
    D = zeros(Float64,N)
    P[k] = norm_B[k]
    D[k] = norm_B[k]
    
    for j in (k-1):-1:i
        P[j] = mu[k,j]*norm_B[j]
        D[j] = D[j+1] + mu[k,j]*P[j]
    end
    
    S = zeros(Float64,N)
    for j in k:-1:(i+1)
        T = mu[k,j-1]/D[j]
        for l in N:-1:(k+1)
            S[l] +=mu[l,j]*P[j]
            mu[l,j] = mu[l,j-1]- T*S[l]
        end
        for l in k:-1:(j+1)
            S[l]+=mu[l-1,j]*P[j]
            mu[l,j] = mu[l-1,j-1] - T*S[l]
        end
    end
    T = 1/D[i]
    
    for l in N:-1:(k+1)
        mu[l,i] = T*(S[l]+mu[l,i]*P[i])
    end
    for l in k:-1:(i+2)
        mu[l,i] = T*(S[l] +mu[l-1,i]*P[i])
    end
    mu[i+1,i] = T*P[i]
    for j in 1:(i-1)
        ips = mu[k,j]
        for l in k:-1:(i+1)
            mu[l,j] = mu[l-1,j]
        end
        mu[i,j] = ips
    end
    
    for j in k:-1:(i+1)
        norm_B[j] = D[j] * norm_B[j-1]/D[j-1]
    end
    norm_B[i] = D[i]
    
    return mu,norm_B
end


function DeepLLL(B::Matrix{Int64},delta::Float64)::Matrix{Int64}
    N,_ = size(B)
    B_p,mu = GSO(B,N)
    norm_B = zeros(Float64,N)
    for p in 1:N
        norm_B[p] = LinearAlgebra.norm(B_p[p,:])^2
    end
    k=2
    while k <= N
        for j in (k-1):-1:1
            #Size-reduce
            if abs(mu[k,j])>0.50
                q = round(mu[k,j])
                B[k,:] -= q*B[j,:]
                for l in 1:j
                    mu[k,l] -=q*mu[j,l]
                end
            end
        end
        C = LinearAlgebra.norm(B[k,:])^2
        i=1
        while i < k
            if C >= delta*norm_B[i]
                C-= (mu[k,i]^2)*norm_B[i]
                i +=1
            else
                #　DeepInsertion
                v = B[k,:]
                for j in k:-1:i+1
                    B[j,:] = B[j-1,:]
                end
                B[i,:] = v
                mu,norm_B = GSOUp_DeepLLL(mu,norm_B,i,k,N)
                k = max(i,2) -1
            end
        end
        k+=1
    end
    return B
end
