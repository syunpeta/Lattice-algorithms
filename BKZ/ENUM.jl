using LinearAlgebra
include("utils.jl")
function return_r(k::Int64,r0::Int64,r::Array{Int64})
    if k == 0
        return r0
    else
        return r[k]
    end
end
        
function ENUM(B::Matrix{Int64})
    N,_ = size(B)
    B_p,mu =GSO(B,N) 
    norm_Bp = zeros(Float64,N)
    for i in 1:N
        norm_Bp[i] = norm(B_p[i,:])^2
    end
    
    ##R設定
    R = 0.99*norm_Bp[1]
    ##
    
    sigma = zeros(Float64,(N+1,N))
    r = [i for i in 1:N]
    r0 = 0
    rho = zeros(Float64,N+1)
    v = zeros(Int64,N)
    v[1] = 1
    c = zeros(Float64,N)
    w = zeros(Int64,N)
    last_nonzero = 1
    k = 1
    while true
        rho[k] = rho[k+1] + ((v[k] - c[k])^2) * norm_Bp[k]
        if rho[k] <= R
            if k == 1
                return v
            else
                k-=1
                ##r
                if k== 1
                    r0 = max(r0,r[k])
                else
                    r[k-1] = max(r[k-1],r[k])
                end
                ##
                for i in return_r(k,r0,r):-1:(k+1)
                    sigma[i,k] = sigma[i+1,k] + v[i]*mu[i,k] 
                end
                c[k] = -sigma[k+1,k]
                v[k] = round(c[k])
                w[k] = 1
            end
        else
            k+=1
            if k == (N+1)
                return false
            end
            ##r
            if k==1
                r0 = 1
            else
                r[k-1] = k
            end
            ##
            if k >= last_nonzero
                last_nonzero = k
                v[k]+=1
            else
                if v[k] > c[k]
                    v[k] -= w[k]
                else
                    v[k] += w[k]
                end
                w[k] += 1
            end
        end
    end
end           
