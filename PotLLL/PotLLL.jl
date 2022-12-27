include("LLL.jl")
include("utils.jl")
using LinearAlgebra

function GS_UP_l(BB,mu,GS,l,N)
    GS[1,:] = BB[1,:]
    for i in 2:l
        tmp = zeros(Float64,N)
        for j in 1:i-1
            tmp += mu[i,j]*GS[j,:]
        end
        GS[i,:] = BB[i,:] - tmp
    end
    return GS
end




function  PotLLL(BB::Matrix{Int64},delta::Float64)
    N,_ = size(BB)
    BB = LLL(BB,delta)
    l = 1
    GS,mu = GSO(BB,N)

    while l <= N
        for j in l-1:-1:1
            BB,mu = part_Size_reduce(BB,mu,l,j)
        end
        GS= GS_UP_l(BB,mu,GS,l,N)
        P = 1
        Pmin = 1
        k = 1

        for j in l-1:-1:1
            sum_tmp = 0
            for i in j:l-1
                sum_tmp += (mu[l,i]^2)*norm(GS[i,:])^2
            end
            
            P = P*((norm(GS[l,:])^2  + sum_tmp)/norm(GS[j,:])^2)
            if P < Pmin
                k = j
                Pmin = P
            end
        end

        if delta > Pmin
            v = BB[l,:]
            for j in l:-1:k+1
                BB[j,:] = BB[j-1,:]
            end
            BB[k,:] = v

            GS,mu = GSO(BB,N)
            l=k
        else
            l += 1
        end
    end

    return BB
    
end

