include("LLL.jl")
include("utils.jl")
using LinearAlgebra

function  PotLLL(BB::Matrix{Int64},delta::Float64)
    N,_ = size(BB)
    BB = LLL(BB,delta)
    l = 1
    GS,mu = GSO(BB,N)

    while l <= N
        for j in l-1:-1:1
            part_Size_reduce(BB,mu,l,j)
        end
        GS,mu = GSO(BB,N)
        P = 1
        Pmin = 1
        k = 1

        for j in l-1:-1:1
            sum_tmp = 0
            for i in j:l-1
                sum_tmp += (mu[l,i]^2)*norm(GS[i,:])^2
            end
            
            P = P*((norm(GS[l,:])^2  + sum_tmp)/norm(GS[j,:])^2)
            println(P)
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
        else
            l += 1
        end
    end

    return BB
    
end

BB = [63 -14 -1 84 61;74 -20 23 -32 -52;93 -46 -19 0 -63;93 11 13 60 52;33 -93 12 57 -2]
PotLLL(BB,0.99)

BB = [63 -14 -1 84 61;74 -20 23 -32 -52;93 -46 -19 0 -63;93 11 13 60 52;33 -93 12 57 -2]
PotLLL(BB,0.99)
LLL(BB,0.99)