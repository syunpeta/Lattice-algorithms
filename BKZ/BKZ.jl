include("utils.jl")
include("LLL.jl")
include("ENUM.jl")
include("MLLL.jl")
using LinearAlgebra

function BKZ(B::Matrix{Int64},beta::Int64,delta::Float64)
    N,_ = size(B)
    B= LLL(B,delta)
    BB,mu = GSO(B,N)
    BB_norm = zeros(Float64,N)
    for i in 1:N
        BB_norm[i] = norm(BB[i,:])^2
    end

    z = 0
    k = 0
    
    while z < (N-1)
        k = mod(k+1,N-1)+1
        l = min(k+beta-1,N)
        h = min(l+1,N)
        R = 0.99*BB_norm[k]
        println("k,l,h,z,R",k,l,h,z,R)

        #μ[k,l],Bk...Bl
        mu_kl = zeros(Float64,l-k+1,l-k+1)
        BB_norm_kl = zeros(Float64,l-k+1)
        for i in k:l
            for j in k:l
                mu_kl[i-k+1,j-k+1] = mu[i,j]
            end
            BB_norm_kl[i-k+1] = BB_norm[i]
        end

        v1 = ENUM(mu_kl,BB_norm_kl,R)
        v = zeros(Int64,N)
        if !iszero(v1)
            for i in k:l
                v += B[i,:]*v1[i-k+1]
            end
            println("v",v)
        end

        if !iszero(v)
            z =0
            C = zeros(Int64,h+1,N)
            for i in 1:k-1
                C[i,:] = B[i,:]
            end
            C[k,:] = v
            for i in (k+1):(h+1)
                C[i,:] = B[i-1,:]
            end
            C= MLLL(C,delta)
            for i in 1:h
                B[i,:] = C[i,:]
            end
            BB ,mu = GSO(B,N)
            BB_norm = zeros(Float64,N)
            for i in 1:N
                BB_norm[i] = norm(BB[i,:])^2
            end
        else
            z += 1
            B = LLL(B,delta)
            BB,mu = GSO(B,N)
            BB_norm = zeros(Float64,N)
            for i in 1:N
                BB_norm[i] = norm(BB[i,:])^2
            end
        end
    end
    return B
end



#BB = [63 -14 -1 84 61;74 -20 23 -32 -52;93 -46 -19 0 -63;93 11 13 60 52;33 -93 12 57 -2]
BB = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 41 62 98 112 90 61 101 42 74 11 40 62 101 84 110 90 ;0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 90 41 62 98 112 90 61 101 42 74 11 40 62 101 84 110 ;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 110 90 41 62 98 112 90 61 101 42 74 11 40 62 101 84 ;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 84 110 90 41 62 98 112 90 61 101 42 74 11 40 62 101 ;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 101 84 110 90 41 62 98 112 90 61 101 42 74 11 40 62 ;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 62 101 84 110 90 41 62 98 112 90 61 101 42 74 11 40 ;0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 40 62 101 84 110 90 41 62 98 112 90 61 101 42 74 11 ;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 11 40 62 101 84 110 90 41 62 98 112 90 61 101 42 74 ;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 74 11 40 62 101 84 110 90 41 62 98 112 90 61 101 42 ;0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 42 74 11 40 62 101 84 110 90 41 62 98 112 90 61 101 ;0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 101 42 74 11 40 62 101 84 110 90 41 62 98 112 90 61 ;0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 61 101 42 74 11 40 62 101 84 110 90 41 62 98 112 90 ;0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 90 61 101 42 74 11 40 62 101 84 110 90 41 62 98 112 ;0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 112 90 61 101 42 74 11 40 62 101 84 110 90 41 62 98 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 98 112 90 61 101 42 74 11 40 62 101 84 110 90 41 62 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 62 98 112 90 61 101 42 74 11 40 62 101 84 110 90 41 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 0 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 131 ]
show(stdout, "text/plain", BB)
BB =BKZ(BB,5,0.99)
show(stdout, "text/plain", BB)

