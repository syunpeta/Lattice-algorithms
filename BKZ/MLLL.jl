using LinearAlgebra
include("utils.jl") 

function MLLL(B::Matrix{Int64},δ::Float64)::Matrix{Int64}
    
    h,n = size(B)
    z,g = h,1
    B_p = zeros(Float64,(h,n))    
    norm_Bp = zeros(Float64,h)
    mu= zeros(Float64,(h,h))
        
    while g <= z
        if iszero(B[g,:])
            if g < z
                B[g,:],B[z,:] = B[z,:],B[g,:]
            end
            z -=1
        end
        B_p[g,:] = B[g,:]
        for j in 1:(g-1)
            mu[g,j] = (B[g,:]⋅B_p[j,:])/norm_Bp[j]
            B_p[g,:] -= mu[g,j]*B_p[j,:]
        end
        norm_Bp[g] = norm(B_p[g,:])^2
        mu[g,g] = 1
        if g == 1
            g = 2
        else
            l,k = g,g
            startgain = false
            while k <= l && !startgain
                #Size-reduce
                B,mu = part_Size_reduce(B,mu,k,k-1)

                ν = mu[k,k-1]
                temp = norm_Bp[k]+(ν^2)*norm_Bp[k-1]
                
                if temp >= δ*norm_Bp[k-1]
                    for j in (k-2):-1:1
                       B,mu = part_Size_reduce(B,mu,k,j)
                    end
                    k+=1
                else
                    if iszero(B[k,:])
                        if k < z
                            B[z,:],B[k,:] = B[k,:],B[z,:]
                        end
                        z -= 1
                        g = k
                        startgain = true
                    else
                        B[k-1,:],B[k,:] = B[k,:],B[k-1,:]
                        for j in 1:k-2
                            mu[k,j],mu[k-1,j] = mu[k-1,j],mu[k,j]
                        end
                        if temp != 0
                            if norm_Bp[k] == 0
                                norm_Bp[k-1] = temp
                                Bp[k-1,:] *= ν
                                mu[k,k-1] = 1/ν
                                for i in (k+1):l
                                    mu[i,k-1] /= ν
                                end
                            else
                                t = norm_Bp[k-1]/temp
                                mu[k,k-1] = ν*t
                                w = B_p[k-1,:]
                                B_p[k-1,:] = B_p[k,:] + ν*w
                                norm_Bp[k-1] = temp
                                if k <= l
                                    B_p[k,:] = -mu[k,k-1]*B_p[k,:] + (norm_Bp[k]/temp)*w
                                    norm_Bp[k] *=t
                                end
                                for i in (k+1):l
                                    t = mu[i,k]
                                    mu[i,k] = mu[i,k-1] - ν*t
                                    mu[i,k-1] = t + mu[k,k-1]*mu[i,k]
                                end
                            end
                        else
                            norm_Bp[k],norm_Bp[k-1] = norm_Bp[k-1],norm_Bp[k]
                            Bp[k,:],Bp[k-1,:] = Bp[k-1,:],Bp[k,:]
                            for i in (k+1):l
                                mu[i,k],mu[i,k-1] = mu[i,k-1],mu[i,k]
                            end
                        end
                        
                        k = max(k-1,2)
                    end
                end
            end
            if !startgain
                g+=1
            end
        end
      
    end
    return B
end

