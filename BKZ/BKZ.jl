include("DeepLLL.jl")
include("utils.jl")
include("LLL.jl")
include("MLLL.jl")
include("ENUM.jl")

using LinearAlgebra

function BKZ(B::Matrix{Int64},beta::Int64,delta::Float64)
    N,_ = size(B)
    B,mu,norm_Bp = LLL(B,delta)
    z = 0
    k = 0
end