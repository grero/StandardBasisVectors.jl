module StandardBasisVectors
using LinearAlgebra

export standardize_basis

struct RotationMatrix{T<:Real}
    R::Matrix{T}
    Rt::Matrix{T}
    R2::Matrix{T}
    didx::AbstractVector{Int64}
    N::Int64
end

function RotationMatrix(N::Integer)
    R = diagm(fill(1.0,N)) 
    Rt = diagm(fill(1.0,N)) 
    R2 = fill(0.0, N,N)
    didx = diagind(Rt)
    RotationMatrix(R,Rt,R2,didx,N)
end

function reset!(M, v, didx)
    fill!(M, 0.0)
    view(M,didx) .= v
end

function rot(N::Integer, i,j,θ::T) where T <: Real
    R = diagm(fill(one(T),N) )
    rot!(R, i,j,θ)
    return R
end

function rot!(R::AbstractMatrix{T}, i,j, θ) where T <: Real
    s = sin(θ)
    c = cos(θ)
    R[i,i] = c
    R[j,j] = c
    R[i,j] = -s
    R[j,i] = s
    nothing
end

function get_rotation_matrix(N::Integer, i::Integer, θ::Vector{T}) where T <:Real
    R = diagm(fill(one(T),N)) 
    Rt = diagm(fill(one(T),N)) 
    didx = diagind(Rt)
    for j in range(N, stop=i+1, step=-1)
        reset!(Rt, 1.0, didx)
        rot!(Rt, i,j,θ[j])
        R .= Rt*R
    end
    R
end

function get_rotation_matrix!(R::RotationMatrix{T}, i::Integer, θ::Vector{T}) where T <: Real
    reset!(R.R, 1.0, R.didx)
    for j in range(R.N, stop=i+1, step=-1)
        reset!(R.Rt, 1.0, R.didx)
        rot!(R.Rt, i,j,θ[j])
        mul!(R.R2, R.Rt, R.R)
        copy!(R.R, R.R2)
    end
end

function get_rotation_angles(N::Integer,i::Integer,v::AbstractVector{T}) where T <: Real
    θ = zeros(T,N)
    r = one(T)
    θ[N] = asin(v[N])
    for j in range(N-1, stop=i+1, step=-1)
        y = v[j]
        r = r*cos(θ[j+1])
        if r != 0
            θ[j] = asin(y/r)
        else
            θ[j] = 0.0
        end
    end
    θ
end

function reduce_dimensions!(Vo::AbstractMatrix{T}, R::RotationMatrix{T}, i::Integer, V::AbstractMatrix{T}) where T <: Real
    if i == R.N
        copy!(Vo, V)
        return zeros(T, R.N)
    end
    v = view(V,:,i)
    θ = get_rotation_angles(R.N,i,v)
    get_rotation_matrix!(R, i, θ)
    mul!(Vo, R.R',V)
    return θ
end

function reduce_dimensions(N::Integer, i::Integer, V::AbstractMatrix{T}) where T <: Real
    if i == N
        return V, zeros(T, N)
    end
    v = view(V,:,i)
    θ = get_rotation_angles(N,i,v)
    R = get_rotation_matrix(N, i, θ)
    return R'*V, θ
end

function orient_vectors(V::AbstractMatrix{T}) where T <: Real
    size(V,1) == size(V,2) || error("Matrix should be square")
    _V = similar(V)
    _Vt = similar(V)
    _V .= V
    N = size(V,1)
    signflip = zeros(T, N)
    θ = zeros(T,size(V)...) 
    R = RotationMatrix(N)
    for i in 1:N
        @inbounds if _V[i,i] >= 0.0
            signflip[i] = 1.0
        else
            signflip[i] = -1.0
        end
        si = signflip[i]
        for j in 1:N 
            @inbounds _V[j,i] *= si
        end

        θ[:,i] = reduce_dimensions!(_Vt, R, i,_V)
        copy!(_V, _Vt)
    end
    return V*diagm(signflip), signflip,θ 
end

function standardize_basis(Q::AbstractMatrix{T}) where T <: Real
    V,signs,θ = orient_vectors(Q)
    return V
end

end # module
