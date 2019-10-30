module Hubbard

using LinearAlgebra: ⋅, norm
using LinearMaps: FunctionMap

squares = (0:8) .^ 2
sizes = filter(n -> any(w -> (n - w) in squares, squares), 1:64)

"""
    hops(n)

Computes neighbor bit masks for the `n` site square lattice with
periodic boundary conditions.

# Examples
```jldoctest
julia> hops(8)
16-element Array{BitArray{1},1}:
 [1, 0, 0, 0, 1, 0, 0, 0]
 [0, 1, 1, 0, 0, 0, 0, 0]
 [0, 0, 1, 1, 0, 0, 0, 0]
 [0, 0, 0, 1, 0, 0, 0, 1]
 [0, 0, 0, 0, 1, 1, 0, 0]
 [0, 0, 0, 0, 0, 1, 1, 0]
 [1, 0, 0, 0, 0, 0, 1, 0]
 [0, 1, 0, 0, 0, 0, 0, 1]
 [1, 0, 1, 0, 0, 0, 0, 0]
 [0, 1, 0, 0, 1, 0, 0, 0]
 [0, 0, 1, 0, 0, 1, 0, 0]
 [0, 0, 0, 1, 0, 0, 1, 0]
 [0, 0, 0, 1, 1, 0, 0, 0]
 [0, 0, 0, 0, 0, 1, 0, 1]
 [0, 1, 0, 0, 0, 0, 1, 0]
 [1, 0, 0, 0, 0, 0, 0, 1]

```
"""
function hops(n::Integer)::Array{BitVector}
    u, v = first([u, v] for u=isqrt(n):-1:1 for v=0:u if sum(abs2, [u, v]) == n)
    # apply periodic boundary conditions
    fold(x) = foldl((x, y) -> x - y * fld(x ⋅ y, n), eachrow([u v; -v u]), init = x)
    sites = [[x, y] for y=0:(u + v) for x=-v:u if fold([x, y]) == [x, y]]
    masks(y) = map(x -> BitArray(z ∈ (x, fold(x + y)) for z ∈ sites), sites)
    return [masks([1, 0]); masks([0, 1])]
end

"""
    index(mask)

Compute the index of the combination of `count_ones(mask)` items
in lexicographic order

# Examples
```jldoctest
julia> map(index ∘ BitVector, eachrow(BitArray([
       1 1 0 0; 1 0 1 0; 0 1 1 0; 1 0 0 1; 0 1 0 1; 0 0 1 1
       ])))
6-element Array{Int64,1}:
 1
 2
 3
 4
 5
 6

```
"""
function index(mask::BitVector)::Integer
    n = k = nmkp1 = nCkm1 = 1
    nCk = index = 0
    for bit ∈ mask
        if bit
            index += nCk
            nCkm1 *= n
            nCkm1 ÷= k
            nCk *= n
            k += 1
            nCk ÷= k
        else
            nCk += nCkm1
            nCkm1 *= n
            nCkm1 ÷= nmkp1
            nmkp1 += 1
        end
        n += 1
    end
    return index + 1
end

"""
    mask(index, n, k)

Compute the binary mask of the combination of `k` items out of `n`
at position `index` in lexicographic order

# Examples
```jldoctest
julia> map(i -> mask(i, 4, 2), 1:6)
6-element Array{BitArray{1},1}:
 [1, 1, 0, 0]
 [1, 0, 1, 0]
 [0, 1, 1, 0]
 [1, 0, 0, 1]
 [0, 1, 0, 1]
 [0, 0, 1, 1]
 
```
"""
function mask(index::Integer, n::Integer, k::Integer)::BitVector
    index -= 1
    nCk = binomial(n, k)
    nmk = n - k
    mask = falses(n)
    while true
        if index >= nCk
            mask[1] = true
            if k == 1
                break
            end
            index -= nCk
            nCk *= k
            nCk ÷= n
            k -= 1
        else
            nCk *= nmk
            nCk ÷= n
            nmk -= 1
        end
        mask >>= 1
        n -= 1
    end
    return mask >> n
end

"""
    hamiltonian(n, k, t, U)

Implements the [Hubbard hamiltonian](https://en.wikipedia.org/wiki/Hubbard_model)
as a [mutating linear map](https://github.com/Jutho/LinearMaps.jl).

# Arguments
- `n::Integer`: sites in the square lattice with periodic boundary conditions
- `k::Integer`: spin up and spin down fermions
- `t::Real`: hopping integral
- `U::Real`: strength of on-site interaction 

# Examples
```jldoctest
julia> h = hamiltonian(8, 4, 1., 4.)
LinearMaps.FunctionMap{Float64}(multiply!, 4900, 4900; ismutating=true, issymmetric=true, ishermitian=true, isposdef=false)

julia> using Arpack: eigs

julia> λ, ϕ = eigs(h, nev=1, which=:SR)

julia> λ
1-element Array{Float64,1}:
 -11.775702792136244

 julia> size(ϕ)
 (4900, 1)

```
"""
function hamiltonian(n::Integer, k::Integer, t::T, U::T) where T <: Real
    nCk = binomial(n, k)
    hop = hops(n)
    function multiply!(y, x)
        fill!(y, 0.0)
        xm = reshape(x, nCk, :)
        ym = reshape(y, nCk, :)
        a = falses(n)
        @inbounds for i ∈ CartesianIndices(xm)
            w = xm[i]
            if w == 0
                continue
            end

            l, h = Tuple(i)
            u, d = mask(h, n, k), mask(l, n, k)
            map!(&, a, u, d)
            ym[i] += U * w * count(a)

            tw = t * w
            for p ∈ hop
                map!(&, a, u, p)
                if any(a) && a ≠ p
                    map!(⊻, a, u, p)
                    ym[l, index(a)] -= tw
                end
                map!(&, a, d, p)
                if any(a) && a ≠ p
                    map!(⊻, a, d, p)
                    ym[index(a), h] -= tw
                end
            end
        end
    end
    dimension = abs2(binomial(n, k))
    return FunctionMap{T}(multiply!, dimension, ismutating=true, issymmetric=true)
end

end