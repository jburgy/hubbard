module Hubbard

using LinearAlgebra: ⋅, norm, SymTridiagonal, Eigen
import LinearAlgebra.eigen
export Model

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
    u, v = first(Iterators.filter(b -> b ⋅ b == n, [u, v] for u=isqrt(n):-1:1 for v=0:u))
    # apply periodic boundary conditions
    fold(x) = foldl((x, y) -> x - y * fld(x ⋅ y, n), eachrow([u v; -v u]), init = x)
    sites = Iterators.filter(x -> fold(x) == x, [x, y] for y=0:(u + v) for x=-v:u)
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
 0
 1
 2
 3
 4
 5

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
    return index
end

"""
    mask(index, n, k)

Compute the binary mask of the combination of `k` items out of `n`
at position `index` in lexicographic order

# Examples
```jldoctest
julia> map(i -> mask(i, 4, 2), 0:5)
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

function hamiltonian(n::Integer, k::Integer, t::Real, U::Real)::Function
    nCk = binomial(n, k)
    hop = hops(n)
    empty = falses(n)
    function multiply!(y, x)
        for (i, w) ∈ enumerate(x)
            if w == 0
                continue
            end

            h, l = divrem(i - 1, nCk)
            u, d = mask(h, n, k), mask(l, n, k)
            y[i] += U * w * count(u .& d)

            tw = t * w
            for p ∈ hop
                a = u .& p
                if a ≠ empty && a ≠ p
                    y[muladd(nCk, index(u .⊻ p), l) + 1] -= tw
                end
                b = d .& p
                if b ≠ empty && b ≠ p
                    y[muladd(nCk, h, index(d .⊻ p)) + 1] -= tw
                end
            end
        end
    end
end

function lanczos(dimension::Integer, vectors::Matrix, multiply!::Function)::Union{SymTridiagonal, Matrix}
    u = zeros(dimension)
    v = zeros(dimension)
    w = zeros(dimension)
    steps, range = size(vectors)
    z = zeros(dimension, range)

    alpha = zeros(steps)
    beta = zeros(steps - 1)

    v[1] = 1
    @. z += vectors[1, :]' * v
    multiply!(w, v)

    alpha[1] = v ⋅ w
    @. w -= alpha[1] * v

    for step = 2:steps
        beta[step - 1] = norm(w)
        @. begin
            u = v
            v = w / beta[step - 1]
            w = 0
            z += vectors[step, :]' * v
        end
        multiply!(w, v)
        alpha[step] = v ⋅ w
        @. w -= alpha[step] * v + beta[step - 1] * u
    end
    return isempty(z) ? SymTridiagonal(alpha, beta) : z
end

struct Model
    # lattice
    n::Integer
    k::Integer
    # interactions
    t::Real
    U::Real
end

function eigen(m::Model, irange::UnitRange=1:1; steps::Integer=100)
    dimension = abs2(binomial(m.n, m.k))
    multiply! = hamiltonian(m.n, m.k, m.t, m.U)

    m = lanczos(dimension, Matrix(undef, steps, 0), multiply!)
    e = eigen(m, irange)
    z = lanczos(dimension, e.vectors, multiply!)
    return Eigen(e.values, z)
end

end