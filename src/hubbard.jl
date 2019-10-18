module Hubbard

using Base: binomial, divrem, muladd, zeros
using LinearAlgebra: norm, SymTridiagonal, ⋅
import LinearAlgebra.eigvals
export Model, eigvals

squares = (0:8) .^ 2
sizes = filter(n -> any(w -> (n - w) in squares, squares), 1:64)

"""
    hops(n)

Computes bit masks neighbors for the square lattice with `n` sites

# Examples
```jldoctest
julia> hops(8)
16-element Array{BitArray{1},1}:
 [1, 0, 0, 0, 1, 0, 0, 0]
 [1, 0, 1, 0, 0, 0, 0, 0]
 [0, 1, 1, 0, 0, 0, 0, 0]
 [0, 1, 0, 0, 1, 0, 0, 0]
 [0, 0, 1, 1, 0, 0, 0, 0]
 [0, 0, 1, 0, 0, 1, 0, 0]
 [0, 0, 0, 1, 0, 0, 0, 1]
 [0, 0, 0, 1, 0, 0, 1, 0]
 [0, 0, 0, 0, 1, 1, 0, 0]
 [0, 0, 0, 1, 1, 0, 0, 0]
 [0, 0, 0, 0, 0, 1, 1, 0]
 [0, 0, 0, 0, 0, 1, 0, 1]
 [1, 0, 0, 0, 0, 0, 1, 0]
 [0, 1, 0, 0, 0, 0, 1, 0]
 [0, 1, 0, 0, 0, 0, 0, 1]
 [1, 0, 0, 0, 0, 0, 0, 1]

```
"""
function hops(n)
    spans((u, v)) = abs2(u) + abs2(v) == n
    u, v = first(filter(spans, collect((u, v) for u=1:isqrt(n) for v=0:u)))
    function fold((x, y))
        z = u * x + v * y
        while 0 > z
            x += u
            y += v
            z += n
        end
        while z ≥ n
            x -= u
            y -= v
            z -= n
        end
        z = u * y - v * x
        while 0 > z
            x -= v
            y += u
            z += n
        end
        while z ≥ n
            x += v
            y -= u
            z -= n
        end
        x, y
    end
    inside(s) = fold(s) == s
    sites = filter(inside, collect((x, y) for y=0:(u + v) for x=-v:u))
    masks = map(sites) do s
        a = fold(s .+ (1, 0))
        b = fold(s .+ (0, 1))
        a = map(r -> r ∈ (s, a), sites)
        b = map(r -> r ∈ (s, b), sites)
        BitVector(a), BitVector(b)
    end
    collect(Iterators.flatten(masks))
end

"""
    index(mask)

Compute the index of the combination of `count_ones(mask)` items
in lexicographic order

# Examples
```jldoctest
julia> map(index, eachrow(BitArray([
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
function index(mask)
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
    index
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
function mask(index, n, k)
    nCk = binomial(n, k)
    nmk = n - k
    mask = BitVector(undef, n)
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
    mask >> n
end

function hamiltonian(n, k, t, U)
    nCk = binomial(n, k)
    hop = hops(n)
    function multiply(x, y)
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
                if a ≠ falses(n) && a ≠ p
                    y[muladd(nCk, index(u .⊻ p), l) + 1] -= tw
                end
                b = d .& p
                if b ≠ falses(n) && b ≠ p
                    y[muladd(nCk, h, index(d .⊻ p)) + 1] -= tw
                end
            end
        end
    end
end

function lanczos(dimension, steps, multiply)
    u = zeros(dimension)
    v = zeros(dimension)
    w = zeros(dimension)

    alpha = zeros(steps)
    beta = zeros(steps - 1)

    v[1] = 1
    multiply(v, w)

    alpha[1] = v ⋅ w
    w .-= alpha[1] .* v

    for step = 2:steps
        beta[step - 1] = norm(w)
        u .= v
        v .= w / beta[step - 1]
        w .= 0
        multiply(v, w)
        alpha[step] = v ⋅ w
        w .-= alpha[step] .* v - beta[step - 1] .* u
    end
    SymTridiagonal(alpha, beta)
end

struct Model
    # lattice
    n::Int
    k::Int
    # interactions
    t::Real
    U::Real
end

function eigvals(m::Model, steps::Int)
    dimension = binomial(m.n, m.k)^2
    multiply = hamiltonian(m.n, m.k, m.t, m.U)
    matrix = lanczos(dimension, steps, multiply)
    eigvals(matrix)
end

end