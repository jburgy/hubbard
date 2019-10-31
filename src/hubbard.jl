module Hubbard

using LinearMaps: FunctionMap

"""
    hops(n)

Computes neighbor bit masks for the `n` site square lattice with
periodic boundary conditions.

# Examples
```jldoctest
julia> map(h -> digits(h, base = 2, pad = 8), hops(8))
16-element Array{Array{Int64,1},1}:
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
function hops(n::Integer)::Vector
    a = first([u, v] for u=isqrt(n):-1:1 for v=0:u if sum(abs2, [u, v]) == n)
    b = [0 -1; 1 0]a  # orthogonal basis vector
    # apply periodic boundary conditions
    fold(x) = isequal(reduce((x, y) -> x - y * fld(x'y, n), [a, b], init = x))
    inside(x) = fold(x)(x)
    sites = [[x, y] for y=0:sum(a) for x=b[1]:b[2] if inside([x, y])]
    masks(y) = [(1<<(i-1))|(1<<(findfirst(fold(x + y), sites)-1)) for (i, x) in enumerate(sites)]
    return [masks([1, 0]); masks([0, 1])]
end

"""
    index(mask)

Compute the index of the combination of `count_ones(mask)` items
in lexicographic order

# Examples
```jldoctest
julia> c = [binomial(i, j) for i=0:3, j=1:2];

julia> [index(parse(Int, b, base = 2), c) for b ∈ ["0011", "0101", "0110", "1001", "1010", "1100"]]
6-element Array{Int64,1}:
 1
 2
 3
 4
 5
 6

```
"""
function index(m::Integer, c::Matrix)::Integer
    i = 1
    @inbounds for k ∈ 1:size(c, 2)
        i += c[trailing_zeros(m) + 1, k]
        m &= m - 1
    end
    return i
end

"""
    mask(index, n, k)

Compute the binary mask of the combination of `k` items out of `n`
at position `index` in lexicographic order

# Examples
```jldoctest
julia> c = [binomial(i, j) for i=0:3, j=1:2];

julia> [digits(mask(i, c), base = 2, pad = 4) for i=1:6]
6-element Array{Array{Int64,1},1}:
 [1, 1, 0, 0]
 [1, 0, 1, 0]
 [0, 1, 1, 0]
 [1, 0, 0, 1]
 [0, 1, 0, 1]
 [0, 0, 1, 1]
 
```
"""
function mask(i::Integer, c::Matrix)::Integer
    i -= 1
    n, k = size(c)
    m = 0
    @inbounds while k > 0
        j = i - c[n, k]
        if j >= 0
            i = j
            k -= 1
            m |= 1
        end
        n -= 1
        m <<= 1
    end
    return m << (n - 1)
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
    c = [binomial(i, j) for i=0:(n-1), j=1:k]
    nCk = binomial(n, k)
    hop = hops(n)
    function mult!(y, x)
        fill!(y, 0.0)
        xm = reshape(x, nCk, :)
        ym = reshape(y, nCk, :)
        a = falses(n)
        @inbounds for i ∈ CartesianIndices(xm)
            w = xm[i]
            w == 0 && continue
 
            l, h = i.I
            u, d = mask(h, c), mask(l, c)
            ym[i] += U * w * count_ones(u & d)
 
            tw = t * w
            for p ∈ hop
                a = u ⊻ p
                if a ≠ u && count_ones(a) == k
                    ym[l, index(a, c)] -= tw
                end
                a = d ⊻ p
                if a ≠ d && count_ones(a) == k
                    ym[index(a, c), h] -= tw
                end
            end
        end
    end
    return FunctionMap{T}(mult!, abs2(nCk), ismutating=true, issymmetric=true)
end

end