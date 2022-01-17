module Hubbard

import Base: findall, isequal, length
using LinearMaps: FunctionMap
export TiltedSquare, hamiltonian, hops, index, mask

"""
    TiltedSquare{T}

represents a compact subset of the Euclidean square lattice with length `N`.
"""
struct TiltedSquare{T<:Integer}
    sites::Vector{Vector{T}}
    a::Vector{T}
    b::Vector{T}

    function TiltedSquare{T}(n::T) where T<:Integer
        a = first([u, v] for u=isqrt(n):-1:1 for v=0:u if sum(abs2, [u, v]) == n)
        b = [0 -1; 1 0]a
        s0 = new(Vector{Vector{T}}(undef, n), a, b)
        return new([[x, y] for y=0:sum(a) for x=b[1]:b[2] if [x, y] ∈ s0], a, b)
    end    
end

"""
    isequal(x, s)

Implement https://en.wikipedia.org/wiki/Periodic_boundary_conditions in tilted
square lattice `s`

# Examples
```jldoctest
julia> s = TiltedSquare{Int64}(8)
TiltedSquare{Int64}([[0, 0], [-1, 1], [0, 1], [1, 1], [-1, 2], [0, 2], [1, 2], [0, 3]], [2, 2], [-2, 2])

julia> isequal(s.a, s)([0, 0])
true

````
"""
isequal(x::Vector{T}, s::TiltedSquare{T}) where T<:Integer = isequal(
    reduce((u, v) -> u - v * fld(u'v, length(s.sites)), [s.a, s.b], init = x)
)

∈(x::Vector{T}, s::TiltedSquare{T}) where T<:Integer = isequal(x, s)(x)

"""
    findall(r::Matrix{T}, s::TiltedSquare{T})

Convert application of rotation matrix `r` to all sites in tilted square `s` to a
permutation `1:length(s.sites)` assuming periodic boundary conditions.

See the application of a rotation by `π/4` to the 8-site tilted square for example.
The sites of `TiltedSquare{Int64}(8)` are labeled as follows:

      8
    5 6 7 
    2 3 4
      1

Rotating these sites by `π/4` results in

      7 4
    8 6 3 1
      5 2

Translating indices which landed outside the original shape (i.e. all but `1` and `4`)
by `a = [2, 2]` gives

      7
    8 6 3
    4 5 2
      1

```jldoctest
julia> findall([0 1; -1 0], TiltedSquare{Int64}(8))
8-element Vector{Int64}:
 1
 4
 5
 2
 8
 6
 3
 7

```
"""
findall(r::Matrix{T}, s::TiltedSquare{T}) where T<:Integer = [findfirst(isequal(r * site, s), s.sites) for site ∈ s.sites]

"""
    findall(d::Vector{T}, s::TiltedSquare{T})

Convert application of displacement `d` to all sites in tilted square `s` to a
permutation of `1:length(s.sites)` assuming periodic boundary conditions.

    8→2'
  5→6→7→1'
  2→3→4→8'
    1→5'

```jldoctest
julia> findall([1, 0], TiltedSquare{Int64}(8))
8-element Vector{Int64}:
 5
 3
 4
 8
 6
 7
 1
 2

"""
findall(d::Vector{T}, s::TiltedSquare{T}) where T<:Integer = [findfirst(isequal(site + d, s), s.sites) for site ∈ s.sites]


"""
    hops(tilted_square)

Computes neighbor bit masks for the `n` site square lattice with
periodic boundary conditions.

# Examples
```jldoctest
julia> [digits(h, base = 2, pad = 8) for h ∈ hops(8)]
16-element Vector{Vector{Int64}}:
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
function hops(s::TiltedSquare)::Vector
    masks(y) = [(1<<(i-1))|(1<<(j-1)) for (i, j) in enumerate(findall(y, s))]
    return [masks([1, 0]); masks([0, 1])]
end

hops(n::T) where T<:Integer = hops(TiltedSquare{T}(n))


"""
    symmetries()
"""
function symmetries(s::TiltedSquare)::Vector
    mask(r) = reduce((a, b) -> a << 1 | b, (a == b for (a, b) ∈ enumerate(findall(r, s))))
    return [mask(r) for r ∈ [[1 0; 0 -1], [-1 0; 0 1], [0 -1; 1 0], [-1 0; 0 -1]]]
end

"""
    index(m, c)

Compute the index of the combination of `count_ones(mask)` items
in lexicographic order

# Examples
```jldoctest
julia> c = [binomial(i, j) for i=0:3, j=1:2];

julia> [index(parse(Int, b, base = 2), c) for b ∈ ["0011", "0101", "0110", "1001", "1010", "1100"]]
6-element Vector{Int64}:
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
    mask(i, c)

Compute the binary mask of the combination of `k` items out of `n`
at position `index` in lexicographic order

# Examples
```jldoctest
julia> c = [binomial(i, j) for i=0:3, j=1:2];

julia> [digits(mask(i, c), base = 2, pad = 4) for i=1:6]
6-element Vector{Vector{Int64}}:
 [1, 1, 0, 0]
 [1, 0, 1, 0]
 [0, 1, 1, 0]
 [1, 0, 0, 1]
 [0, 1, 0, 1]
 [0, 0, 1, 1]
 
```
"""
function mask(i::Integer, c::Matrix)::Integer
    k = size(c, 2)
    m = 0
    @inbounds for n ∈ size(c, 1):-1:1
        j = i - c[n, k]
        if j ≥ 1
            m |= 1
            k == 1 && return m << (n - 1)
            i = j
            k -= 1
        end
        m <<= 1
    end
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
4900×4900 LinearMaps.FunctionMap{Float64}(mult!; ismutating=true, issymmetric=true, ishermitian=true, isposdef=false)

julia> using Arpack: eigs

julia> λ, ϕ = eigs(h, nev=1, which=:SR);

julia> λ
1-element Vector{Float64}:
 -11.775702792136236

julia> size(ϕ)
(4900, 1)

```
"""
function hamiltonian(n::Integer, k::Integer, t::T, U::T) where T <: Real
    c = [binomial(i, j) for i=0:(n-1), j=1:k]
    nCk = binomial(n, k)
    hop = hops(n)
    ok = isequal(k) ∘ count_ones
    function mult!(y, x)
        fill!(y, 0.0)
        xm = reshape(x, nCk, :)
        ym = reshape(y, nCk, :)
        @inbounds for i ∈ CartesianIndices(xm)
            w = xm[i]
            w == 0 && continue
 
            l, h = i.I
            u, d = mask(h, c), mask(l, c)
            ym[i] += U * w * count_ones(u & d)
 
            tw = t * w
            for p ∈ hop
                a = u ⊻ p
                if a ≠ u && ok(a)
                    ym[l, index(a, c)] -= tw
                end
                a = d ⊻ p
                if a ≠ d && ok(a)
                    ym[index(a, c), h] -= tw
                end
            end
        end
    end
    return FunctionMap{T}(mult!, abs2(nCk), ismutating=true, issymmetric=true)
end

end