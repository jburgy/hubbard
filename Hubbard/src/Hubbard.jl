module Hubbard

using LinearMaps: FunctionMap
using StaticArrays
export TiltedSquare, hamiltonian, hops, index, mask, symmetries, ρ

const ρ = SA[0 -1; 1 0]  # Rotation by `½π`
const σ = SA[-1 0; 0 1]  # Reflection about `x` axis

"""
    restrict(x, tilt)

Apply https://en.wikipedia.org/wiki/Periodic_boundary_conditions to restrict point `x`
within the `tilt × ρ * tilt` titled square
"""
function restrict(x::SVector{2,Int}, tilt::SVector{2,Int})::SVector{2,Int}
    reduce((u, v) -> u - v * fld(u'v, sum(abs2, tilt)), [tilt, ρ * tilt], init = x)
end

restrict(tilt::SVector{2,Int}) = Base.Fix2(restrict, tilt)

"""
An integer greater than one can be written as a sum of two squares if and only if
its prime decomposition contains no factor pᵏ, where prime p ≡ 3 (mod 4) and k is odd.
"""
function istwosquares(n::Int)::Bool
    p = 2
    k = 0
    i = false
    while p ≤ n
        div, rem = divrem(n, p)
        if rem ≠ 0
            isodd(k) && return false
            p += 1
            k = 0
            i = p % 4 == 3
        else
            n = div
            k += i
        end
    end
    return iseven(k)
end

"""
    TiltedSquare

represents a compact subset of the Euclidean square lattice with `n` sites.
"""
struct TiltedSquare{N}
    tilt::SVector{2,Int}
    sites::SVector{N,SVector{2,Int}}

    function TiltedSquare{N}() where N
        istwosquares(N) || throw(DomainError("$n is not a sum of two squares"))
        pred = ==(N) ∘ Base.Fix1(sum, abs2)
        tilt = first(t for t ∈ (SVector(isqrt(N - v^2), v) for v ∈ 0:isqrt(N÷2)) if pred(t))
        pred = isequal ∘ Base.Fix2(restrict, tilt)
        sites = [site for site ∈ (SVector(x, y) for y ∈ 0:sum(tilt) for x ∈ -tilt[2]:tilt[1]) if pred(site)(site)]
        return new(tilt, SVector{N,SVector{2,Int}}(sites))
    end
end

Base.isequal(s::TiltedSquare) = isequal ∘ restrict(s.tilt)
Base.iterate(s::TiltedSquare, args...) = Base.iterate(s.sites, args...)
Base.length(s::TiltedSquare) = Base.length(s.sites)

"""
    Base.findall(r::Matrix{Int}, s::TiltedSquare)

Convert application of rotation matrix `r` to all sites in tilted square `s` to a
permutation `1:length(s.sites)` assuming periodic boundary conditions.

See the application of a rotation by `½π` to the 8-site tilted square for example.
The sites of `TiltedSquare{8}()` are labeled as follows:

      8
    5 6 7 
    2 3 4
      1

Rotating these sites by `½π` results in

      7 4
    8 6 3 1
      5 2

Translating indices which landed outside the original shape (i.e. all but `1` and `4`)
by `a = [2, 2]` gives

      7
    8 6 3
    4 5 2
      1

Interpreting this new arrangement as a permutation, site `1` stayed in place, `2` moved to
position `4`, `3` moved to position `7`, ... → `(1 4 7 2 3 6 8 5)`

```jldoctest
julia> findall(ρ, TiltedSquare{8}())
8-element Vector{Int64}:
 1
 4
 7
 2
 3
 6
 8
 5

```
"""
function Base.findall(r::SMatrix{2,2,Int,4}, s::TiltedSquare)::Vector{Int}
    [findfirst(isequal(s)(r * site), s.sites) for site ∈ s]
end

"""
    findall(d::Vector{T}, s::TiltedSquare)

Convert application of displacement `d` to all sites in tilted square `s` to a
permutation of `1:length(s.sites)` assuming periodic boundary conditions.

      8→2'
    5→6→7→1'
    2→3→4→8'
      1→5'

```jldoctest
julia> findall([1, 0], TiltedSquare{8}())
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
function Base.findall(d::SVector{2,Int}, s::TiltedSquare)::Vector{Int}
    [findfirst(isequal(s)(site + d), s.sites) for site ∈ s]
end

"""
    hops(tilted_square)

Computes neighbor bit masks for the `n` site square lattice with
periodic boundary conditions.

# Examples
```jldoctest
julia> [digits(h, base = 2, pad = 8) for h ∈ hops(TiltedSquare{8}())]
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
function hops(s::TiltedSquare)::Vector{Int}
    masks(y) = [(1<<(i-1))|(1<<(j-1)) for (i, j) in enumerate(findall(y, s))]
    return [masks(SA[1, 0]); masks(SA[0, 1])]
end


"""
    symmetries(s)

Elements of dihedral group D₄ specialized to the `n` site tilted square (minus
translations by `tilt` and `ρ * tilt` as those are already taken care of by `restrict`)
as binary masks.  The group is generated by 2 elements:

1. Rotation by `½π`: `ρ`
1. Reflection about `x` axis: ``

Remaining elements can be represented as products of these generators:

1. Inversion: `ρ²`
1. Rotation by `-½π`: `ρ³ ≡ ρ⁻¹`
1. Reflection about anti-diagonal: `ρ σ`
1. Reflection about `y` axis: `ρ² σ ≡ ρ σ ρ⁻¹ ≡ τ`
1. Reflection about diagonal: `ρ³ σ ≡ ρ τ`

We already know that rotation by `½π` is equivalent to the `(1 4 7 2 3 6 8 5)` permutation.
Because fermions are indistinghishable, we convert this permutation to a bit mask where 1
indicates the site moves.  `(1 4 7 2 3 6 8 5)` becomes `0b011111011`

Reflection about the `x` axis does

      8
    7 6 5
    4 3 2
      1

or `(1 4 3 2 7 6 5 8)` as a permutation and `0b01011010` as a bitmask.

Let us confirm D₄ by working out one repeated application: `ρ σ`

                 5
      5 2      8 6 3
    8 6 3 1    2 7 4
      7 4        1

or `(1 2 7 4 8 6 3 5)` as a permutation.  We can easily verify that this equals the product of
`(1 4 7 2 3 6 8 5)(1 4 3 2 7 6 5 8)` by following the orbit of each site:

    1 → 1 → 1
    2 → 4 → 2
    3 → 3 → 7
    4 → 2 → 4
    5 → 7 → 8
    6 → 6 → 6
    7 → 5 → 3
    8 → 8 → 5

```jldoctest
julia> [digits(s, base = 2, pad = 8) for s ∈ symmetries(TiltedSquare{8}())]
7-element Vector{Vector{Int64}}:
 [0, 0, 0, 0, 0, 0, 0, 0]
 [0, 1, 1, 1, 1, 0, 1, 1]
 [0, 0, 1, 0, 1, 0, 1, 1]
 [0, 1, 1, 1, 1, 0, 1, 1]
 [0, 0, 1, 0, 1, 0, 1, 1]
 [0, 1, 1, 1, 0, 0, 0, 1]
 [0, 0, 1, 0, 1, 0, 1, 1]

```
"""
function symmetries(s::TiltedSquare)::Vector{Int}
    mask(r) = reduce((a, b) -> a << 1 | b, (a ≠ b for (a, b) ∈ Iterators.reverse(enumerate(findall(r, s)))))
    return [mask(r) for r ∈ [SA[1 0; 0 1], ρ, ρ^2, ρ^3, ρ * σ, ρ^2 * σ, ρ^3 * σ]]
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

index(n::Tuple{Int, Int}, c::Matrix) = CartesianIndex(Tuple(index(m, c) for m ∈ n))

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
3414×3414 LinearMaps.FunctionMap{Float64}(mult!; ismutating=true, issymmetric=true, ishermitian=true, isposdef=false)

julia> using Arpack: eigs

julia> λ, ϕ = eigs(h, nev=1, which=:SR);

julia> λ
1-element Vector{Float64}:
 -11.859367228698792

julia> size(ϕ)
(3414, 1)

```
"""
function hamiltonian(n::Integer, k::Integer, t::T, U::T) where T <: Real
    c = [binomial(i, j) for i=0:(n-1), j=1:k]
    nCk = binomial(n, k)
    s = TiltedSquare{n}()
    hop = hops(s)
    ok = isequal(k) ∘ count_ones

    # Apply all symmetries to all combinations
    # TODO: separate nup and ndn (now k = nup = ndn)
    masks = Base.Fix2(mask, c).(1:nCk) .⊻ unique(symmetries(s))'

    minpairs = [
        minimum(t for t ∈ zip(m, n) if all(ok, t))
        for m ∈ eachrow(masks), n ∈ eachrow(masks)
    ]

    uniquepairs = sort(vec(minpairs))
    unique!(uniquepairs)

    minindices = Base.Fix1(searchsortedfirst, uniquepairs).(minpairs)
    function mult!(y, x)
        fill!(y, 0.0)
        @inbounds for ((i, w), (u, d)) ∈ zip(enumerate(x), uniquepairs)
            w == 0 && continue

            y[i] += U * w * count_ones(u & d)

            iu = index(u, c)
            id = index(d, c)

            tw = t * w
            for p ∈ hop
                a = u ⊻ p
                if a ≠ u && ok(a)
                    y[minindices[index(a, c), id]] -= tw
                end
                a = d ⊻ p
                if a ≠ d && ok(a)
                    y[minindices[iu, index(a, c)]] -= tw
                end
            end
        end
    end
    return FunctionMap{T}(mult!, length(uniquepairs), ismutating=true, issymmetric=true)
end

end