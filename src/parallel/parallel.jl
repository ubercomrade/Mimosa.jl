# Execution control for computational kernels.
#
# Targets and comparison pairs remain serial. Parallelism is applied inside
# independent computational kernels such as sequence scanning, normalization,
# anchor collection, and profile alignment.

"""
    Execution(ntasks=1)

Control the maximum number of Julia tasks used by computational kernels.

`Execution(1)` uses the allocation-free sequential path. `Execution(n)` with
`n > 1` schedules at most `n` worker tasks on Julia's thread pool, capped by
`Threads.nthreads()` and the amount of available work.

Results are written into pre-allocated slots indexed by their original
position, preserving deterministic output order.
"""
struct Execution{N}
    ntasks::Int

    function Execution{N}() where {N}
        N isa Int || throw(ArgumentError("ntasks must be an Int, got $N."))
        N < 1 && throw(ArgumentError("ntasks must be ≥ 1, got $N."))
        return new{N}(N)
    end
end

Execution(ntasks::Integer=1) = Execution{Int(ntasks)}()

function Base.show(io::IO, execution::Execution)
    return print(io, "Execution(ntasks=$(execution.ntasks))")
end

function _effective_ntasks(::Execution{N}, n::Int) where {N}
    return min(N, n, Threads.nthreads())
end

@inline function _chunk_bounds(n::Int, chunk::Int, nchunks::Int)
    return (fld((chunk - 1) * n, nchunks) + 1, fld(chunk * n, nchunks))
end

"""
    _parallel_for(f!, execution, n)

Execute `f!(i)` for every `i in 1:n`. Work is claimed dynamically to balance
independent items with unknown or moderately varying cost. The callback must
not mutate shared state except for pre-allocated slots owned by `i`.
"""
function _parallel_for(f!, ::Execution{1}, n::Int)
    for i in 1:n
        f!(i)
    end
    return nothing
end

function _parallel_for(f!, execution::Execution, n::Int)
    n <= 0 && return nothing
    ntasks = _effective_ntasks(execution, n)

    if ntasks <= 1
        for i in 1:n
            f!(i)
        end
        return nothing
    end

    next_index = Threads.Atomic{Int}(1)
    @sync for _ in 1:ntasks
        Threads.@spawn begin
            while true
                i = Threads.atomic_add!(next_index, 1)
                i > n && break
                f!(i)
            end
        end
    end
    return nothing
end

"""
    _parallel_chunks(f!, execution, n)

Divide `1:n` into at most `execution.ntasks` contiguous chunks and call
`f!(first, last, chunk)` once for each chunk. This is intended for dense loops
that need per-worker scratch space or deterministic local reductions.
"""
function _parallel_chunks(f!, ::Execution{1}, n::Int)
    n > 0 && f!(1, n, 1)
    return nothing
end

function _parallel_chunks(f!, execution::Execution, n::Int)
    n <= 0 && return nothing
    nchunks = _effective_ntasks(execution, n)

    if nchunks <= 1
        f!(1, n, 1)
        return nothing
    end

    @sync for chunk in 1:nchunks
        first, last = _chunk_bounds(n, chunk, nchunks)
        Threads.@spawn f!(first, last, chunk)
    end
    return nothing
end

"""
    _parallel_for_weighted(f!, execution, costs)

Execute indices in contiguous, approximately equal-cost blocks. Blocks are
smaller than a worker's full share so ragged heavy items remain dynamically
balanced while each atomic queue operation claims several adjacent items.
"""
function _parallel_for_weighted(
    f!, execution::Execution{1}, costs::AbstractVector{<:Integer}
)
    return _parallel_for(f!, execution, length(costs))
end

function _parallel_for_weighted(f!, execution::Execution, costs::AbstractVector{<:Integer})
    n = length(costs)
    n == 0 && return nothing
    ntasks = _effective_ntasks(execution, n)
    ntasks <= 1 && return _parallel_for(f!, execution, n)

    total_cost = 0
    for cost in costs
        total_cost += max(Int(cost), 1)
    end
    target_cost = max(1, cld(total_cost, ntasks * 8))

    ranges = UnitRange{Int}[]
    sizehint!(ranges, min(n, ntasks * 8))
    start = 1
    accumulated = 0
    for i in 1:n
        accumulated += max(Int(costs[i]), 1)
        if accumulated >= target_cost
            push!(ranges, start:i)
            start = i + 1
            accumulated = 0
        end
    end
    start <= n && push!(ranges, start:n)

    return _parallel_for(execution, length(ranges)) do block_index
        for i in ranges[block_index]
            f!(i)
        end
    end
end
