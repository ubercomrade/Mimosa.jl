# Execution policies for top-level parallelism.
#
# Parallelism in Mimosa is applied at the *top independent level* — sequences,
# target models, or comparison pairs — never inside inner scanning kernels.
# The kernels themselves are serial and composable. This keeps results
# deterministic regardless of thread count.
#
# Design per ADR 0004 (parallelism-and-rng):
#   - `SerialExecution` processes items in order.
#   - `ThreadedExecution(ntasks)` uses at most `ntasks` worker tasks and never
#     exceeds the number of Julia threads available to the process.
#   - Nested calls execute serially inside an existing parallel worker.

"""
    ExecutionPolicy

Abstract supertype for parallel execution policies.

Concrete policies:
- [`SerialExecution`](@ref): process items sequentially (default).
- [`ThreadedExecution`](@ref): process items using multiple Julia threads.
"""
abstract type ExecutionPolicy end

"""
    SerialExecution

Execute items in sequential order. This is the default policy and guarantees
deterministic, single-threaded execution.
"""
struct SerialExecution <: ExecutionPolicy end

"""
    ThreadedExecution

Execute items using multiple Julia tasks running on separate threads.

Fields:
- `ntasks::Int`: maximum number of concurrent tasks. Defaults to
  `Threads.nthreads()` when constructed with no arguments.

Results are written into pre-allocated slots indexed by original position,
so the output order and values are identical to `SerialExecution`.
"""
struct ThreadedExecution <: ExecutionPolicy
    ntasks::Int

    function ThreadedExecution(ntasks::Integer)
        ntasks < 1 && throw(ArgumentError("ntasks must be ≥ 1, got $ntasks."))
        return new(Int(ntasks))
    end
end

ThreadedExecution() = ThreadedExecution(max(1, Threads.nthreads()))

function Base.show(io::IO, ::SerialExecution)
    return print(io, "SerialExecution()")
end

function Base.show(io::IO, pol::ThreadedExecution)
    return print(io, "ThreadedExecution(ntasks=$(pol.ntasks))")
end

# ── Parallel map helper ───────────────────────────────────────────────────
#
# `_parallel_for(policy, n, f!)` — iterate `f!(i)` for `i in 1:n`.
# Under `SerialExecution` this is a simple loop. Under `ThreadedExecution`
# it uses a bounded dynamic queue with at most `ntasks` workers.
#
# The caller is responsible for pre-allocating result slots and ensuring
# `f!` is thread-safe (no shared mutable state, no `push!` to shared vectors).

"""
    _parallel_for(f!, policy::ExecutionPolicy, n::Int)

Execute `f!(i)` for each `i in 1:n` according to `policy`.

Under `SerialExecution`, this is a simple `for` loop. Under
`ThreadedExecution`, indices are claimed from a bounded dynamic queue. The
function `f!` must be thread-safe: no shared mutable
state, no `push!` to shared vectors, results written only to pre-allocated
slots indexed by `i`.
"""
function _parallel_for end

const _PARALLEL_DEPTH_KEY = gensym(:mimosa_parallel_depth)

function _effective_ntasks(pol::ThreadedExecution, n::Int)
    return min(pol.ntasks, n, max(1, Threads.nthreads()))
end

_is_effectively_threaded(::SerialExecution) = false
_is_effectively_threaded(pol::ThreadedExecution) = min(pol.ntasks, Threads.nthreads()) > 1

function _validate_execution_levels(
    outer_execution::ExecutionPolicy, scan_execution::ExecutionPolicy
)
    _is_effectively_threaded(outer_execution) && _is_effectively_threaded(scan_execution) &&
        throw(ArgumentError("outer_execution and scan_execution cannot both use multiple threads."))
    return nothing
end

function _in_parallel_region()
    return get(task_local_storage(), _PARALLEL_DEPTH_KEY, 0) > 0
end

function _with_parallel_region(f)
    tls = task_local_storage()
    previous = get(tls, _PARALLEL_DEPTH_KEY, 0)
    tls[_PARALLEL_DEPTH_KEY] = previous + 1
    try
        return f()
    finally
        if previous == 0
            delete!(tls, _PARALLEL_DEPTH_KEY)
        else
            tls[_PARALLEL_DEPTH_KEY] = previous
        end
    end
end

function _parallel_for(f!, ::SerialExecution, n::Int)
    @inbounds for i in 1:n
        f!(i)
    end
    return nothing
end

function _parallel_for(f!, pol::ThreadedExecution, n::Int)
    n <= 0 && return nothing
    _in_parallel_region() && return _parallel_for(f!, SerialExecution(), n)

    ntasks = _effective_ntasks(pol, n)
    ntasks <= 1 && return _parallel_for(f!, SerialExecution(), n)

    # A bounded queue avoids stranding a worker behind a long ragged item.
    next_index = Threads.Atomic{Int}(1)
    @sync for t in 1:ntasks
        Threads.@spawn begin
            _with_parallel_region() do
                while true
                    i = Threads.atomic_add!(next_index, 1)
                    i > n && break
                    f!(i)
                end
            end
        end
    end

    return nothing
end

"""
    _parallel_for_weighted(f!, policy, costs)

Execute indices using contiguous, approximately equal-cost blocks. Blocks are
smaller than a worker's full share so ragged heavy items can still be balanced
dynamically, while each atomic queue operation claims several adjacent items.
"""
function _parallel_for_weighted(f!, ::SerialExecution, costs::AbstractVector{<:Integer})
    return _parallel_for(f!, SerialExecution(), length(costs))
end

function _parallel_for_weighted(
    f!, pol::ThreadedExecution, costs::AbstractVector{<:Integer}
)
    n = length(costs)
    n == 0 && return nothing
    _in_parallel_region() && return _parallel_for(f!, SerialExecution(), n)

    ntasks = _effective_ntasks(pol, n)
    ntasks <= 1 && return _parallel_for(f!, SerialExecution(), n)

    total_cost = 0
    @inbounds for cost in costs
        total_cost += max(Int(cost), 1)
    end
    target_cost = max(1, cld(total_cost, ntasks * 8))

    ranges = UnitRange{Int}[]
    start = 1
    accumulated = 0
    @inbounds for i in 1:n
        accumulated += max(Int(costs[i]), 1)
        if accumulated >= target_cost
            push!(ranges, start:i)
            start = i + 1
            accumulated = 0
        end
    end
    start <= n && push!(ranges, start:n)

    return _parallel_for(pol, length(ranges)) do block_index
        @inbounds for i in ranges[block_index]
            f!(i)
        end
    end
end
