using Printf

"""
    ProgressBar(io=stderr; width=28, refresh_seconds=0.1)

A throttled terminal progress renderer for one-to-many comparisons and
`build_null`. Pass an instance as `on_progress=ProgressBar()`. Output is sent
to `stderr` by default, keeping machine-readable results on `stdout`.
"""
mutable struct ProgressBar{I<:IO}
    io::I
    width::Int
    refresh_ns::UInt64
    stage::Symbol
    started_ns::UInt64
    updated_ns::UInt64
    active::Bool
end

function ProgressBar(io::IO=stderr; width::Int=28, refresh_seconds::Real=0.1)
    width > 0 || throw(ArgumentError("progress width must be positive."))
    isfinite(refresh_seconds) && refresh_seconds >= 0 ||
        throw(ArgumentError("progress refresh interval must be non-negative."))
    refresh_ns = UInt64(round(refresh_seconds * 1.0e9))
    return ProgressBar(io, width, refresh_ns, :none, UInt64(0), UInt64(0), false)
end

function _progress_stage_label(stage::Symbol)
    stage === :compare && return "Comparing targets"
    stage === :prepare && return "Preparing models"
    stage === :null && return "Building null"
    return String(stage)
end

function _format_progress_duration(seconds::Real)
    isfinite(seconds) || return "--:--"
    total_seconds = max(0, round(Int, seconds))
    hours, remainder = divrem(total_seconds, 3600)
    minutes, seconds_part = divrem(remainder, 60)
    hours > 0 && return @sprintf("%d:%02d:%02d", hours, minutes, seconds_part)
    return @sprintf("%02d:%02d", minutes, seconds_part)
end

function (progress::ProgressBar)(event)
    now = time_ns()
    stage_changed = event.stage !== progress.stage
    complete = event.total == 0 || event.current >= event.total
    if !stage_changed && !complete && now - progress.updated_ns < progress.refresh_ns
        return nothing
    end

    if stage_changed
        progress.active && println(progress.io)
        progress.stage = event.stage
        progress.started_ns = now
    end

    total = max(event.total, 0)
    current = clamp(event.current, 0, total)
    fraction = total == 0 ? 1.0 : current / total
    filled = clamp(floor(Int, progress.width * fraction), 0, progress.width)
    elapsed = Float64(now - progress.started_ns) / 1.0e9
    eta = current > 0 && current < total ? elapsed * (total - current) / current : 0.0
    timing = if complete
        "elapsed $(_format_progress_duration(elapsed))"
    elseif current == 0
        "ETA --:--"
    else
        "ETA $(_format_progress_duration(eta))"
    end
    percent = round(Int, 100 * fraction)
    label = rpad(_progress_stage_label(event.stage), 17)
    bar = repeat("█", filled) * repeat("░", progress.width - filled)

    print(
        progress.io,
        '\r',
        label,
        " [",
        bar,
        "] ",
        current,
        '/',
        total,
        ' ',
        lpad(percent, 3),
        "%  ",
        timing,
    )
    flush(progress.io)
    progress.updated_ns = now
    progress.active = !complete
    complete && println(progress.io)
    return nothing
end

function _finish_progress!(::Nothing)
    return nothing
end

function _finish_progress!(progress::ProgressBar)
    progress.active && println(progress.io)
    progress.active = false
    return nothing
end
