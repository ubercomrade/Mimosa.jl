# Adapters from raw PFM files to ready-to-scan PWM models.

function _read_meme_pwm(
    path::AbstractString; index::Integer=0, background::AbstractFloat=0.25f0
)
    pfm = read_meme(path; index=index)
    return pwm_from_pfm(pfm.frequencies; background=background, name=pfm.name)
end

function _read_pfm_pwm(path::AbstractString; background::AbstractFloat=0.25f0)
    pfm = read_pfm(path)
    return pwm_from_pfm(pfm.frequencies; background=background, name=pfm.name)
end
