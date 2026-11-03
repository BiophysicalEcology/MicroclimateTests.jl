# data_corrections.jl — the kinds of correction tables/data_corrections.toml applies to mistakes
# in the OzFlux files, and how each changes the data.

"""
A change to one site's OzFlux data. Each carries the `reason` for it.
"""
abstract type DataCorrection end

"Values of `variable` multiplied by `factor`: units mislabelled in the file"
@kwdef struct Rescale <: DataCorrection
    variable::String
    factor::Float64
    reason::String
end

"Values of `variable` at or below `threshold` multiplied by `factor`: units switch within the record"
@kwdef struct RescaleBelow <: DataCorrection
    variable::String
    threshold::Float64
    factor::Float64
    reason::String
end

"Values of `variable` outside `low`–`high` dropped: physically implausible, though QC passes them"
@kwdef struct ImplausibleRange <: DataCorrection
    variable::String
    low::Float64
    high::Float64
    reason::String
end

"Runs of identical non-zero readings at least `duration_minutes` long dropped, in variables not starting with one of `exempt_prefixes`: a stuck sensor"
@kwdef struct Flatline <: DataCorrection
    duration_minutes::Int
    exempt_prefixes::Vector{String}
    reason::String
end

"`variable` is soil moisture as a fraction of saturation, so is multiplied by the porosity"
@kwdef struct SaturationFraction <: DataCorrection
    variable::String
    reason::String
end

"The file's global attribute `attribute` replaced by `value`"
@kwdef struct AttributeOverride <: DataCorrection
    attribute::String
    value::String
    reason::String
end

"`variable`'s depth range (m), in place of its file's height"
@kwdef struct SensorDepth <: DataCorrection
    variable::String
    top_m::Float64
    bottom_m::Float64
    reason::String
end

"`variable` dropped wherever it is read, forcing or comparison: a faulty sensor"
@kwdef struct FaultySensor <: DataCorrection
    variable::String
    reason::String
end

"The site is not run"
@kwdef struct ExcludedSite <: DataCorrection
    reason::String
end

const CORRECTION_KINDS = Dict(string(nameof(T)) => T for T in
    (Rescale, RescaleBelow, ImplausibleRange, Flatline, SaturationFraction, AttributeOverride, SensorDepth, FaultySensor, ExcludedSite))

# A table entry as its correction
correction(entry) = CORRECTION_KINDS[entry["type"]](; (Symbol(k) => v for (k, v) in entry if k != "type")...)

const DATA_CORRECTIONS_TABLE = read_table("data_corrections")
const ALL_SITE_CORRECTIONS = DataCorrection[correction(entry) for entry in DATA_CORRECTIONS_TABLE["all_sites"]]
const SITE_DATA_CORRECTIONS = Dict{String,Vector{DataCorrection}}(site_name => correction.(entries)
    for (site_name, entries) in DATA_CORRECTIONS_TABLE["sites"])

site_corrections(site_name) = [ALL_SITE_CORRECTIONS; get(SITE_DATA_CORRECTIONS, site_name, DataCorrection[])]

has_correction(corrections, type, variable) = any(c -> c isa type && c.variable == variable, corrections)
excluded(site_name) = any(c -> c isa ExcludedSite, site_corrections(site_name))

# One line per correction of a site's own, with its reason, for the run log
function print_corrections(site_name)
    for c in get(SITE_DATA_CORRECTIONS, site_name, DataCorrection[])
        fields = join((repr(getfield(c, f)) for f in fieldnames(typeof(c)) if f != :reason), ", ")
        println("  Correction $(nameof(typeof(c)))($fields): $(c.reason)")
    end
end

# `vals` of `variable` after every value correction in `corrections`; `time_step` is the
# file's record interval, for Flatline
correct_values(corrections, variable, vals, time_step) =
    foldl((v, c) -> correct(c, variable, v, time_step), corrections; init=vals)

correct(::DataCorrection, variable, vals, time_step) = vals
correct(c::FaultySensor, variable, vals, time_step) = variable == c.variable ? fill(missing, length(vals)) : vals
correct(c::Rescale, variable, vals, time_step) = variable == c.variable ? vals .* c.factor : vals
correct(c::RescaleBelow, variable, vals, time_step) = variable == c.variable ?
    [ismissing(v) || v > c.threshold ? v : v * c.factor for v in vals] : vals
correct(c::ImplausibleRange, variable, vals, time_step) = variable == c.variable ?
    [ismissing(v) || c.low <= v <= c.high ? v : missing for v in vals] : vals
function correct(c::Flatline, variable, vals, time_step)
    any(prefix -> startswith(variable, prefix), c.exempt_prefixes) && return vals
    min_run = Minute(c.duration_minutes) ÷ time_step
    out = Vector{Union{Missing,Float64}}(vals)
    i = 1
    while i <= length(out)
        j = i
        while j < length(out) && !ismissing(out[i]) && !ismissing(out[j+1]) && out[j+1] == out[i]
            j += 1
        end
        !ismissing(out[i]) && out[i] != 0.0 && j - i + 1 >= min_run && (out[i:j] .= missing)
        i = j + 1
    end
    return out
end
