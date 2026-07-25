# Stage 7: local ACCESS-S2 seasonal-forecast reader.
#
# <YYYYMMDD>_<var>.nc: one issue date and variable per file, dims (lon, lat,
#  time, ensemble), time already CF-decoded to real `DateTime`s, ensemble = 1:99, 
# variables rain/tmax/tmin/radn/evap/vapr.
#
# Implemented as a MicroclimateMapper weather-source + Loader
# extension based on `MicroclimateMapper.loader(::Type{<:SILO}) = 
# PointQuery()` in points_australia.jl).
#
# Grid-mode `_load_layers` here -- points-mode
# (`_load_layers_at_points`) available via MicroclimateMapper's generic
# `Loader` fallback (`_load_layers_at_points(loader::Loader, source, ...) =
# ... _load_layers(loader, source, ...) ...`).

using MicroclimateMapper
using Rasters, Rasters.Extents
using Dates

const ACCESS_S2_DIR = joinpath(ENV["RASTERDATASOURCES_PATH"], "ACCESS-S2")

struct AccessS2Loader <: MicroclimateMapper.Loader end

# Issue date + ensemble member are encoded as type parameters (Rata Die day
# number + plain Int)
struct AccessS2{IssueDateDays,Member} end

AccessS2(issue_date::Date, member::Integer) = AccessS2{Dates.value(issue_date),Int(member)}

_issue_date(::Type{AccessS2{D,M}}) where {D,M} = Date(Dates.UTD(D))
_ensemble_member(::Type{AccessS2{D,M}}) where {D,M} = M

MicroclimateMapper.weather_calendar(::Type{<:AccessS2}) = MicroclimateMapper.Daily()
MicroclimateMapper.loader(::Type{<:AccessS2}) = AccessS2Loader()

function MicroclimateMapper.variables(::Type{<:AccessS2})
    (
        MicroclimateMapper.Variable(MicroclimateMapper.Temperature(MicroclimateMapper.Maximum()), :tmax, u"°C"),
        MicroclimateMapper.Variable(MicroclimateMapper.Temperature(MicroclimateMapper.Minimum()), :tmin, u"°C"),
        MicroclimateMapper.Variable(MicroclimateMapper.Rainfall(), :rain, u"kg/m^2"),
        MicroclimateMapper.Variable(MicroclimateMapper.GlobalRadiation(), :radn, u"MJ/m^2/d"),
        MicroclimateMapper.Variable(MicroclimateMapper.ActualVapourPressure(), :vapr, u"hPa"),
    )
end

_access_s2_path(issue_date::Date, name::Symbol) =
    joinpath(ACCESS_S2_DIR, "$(Dates.format(issue_date, "yyyymmdd"))_$(name).nc")

# MicroclimateMapper's `_ti_range_for_dates` assumes the loaded Ti axis is
# a complete year. ACCESS-S2's real data starts on an arbitrary issue date,
# so need to pad it out.
function _pad_to_years(real_layer::Raster, issue_date::Date, years)
    years_v = collect(years)
    total_days = sum(Dates.daysinyear, years_v)
    yi = findfirst(==(year(issue_date)), years_v)
    yi === nothing && error(
        "ACCESS-S2 issue date $issue_date's year is not in the requested `years` " *
        "range $years -- request simulation `dates` within the issue date's forecast window.",
    )
    offset = sum(Dates.daysinyear, years_v[1:yi-1]; init=0) + (dayofyear(issue_date) - 1)
    nsteps = min(size(real_layer, Ti), total_days - offset)
    filled = Rasters.replace_missing(real_layer, NaN)
    spatial_dims = dims(filled, (X, Y))
    padded = fill(NaN, size(filled, X), size(filled, Y), total_days)
    padded[:, :, offset+1:offset+nsteps] .= parent(filled)[:, :, 1:nsteps]
    return Raster(padded, (spatial_dims..., Ti(1:total_days)); crs=crs(filled))
end

# All 99 ensemble members open + crop the same file over the same area,
# differing only in which member gets read afterward -- cache the lazy
# cropped raster per (file, area) so repeated members skip the file-open and
# crop planning, while still only reading one member's data (not all 99) at
# a time to stay within memory.
const _ACCESS_S2_CROP_CACHE = Dict{Tuple{String,Extents.Extent}, Raster}()

function _cropped_lazy(path::String, name::Symbol, area::Extents.Extent)
    get!(_ACCESS_S2_CROP_CACHE, (path, area)) do
        full = Raster(path; name, lazy=true)
        crop(full; to=area, touches=true)
    end
end

function MicroclimateMapper._load_layers(::AccessS2Loader, source::Type{<:AccessS2},
                                          fields::Tuple, area::Extents.Extent, years)
    issue_date = _issue_date(source)
    member = _ensemble_member(source)
    layers = map(fields) do name
        @info "  loading ACCESS-S2 $name (issue $issue_date, member $member)..."
        path = _access_s2_path(issue_date, name)
        cropped = _cropped_lazy(path, name, area)
        real_layer = read(cropped[Dim{:ensemble}(Rasters.At(member))])
        _pad_to_years(real_layer, issue_date, years)
    end
    return NamedTuple{fields}(layers)
end
