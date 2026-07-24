# Stage 7: local ACCESS-S2 seasonal-forecast reader.
#
# Real local files (supplied by a third party, already on disk -- no network
# access needed) live at z:/ACCESS-S2/<YYYYMMDD>_<var>.nc: one issue date and
# variable per file, dims (lon, lat, time, ensemble), time already CF-decoded
# to real `DateTime`s, ensemble = 1:99, variables rain/tmax/tmin/radn/evap/vapr.
# This does *not* match RasterDataSources.jl's `ACCESSS` type (which assumes a
# live THREDDS download into a per-layer-subfolder layout), so this reader is
# kept entirely project-local rather than upstreamed -- see Stage 7 of the plan.
#
# Implemented as a genuine MicroclimateMapper weather-source + Loader
# extension (the same multiple-dispatch extension idiom already used for
# `MicroclimateMapper.loader(::Type{<:SILO}) = PointQuery()` in
# points_australia.jl) so the *existing* mean-temperature/subdaily-
# disaggregation/derive! pipeline is reused as-is: the variable set here
# (tmax/tmin/rain/radn/vapr) mirrors WorldClim's declaration almost exactly
# (`MicroclimateMapper.jl/src/climate/worldclim.jl`), and `vapr` (actual
# vapour pressure) feeds the same vapour-pressure-deficit/relative-humidity
# derivation SILO drives from RH instead -- no cloud-cover variable is needed
# either, since SILO also has none and MicroclimateMapper derives cloud
# fraction from global_radiation vs. a computed clear-sky reference
# internally (confirmed against NicheMapR's own R implementation,
# `NicheMapR/R/micro_access_s2.R`, which does this ratio explicitly).
#
# Only grid-mode `_load_layers` needs writing -- points-mode
# (`_load_layers_at_points`) comes for free via MicroclimateMapper's generic
# `Loader` fallback (`_load_layers_at_points(loader::Loader, source, ...) =
# ... _load_layers(loader, source, ...) ...`).

using MicroclimateMapper
using Rasters, Rasters.Extents
using Dates

const ACCESS_S2_DIR = "z:/ACCESS-S2"

struct AccessS2Loader <: MicroclimateMapper.Loader end

# Issue date + ensemble member are encoded as type parameters (Rata Die day
# number + plain Int) rather than struct fields: every MicroclimateMapper
# weather source is dispatched on as `::Type{<:Source}` (e.g. `SILO`, not an
# instance of it), so `model.weather_source` must stay a `Type` -- the only
# way to parametrise "which issue date / member" while keeping that
# convention is to bake the value into the type itself.
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
        # `:rain` is mm depth; numerically equal to kg/m^2 (matches SILO's :daily_rain).
        MicroclimateMapper.Variable(MicroclimateMapper.Rainfall(), :rain, u"kg/m^2"),
        MicroclimateMapper.Variable(MicroclimateMapper.GlobalRadiation(), :radn, u"MJ/m^2/d"),
        # Actual vapour pressure direct from the forecast -- no RH-derivation
        # needed (unlike SILO, which only has rh_tmax/rh_tmin).
        MicroclimateMapper.Variable(MicroclimateMapper.ActualVapourPressure(), :vapr, u"hPa"),
    )
end

_access_s2_path(issue_date::Date, name::Symbol) =
    joinpath(ACCESS_S2_DIR, "$(Dates.format(issue_date, "yyyymmdd"))_$(name).nc")

# MicroclimateMapper's Ti-range arithmetic (`_ti_range_for_dates`) assumes the
# loaded Ti axis is a contiguous daily sequence starting Jan 1 of `years[1]`
# and running through Dec 31 of `years[end]` (matching how every other loader
# builds `Ti(1:nsteps)` from a full-year span) -- it maps requested calendar
# dates to Ti positions by day-of-year arithmetic, not by inspecting the
# data's real dates. ACCESS-S2's real data starts on an arbitrary issue date,
# not Jan 1, so the real (issue_date:issue_date+549) slice is padded out to a
# full `years`-spanning array at the matching calendar offset; the padding is
# never touched as long as the caller only requests `dates` within the real
# forecast window (issue_date to issue_date+549 days).
function _pad_to_years(real_layer::Raster, issue_date::Date, years)
    years_v = collect(years)
    total_days = sum(Dates.daysinyear, years_v)
    yi = findfirst(==(year(issue_date)), years_v)
    yi === nothing && error(
        "ACCESS-S2 issue date $issue_date's year is not in the requested `years` " *
        "range $years -- request simulation `dates` within the issue date's forecast window.",
    )
    offset = sum(Dates.daysinyear, years_v[1:yi-1]; init=0) + (dayofyear(issue_date) - 1)
    # The real file always has its full ~550-day forecast horizon regardless
    # of how much of it the caller's `years` span actually covers -- truncate
    # to whatever fits (the caller only slices `[Ti(ti_start:ti_end)]` within
    # its own requested dates afterwards, so days beyond `years` are unneeded).
    nsteps = min(size(real_layer, Ti), total_days - offset)
    filled = Rasters.replace_missing(real_layer, NaN)
    spatial_dims = dims(filled, (X, Y))
    padded = fill(NaN, size(filled, X), size(filled, Y), total_days)
    padded[:, :, offset+1:offset+nsteps] .= parent(filled)[:, :, 1:nsteps]
    return Raster(padded, (spatial_dims..., Ti(1:total_days)); crs=crs(filled))
end

function MicroclimateMapper._load_layers(::AccessS2Loader, source::Type{<:AccessS2},
                                          fields::Tuple, area::Extents.Extent, years)
    issue_date = _issue_date(source)
    member = _ensemble_member(source)
    layers = map(fields) do name
        @info "  loading ACCESS-S2 $name (issue $issue_date, member $member)..."
        path = _access_s2_path(issue_date, name)
        full = Raster(path; name, lazy=true)
        slice = full[Dim{:ensemble}(Rasters.At(member))]
        real_layer = read(crop(slice; to=area, touches=true))
        _pad_to_years(real_layer, issue_date, years)
    end
    return NamedTuple{fields}(layers)
end
