# sites.jl — per-site inputs derived from each site's OzFlux file and MODIS, so every site on
# disk runs without hand-written entries. tables/sites.toml overrides these only where a site
# has measurements (Whroo's soil pit, leaf area index).

using Downloads

# A site's global attributes, read once, with data_corrections.jl's AttributeOverrides
const _SITE_ATTRIBUTES = Dict{String,Dict{String,Any}}()
site_attributes(site_name) = get!(_SITE_ATTRIBUTES, site_name) do
    attributes = NCDataset(ds -> Dict{String,Any}(ds.attrib), ozflux_site_path(OZFLUX_DATA_DIR, site_name))
    for c in site_corrections(site_name)
        c isa AttributeOverride && (attributes[c.attribute] = c.value)
    end
    attributes
end

site_canopy_height(site_name) = _parse_length_attr(site_attributes(site_name)["canopy_height"])

# (height_m, name) of the air temperature, wind speed and absolute humidity sensors the model is
# compared against: one per quantity and height (forcing_selection.jl's preference), each with a
# single height below the reference height and data in `years`, less the forcing and its backups, sonic
# temperatures (virtual, not air, temperature) and data_corrections.jl's FaultySensor ones
function discover_height_series(site_name, years, forcing, backups)
    corrections = site_corrections(site_name)
    codes = OZFLUX_QUANTITY_CODES
    forcing_names = Set(String[forcing.air_temperature, forcing.humidity, forcing.wind_speed,
        (String(n) for names in values(backups) for n in names)...])
    NCDataset(ozflux_site_path(OZFLUX_DATA_DIR, site_name)) do ds
        rows = findall(t -> Dates.year(t) in years, ds["time"][:])
        series = Tuple{Float64,String}[]
        for code in (codes.air_temperature, codes.wind_speed, codes.absolute_humidity)
            candidates = sort(filter(name -> measures(name, code) && !(name in forcing_names) &&
                !measures_sonic_temperature(name) && !has_correction(corrections, FaultySensor, name) && has_data(ds, name, rows, corrections),
                collect(keys(ds))); by=preference)
            taken = Float64[]
            for name in candidates
                heights = sensor_heights(ds, name)
                length(heights) == 1 && only(heights) * u"m" < forcing.reference_height &&
                    !(only(heights) in taken) || continue
                push!(taken, only(heights))
                push!(series, (only(heights), name))
            end
        end
        sort(series)
    end
end

# Whether a site's `soil` description names a peat: "peat", or the Australian Soil
# Classification's Organosol order
file_names_peat(site_name) = occursin(r"\b(peat|organosol)"i, string(get(site_attributes(site_name), "soil", "")))

# Where a site's soil comes from: its soil_source setting (:horizons), else :peat where its
# file names a peat, which SLGA's 90 m grid misses, else :slga
function soil_source(site_name)
    source = site_setting(site_name, "soil_source")
    source === nothing || return Symbol(source)
    return file_names_peat(site_name) ? :peat : :slga
end

# Minutes the site's local standard time is ahead of UTC, from its file's time_zone
site_utc_offset(site_name) =
    Minute(OZFLUX_CONVENTIONS["time_zone_utc_offset_minutes"][site_attributes(site_name)["time_zone"]])

# Calendar years a site runs over, from FIRST_MODIS_YEAR so MODIS gives its leaf area index:
# the whole year with the most half-hours of good tower forcing and soil temperature (the
# latest on ties), else every year the record touches, for records shorter than any whole year.
# Times label each half-hour's end, so a year's last record is 00:00 on 1 January.
function run_years(site_name)
    NCDataset(ozflux_site_path(OZFLUX_DATA_DIR, site_name)) do ds
        times = ds["time"][:]
        years = Dates.year.(times .- Minute(30))
        good = trues(length(years))
        for name in BEST_YEAR_VARIABLES
            good .&= haskey(ds, name * "_QCFlag") ? in.(ds[name * "_QCFlag"][:], Ref(OZFLUX_GOOD_QC_FLAGS)) : false
        end
        whole(y) = first(times) <= DateTime(y, 1, 1, 0, 30) && last(times) >= DateTime(y + 1, 1, 1)
        recorded = filter(>=(FIRST_MODIS_YEAR), unique(years))
        isempty(recorded) && error("$site_name has no record from $FIRST_MODIS_YEAR")
        candidates = filter(whole, recorded)
        isempty(candidates) && return recorded
        counts = [count(good[years .== y]) for y in candidates]
        [candidates[findlast(==(maximum(counts)), counts)]]
    end
end

# First and last day of the site's record within `years`
function record_period(site_name, years)
    NCDataset(ozflux_site_path(OZFLUX_DATA_DIR, site_name)) do ds
        days = filter(d -> Dates.year(d) in years, Date.(ds["time"][:] .- Minute(30)))
        return first(days), last(days)
    end
end

# MODIS date code for the subset service: A, the year, then the day of year
_modis_date(date) = "A$(Dates.year(date))$(lpad(dayofyear(date), 3, '0'))"

# Median MODIS MCD15A3H leaf area index (500 m, 4-day composites) over the site's record within
# `years`, from the ORNL DAAC subset service, cached in MODIS_LAI_DIR. The service returns at
# most 10 composites (40 days) per request, so the requests run concurrently, retried with
# backoff since the service stalls at times; values above 10 m²/m² are fill flags.
function modis_leaf_area_index(site_name, latitude, longitude, years)
    first_day, last_day = record_period(site_name, years)
    path = joinpath(MODIS_LAI_DIR, "$(site_name)_$(first_day)_$(last_day).csv")
    if !isfile(path)
        mkpath(MODIS_LAI_DIR)
        lat, lon = ustrip(u"°", latitude), ustrip(u"°", longitude)
        # A request can't span two years
        windows = [(start, min(start + Day(39), year_end))
            for y in years
            for (year_start, year_end) in ((max(first_day, Date(y, 1, 1)), min(last_day, Date(y, 12, 31))),)
            for start in year_start:Day(40):year_end]
        rows = reduce(vcat, asyncmap(windows; ntasks=length(windows)) do (start, stop)
            window = "startDate=$(_modis_date(start))&endDate=$(_modis_date(stop))"
            url = "https://modis.ornl.gov/rst/api/v1/MCD15A3H/subset?latitude=$lat&longitude=$lon&$window&kmAboveBelow=0&kmLeftRight=0&band=Lai_500m"
            fetch() = String(take!(Downloads.download(url, IOBuffer(); headers=["Accept" => "text/csv"])))
            body = retry(fetch; delays=ExponentialBackOff(n=4, first_delay=5.0))()
            map(split(strip(body), '\n')) do line
                fields = split(line, ',')
                (date=Date(parse(Int, fields[3][2:5])) + Day(parse(Int, fields[3][6:8]) - 1),
                 leaf_area_index=parse(Int, fields[end]) / 10)
            end
        end)
        CSV.write(path, DataFrame(sort(rows; by=r -> r.date)))
    end
    return median(filter(<=(10.0), CSV.read(path, DataFrame).leaf_area_index))
end

# Site elevation from the file's `elevation` or `altitude` attribute, else SRTM (90 m) from the
# OpenTopoData service, cached in SRTM_ELEVATION_DIR (8 sites name neither)
function site_elevation(site_name, latitude, longitude)
    attributes = site_attributes(site_name)
    for key in ("elevation", "altitude")
        haskey(attributes, key) && return _parse_length_attr(attributes[key])
    end
    path = joinpath(SRTM_ELEVATION_DIR, "$site_name.txt")
    if !isfile(path)
        mkpath(SRTM_ELEVATION_DIR)
        url = "https://api.opentopodata.org/v1/srtm90m?locations=$(ustrip(u"°", latitude)),$(ustrip(u"°", longitude))"
        body = String(take!(Downloads.download(url, IOBuffer())))
        write(path, match(r"\"elevation\":\s*([-\d.]+)", body).captures[1])
    end
    return parse(Float64, read(path, String)) * u"m"
end

# Hourly air pressure: the tower's, gap-filled; the standard atmosphere at the site's elevation
# where the file has none (Ridgefield)
function site_pressure(ps, elevation)
    any(!ismissing, ps) && return uconvert.(u"Pa", fill_gaps(ps) .* u"kPa")
    return fill(uconvert(u"Pa", atmospheric_pressure(elevation)), length(ps))
end

# Measured leaf area index where a site has it, MODIS otherwise
function site_leaf_area_index(site_name, latitude, longitude, years)
    measured = site_setting(site_name, "measured_leaf_area_index")
    return measured === nothing ? modis_leaf_area_index(site_name, latitude, longitude, years) : measured
end

# Every site has a canopy and litter; woody sites' differ from those of grass, crops and pasture
is_woody(site_name) = site_canopy_height(site_name) >= site_setting(site_name, "woody_canopy_height")
# Stack LitterLayer takes cover_fraction, not area_per_mass/unevenness: the cover is the reference's
# 1 - gap fraction, gap = exp(-τ) (unevenness 0) or (1 + unevenness·τ)^(-1/unevenness), τ = area_per_mass·dry_mass
function site_litter(site_name)
    litter_choice == "none" && return isdefined(Microclimate, :NoLitter) ? Microclimate.NoLitter() : nothing
    settings = site_setting(site_name, is_woody(site_name) ? "woody_litter" : "pasture_litter")
    (; dry_mass, thickness, area_per_mass, albedo) = settings
    unevenness = get(settings, :unevenness, 0.0)
    τ = uconvert(NoUnits, area_per_mass * dry_mass)
    gap = unevenness == 0 ? exp(-τ) : (1.0 + unevenness * τ)^(-1 / unevenness)
    return Microclimate.LitterLayer(; dry_mass, thickness, cover_fraction = 1.0 - gap, albedo)
end
site_ground_roughness_height(site_name) = site_setting(site_name, "woody_ground_roughness_height")
