# debug_site.jl
#
# Flat, single-site debugging entry point — calls `run_site` directly with
# no try/catch, so errors propagate normally with a full stack trace pointing
# at the real failure line. Step through top-to-bottom (Ctrl+Enter /
# Shift+Enter per block in the VS Code Julia extension), or
# `include("debug_site.jl")` from the REPL to run it all at once.
#
# All model setup and per-site logic lives in pipeline.jl (shared with
# comparison.jl) and config.jl (shared model/soil-source/weather-source
# parameters — change `weather_source_choice`, `pedotransfer_model_choice`,
# etc. there to compare variants).

using Microclimate, MicroclimateMapper, Unitful
using CSV, DataFrames, Dates, Statistics, Printf, Plots
using Rasters, ArchGDAL, NCDatasets, RasterDataSources, PointDataSources
using Rasters.Extents: Extent
using GeoInterface: Wrappers as GIW
using Serialization

include("config.jl")
include("utils.jl")
include("pipeline.jl")

# ── Pick the site to debug here ───────────────────────────────────────────────
site_name = "a1"

sim_start = Date(2007, 1, 1)
sim_end   = Date(2010, 12, 31)
auto_date_range = true

plot_start = nothing
plot_end   = nothing

siteinfo = DataFrame(CSV.File(OZNET_SITEINFO))
row = siteinfo[findfirst(==(site_name), siteinfo.name), :]

weather_cache = Dict{Tuple{Symbol,String,Date,Date}, Any}()

@time site_rows = run_site(row, weather_cache;
    sim_start, sim_end, auto_date_range, max_sim_years, plot_start, plot_end)

using MySQL
using DBInterface
using DataFrames

host = "115.146.93.180"
uid = "general"
pwd = "predecol"


conn = DBInterface.connect(
    MySQL.Connection,
    host,
    uid,
    pwd;
    db = "AWAPMoist",
    port = 3306,
    ssl_mode = MySQL.API.SSL_MODE_DISABLED,
    ssl_enforce = false,
    ssl_verify_server_cert = false,
)
#databases = DataFrame(DBInterface.execute(conn, "SHOW DATABASES"))
tables = DataFrame(DBInterface.execute(conn, "SHOW TABLES FROM AWAPMoist"))
display(tables)
columns = DataFrame(DBInterface.execute(
    conn,
    "SHOW COLUMNS FROM AWAPMoist.1960"
))
display(columns)
columns = DataFrame(DBInterface.execute(
    conn,
    "SHOW COLUMNS FROM AWAPMoist.latlon"
))
display(columns)
columns = DataFrame(DBInterface.execute(
    conn,
    "SHOW COLUMNS FROM AWAPMoist.latlon_char"
))
display(columns)


function query_awap_moisture(
    conn,
    start_date::Date,
    end_date::Date;
    id::Union{Nothing,Integer}=nothing,
)
    start_date <= end_date ||
        error("start_date must not be after end_date")

    table_df = DataFrame(
        DBInterface.execute(conn, "SHOW TABLES FROM AWAPMoist")
    )

    available_years = Set(
        parse(Int, name)
        for name in string.(table_df[:, 1])
        if occursin(r"^\d{4}$", name)
    )

    requested_years = collect(year(start_date):year(end_date))
    years = filter(in(available_years), requested_years)

    isempty(years) &&
        error("No AWAPMoist tables exist between $start_date and $end_date")

    missing_years = setdiff(requested_years, years)
    if !isempty(missing_years)
        @warn "No AWAPMoist tables found for these years" missing_years
    end

    clauses = String[]
    parameters = Any[]

    for y in years
        first_day = y == year(start_date) ? dayofyear(start_date) : 1
        last_day  = y == year(end_date)   ? dayofyear(end_date)   :
                    daysinyear(y)

        if id === nothing
            push!(clauses, """
                SELECT
                    $y AS yr,
                    id,
                    day,
                    WRel1,
                    WRel2
                FROM AWAPMoist.`$y`
                WHERE day BETWEEN ? AND ?
            """)

            append!(parameters, (first_day, last_day))
        else
            push!(clauses, """
                SELECT
                    $y AS yr,
                    id,
                    day,
                    WRel1,
                    WRel2
                FROM AWAPMoist.`$y`
                WHERE day BETWEEN ? AND ?
                  AND id = ?
            """)

            append!(parameters, (first_day, last_day, id))
        end
    end

    sql = join(clauses, "\nUNION ALL\n") *
          "\nORDER BY yr, day, id"

    # Parameter binding in MySQL.jl requires a prepared statement.
    stmt = DBInterface.prepare(conn, sql)

    result = try
        DataFrame(DBInterface.execute(stmt, Tuple(parameters)))
    finally
        DBInterface.close!(stmt)
    end

    if !isempty(result)
        result.date = [
            Date(y, 1, 1) + Day(d - 1)
            for (y, d) in zip(result.yr, result.day)
        ]

        select!(result, :id, :date, :WRel1, :WRel2)
    end

    return result
end

function nearest_awap_cell(conn, lon, lat)
    sql = """
        SELECT id, latitude, longitude
        FROM AWAPMoist.latlon
        ORDER BY
            POWER(latitude - ?, 2) +
            POWER(longitude - ?, 2)
        LIMIT 1
    """

    stmt = DBInterface.prepare(conn, sql)

    try
        DataFrame(DBInterface.execute(stmt, (lat, lon)))
    finally
        DBInterface.close!(stmt)
    end
end

lon = siteinfo.lon[row]
lat = siteinfo.lat[row]

cell = nearest_awap_cell(conn, lon, lat)
display(cell)

cell_id = cell.id[1]

moisture = query_awap_moisture(
    conn,
    Date(2010, 1, 1),
    Date(2010, 12, 31);
    id=cell_id,
)
display(moisture)