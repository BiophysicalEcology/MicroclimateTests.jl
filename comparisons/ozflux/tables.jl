# tables.jl — the data tables in tables/ (OzFlux naming conventions, data corrections, site
# settings, soil textures), read into Julia values. Every fact about the datasets and sites
# lives in those tables, none in code.

using TOML

const TABLES_DIR = joinpath(@__DIR__, "tables")
read_table(name) = TOML.parsefile(joinpath(TABLES_DIR, "$name.toml"))

# A table value as Julia: "<number> <unit>" strings as quantities, tables as NamedTuples
const QUANTITY_PATTERN = r"^(-?[\d.]+(?:e-?\d+)?) (\S+)$"
function table_value(value::AbstractString)
    m = match(QUANTITY_PATTERN, value)
    return m === nothing ? String(value) : parse(Float64, m[1]) * uparse(m[2])
end
table_value(value::AbstractDict) = NamedTuple(Symbol(k) => table_value(v) for (k, v) in value)
table_value(value::AbstractVector) = map(table_value, value)
table_value(value) = value

# ── OzFlux naming conventions ─────────────────────────────────────────────────
const OZFLUX_CONVENTIONS = read_table("ozflux_conventions")
const OZFLUX_GOOD_QC_FLAGS = OZFLUX_CONVENTIONS["good_qc_flags"]
const OZFLUX_STATISTIC_SUFFIXES = OZFLUX_CONVENTIONS["statistic_suffixes"]
const BEST_YEAR_VARIABLES = OZFLUX_CONVENTIONS["best_year_variables"]
const OZFLUX_QUANTITY_CODES = table_value(OZFLUX_CONVENTIONS["quantity_codes"])
const OZFLUX_INSTRUMENT_CODES = Dict(code => table_value(instrument) for (code, instrument) in OZFLUX_CONVENTIONS["instrument_codes"])

# ── Soil textures ─────────────────────────────────────────────────────────────
const CAMPBELL_NORMAN_TEXTURES = Dict(Symbol(name) => table_value(texture) for (name, texture) in read_table("soil_textures"))
const CAMPBELL_NORMAN_MICROPOINT_NAME = Dict(name => texture.micropoint_name for (name, texture) in CAMPBELL_NORMAN_TEXTURES)
const PEAT = table_value(read_table("peat"))

# ── Site settings ─────────────────────────────────────────────────────────────
const SITE_SETTINGS = read_table("sites")

# A site's setting `key`: its own, else the default, merged key by key for sub-tables; nothing
# when neither sets it
function site_setting(site_name, key)
    default = get(SITE_SETTINGS["defaults"], key, nothing)
    own = get(get(SITE_SETTINGS["sites"], site_name, Dict{String,Any}()), key, nothing)
    value = own isa AbstractDict && default isa AbstractDict ? merge(default, own) : something(own, default, Some(nothing))
    return value === nothing ? nothing : table_value(value)
end

site_albedo(site_name) = site_setting(site_name, "albedo")
leaf_parameters(site_name) = LeafParameters(; site_setting(site_name, "leaf")...)
pai_shape(site_name) = Symbol(site_setting(site_name, "pai_shape"))
forcing_backups(site_name) = something(site_setting(site_name, "forcing_backups"), (;))
soil_horizons(site_name) = site_setting(site_name, "soil_horizons")
stomatal_closure_potential(site_name) = site_setting(site_name, "stomatal_closure_potential")
soil_hydraulic_model(site_name) = Microclimate.example_soil_hydraulic_model(;
    stomatal_closure_potential = stomatal_closure_potential(site_name),
    rainfall_entry_mode = rainfall_entry_mode_choice,
)
