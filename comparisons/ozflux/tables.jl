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
const LAND_COVER_CLASSES = read_table("land_cover_classes")["classes"]

# A site's own sites.toml table, empty when it has none
own_site_settings(site_name) = get(SITE_SETTINGS["sites"], site_name, Dict{String,Any}())

# A site's land-cover class: its own land_cover_class, else land_cover_classes.toml's for its FLUXNET ID
function land_cover_class(site_name)
    own = get(own_site_settings(site_name), "land_cover_class", nothing)
    own === nothing || return own
    id = get(site_attributes(site_name), "fluxnet_id", "")
    haskey(LAND_COVER_CLASSES, id) || error("$site_name: FLUXNET ID \"$id\" not in land_cover_classes.toml; set its land_cover_class in sites.toml")
    return LAND_COVER_CLASSES[id]
end

# A site's vegetation type (sites.toml's [vegetation.<type>]), from its land-cover class
function vegetation_type(site_name)
    class = land_cover_class(site_name)
    haskey(SITE_SETTINGS["land_cover_vegetation"], class) || error("$site_name: land-cover class $class has no entry in sites.toml's [land_cover_vegetation]")
    return SITE_SETTINGS["land_cover_vegetation"][class]
end

# Merges key by key where both are tables, else the later replaces the earlier
layer_setting(earlier, later) = earlier isa AbstractDict && later isa AbstractDict ? merge(earlier, later) : something(later, Some(earlier))

# A site's setting `key`: the default, then its vegetation type's, then its own; nothing when none sets it
function site_setting(site_name, key)
    default = get(SITE_SETTINGS["defaults"], key, nothing)
    vegetation = get(SITE_SETTINGS["vegetation"][vegetation_type(site_name)], key, nothing)
    own = get(own_site_settings(site_name), key, nothing)
    value = layer_setting(layer_setting(default, vegetation), own)
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
