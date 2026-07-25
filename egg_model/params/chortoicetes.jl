# egg traits

# morphometric traits
initial_egg_mass        = 3.6u"mg"
minimum_egg_mass        = 2.6u"mg"
maximum_egg_mass        = 2.3 * initial_egg_mass
dry_mass                = 0.1 * initial_egg_mass
axis_ratio              = 1.82 / 0.69
egg_density             = 1000.0u"kg/m^3"

# hydric traits
base_K_e                = 2.347802e-9 * u"kg/m^2/s/(J/kg)"
hydraulic_conductance   = 2.347802e-9 * u"kg/m^2/s/(J/kg)"
specific_hydration      = 3.04e-4 * u"m^3/m^3/(J/kg)"
conduction_fraction     = 0.5
skin_wetness            = 0.0035

# arrest traits
cold_temperature        = u"K"(12.5u"°C")
diapause_window         = (0.45, 0.50) # development thresholds
quiescence_windows      = ((0.25, 0.30), (0.45, 0.50)) # development thresholds
cold_hour_threshold     = 0.0u"d" # 30u"d"
diapause_hour_threshold = 1240.0u"hr"

# hydric limits
desiccation_tolerance   = 0.6
critical_water_ratio    = 0.6 # TODO resolve difference with desiccation_tolerance

conductance_threshold   = 0.25
wetness_threshold       = 0.45
dormant_conductance     = 0.0u"kg/m^2/s/(J/kg)"
active_conductance      = base_K_e * 3
dormant_wetness         = skin_wetness
active_wetness          = skin_wetness * 100

# Arrhenius thermal response (TODO check if can/need to add Kelvin)
T_A                     = 6641.6175
T_AL                    = 33600.0
T_AH                    = 48000.0
T_L                     = 289.15
T_H                     = 314.65
T_ref                   = 301.65
rate_at_reference       = 1 / 17.4
rate_unit               = 1.0u"d^-1"
lower_lethal_temperature = u"K"(-5.0u"°C")
upper_lethal_temperature = u"K"(52.0u"°C")