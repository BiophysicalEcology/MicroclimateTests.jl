using ThermalPhysiology
using Unitful

# Always pass Unitful temperatures to ThermalPhysiology models — a bare Float64
# is interpreted as °C by their `_K` conversion, not K, even though the models'
# own fields store bare-Float64 Kelvin internally.

# builds a RateModel from the 5-parameter Sharpe-DeMichele/Schoolfield form,
# matching the DEBtool-normalised (rate_at_reference exact at T_ref) variant.
# `rate_at_reference` is unitless by ThermalPhysiology's own convention (its
# fields are fixed Float64); `rate_unit` attaches the actual physical unit.
function arrhenius_development_model(;
    T_A, T_AL, T_AH, T_L, T_H, T_ref, rate_at_reference, rate_unit=1.0u"hr^-1",
)
    tpc = SharpSchoolDEBModel(; T_A, T_AL, T_AH, T_L, T_H, T_ref, rate_at_reference)
    RateModel(tpc, rate_unit)
end
