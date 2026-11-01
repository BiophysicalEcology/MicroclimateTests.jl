# Whroo shortwave albedo: findings and plan

## Problem
Predicted albedo (upwelling / downwelling shortwave, Fsu/Fsd) is too high at Whroo 2016:
midday predicted 0.165 vs measured 0.107; 7am/6pm predicted ~0.19 vs measured ~0.15.

## Evidence
- Tower is right: MODIS MCD43A black-sky albedo for the Whroo pixel (ORNL site `au_victoria_whroo`)
  is 0.110–0.134 shortwave, 0.024–0.049 visible, 0.165–0.204 near-infrared (2016).
- The two-stream code is correct: matches an independent numerical two-stream solution to 4 decimals,
  diffuse and direct. (Exactly black leaves give 0/0 → 0; irrelevant in practice.)
- Leaves are bright, not dark: Hume, McVicar & Roderick (2002, CRC Catchment Hydrology 02/3), 8 eucalypt species
  incl. E. polyanthemos and E. pauciflora: red reflectance 0.11–0.14, transmittance 0.00–0.03;
  near-infrared (725–1100 nm) reflectance 0.49–0.52, transmittance 0.45 beam / 0.25 diffuse.
- Ground is bright: pale grey litter (TERN panorama photos) and topsoil 7.5YR 5/3 dry
  (Drake 2015 pit) → ~0.23 by Post et al. (2000), albedo = 0.069 × Munsell value − 0.114.
- Shrubs are sparse: TERN NVIS structure survey 2014, core 1 ha: trees 5–30 m, 20 % crown cover;
  shrubs 0–2 m, 5 % foliage cover. The hemispherical-photo LAI (0.8) is taken at 1.3 m, so excludes shrubs.
- Leaf angle, plant area index (0.8 → 1.5) and a visible/near-infrared band split each change
  predicted albedo by < 0.02. With measured leaf and ground optics, the evenly mixed canopy
  cannot reach the measured albedo.
- MODIS visible albedo is lowest when the midday sun is lowest (Jul/Sep), consistent with tall
  sparse crowns shadowing pale ground.
- CABLE has no structural fix: same evenly mixed albedo (`cbl_albedo.F90`), its clumping factor (`Kclump`)
  lives only in POP. It matches eucalypt albedo by making leaves too dark (near-infrared reflectance = transmittance = 0.29).

## Plan
1. Add crown structure to canopy shortwave: a clumping index (Nilson/Chen) as the minimal version,
   or a geometric-optical crown-shadowing term (Ni-Meister et al. 2010; Kobayashi et al. 2012, CLM 3D canopy).
   Inputs from TERN data: crown cover 20 %, tree height 5–30 m, stem density.
2. Split shortwave into visible and near-infrared bands, so each can be checked against MODIS band albedos.
3. Then set leaf optics from Hume et al. (2002) and ground albedo from the soil colour/litter evidence,
   without tuning, and compare the daily albedo curve and MODIS seasonal bands.
