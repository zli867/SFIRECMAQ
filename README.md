# SFIRE-F0AM-CMAQ v1.0

Offline coupling of a fire behavior model (WRF-SFIRE) with a chemical transport model (CMAQ) for prescribed-burn smoke. The code turns the inert smoke tracers simulated by WRF-SFIRE into reactive emission inputs for CMAQ, keeping the plume structure resolved by the fire model.

## Coupling algorithms

| Algorithm | Script | Description |
|---|---|---|
| **SFIRE-F0AM-CMAQ** | `ChemEmissionGenerator.py` | Mass flux of WRF-SFIRE smoke across a coupling surface, with near-source chemistry simulated by F0AM smoke and no-smoke box models. |
| **SFIRE-CMAQ** | `VanillaEmissionGenerator.py` | Mass flux of WRF-SFIRE smoke across the coupling surface, without near-source chemistry. |
| **Freitas-SFIRE-CMAQ** | `FreitasEmissionGenerator.py` | WRF-SFIRE emissions and heat flux passed to the Freitas plume rise scheme; emissions injected at the burn-unit centroid. |

Each script writes a 3-D, time-varying CMAQ emission file (IOAPI netCDF).

## Repository structure

```
ChemEmissionGenerator.py     SFIRE-F0AM-CMAQ driver
VanillaEmissionGenerator.py  SFIRE-CMAQ driver
FreitasEmissionGenerator.py  Freitas-SFIRE-CMAQ driver
CouplingHelper.py            Coupling distance, coupling surface, plume mask, mass flux
utils.py                     WRF/CMAQ grid handling, emission file writing
coupling_utils/              Emission, entrainment, photolysis, and transport terms for F0AM
interfaces/                  Converters: SFIRE->CMAQ, SFIRE->F0AM, F0AM->CMAQ
MECHS/                       Species mapping between WRF-SFIRE tracers, F0AM (CB6), and CMAQ (CB6-AERO7)
plumerise/                   Freitas plume rise model (Fortran source + Python wrapper)
F0AM/                        F0AM v4.3 and the box-model driver SFIRECMAQ.m
CCTM/                        Modified CMAQ v5.3.2 CCTM source (see below)
```

## Requirements

- **Python 3.10** with `numpy`, `scipy`, `netCDF4`, `shapely`, `pyproj`, `kneed`, `matplotlib`
- **MATLAB** (proprietary; required by F0AM, used only by SFIRE-F0AM-CMAQ)
- **F0AM v4.3** (Wolfe et al., 2016), included in `F0AM/Model`
- **CMAQ v5.3.2**, compiled from the modified source in `CCTM/`
- A Fortran compiler, to build the Freitas plume rise module

The included `plumerise/plumeRiseFreitas.cpython-310-darwin.so` was built for macOS with Python 3.10. On other systems, rebuild it:

```bash
cd plumerise
python -m numpy.f2py -c plumeRiseFreitas.f -m plumeRiseFreitas
```

## Inputs

- WRF-SFIRE output (`wrfout`) with inert smoke tracers
- Burn-unit perimeters (GeoJSON-style JSON)
- CMAQ MCIP files (`METCRO3D`, and `METDOT3D` for Freitas-SFIRE-CMAQ)
- For SFIRE-F0AM-CMAQ only: CMAQ background run without the studied fire (`ACONC` concentrations and `PHOTDIAG2` photolysis rates)

## Usage

1. Edit the file paths, burn dates, and tracer names at the top of the driver script. For SFIRE-F0AM-CMAQ, also set `matlab_bin` and `f0am_dir`.
2. Run the driver, e.g.
   ```bash
   python ChemEmissionGenerator.py
   ```
3. Use the output file as a gridded emission stream in the modified CMAQ, with a biomass-burning source tag in the emission control namelist.

## Modifications to CMAQ v5.3.2

The CCTM source in `CCTM/` is modified to:
- turn off linear temporal interpolation of emissions, so the time profile of the emission input is kept;
- accept secondary species (e.g., ozone) in the emission input;
- allow negative emission values.