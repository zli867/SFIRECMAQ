import os.path
import numpy as np                                  # <<< NEW: to zero the emissions dict
from interfaces.F0AM2CMAQ import calculate_chem_emission_tensor
import netCDF4 as nc
from utils import CMAQGridInfo, write_fire_emissions, WRFGridInfo
from shapely import geometry
from datetime import datetime
import json
from CouplingHelper import calculate_coupling_distance, coupling_surface_buffer, vanilla_emission_tensor, \
    coarse_emission_tracer, inside_mass, mask_volume, calculate_plume_mask, run_matlab
from shapely.geometry import shape
from interfaces.SFIRE2F0AM import createF0AMInputs
from CouplingHelper import get_plume_met, inside_smoke_mask
from coupling_utils.Transport import calculate_transport_ratio
from coupling_utils.Emission import calculate_mapped_emission
from coupling_utils.Entrainment import calculate_entrain_mass
from coupling_utils.Photolysis import calculate_photo_rate
from scipy.io import savemat

# select_dates = [datetime(2021, 3, 22), datetime(2021, 4, 7), datetime(2021, 4, 20)]
select_date = datetime(2021, 4, 20)
filename = "/Volumes/T9/Coupling/WRF-SFIRE/wrfout_d01_%s_00:00:00" % (select_date.strftime("%Y-%m-%d"))
fire_file = "/Volumes/Shield/SFIRECMAQ/data/FtBn_BurnInfo.json"
bkg_CMAQ_file = "/Volumes/T9/Coupling/CMAQ/bkg/CCTM_ACONC_v532_intel_FtBnd03_%s.nc" % select_date.strftime("%Y%m%d")
photo_diag_file = "/Volumes/T9/Coupling/CMAQ/bkg/CCTM_PHOTDIAG2_v532_intel_FtBnd03_%s.nc" % select_date.strftime("%Y%m%d")
metcros_file = "/Volumes/T9/Coupling/CMAQ/mcip/FtBnd03/METCRO3D_FtBnd03_%s.nc" % select_date.strftime("%Y%m%d")
work_dir = "/Users/zongrunli/Desktop/SFIRECMAQV2/work_dir/%s" % select_date.strftime("%Y%m%d")
matlab_bin = "/Applications/MATLAB_R2025b.app/bin/matlab"
f0am_dir = "/Users/zongrunli/Desktop/SFIRECMAQV2/F0AM"
emission_flux_output = "/Users/zongrunli/Desktop/sfirecmaq_emiss/emis_mole_surface_%s_hires1_cmaq_cb6_2016fh_16j_fire_CP_CPChem_cln2.ncf" % select_date.strftime("%Y%m%d")

os.makedirs(work_dir, exist_ok=True)
with open(fire_file) as json_file:
    fire_events = json.load(json_file)
rx_polygons = []
for fire_event in fire_events["fires"]:
    fire_date = datetime.strptime(fire_event["date"], "%Y-%m-%d")
    if fire_date == select_date:
        rx_polygons.append(shape(fire_event["perimeter"]))
rx_polygon = geometry.MultiPolygon(rx_polygons)
# Some maybe not good for merge, add a check better
unit_centroid = rx_polygon.centroid
ignition_lon = unit_centroid.x
ignition_lat = unit_centroid.y

tracer_names = ["tr17_1", "tr17_2", "tr17_3", "tr17_4", "tr17_5", "tr17_6", "tr17_7", "tr17_8"]
smoke_tracer = "tr17_2"
quantile_vertical = 0.1
quantile_vertical_low = 0.9
quantile_horizontal_high = 0.1
# unit: km
bin_step = 1
wrf_fire_data = nc.Dataset(filename)
wrf_info = WRFGridInfo(wrf_fire_data)
cmaq_data = nc.Dataset(bkg_CMAQ_file)
cmaq_info = CMAQGridInfo(cmaq_data)
coupling_distance = calculate_coupling_distance(smoke_tracer, wrf_fire_data, ignition_lon, ignition_lat)

print("Coupling distance is %.2f" % coupling_distance)
smoke_mask = calculate_plume_mask(smoke_tracer, wrf_fire_data, quantile_vertical_high=0.1, quantile_vertical_low=0.9, quantile_horizontal=0.1)
inside_mask = inside_smoke_mask(smoke_mask, wrf_fire_data, ignition_lon, ignition_lat, coupling_distance)
smoke_volume = mask_volume(wrf_fire_data, inside_mask)
smoke_met = get_plume_met(wrf_fire_data, inside_mask)
total_inside_mass = inside_mass(smoke_tracer, wrf_fire_data, ignition_lon, ignition_lat, coupling_distance)  # ug
flux_surface_dict = coupling_surface_buffer(wrf_fire_data, ignition_lon, ignition_lat, coupling_distance)
sfire_emission_tensor = vanilla_emission_tensor(smoke_tracer, wrf_fire_data, flux_surface_dict)
coarse_emission_tensor = coarse_emission_tracer(wrf_fire_data, cmaq_info, sfire_emission_tensor)

# Smoke Volume
print("Smoke Property Calculating...")
smoke_volume_data = {"time": wrf_info["time"], "V": smoke_volume}

# Transport Ratio
print("Smoke Transport Calculating...")
transport_ratio = calculate_transport_ratio(total_inside_mass, sfire_emission_tensor, wrf_info)

# Emissions
print("Emission Calculating...")
emissions = calculate_mapped_emission(wrf_fire_data, tracer_names)
print(emissions.keys())
print(emissions["time"])
# Background CMAQ
print("Entrainment Calculating...")
bkg_conc = calculate_entrain_mass(bkg_CMAQ_file, wrf_info, smoke_volume)

# ============================ BOX 1: FIRE ============================
print("F0AM Calculating (FIRE box)")
foam_input = createF0AMInputs(smoke_met, smoke_volume_data, emissions, transport_ratio, bkg_conc)
foam_input["J"] = calculate_photo_rate(photo_diag_file, foam_input["time"])
time_vals = foam_input['time']
foam_input['time'] = [dt.isoformat() for dt in time_vals]
savemat(os.path.join(work_dir, "foam_input.mat"), foam_input)
run_matlab("SFIRECMAQ", f"{work_dir}/foam_input.mat",     f"{work_dir}/sfire_res.mat", matlab_bin, f0am_dir)

print("F0AM Calculating (NO-FIRE background box)")
foam_input_bkg = createF0AMInputs(smoke_met, smoke_volume_data, emissions, transport_ratio, bkg_conc, no_fire=True)
foam_input_bkg["J"] = calculate_photo_rate(photo_diag_file, foam_input_bkg["time"])
foam_input_bkg["time"] = [dt.isoformat() for dt in foam_input_bkg["time"]]
savemat(os.path.join(work_dir, "foam_input_bkg.mat"), foam_input_bkg)
run_matlab("SFIRECMAQ", f"{work_dir}/foam_input_bkg.mat", f"{work_dir}/sfire_res_bkg.mat", matlab_bin, f0am_dir)

# ================ CONVERT: IMPACT = FIRE - NO-FIRE ================
foam_output = os.path.join(work_dir, "sfire_res.mat")
foam_output_bkg = os.path.join(work_dir, "sfire_res_bkg.mat")
foam_input["time"] = time_vals
CMAQ_emissions = calculate_chem_emission_tensor(foam_input, foam_output, foam_output_bkg,
                                                coarse_emission_tensor, cmaq_info)

write_fire_emissions(CMAQ_emissions, metcros_file, emission_flux_output)