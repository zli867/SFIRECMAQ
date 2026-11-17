import numpy as np
import netCDF4 as nc
from utils import write_fire_emissions, WRFGridInfo, get_sfire_emissions, dailyVerticalHourlyEmission, calculate_Area_FRP
from shapely import geometry
from datetime import datetime
import json
from MECHS.Mechanisms import tracer_to_species, SFIRE_to_CMAQ
from shapely.geometry import shape
import pickle


select_dates = [datetime(2021, 3, 22), datetime(2021, 4, 7),
                datetime(2021, 4, 20), datetime(2022, 4, 25)]
# select_dates = [datetime(2021, 3, 20)]
for select_date in select_dates:
    filename = "/Volumes/T9/Coupling/WRF-SFIRE/wrfout_d01_%s_00:00:00" % (select_date.strftime("%Y-%m-%d"))
    fire_file = "/Volumes/Shield/SFIRECMAQ/data/FtBn_BurnInfo.json"
    metcros_file = "/Volumes/T9/Coupling/CMAQ/mcip/FtBnd03/METCRO3D_FtBnd03_%s.nc" % select_date.strftime("%Y%m%d")
    metdot_file = "/Volumes/T9/Coupling/CMAQ/mcip/FtBnd03/METDOT3D_FtBnd03_%s.nc" % select_date.strftime("%Y%m%d")
    emission_flux_output = "/Users/zongrunli/Desktop/sfirecmaq_emiss/emis_mole_surface_%s_hires1_cmaq_cb6_2016fh_16j_fire_CP_freitas.ncf" % select_date.strftime("%Y%m%d")
    tracer_names = ["tr17_1", "tr17_2", "tr17_3", "tr17_4", "tr17_5", "tr17_6", "tr17_7", "tr17_8"]

    with open(fire_file) as json_file:
        fire_events = json.load(json_file)
    rx_polygons = []
    for fire_event in fire_events["fires"]:
        fire_date = datetime.strptime(fire_event["date"], "%Y-%m-%d")
        if fire_date == select_date:
            rx_polygons.append(shape(fire_event["perimeter"]))
    rx_polygon = geometry.MultiPolygon(rx_polygons)
    unit_centroid = rx_polygon.centroid
    ignition_lon = unit_centroid.x
    ignition_lat = unit_centroid.y
    wrf_fire_data = nc.Dataset(filename)
    wrf_info = WRFGridInfo(wrf_fire_data)
    burned_area, frp_values = calculate_Area_FRP(wrf_fire_data)
    # Get Hourly Emissions
    emissions_dict = get_sfire_emissions(wrf_fire_data, tracer_names)
    start_time = None
    for i in range(0, len(emissions_dict["time"])):
        if emissions_dict[tracer_names[0]][i] > 0:
            start_time = emissions_dict["time"][i]
            break
    # Mapping to CMAQ
    CMAQ_emissions = {}
    for tracer_name in tracer_names:
        cur_specie_name = tracer_to_species(tracer_name)
        cur_mapping = SFIRE_to_CMAQ(cur_specie_name)
        cur_emission = emissions_dict[tracer_name]
        for key in cur_mapping:
            CMAQ_emissions[key] = cur_mapping[key][0] * cur_emission
    # Get Time Averaged FRP and Burned Area
    burned_area = (burned_area[1:] + burned_area[:-1]) / 2
    frp_values = (frp_values[1:] + frp_values[:-1]) / 2
    CMAQ_emissions["time"] = emissions_dict["time"]
    CMAQ_emissions["FRP"] = frp_values
    CMAQ_emissions["area"] = burned_area
    CMAQ_emissions["start_time"] = start_time

    emission_tensor, plume_info = dailyVerticalHourlyEmission(metcros_file, metdot_file, ignition_lon, ignition_lat, CMAQ_emissions)
    write_fire_emissions(emission_tensor, metcros_file, emission_flux_output)
    with open('plume_data_%s.pkl' % datetime.strftime(select_date, "%Y-%m-%d"), 'wb') as f:
        pickle.dump(plume_info, f, protocol=pickle.HIGHEST_PROTOCOL)