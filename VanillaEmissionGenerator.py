from interfaces.SFIRE2CMAQ import calculate_vanilla_emission_tensor
import netCDF4 as nc
from utils import CMAQGridInfo, write_fire_emissions, WRFGridInfo
from shapely import geometry
from datetime import datetime
import json
from CouplingHelper import calculate_coupling_distance
from shapely.geometry import shape

select_dates = [datetime(2021, 3, 20), datetime(2021, 3, 22), datetime(2021, 4, 7),
                datetime(2021, 4, 20), datetime(2022, 4, 25)]
for select_date in select_dates:
    filename = "/Volumes/T9/Coupling/WRF-SFIRE/wrfout_d01_%s_00:00:00" % (select_date.strftime("%Y-%m-%d"))
    fire_file = "/Volumes/Shield/SFIRECMAQ/data/FtBn_BurnInfo.json"
    metcros_file = "/Volumes/T9/Coupling/CMAQ/mcip/FtBnd03/METCRO3D_FtBnd03_%s.nc" % select_date.strftime("%Y%m%d")
    emission_flux_output = "/Users/zongrunli/Desktop/sfirecmaq_emiss/emis_mole_surface_%s_hires1_cmaq_cb6_2016fh_16j_fire_CP_no_chem.ncf" % select_date.strftime("%Y%m%d")

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

    tracer_names = ["tr17_1", "tr17_2", "tr17_3", "tr17_4", "tr17_5", "tr17_6", "tr17_7", "tr17_8"]
    smoke_tracer = "tr17_2"
    quantile_vertical = 0.1
    quantile_vertical_low = 0.9
    quantile_horizontal_high = 0.1
    # unit: km
    bin_step = 1
    wrf_fire_data = nc.Dataset(filename)
    wrf_info = WRFGridInfo(wrf_fire_data)
    cmaq_info = CMAQGridInfo(nc.Dataset(metcros_file))
    coupling_distance = calculate_coupling_distance(smoke_tracer, wrf_fire_data, ignition_lon, ignition_lat)

    print("Coupling distance is %.2f" % coupling_distance)

    CMAQ_emissions = calculate_vanilla_emission_tensor(smoke_tracer, tracer_names, wrf_fire_data,
                                                       ignition_lon, ignition_lat, coupling_distance, cmaq_info)
    write_fire_emissions(CMAQ_emissions, metcros_file, emission_flux_output)