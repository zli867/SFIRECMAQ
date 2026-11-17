import numpy as np
import pyproj
from datetime import datetime
import netCDF4 as nc
from MECHS.HeaderInfo import HEADERINFO
from MECHS.Mechanisms import get_molecule_weight
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy import stats
from plumerise import plumeRiseFreitas


def heightAboveGround(ds):
    PH = ds["PH"][:]
    PHB = ds["PHB"][:]
    HGT = ds["HGT"][:]
    HGT = HGT[:, np.newaxis, :, :]
    height_above_sea_level = (PH + PHB) / 9.8
    return height_above_sea_level - HGT


def grid_cell_plume_height(vertical_column_weighted_tracer, vertical_column_thickness, quantiles):
    # This calculate the plume height above ground
    if np.sum(vertical_column_weighted_tracer) == 0:
        return 0
    vertical_column_tracer_normalized = vertical_column_weighted_tracer / np.sum(vertical_column_weighted_tracer)
    target_ratio = 1 - quantiles
    layer_num = len(vertical_column_tracer_normalized)
    current_sum = 0
    result_idx = None
    for layer_idx in range(0, layer_num):
        current_sum += vertical_column_tracer_normalized[layer_idx]
        if current_sum >= target_ratio:
            result_idx = layer_idx
            break
    # current plume height is sum of thickness from 0, 1, ..., result_idx
    # search index will be [result_idx - 1, result_idx] -> thickness[result_idx]
    residual_ratio = current_sum - target_ratio
    overflow_mass_ratio = (residual_ratio * np.sum(vertical_column_weighted_tracer)) / vertical_column_weighted_tracer[
        result_idx]
    overflow_height = overflow_mass_ratio * vertical_column_thickness[result_idx]
    plume_height = 0
    for idx in range(0, result_idx + 1):
        plume_height += vertical_column_thickness[idx]
    plume_height -= overflow_height
    return plume_height


def convertByteToStr(byte_list):
    str1 = ""
    for i in range(0, len(byte_list)):
        str1 = str1 + byte_list[i].decode('UTF-8')
    return str1


def CMAQGridInfo(met_ds):
    time_data = met_ds['TFLAG'][:]
    lat_1 = met_ds.getncattr('P_ALP')
    lat_2 = met_ds.getncattr('P_BET')
    lat_0 = met_ds.getncattr('YCENT')
    lon_0 = met_ds.getncattr('XCENT')
    crs = pyproj.Proj("+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=" + str(lat_1)
                      + " +lat_2=" + str(lat_2) + " +lat_0=" + str(lat_0) +
                      " +lon_0=" + str(lon_0))
    xcell = met_ds.getncattr('XCELL')
    ycell = met_ds.getncattr('YCELL')
    xorig = met_ds.getncattr('XORIG')
    yorig = met_ds.getncattr('YORIG')

    ncols = met_ds.getncattr('NCOLS')
    nrows = met_ds.getncattr('NROWS')

    # > for X, Y cell centers
    x_center_range = np.linspace(xorig + xcell / 2, (xorig + xcell / 2) + xcell * (ncols - 1), ncols)
    y_center_range = np.linspace(yorig + ycell / 2, (yorig + ycell / 2) + ycell * (nrows - 1), nrows)

    Xcenters, Ycenters = np.meshgrid(x_center_range, y_center_range)

    # > for X, Y cell boundaries (i.e., cell corners)
    x_bound_range = np.linspace(xorig, xorig + xcell * ncols, ncols + 1)
    y_bound_range = np.linspace(yorig, yorig + ycell * nrows, nrows + 1)

    Xbounds, Ybounds = np.meshgrid(x_bound_range, y_bound_range)

    x_max = np.max(Xbounds)
    x_min = np.min(Xbounds)
    y_max = np.max(Ybounds)
    y_min = np.min(Ybounds)

    cmaq_time_array = []
    for i in range(0, time_data.shape[0]):
        time_data_tmp = time_data[i, 0, :]
        time_str = str(time_data_tmp[0]) + str(time_data_tmp[1]).rjust(6, '0')
        parsed = datetime.strptime(time_str, '%Y%j%H%M%S')
        cmaq_time_array.append(parsed)

    lon_ctr, lat_ctr = crs(Xcenters, Ycenters, inverse=True)
    lat_min = np.min(lat_ctr)
    lat_max = np.max(lat_ctr)
    lon_min = np.min(lon_ctr)
    lon_max = np.max(lon_ctr)
    res_dict = {"crs": crs, "X_ctr": Xcenters, "Y_ctr": Ycenters, "X_bdry": [x_min, x_max], "Y_bdry": [y_min, y_max],
                "time": cmaq_time_array, "Latitude": lat_ctr, "Longitude": lon_ctr, "Lat_bdry": [lat_min, lat_max],
                "Lon_bdry": [lon_min, lon_max], "XCELL": xcell, "YCELL": ycell}
    return res_dict


def generate_TFLAG(ds, var_len):
    tflag_src = ds["TFLAG"][:]
    tflag_meta = tflag_src[:, 0, :]
    tflag_meta = tflag_meta[:, np.newaxis, :]
    tflag_new = np.repeat(tflag_meta, var_len, axis=1)
    return tflag_new


def plotPolygons(polygon_list, ax, color):
    for current_polygon in polygon_list:
        if current_polygon.geom_type == "MultiPolygon":
            for geom in current_polygon.geoms:
                xs, ys = geom.exterior.xy
                ax.plot(xs, ys, color)
        else:
            xs, ys = current_polygon.exterior.xy
            ax.plot(xs, ys, color)


def WRFGridInfo(ds):
    crs = pyproj.Proj(proj='lcc',  # projection type: Lambert Conformal Conic
                      lat_1=ds.TRUELAT1, lat_2=ds.TRUELAT2,  # Cone intersects with the sphere
                      lat_0=ds.MOAD_CEN_LAT, lon_0=ds.STAND_LON,  # Center point
                      a=6370000, b=6370000)  # This is it! The Earth is a perfect sphere
    # Grid parameters
    dx, dy = ds.DX, ds.DY
    nx, ny = len(ds.dimensions["west_east"]), len(ds.dimensions["south_north"])
    # Down left corner of the domain
    e, n = crs(ds.CEN_LON, ds.CEN_LAT)
    x0 = -(nx - 1) / 2. * dx + e
    y0 = -(ny - 1) / 2. * dy + n
    # 2d grid
    xx, yy = np.meshgrid(np.arange(nx) * dx + x0, np.arange(ny) * dy + y0)
    Xcenters = xx
    Ycenters = yy
    lat_ctr = ds["XLAT"][0, :, :]
    lon_ctr = ds["XLONG"][0, :, :]
    xcell = dx
    ycell = dy
    wrf_time_array = []
    wrf_time = ds["Times"][:]
    for i in range(0, wrf_time.shape[0]):
        current_time_str = ""
        for j in range(0, wrf_time.shape[1]):
            current_time_str = current_time_str + wrf_time[i][j].decode()
        current_time_obj = datetime.strptime(current_time_str, '%Y-%m-%d_%H:%M:%S')
        wrf_time_array.append(current_time_obj)
    res_dict = {"crs": crs, "X": Xcenters, "Y": Ycenters,
                "time": wrf_time_array,
                "Lat": lat_ctr, "Lon": lon_ctr,
                "XCELL": xcell, "YCELL": ycell}
    return res_dict


def plume_grid_cell(vertical_column_weighted_tracer, quantiles):
    if np.sum(vertical_column_weighted_tracer) == 0:
        return 0
    vertical_column_tracer_normalized = vertical_column_weighted_tracer / np.sum(vertical_column_weighted_tracer)
    target_ratio = 1 - quantiles
    layer_num = len(vertical_column_tracer_normalized)
    current_sum = 0
    result_idx = None
    for layer_idx in range(0, layer_num):
        current_sum += vertical_column_tracer_normalized[layer_idx]
        if current_sum >= target_ratio:
            result_idx = layer_idx
            break
    return result_idx


def extract_wrf_time_array(ds):
    wrf_time_array = []
    wrf_time = ds["Times"][:]
    for i in range(0, wrf_time.shape[0]):
        current_time_str = ""
        for j in range(0, wrf_time.shape[1]):
            current_time_str = current_time_str + wrf_time[i][j].decode()
        current_time_obj = datetime.strptime(current_time_str, '%Y-%m-%d_%H:%M:%S')
        wrf_time_array.append(current_time_obj)
    return wrf_time_array


def findSpatialIndex(fire_x, fire_y, X_ctr, Y_ctr):
    """

    :param fire_x: X of fire location in CMAQ projection
    :param fire_y: Y of fire location in CMAQ projection
    :param X_ctr: CMAQ grid X center
    :param Y_ctr: CMAQ grid Y center
    :return: x_idx, y_idx which are the fire location in CMAQ grid
    """
    dist = np.sqrt((X_ctr - fire_x) ** 2 + (Y_ctr - fire_y) ** 2)
    x_idx, y_idx = np.unravel_index(np.argmin(dist, axis=None), dist.shape)
    return x_idx, y_idx


def get_sfire_emissions(wrf_data, select_species):
    # WRF-SFIRE Emission calculation:
    # fgip, &            ! initial fuel load kg/m^2
    # fuel_frac_burnt, &        ! fuel fraction burned this step kg/kg
    # https://github.com/openwfm/WRF-SFIRE/blob/e2bb84f1a1c9f2c0a77986e4b33141d3368f5ef9/phys/module_fr_sfire_atm.F#L724C16-L724C26
    emission_file = "/Volumes/Shield/SFIRECMAQ/data/namelist.fire_emissions"
    dx = wrf_data.getncattr('DX')
    dy = wrf_data.getncattr('DY')
    sub_grid_ratio = int(wrf_data["FXLONG"].shape[1] / (wrf_data["XLONG"].shape[1] + 1))
    grid_area = (dx / sub_grid_ratio) * (dy / sub_grid_ratio)  # m^2
    wrf_time_array = extract_wrf_time_array(wrf_data)
    # Fuel Loading (kg/m^2)
    fgi = np.array([0.166, 0.897, 0.675, 2.468, 0.785, 1.345, 1.092, 1.121, 0.780, 2.694, 2.582, 7.749, 13.024])
    species_mapping = {
        'tr17_1': {'unit': 'g/kg'},
        'tr17_2': {'unit': 'g/kg'},
        'tr17_3': {'unit': 'g/kg'},
        'tr17_4': {'unit': 'g/kg'},
        'tr17_5': {'unit': 'g/kg'},
        'tr17_6': {'unit': 'g/kg'},
        'tr17_7': {'unit': 'g/kg'},
        'tr17_8': {'unit': 'g/kg'},
        'area': {'unit': 'm2'}
    }
    emiss_mapping = {}
    species_name = []
    # Using readlines()
    file1 = open(emission_file, 'r')
    lines = file1.readlines()
    for line in lines:
        if line[0] == '!' or line[0] == '&' or line[0] == '/':
            continue
        elif 'compatible_chem_opt' in line or 'printsums' in line:
            continue
        else:
            line = line.replace("\n", "")
            line = line.replace(' ', '')
            split_line = line.split('=')
            species = split_line[0]
            emiss_factor = split_line[1].split(',')
            emiss_factor_arry = np.zeros(13)
            for i in range(0, len(emiss_factor_arry)):
                emiss_factor_arry[i] = float(emiss_factor[i])
            emiss_mapping[species] = emiss_factor_arry
            species_name.append(species)
    print(species_name)
    print("Remeber to check emission table file: %s" % emission_file)
    # Generate factor matrix (species * category)
    emiss_factor_matrix = np.zeros((len(species_name), 13))
    for i in range(0, len(species_name)):
        emiss_factor_matrix[i, :] = emiss_mapping[species_name[i]]

    # Read NFUEL_CAT, FUEL_FRAC
    fuel_cat = wrf_data["NFUEL_CAT"][:][0, :-sub_grid_ratio, :-sub_grid_ratio]
    fuel_frac = wrf_data["FUEL_FRAC"][:][:, :-sub_grid_ratio, :-sub_grid_ratio]

    # DO NOT APPLY AREA FILTER
    # UNCOMMENT CODE BELOW
    mask_matrix = np.zeros((13, fuel_cat.shape[0], fuel_cat.shape[1]))
    for i in range(0, 13):
        fuel_cat_num = 1 + i
        mask_matrix_tmp = np.zeros(fuel_cat.shape)
        mask_matrix_tmp[fuel_cat == fuel_cat_num] = 1
        mask_matrix[i, :, :] = mask_matrix_tmp

    total_emiss = np.zeros((fuel_frac.shape[0] - 1, len(species_name)))
    total_consume = np.zeros((fuel_frac.shape[0] - 1, 13))
    for t in range(0, fuel_frac.shape[0] - 1):
        fuel_frac_curr = fuel_frac[t, :, :]
        fuel_frac_nxt = fuel_frac[t + 1, :, :]
        fuel_frac_diff = fuel_frac_curr - fuel_frac_nxt
        consume_area = mask_matrix * fuel_frac_diff * grid_area  # m^2
        consume_fuel = consume_area * fgi[:, None, None]
        consume_fuel_total = np.sum(consume_fuel, axis=(1, 2))
        consume_area_total = np.sum(consume_area, axis=(1, 2))
        # emiss_matrix = (species * category)
        emiss_matrix = consume_fuel_total * emiss_factor_matrix
        total_emiss_species = np.sum(emiss_matrix, axis=1)
        total_emiss[t, :] = total_emiss_species
        total_consume[t, :] = consume_area_total

    # Original Emiss (total_emiss, total_consume)
    total_consume = np.reshape(np.sum(total_consume, axis=1), (-1, 1))
    emiss_time = wrf_time_array[0: -1]
    data_matrix = np.hstack((total_emiss, total_consume))
    # convert unit
    species_name.append('area')
    for i in range(0, len(species_name)):
        if species_mapping[species_name[i]]['unit'] == 'g/kg':
            data_matrix[:, i] = data_matrix[:, i]
        if species_mapping[species_name[i]]['unit'] == 'm2':
            data_matrix[:, i] = data_matrix[:, i]
        else:
            # TODO: Is this true?
            data_matrix[:, i] = data_matrix[:, i]

    res = {}
    res["time"] = emiss_time
    for current_specie in select_species:
        specie_idx = species_name.index(current_specie)
        emission_series = data_matrix[:, specie_idx]
        res[current_specie] = emission_series
    return res


def smoke_concentration(tracer_name, d):
    # https://github.com/openwfm/wrfxpy/blob/335669bab63cbca88e827f288c11d875124559cc/src/vis/vis_utils.py#L149
    P = d.variables['P'][:] + d.variables['PB'][:]
    T00 = d.variables['T00'][:]  # Load the full array first
    T00 = T00[:, np.newaxis, np.newaxis, np.newaxis]  # Then reshape
    T = d.variables['T'][:] + T00  # temperature (K)
    r_d = 287  # specific gas constant (J/kg/K)
    rho = P / (r_d * T)  # dry air density  (kg/m^3)
    s = d.variables[tracer_name][:, :, :, :] * rho
    return s


def CMAQGrid3D(metcros3d):
    """

    :param metcros3d: METCROS3D filename
    :return: a dictionary which includes CMAQ grid definitions (X, Y, Z) and projection definition
    """
    met_ds = nc.Dataset(metcros3d)
    cmaq_height = met_ds['ZF'][:]
    ZH = met_ds['ZH'][:]
    time_data = met_ds['TFLAG'][:]
    lat_1 = met_ds.getncattr('P_ALP')
    lat_2 = met_ds.getncattr('P_BET')
    lat_0 = met_ds.getncattr('YCENT')
    lon_0 = met_ds.getncattr('XCENT')
    crs = pyproj.Proj("+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=" + str(lat_1)
                      + " +lat_2=" + str(lat_2) + " +lat_0=" + str(lat_0) +
                      " +lon_0=" + str(lon_0))
    xcell = met_ds.getncattr('XCELL')
    ycell = met_ds.getncattr('YCELL')
    xorig = met_ds.getncattr('XORIG')
    yorig = met_ds.getncattr('YORIG')

    ncols = met_ds.getncattr('NCOLS')
    nrows = met_ds.getncattr('NROWS')

    #> for X, Y cell centers
    x_center_range = np.linspace(xorig+xcell/2, (xorig+xcell/2)+xcell*(ncols-1), ncols)
    y_center_range = np.linspace(yorig+ycell/2, (yorig+ycell/2)+ycell*(nrows-1), nrows)

    Xcenters, Ycenters = np.meshgrid(x_center_range, y_center_range)

    #> for X, Y cell boundaries (i.e., cell corners)
    x_bound_range = np.linspace(xorig, xorig+xcell*ncols, ncols+1)
    y_bound_range = np.linspace(yorig, yorig+ycell*nrows, nrows+1)

    Xbounds, Ybounds = np.meshgrid(x_bound_range, y_bound_range)

    x_max = np.max(Xbounds)
    x_min = np.min(Xbounds)
    y_max = np.max(Ybounds)
    y_min = np.min(Ybounds)
    lon_ctr, lat_ctr = crs(Xcenters, Ycenters, inverse=True)
    cmaq_time_array = []
    for i in range(0, time_data.shape[0]):
        time_data_tmp = time_data[i, 0, :]
        time_str = str(time_data_tmp[0]) + str(time_data_tmp[1]).rjust(6, '0')
        parsed = datetime.strptime(time_str, '%Y%j%H%M%S')
        cmaq_time_array.append(parsed)

    res_dict = {"crs": crs, "X_ctr": Xcenters, "Y_ctr": Ycenters, "X_bdry": [x_min, x_max], "Y_bdry": [y_min, y_max],
                "height": cmaq_height, "time": cmaq_time_array, "Latitude": lat_ctr, "Longitude": lon_ctr, "mid_height": ZH}
    return res_dict


def write_fire_emissions(emission_tensor_dict, metcoros3d_filename, output_filename):
    var_lists = HEADERINFO['VAR-LIST'].split()
    # generate emission file
    with nc.Dataset(metcoros3d_filename) as src, nc.Dataset(output_filename, "w",
                                                            format='NETCDF3_64BIT_OFFSET') as dst:
        # copy global attributes all at once via dictionary
        globle_dict = src.__dict__
        globle_dict['EXEC_ID'] = '3D fire emission process'
        globle_dict['GDNAM'] = 'FIRE_CMAQ'
        globle_dict['UPNAM'] = 'FIRE_EMISS'
        globle_dict['VAR-LIST'] = HEADERINFO['VAR-LIST']
        globle_dict['NVARS'] = len(var_lists)
        globle_dict['FILEDESC'] = '3D emissions output file from Fire Emission Process '
        dst.setncatts(globle_dict)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            if name == 'VAR':
                dst.createDimension(name, len(var_lists))
            else:
                dst.createDimension(name, len(dimension))
        # copy all file data except for the excluded
        for name in var_lists:
            var_tmp = dst.createVariable(name, 'f4', ('TSTEP', 'LAY', 'ROW', 'COL'))
            # REMEMBER DO NOT CHANGE ALL VARIABLES
            if name in emission_tensor_dict:
                dst[name][:] = np.zeros(var_tmp.shape)
                # generate emission data
                # Daily emission data, hour is 25
                print("Mapping " + name + " into Emission File")
                # Generate emiss based on mapping
                dst[name][:] = emission_tensor_dict[name]
            else:
                print("Set the emission variable to zeros: " + name)
                dst[name][:] = np.zeros(var_tmp.shape)
            # # copy variable attributes all at once via dictionary
            dst[name].setncatts(HEADERINFO[name])
        # TFLAG
        var_tmp = dst.createVariable("TFLAG", src.variables["TFLAG"].datatype, src.variables["TFLAG"].dimensions)
        dst["TFLAG"][:] = generate_TFLAG(src, len(var_lists))
        dst["TFLAG"].setncatts(HEADERINFO["TFLAG"])


def discrete_conc_cmap(invertals, base_color_scheme="Spectral_r", low_value_color='white', high_value_color='purple'):
    cmap = plt.get_cmap(base_color_scheme, invertals)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = colors.to_rgba(low_value_color)
    cmap = colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap.set_over(color=high_value_color, alpha=1.0)
    # cmap.set_under(color=low_value_color, alpha=1.0)
    return cmap

def discrete_pos_neg_cmap(invertals, base_color_scheme="seismic", neutral_color='white', low_value_color='black', high_value_color='purple'):
    invertals = (invertals // 2) * 2 + 1
    cmap = plt.get_cmap(base_color_scheme, invertals)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    cmaplist[invertals // 2] = colors.to_rgba(neutral_color)
    cmap.set_over(color=high_value_color, alpha=1.0)
    cmap.set_under(color=low_value_color, alpha=1.0)
    return cmap

def findDomainIndex(coarse_X, coarse_Y, fine_X, fine_Y):
    m, n = fine_X.shape
    corners = [[0, 0], [0, n - 1]]
    y_range_idx = []
    for corner in corners:
        corner_x, corner_y = findSpatialIndex(fine_X[corner[0], corner[1]],
                                              fine_Y[corner[0], corner[1]], coarse_X,
                                              coarse_Y)
        y_range_idx.append(corner_y)
    start_y_idx = min(y_range_idx)
    end_y_idx = max(y_range_idx)
    corners = [[0, 0], [m - 1, 0]]
    x_range_idx = []
    for corner in corners:
        corner_x, corner_y = findSpatialIndex(fine_X[corner[0], corner[1]],
                                              fine_Y[corner[0], corner[1]], coarse_X,
                                              coarse_Y)
        x_range_idx.append(corner_x)
    start_x_idx = min(x_range_idx)
    end_x_idx = max(x_range_idx)

    return start_x_idx, end_x_idx, start_y_idx, end_y_idx


def upsample(conc_field, ratio):
    upsampled = np.repeat(np.repeat(conc_field, ratio, axis=2), ratio, axis=3)
    return upsampled


def WRFZF(ds):
    H = heightAboveGround(ds)
    z_mass_agl = 0.5 * (H[:, :-1, :, :] + H[:, 1:, :, :])
    return z_mass_agl

def discrete_cmap(invertals, base_color_scheme="viridis"):
    cmap = plt.get_cmap(base_color_scheme, invertals)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    # cmap.set_over(color=high_value_color, alpha=1.0)
    return cmap

def conc_at_obs_CMAQ(cmaq_ds, pollutant_name, obs_lon, obs_lat):
    cmaq_info = CMAQGridInfo(cmaq_ds)
    lat = cmaq_info["Latitude"]
    lon = cmaq_info["Longitude"]
    pollutant = cmaq_ds[pollutant_name]
    pollutant_surface = pollutant[:][:, 0, :, :]
    # Nearest grid
    distance = (lat - obs_lat) ** 2 + (lon - obs_lon) ** 2
    x, y = np.where(distance == np.min(distance))
    pollutant_surface_at_obs = pollutant_surface[:, x, y]
    pollutant_surface_at_obs = pollutant_surface_at_obs.flatten()
    return {
        "time": cmaq_info["time"],
        "conc": pollutant_surface_at_obs
    }


def conc_at_obs_SFIRE(wrf_ds, pollutant_name, obs_lon, obs_lat):
    wrf_info = WRFGridInfo(wrf_ds)
    lat = wrf_info["Lat"]
    lon = wrf_info["Lon"]
    # pollutant = wrf_ds[pollutant_name]
    pollutant = smoke_concentration(pollutant_name, wrf_ds)
    pollutant_surface = pollutant[:][:, 0, :, :]
    # Nearest grid
    distance = (lat - obs_lat) ** 2 + (lon - obs_lon) ** 2
    x, y = np.where(distance == np.min(distance))
    pollutant_surface_at_obs = pollutant_surface[:, x, y]
    pollutant_surface_at_obs = pollutant_surface_at_obs.flatten()
    return {
        "time": wrf_info["time"],
        "conc": pollutant_surface_at_obs
    }


def conc_at_obs_WRFChem(wrf_ds, pollutant_name, obs_lon, obs_lat):
    wrf_info = WRFGridInfo(wrf_ds)
    lat = wrf_info["Lat"]
    lon = wrf_info["Lon"]
    # pollutant = wrf_ds[pollutant_name]
    pollutant = wrf_ds[pollutant_name][:]
    pollutant_surface = pollutant[:][:, 0, :, :]
    # Nearest grid
    distance = (lat - obs_lat) ** 2 + (lon - obs_lon) ** 2
    x, y = np.where(distance == np.min(distance))
    pollutant_surface_at_obs = pollutant_surface[:, x, y]
    pollutant_surface_at_obs = pollutant_surface_at_obs.flatten()
    return {
        "time": wrf_info["time"],
        "conc": pollutant_surface_at_obs
    }


def get_conf_intercept(alpha, lr, X, y):
    coefs = np.r_[[lr.intercept_], lr.coef_]
    X_aux = np.zeros((X.shape[0], X.shape[1] + 1))
    X_aux[:, 1:] = X
    X_aux[:, 0] = 1
    dof = -np.diff(X_aux.shape)[0]
    mse = np.sum((y - lr.predict(X)) ** 2) / dof
    var_params = np.diag(np.linalg.inv(X_aux.T.dot(X_aux)))
    t_val = stats.t.isf(alpha / 2, dof)
    gap = t_val * np.sqrt(mse * var_params)
    return {
        "coeffs": coefs,
        "lower": coefs - gap,
        "upper": coefs + gap
    }

def calculate_Area_FRP(ds):
    # calculate heat flux -> FRP
    # How freitas use FRP
    # ! "beta" factor converts FRP to convective energy
    #       heat_fluxW = beta*(FRP/area)/0.55 ! in W/m^2
    # HFX in WRF-SFIRE is total: FIRE_HFX
    # How to calculate FRP and area for each time? FRP = \sum(area_i * HFX_i) -> Watts
    # Area: area = \sum(area_i)
    dx, dy = ds.DX, ds.DY
    cell_area = dx * dy # unit m2
    heat_flux = ds["GRNHFX"][:] + ds["CANHFX"][:]
    counts = np.sum(heat_flux > 0, axis=(1, 2))
    burned_area = counts * cell_area
    frp_value = np.sum(heat_flux * cell_area, axis=(1, 2)) # units: watts
    burned_area = burned_area / 4046.86 # m2 to acres
    return burned_area, frp_value


def verticalFraction(emission_info, mcip_time, met_info):
    plume_info = {"time": [], "bottom": [], "top": []}
    # temporal fraction * height fraction
    # generate hourly fraction
    mcip_height = met_info["height"]
    lay = mcip_height.shape[1]
    emiss_frac = np.zeros((len(mcip_time), lay))
    for time_idx in range(0, len(mcip_time)):
        if mcip_time[time_idx] < emission_info["start_time"] or mcip_time[time_idx] not in emission_info["time"]:
            continue
        cur_idx = emission_info["time"].index(mcip_time[time_idx])
        cur_frp = emission_info["FRP"][cur_idx]
        cur_area = emission_info["area"][cur_idx]
        # calculate plumerise, frp use watts, area use acres
        _, min_plumerise, max_plumerise, VMD = plumeRiseFreitas.plumerisefreitas(cur_area, met_info["up"][time_idx, :],
                                                                               met_info["vp"][time_idx, :],
                                                                               met_info["ta"][time_idx, :],
                                                                               met_info["pres"][time_idx, :],
                                                                               met_info["rv"][time_idx, :],
                                                                               met_info["zt_rams"][time_idx, :],
                                                                               met_info["zm_rams"][time_idx, :], cur_frp)
        # No injection heights, put all emissions as smoldering
        if np.sum(VMD[:, 0]) < 1e-12:
            if max_plumerise > 1.5 * mcip_height[time_idx, -1] or max_plumerise > 100.0 * lay:
                print(f"WARN {mcip_time[time_idx]}: plume top {max_plumerise:.0f} m "
                      f"but VMD empty -> truncation, not smoldering")
            VMD[:, 0] = 0.0
            VMD[0, 0] = 1.0
        plume_info["time"].append(mcip_time[time_idx])
        plume_info["top"].append(max_plumerise)
        plume_info["bottom"].append(min_plumerise)
        emiss_frac[time_idx, :] = VMD[:, 0]
    plume_info["frac"] = emiss_frac
    return emiss_frac, plume_info


def dailyVerticalHourlyEmission(metcros3d_file, metdot3d_file, lon, lat, emission_info):
    """
    :param metcros3d_file: MCIP 3D Grid Information
    :param finn_df: Processed FINN dataframe
    :return: a emission tensor (species, MCIP time length, LAY, MCIP X length, MCIP Y length)
    """
    gas_species = {
        'CO', 'SO2', 'NH3', 'NO2', 'NO', 'CH4', 'ALD2', 'FORM', 'ACET', 'ALDX',
        'BENZ', 'ETH', 'ETHA', 'ETHY', 'ETOH', 'IOLE', 'ISOP', 'KET', 'MEOH',
        'OLE', 'PAR', 'PRPA', 'TERP', 'TOL', 'XYLMN'
    }
    metcro3d = nc.Dataset(metcros3d_file)
    metdot3d = nc.Dataset(metdot3d_file)
    cmaq_3d_info = CMAQGrid3D(metcros3d_file)
    # CMAQ grid definition
    crs = cmaq_3d_info["crs"]
    Xcenters, Ycenters = cmaq_3d_info["X_ctr"], cmaq_3d_info["Y_ctr"]
    mcip_time = cmaq_3d_info["time"]
    mcip_height = cmaq_3d_info["height"]
    LAY = mcip_height.shape[1]

    res = {}
    for species in emission_info:
        if species in {"time", "FRP", "area", "start_time"}:
            continue

        # species, time, layer, x, y
        emission_tensor = np.zeros((len(mcip_time), LAY, Xcenters.shape[0], Xcenters.shape[1]))
        fire_x, fire_y = crs(lon, lat)
        x_idx, y_idx = findSpatialIndex(fire_x, fire_y, Xcenters, Ycenters)
        # get met info
        up = metdot3d["UWIND"][:, :, x_idx, y_idx]
        vp = metdot3d["VWIND"][:, :, x_idx, y_idx]
        ta = metcro3d["TA"][:, :, x_idx, y_idx]
        pres = metcro3d["PRES"][:, :, x_idx, y_idx]
        rv = metcro3d["QV"][:, :, x_idx, y_idx]
        zt_rams = metcro3d["ZH"][:, :, x_idx, y_idx]
        zm_rams = metcro3d["ZF"][:, :, x_idx, y_idx]
        mcip_loc_height = mcip_height[:, :, x_idx, y_idx]
        met_info = {
            "up": up,
            "vp": vp,
            "ta": ta,
            "pres": pres,
            "rv": rv,
            "zt_rams": zt_rams,
            "zm_rams": zm_rams,
            "height": mcip_loc_height
        }
        # get fraction of each phase, fraction definition : temporal fraction * height fraction
        emission_frac, plume_info = verticalFraction(emission_info, mcip_time, met_info)
        for idx in range(0, len(mcip_time)):
            if mcip_time[idx] < emission_info["start_time"] or mcip_time[idx] not in emission_info["time"]:
                continue
            cur_idx = emission_info["time"].index(mcip_time[idx])
            cur_emiss = emission_info[species][cur_idx]
            emission_tensor[idx, :, x_idx, y_idx] = cur_emiss * emission_frac[idx, :]
        # convert unit here: from gram to g/s or mole/s
        if species in gas_species:
            mw = get_molecule_weight(species)
            res[species] = emission_tensor / (mw * 60 * 20)
        else:
            res[species] = emission_tensor / (60 * 20)
    return res, plume_info