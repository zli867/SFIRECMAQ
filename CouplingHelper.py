import netCDF4
from utils import heightAboveGround, grid_cell_plume_height, convertByteToStr, plume_grid_cell
import numpy as np
from kneed import KneeLocator
from utils import smoke_concentration
from MECHS.Mechanisms import CMAQ_to_F0AM
from utils import CMAQGridInfo, WRFGridInfo
import os

def loss_value(current_distance, current_height, knee_distance):
    current_distance_np = np.array(current_distance)
    current_height_np = np.array(current_height)
    if np.max(current_distance_np) < knee_distance:
        return 0
    else:
        max_height = np.max(current_height_np)
        # why I can always find a value in the array? Because I use bin to calclulate the height,
        # each bin I have a value
        knee_height = current_height_np[current_distance_np == knee_distance][0]
        return (knee_height - max_height) ** 2


def calculate_coupling_distance(tracer_name, wrf_fire_data, ignition_lon, ignition_lat, bin_step=1, quantile_vertical=0.1,
                                quantile_horizontal=0.1):
    # variables extraction
    layer_height = heightAboveGround(wrf_fire_data)
    layer_mid_height = 0.5 * (layer_height[:, 1:, :, :] + layer_height[:, :-1, :, :])
    wrf_tracer = smoke_concentration(tracer_name, wrf_fire_data)
    thickness = layer_height[:, 1:, :, :] - layer_height[:, :-1, :, :]
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    area = dx * dy
    wrf_tracer_mass = wrf_tracer * thickness * area
    xlat = wrf_fire_data["XLAT"][:][0, :, :]
    xlon = wrf_fire_data["XLONG"][:][0, :, :]
    time_str = wrf_fire_data["Times"][:]
    time_frames, layers, m, n = wrf_tracer.shape
    start_time_idx = 0
    end_time_idx = time_frames - 1
    for t_step in range(0, time_frames):
        current_wrf_tracer_mass = wrf_tracer_mass[t_step, :, :, :]
        if np.max(current_wrf_tracer_mass) > 0:
            start_time_idx = t_step
            break

    # store distance of smoke centerline to start point
    distance_res = []
    time_str_res = []
    # store height of smoke height of centerline
    height_res = []
    for t_step in range(start_time_idx, end_time_idx + 1):
        current_wrf_tracer_mass = wrf_tracer_mass[t_step, :, :, :]
        current_thickness = thickness[t_step, :, :, :]
        current_plume_height = np.zeros((m, n))
        current_total_tracer_mass = np.sum(current_wrf_tracer_mass)

        vertical_column_tracer_mass = np.sum(current_wrf_tracer_mass, axis=0)
        tracer_path = np.zeros(vertical_column_tracer_mass.shape)
        tracer_list = []

        for column in range(0, m):
            for row in range(0, n):
                tracer_list.append((vertical_column_tracer_mass[column, row], column, row))
                current_cell_height = grid_cell_plume_height(current_wrf_tracer_mass[:, column, row],
                                                             current_thickness[:, column, row],
                                                             quantile_vertical)
                current_plume_height[column, row] = current_cell_height
        # sort tracer list by mass from largest to smallest
        tracer_list.sort(reverse=True)

        select_total_mass = 0
        target_horizontal_ratio = 1 - quantile_horizontal
        for select_idx in range(0, len(tracer_list)):
            select_mass, max_col, max_row = tracer_list[select_idx]
            select_total_mass += select_mass
            tracer_path[max_col, max_row] = 1
            if select_total_mass >= target_horizontal_ratio * current_total_tracer_mass:
                break

        current_plume_height_selected = current_plume_height * tracer_path

        # Search start grid index
        distance = np.sqrt((xlat - ignition_lat) ** 2 + (xlon - ignition_lon) ** 2)
        idx = np.argwhere(distance == np.min(distance))
        ignition_x, ignition_y = idx[0]

        # distance vs height
        m, n = current_plume_height_selected.shape
        distance = []
        plume_height = []
        for i in range(0, m):
            for j in range(0, n):
                # only include pixels when plume height is greater than 0
                if current_plume_height_selected[i, j] > 0:
                    delta_x = np.abs((i - ignition_x) * (dx / 1000))
                    delta_y = np.abs((j - ignition_y) * (dy / 1000))
                    pixel_distance = np.sqrt(delta_x ** 2 + delta_y ** 2)  # unit km
                    pixel_height = current_plume_height_selected[i, j]
                    distance.append(pixel_distance)
                    plume_height.append(pixel_height)

        distance = np.array(distance)
        plume_height = np.array(plume_height)
        binned_distance = []
        binned_height = []
        start_distance = 0
        current_distance = start_distance
        while current_distance <= np.max(distance):
            lower_bound = current_distance
            upper_bound = current_distance + bin_step
            bin_height = plume_height[(distance >= lower_bound) & (distance < upper_bound)]
            if len(bin_height) == 0:
                mean_bin_height = 0
            else:
                mean_bin_height = np.mean(bin_height)
            mean_distance = (lower_bound + upper_bound) / 2
            binned_distance.append(mean_distance)
            binned_height.append(mean_bin_height)
            current_distance += bin_step
        distance_res.append(binned_distance)
        height_res.append(binned_height)
        current_time_str = convertByteToStr(time_str[t_step])
        time_str_res.append(current_time_str)

    loss_array = []
    distance_array = []
    # calculate start dist and end dist in distance res
    min_bin_dist = None
    max_bin_dist = None
    for i in range(0, len(distance_res)):
        if i == 0:
            min_bin_dist = min(distance_res[i])
            max_bin_dist = max(distance_res[i])
        else:
            current_min = min(distance_res[i])
            current_max = max(distance_res[i])
            min_bin_dist = min([min_bin_dist, current_min])
            max_bin_dist = max([max_bin_dist, current_max])
    distance_values = np.arange(min_bin_dist, max_bin_dist + bin_step, bin_step)
    for distance_value in distance_values:
        total_loss = 0
        for time_idx in range(0, len(time_str_res)):
            current_distance = distance_res[time_idx]
            current_height = height_res[time_idx]
            current_loss = loss_value(current_distance, current_height, distance_value)
            total_loss += current_loss
        loss_array.append(total_loss)
        distance_array.append(distance_value)

    kn = KneeLocator(
        distance_array,
        loss_array,
        curve='convex',
        direction='decreasing',
        interp_method='interp1d',
    )
    global_knee_point = kn.knee
    coupling_distance = global_knee_point
    # I decide to ceil to kilometer to match the WRF-SFIRE grid surface (coupling surface) as same as CMAQ gridsurface
    return np.ceil(coupling_distance)


def coupling_surface_buffer(wrf_fire_data, ignition_lon, ignition_lat, coupling_distance):
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    xlat = wrf_fire_data["XLAT"][:][0, :, :]
    xlon = wrf_fire_data["XLONG"][:][0, :, :]
    # TODO: I would like to not add 1 after buffer_dx, buffer_dy, although it passed the mass conservation
    #  and better for add the tracer mass inside coupling region back, adding back conc field is not realistic
    buffer_dx = int(np.round(coupling_distance * 1000 / dx))
    buffer_dy = int(np.round(coupling_distance * 1000 / dy))

    flux_dict = {(0, 1): [],
                 (1, 0): [],
                 (0, -1): [],
                 (-1, 0): []}
    # Search start grid index
    distance = np.sqrt((xlat - ignition_lat) ** 2 + (xlon - ignition_lon) ** 2)
    idx = np.argwhere(distance == np.min(distance))
    ignition_x, ignition_y = idx[0]

    isx, iex = ignition_x - buffer_dx, ignition_x + buffer_dx
    isy, iey = ignition_y - buffer_dy, ignition_y + buffer_dy
    isx, isy = max(0, isx), max(0, isy)
    x_max, y_max = xlat.shape
    iex, iey = min(iex, x_max - 1), min(iey, y_max - 1)
    # North
    for idx in range(isy, iey + 1):
        flux_dict[(0, 1)].append((iex, idx))
    # South
    for idx in range(isy, iey + 1):
        flux_dict[(0, -1)].append((isx, idx))
    # East
    for idx in range(isx, iex + 1):
        flux_dict[(1, 0)].append((idx, iey))
    # West
    for idx in range(isx, iex + 1):
        flux_dict[(-1, 0)].append((idx, isy))
    return flux_dict


def vanilla_emission_tensor(tracer_name, wrf_fire_data, flux_dict):
    u_wind = wrf_fire_data["U"][:]
    v_wind = wrf_fire_data["V"][:]
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    # variables extraction
    layer_height = heightAboveGround(wrf_fire_data)
    wrf_tracer = smoke_concentration(tracer_name, wrf_fire_data)
    thickness = layer_height[:, 1:, :, :] - layer_height[:, :-1, :, :]
    emission_tracer = np.zeros(wrf_tracer.shape)
    _, _, x_max, y_max = thickness.shape
    # wrf_tracer unit: g/m3
    # wrf_tracer * (u, v) * surface_norm * surface area = ug/s
    for surface_norm in flux_dict.keys():
        if surface_norm[0] == 0:
            # south north surface
            for current_idx in flux_dict[surface_norm]:
                current_x_idx = current_idx[0]
                current_y_idx = current_idx[1]
                # TODO: ADD some error message to remind user to enlarge WRF-SFIRE domain
                current_x_idx = min(x_max - 1, current_x_idx)
                current_y_idx = min(y_max - 1, current_y_idx)
                current_area = thickness[:, :, current_x_idx, current_y_idx] * dx
                v_wind_x_idx_shift = 0
                if surface_norm[1] > 0:
                    v_wind_x_idx_shift = 1
                current_v_wind = v_wind[:, :, v_wind_x_idx_shift + current_x_idx, current_y_idx]
                current_wrf_tracer = wrf_tracer[:, :, current_x_idx, current_y_idx]
                current_flux_value = current_wrf_tracer * surface_norm[1] * current_v_wind * current_area
                # use "+=" since the corner grid cells should be accumulate values
                emission_tracer[:, :, current_x_idx, current_y_idx] += current_flux_value
        else:
            # east west surface
            for current_idx in flux_dict[surface_norm]:
                current_x_idx = current_idx[0]
                current_y_idx = current_idx[1]
                # TODO: ADD some error message to remind user to enlarge WRF-SFIRE domain
                current_x_idx = min(x_max - 1, current_x_idx)
                current_y_idx = min(y_max - 1, current_y_idx)
                current_area = thickness[:, :, current_x_idx, current_y_idx] * dy
                u_wind_y_idx_shift = 0
                if surface_norm[0] < 0:
                    u_wind_y_idx_shift = 1
                current_u_wind = u_wind[:, :, current_x_idx, current_y_idx + u_wind_y_idx_shift]
                current_wrf_tracer = wrf_tracer[:, :, current_x_idx, current_y_idx]
                current_flux_value = current_wrf_tracer * surface_norm[0] * current_u_wind * current_area
                # use "+=" since the corner grid cells should be accumulate values
                emission_tracer[:, :, current_x_idx, current_y_idx] += current_flux_value
    # ug/s to g/s
    emission_tracer = emission_tracer * 1e-6
    return emission_tracer


def coarse_emission_tracer(wrf_fire_data, cmaq_info, emission_tracer):
    cmaq_X = cmaq_info["X_ctr"]
    cmaq_Y = cmaq_info["Y_ctr"]
    xlat = wrf_fire_data["XLAT"][:][0, :, :]
    xlon = wrf_fire_data["XLONG"][:][0, :, :]
    wrf_x, wrf_y = cmaq_info["crs"](xlon, xlat)
    x_size, y_size = wrf_x.shape
    mapping_matrix_X = np.zeros((x_size, y_size))
    mapping_matrix_Y = np.zeros((x_size, y_size))
    nrows_coarse, ncols_coarse = cmaq_X.shape
    time_frames, layers, _, _ = emission_tracer.shape
    for i in range(0, x_size):
        for j in range(0, y_size):
            distance = (cmaq_X - wrf_x[i, j]) ** 2 + (cmaq_Y - wrf_y[i, j]) ** 2
            x_idx, y_idx = np.where(distance == np.min(distance))
            # TODO: Does the vertical height matched???
            mapping_matrix_X[i, j] = x_idx
            mapping_matrix_Y[i, j] = y_idx
    coarse_emiss_tracer = np.zeros((time_frames, layers, nrows_coarse, ncols_coarse))
    for i in range(0, x_size):
        for j in range(0, y_size):
            coarse_x_idx = int(mapping_matrix_X[i, j])
            coarse_y_idx = int(mapping_matrix_Y[i, j])
            coarse_emiss_tracer[:, :, coarse_x_idx, coarse_y_idx] += emission_tracer[:, :, i, j]
    return coarse_emiss_tracer


def calculate_2D_plume_height(tracer_name, wrf_fire_data, t_step, quantile_vertical=0.1, quantile_horizontal=0.1):
    # variables extraction
    layer_height = heightAboveGround(wrf_fire_data)
    wrf_tracer = smoke_concentration(tracer_name, wrf_fire_data)
    thickness = layer_height[:, 1:, :, :] - layer_height[:, :-1, :, :]
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    area = dx * dy
    wrf_tracer_mass = wrf_tracer * thickness * area
    xlat = wrf_fire_data["XLAT"][:][0, :, :]
    xlon = wrf_fire_data["XLONG"][:][0, :, :]
    time_str = wrf_fire_data["Times"][:]
    time_frames, layers, m, n = wrf_tracer.shape
    current_wrf_tracer_mass = wrf_tracer_mass[t_step, :, :, :]
    current_thickness = thickness[t_step, :, :, :]
    current_plume_height = np.zeros((m, n))
    current_total_tracer_mass = np.sum(current_wrf_tracer_mass)

    vertical_column_tracer_mass = np.sum(current_wrf_tracer_mass, axis=0)
    tracer_path = np.zeros(vertical_column_tracer_mass.shape)
    tracer_list = []

    for column in range(0, m):
        for row in range(0, n):
            tracer_list.append((vertical_column_tracer_mass[column, row], column, row))
            current_cell_height = grid_cell_plume_height(current_wrf_tracer_mass[:, column, row],
                                                         current_thickness[:, column, row],
                                                         quantile_vertical)
            current_plume_height[column, row] = current_cell_height
    # sort tracer list by mass from largest to smallest
    tracer_list.sort(reverse=True)

    select_total_mass = 0
    target_horizontal_ratio = 1 - quantile_horizontal
    for select_idx in range(0, len(tracer_list)):
        select_mass, max_col, max_row = tracer_list[select_idx]
        select_total_mass += select_mass
        tracer_path[max_col, max_row] = 1
        if select_total_mass >= target_horizontal_ratio * current_total_tracer_mass:
            break

    current_plume_height_selected = current_plume_height * tracer_path
    current_time_str = convertByteToStr(time_str[t_step])
    return current_time_str, current_plume_height_selected


def calculate_plume_mask(tracer_name, wrf_fire_data, quantile_vertical_high=0.1, quantile_vertical_low=0.9, quantile_horizontal=0.1):
    # variables extraction
    layer_height = heightAboveGround(wrf_fire_data)
    wrf_tracer = smoke_concentration(tracer_name, wrf_fire_data)
    thickness = layer_height[:, 1:, :, :] - layer_height[:, :-1, :, :]
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    area = dx * dy
    wrf_tracer_mass = wrf_tracer * thickness * area
    time_frames, layers, m, n = wrf_tracer.shape
    global_smoke_mask = []
    for t_step in range(0, time_frames):
        current_wrf_tracer_mass = wrf_tracer_mass[t_step, :, :, :]
        current_total_tracer_mass = np.sum(current_wrf_tracer_mass)
        smoke_mask = np.zeros(current_wrf_tracer_mass.shape, dtype=bool)
        # If there is no mass, set the smoke mask as False
        if current_total_tracer_mass == 0:
            global_smoke_mask.append(smoke_mask[np.newaxis, :, :])
            continue

        vertical_column_tracer_mass = np.sum(current_wrf_tracer_mass, axis=0)
        tracer_path = []
        tracer_list = []

        for column in range(0, m):
            for row in range(0, n):
                tracer_list.append((vertical_column_tracer_mass[column, row], column, row))
        # sort tracer list by mass from largest to smallest
        tracer_list.sort(reverse=True)

        select_total_mass = 0
        target_horizontal_ratio = 1 - quantile_horizontal
        for select_idx in range(0, len(tracer_list)):
            select_mass, max_col, max_row = tracer_list[select_idx]
            select_total_mass += select_mass
            tracer_path.append((max_col, max_row))
            if select_total_mass >= target_horizontal_ratio * current_total_tracer_mass:
                break
        for select_idx in range(0, len(tracer_path)):
            max_col, max_row = tracer_path[select_idx]
            upper_idx, lower_idx = plume_grid_cell(current_wrf_tracer_mass[:, max_col, max_row], quantile_vertical_high), plume_grid_cell(current_wrf_tracer_mass[:, max_col, max_row], quantile_vertical_low)
            smoke_mask[lower_idx: upper_idx + 1, max_col, max_row] = True
        global_smoke_mask.append(smoke_mask[np.newaxis, :, :])
    return np.concatenate(global_smoke_mask)


def inside_smoke_mask(smoke_mask, wrf_fire_data, ignition_lon, ignition_lat, coupling_distance):
    # Search start grid index
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    xlat = wrf_fire_data["XLAT"][:][0, :, :]
    xlon = wrf_fire_data["XLONG"][:][0, :, :]
    distance = np.sqrt((xlat - ignition_lat) ** 2 + (xlon - ignition_lon) ** 2)
    idx = np.argwhere(distance == np.min(distance))
    ignition_x, ignition_y = idx[0]
    buffer_dx = int(np.round(coupling_distance * 1000 / dx))
    buffer_dy = int(np.round(coupling_distance * 1000 / dy))
    _, _, x_max, y_max = smoke_mask.shape
    isx, iex = ignition_x - buffer_dx, ignition_x + buffer_dx
    isy, iey = ignition_y - buffer_dy, ignition_y + buffer_dy
    isx, isy = max(0, isx), max(0, isy)
    iex, iey = min(iex, x_max - 1), min(iey, y_max - 1)
    mask = np.zeros(smoke_mask.shape, dtype=bool)
    mask[:, :, isx: iex + 1, isy: iey + 1] = True
    inside_mask = (smoke_mask & mask)
    return inside_mask


def RELHUM(T, P, QV):
    # https://github.com/USEPA/CMAQ/blob/main/CCTM/src/emis/emis/SSEMIS.F
    # https://github.com/NCAR/wrf-python/blob/b439dce6953fec8ee97ef44bb6610edb2f011eef/fortran/wrf_user.f90
    # DCOMPUTERH
    EZERO = 6.112
    ESLCON1 = 17.67
    ESLCON2 = 29.65
    CELKEL = 273.15
    EPS = 0.622
    es = EZERO * np.exp(ESLCON1 * (T - CELKEL)/(T - ESLCON2))
    qvs = EPS * es / (0.01 * P - (1.0 - EPS) * es)
    rh = 100.0 * np.clip(QV / qvs, 0.0, 1.0)
    return rh


def temperature(T, P, PB):
    # reference: https://mailman.ucar.edu/pipermail/wrf-users/2010/001896.html
    # reference: https://github.com/NCAR/wrf-python/blob/b439dce6953fec8ee97ef44bb6610edb2f011eef/fortran/wrf_user.f90#L33
    RD = 287
    CP = 7 * RD / 2
    kappa = RD / CP
    total_potential = T + 300
    total_pressure = pressure(P, PB)
    normal_T = total_potential * ((total_pressure / 100000.0) ** kappa)
    return normal_T


def pressure(P, PB):
    return PB + P


def get_plume_met(wrf_fire_data, smoke_mask):
    T = wrf_fire_data["T"][:]
    QV = wrf_fire_data["QVAPOR"][:]
    P = wrf_fire_data["P"][:]
    PB = wrf_fire_data["PB"][:]
    temp = temperature(T, P, PB)
    press = pressure(P, PB)
    RH = RELHUM(temp, press, QV)
    grid_info = WRFGridInfo(wrf_fire_data)
    res = {
        "time": [], "T": [], "P": [], "RH": []
    }
    for t_idx in range(smoke_mask.shape[0]):
        cur_mask = smoke_mask[t_idx]
        if not np.any(cur_mask):
            continue
        res["time"].append(grid_info["time"][t_idx])
        res["T"].append(np.mean(temp[t_idx][cur_mask]))
        res["P"].append(np.mean(press[t_idx][cur_mask]))
        res["RH"].append(np.mean(RH[t_idx][cur_mask]))
    return res


def inside_mass(tracer_name, wrf_fire_data, ignition_lon, ignition_lat, coupling_distance):
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    xlat = wrf_fire_data["XLAT"][:][0, :, :]
    xlon = wrf_fire_data["XLONG"][:][0, :, :]
    buffer_dx = int(np.round(coupling_distance * 1000 / dx))
    buffer_dy = int(np.round(coupling_distance * 1000 / dy))
    layer_height = heightAboveGround(wrf_fire_data)
    wrf_tracer = smoke_concentration(tracer_name, wrf_fire_data)
    thickness = layer_height[:, 1:, :, :] - layer_height[:, :-1, :, :]
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    area = dx * dy
    wrf_tracer_mass = wrf_tracer * thickness * area
    # Search start grid index
    distance = np.sqrt((xlat - ignition_lat) ** 2 + (xlon - ignition_lon) ** 2)
    idx = np.argwhere(distance == np.min(distance))
    ignition_x, ignition_y = idx[0]
    _, _, x_max, y_max = wrf_tracer_mass.shape
    isx, iex = ignition_x - buffer_dx, ignition_x + buffer_dx
    isy, iey = ignition_y - buffer_dy, ignition_y + buffer_dy
    isx, isy = max(0, isx), max(0, isy)
    iex, iey = min(iex, x_max - 1), min(iey, y_max - 1)
    outside_mass = wrf_tracer_mass.copy()
    outside_mass[:, :, isx: iex + 1, isy: iey + 1] = 0
    in_mass = wrf_tracer_mass - outside_mass
    return in_mass


def mask_volume(wrf_fire_data, grid_mask):
    dx = wrf_fire_data.getncattr('DX')
    dy = wrf_fire_data.getncattr('DY')
    layer_height = heightAboveGround(wrf_fire_data)
    thickness = layer_height[:, 1:, :, :] - layer_height[:, :-1, :, :]
    area = dx * dy
    volumes = area * thickness
    return volumes * grid_mask


def get_emission_prob(coarse_emission_tensor):
    total_emiss = np.sum(coarse_emission_tensor, axis=(1, 2, 3))
    coarse_emission_tensor_prob = np.zeros(coarse_emission_tensor.shape)
    for i in range(0, len(total_emiss)):
        if total_emiss[i] == 0:
            continue
        else:
            coarse_emission_tensor_prob[i, :, :, :] = coarse_emission_tensor[i, :, :, :] / total_emiss[i]
    return coarse_emission_tensor_prob


def run_matlab(script_name, in_path, out_path, matlab_bin, f0am_dir):
    if os.path.exists(out_path):
        os.remove(out_path)                       # no stale reuse
    cmd = (f'{matlab_bin} -batch '
           f'"addpath(genpath(\'{f0am_dir}\')); '
           f'INPUT_MAT=\'{in_path}\'; OUTPUT_MAT=\'{out_path}\'; {script_name}"')
    print("Running MATLAB:", script_name)
    rc = os.system(cmd)
    if not os.path.exists(out_path):              # loud check
        raise RuntimeError(f"{script_name} did not write {out_path} (rc={rc})")
    return rc