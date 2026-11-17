from utils import findDomainIndex, CMAQGridInfo, upsample
import netCDF4
from MECHS.Mechanisms import CMAQ_to_F0AM
import numpy as np


def calculate_entrain_mass(bkg_CMAQ_file, wrf_info, smoke_volume):
    # return unit is ppb * m^3
    CMAQ_ds = netCDF4.Dataset(bkg_CMAQ_file)
    cmaq_info = CMAQGridInfo(CMAQ_ds)
    sfire_dx = wrf_info["X"][0, 1] - wrf_info["X"][0, 0]
    cmaq_dx = cmaq_info["X_ctr"][0, 1] - cmaq_info["X_ctr"][0, 0]
    grid_ratio = int(cmaq_dx / sfire_dx)
    x_s, x_e, y_s, y_e = findDomainIndex(cmaq_info["X_ctr"], cmaq_info["Y_ctr"], wrf_info["X"], wrf_info["Y"])
    diff_v = smoke_volume[1:, :, :, :] - smoke_volume[:-1, :, :, :]
    diff_v[diff_v < 0] = 0

    res = {"time": cmaq_info["time"]}
    for var_name in CMAQ_ds.variables:
        var = CMAQ_ds.variables[var_name]
        cur_species_name = CMAQ_to_F0AM(var_name)
        if cur_species_name is not None:
            cur_conc = var[:, :, x_s: x_e + 1, y_s: y_e + 1] * 1e3
            upsample_conc = upsample(cur_conc, grid_ratio)
            cur_mass = np.sum(upsample_conc * diff_v, axis=(1, 2, 3))
            res[cur_species_name] = cur_mass
    return res