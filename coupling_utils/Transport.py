from CouplingHelper import inside_mass, coupling_surface_buffer, vanilla_emission_tensor
from utils import WRFGridInfo
import numpy as np


def calculate_transport_ratio(total_inside_mass, transport_mass, wrf_info):
    t_inside_mass = np.sum(total_inside_mass, axis=(1, 2, 3))
    # from g/s to ug
    transport_mass = (60 * 20) * 1e6 * transport_mass
    t_transport_mass = np.sum(transport_mass, axis=(1, 2, 3))
    c_transport_mass = (t_transport_mass[:-1] + t_transport_mass[1:]) / 2
    # calculate transport ratio, T_tr^t / M_tr^{t + 1}
    transport_ratio = {"time": [], "ratio": []}
    for i in range(1, len(t_inside_mass)):
        if t_inside_mass[i] == 0:
            continue
        transport_ratio["ratio"].append(c_transport_mass[i - 1] / t_inside_mass[i])
        transport_ratio["time"].append(wrf_info["time"][i - 1])
    return transport_ratio
