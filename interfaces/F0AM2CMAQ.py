import numpy as np
import scipy.io
from datetime import datetime
from MECHS.Mechanisms import SFIRE_to_F0AM
from MECHS.Mechanisms import F0AM_to_CMAQ
from CouplingHelper import get_emission_prob


def convert_mat2dict(filename):
    mat = scipy.io.loadmat(filename)
    conc_vals = mat['s'][0][0]
    var_names = conc_vals.dtype.names
    res = {}
    for i in range(0, len(var_names)):
        res[var_names[i]] = conc_vals[i].flatten()
    res["time"] = [datetime.utcfromtimestamp(ele) for ele in res["time"]]
    return res


def calculate_chem_emission_tensor(f0am_input, f0am_output_file, f0am_output_bkg_file,
                                   coarse_emission_tensor, cmaq_info):
    """
    TWO-BOX version.
    The injected flux is the SIGNED SMOKE IMPACT = (fire box) - (no-fire box):
        inject(L,r,c) = prob(L,r,c) * ( M_fire - M_bkg ) * conv

    M_fire / M_bkg are the transported-out masses from two identically driven boxes
    (same entrainment, met, photolysis, transport); they differ ONLY in fire emissions.

    Why this replaces the ACONC background subtraction:
      - the no-fire box IS the reference state, evolved through the SAME nonlinear
        chemistry, so the difference is the true fire perturbation (not a linear
        approximation);
      - no need to choose "which background to subtract";
      - entrainment errors largely cancel, since both boxes entrain identically.

    Result is signed: negative near the source (fire NOx titrating O3), positive
    downwind. Downstream requirements:
      - write_fire_emissions must NOT clip negatives
      - CMAQ: EMISS_NEG_CHECK relaxed + MAX(CONC,0) after the vdiff solve
    """
    foam_res = convert_mat2dict(f0am_output_file)        # fire box
    bkg_res = convert_mat2dict(f0am_output_bkg_file)     # <<< NEW: no-fire box
    prob = get_emission_prob(coarse_emission_tensor)     # (T, LAY_e, ROW, COL)
    emiss_times = list(f0am_input["time"])
    n = len(emiss_times)
    dt_sec = int((emiss_times[1] - emiss_times[0]).total_seconds())

    # per-step conversion factor (aligned to f0am_input time)
    conv = np.zeros(n)  # ppb*m^3 -> mol/s
    ti_list = []        # index of each F0AM step on the CMAQ time axis
    for i in range(n):
        cur_P = f0am_input["P"][i] * 100  # mb -> Pa
        cur_T = f0am_input["T"][i]
        cur_ratio = f0am_input["ratio"][i]
        conv[i] = cur_ratio * cur_P / (8314.0 * cur_T) * 1e-6 / dt_sec
        ti_list.append(cmaq_info["time"].index(emiss_times[i]))

    CMAQ_emissions = {}
    missing = []
    for species in foam_res:
        if species == "time":
            continue
        cmaq_name = F0AM_to_CMAQ(species)

        cur_tensor = np.zeros(coarse_emission_tensor.shape)
        for i in range(n):
            cur_time = emiss_times[i]
            t_idx = foam_res["time"].index(cur_time)
            ti = ti_list[i]

            m_fire = foam_res[species][t_idx + 1]
            if species in bkg_res:
                tb_idx = bkg_res["time"].index(cur_time)
                m_bkg = bkg_res[species][tb_idx + 1]
            else:
                m_bkg = 0.0
                if species not in missing:
                    missing.append(species)
            nxt_mass = m_fire - m_bkg
            total_emiss = nxt_mass * conv[i]
            cur_tensor[ti] = prob[ti] * total_emiss
        CMAQ_emissions[cmaq_name] = cur_tensor

    if missing:
        print("WARNING: species missing from the no-fire box (treated as 0):", missing)
    print("TWO-BOX impact computed for %d species" % (len(foam_res) - 1))

    # PM species
    pm_mapping = SFIRE_to_F0AM("p25")
    for key in pm_mapping:
        CMAQ_emissions[key] = coarse_emission_tensor * pm_mapping[key][0]
    # TODO: Add NH3
    return CMAQ_emissions
