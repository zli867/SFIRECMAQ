import numpy as np
from datetime import datetime
from MECHS.Mechanisms import get_molecule_weight
import pickle


def createF0AMInputs(smoke_met, smoke_volume, emissions, transport_ratio, bkg_conc,
                     no_fire=False):
    """
     Added no_fire flag for the two-box method.
     no_fire=True -> 'conc' contains ONLY the entrainment/background term (EB_conc),
     i.e. the fire emission contribution is not added. The REAL emissions are still
     used to pick the start time and time window, so the fire box and the no-fire box
     end up on IDENTICAL time axes (required for the fire-minus-no-fire differencing).
    """
    # Create Input file for F0AM
    simulation_start_time = None
    for i in range(0, len(emissions["time"])):
        if emissions["PMC"][i] > 0:
            simulation_start_time = emissions["time"][i]
            break

    start_idx = emissions["time"].index(simulation_start_time)
    foam_input = {"time": emissions["time"][start_idx:], "E": {}}
    for key in emissions.keys():
        if key != "time":
            foam_input["E"][key] = emissions[key][start_idx:]

    # get metereological conditions
    delta_t = emissions["time"][1] - emissions["time"][0]
    start_idx = smoke_met["time"].index((foam_input["time"][0] + delta_t))
    end_idx = smoke_met["time"].index((foam_input["time"][-1] + delta_t))
    foam_input["T"] = smoke_met["T"][start_idx: end_idx + 1]
    foam_input["RH"] = smoke_met["RH"][start_idx: end_idx + 1]
    foam_input["P"] = smoke_met["P"][start_idx: end_idx + 1]
    foam_input["V"] = np.zeros(len(foam_input["time"]))
    foam_input["EB"] = {}
    for i, t in enumerate(foam_input["time"]):
        cur_idx = smoke_volume["time"].index(t)
        nxt_idx = cur_idx + 1
        cur_V, nxt_V = np.sum(smoke_volume["V"][cur_idx, :, :, :]), np.sum(smoke_volume["V"][nxt_idx, :, :, :])
        avg_V = (cur_V + nxt_V) / 2
        foam_input["V"][i] = avg_V
        cmaq_idx = bkg_conc["time"].index(t)
        for key in bkg_conc:
            if key == "time":
                continue
            if key not in foam_input["EB"]:
                foam_input["EB"][key] = [bkg_conc[key][cmaq_idx] / avg_V]
            else:
                foam_input["EB"][key].append(bkg_conc[key][cmaq_idx] / avg_V)
    foam_input["ratio"] = []
    for t in foam_input["time"]:
        t_idx = transport_ratio["time"].index(t)
        cur_ratio = transport_ratio["ratio"][t_idx]
        foam_input["ratio"].append(cur_ratio)

    foam_input["conc"] = {}
    foam_input["E_conc"] = {}
    foam_input["EB_conc"] = {}
    for specie in foam_input["EB"].keys():
        cur_conc = foam_input["EB"][specie]
        foam_input["EB_conc"][specie] = foam_input["EB"][specie]
        if specie in foam_input["E"].keys():
            M = get_molecule_weight(specie)
            # unit: ug/m3
            m_conc = foam_input["E"][specie] / foam_input["V"]
            # ppm
            emit_conc = (8.314 * np.array(foam_input["T"])) / (np.array(foam_input["P"]) * M) * m_conc
            # F0AM needs ppb
            emit_conc = 1e3 * emit_conc
            # <<< CHANGED: skip the fire emission term for the no-fire (background) box
            if not no_fire:
                cur_conc = cur_conc + emit_conc
                foam_input["E_conc"][specie] = emit_conc
            else:
                foam_input["E_conc"][specie] = 0.0
        foam_input["conc"][specie] = cur_conc
    # DEBUG LEVEL
    with open('f0am_input_bkg.pickle' if no_fire else 'f0am_input.pickle', 'wb') as handle:
        pickle.dump(foam_input, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Finalize this input data
    select_keys = ["P", "T", "RH", "conc", "time", "ratio", "V"]
    # Pa to milibar
    foam_input["P"] = np.array(foam_input["P"]) * 0.01
    foam_input = {k: foam_input[k] for k in select_keys}
    return foam_input