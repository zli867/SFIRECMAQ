import numpy as np
from CouplingHelper import coarse_emission_tracer, coupling_surface_buffer, vanilla_emission_tensor, get_emission_prob
from MECHS.Mechanisms import SFIRE_to_CMAQ, tracer_to_species


def calculate_vanilla_emission_tensor(smoke_tracer, tracer_names, wrf_fire_data, ignition_lon, ignition_lat,
                                      coupling_distance, cmaq_info):
    flux_surface_dict = coupling_surface_buffer(wrf_fire_data, ignition_lon, ignition_lat, coupling_distance)
    coarse_emission_tensor_prob = None
    total_mass = {}
    for tracer_name in tracer_names:
        sfire_emission_tensor = vanilla_emission_tensor(tracer_name, wrf_fire_data, flux_surface_dict)
        coarse_emission_tensor = coarse_emission_tracer(wrf_fire_data, cmaq_info, sfire_emission_tensor)
        total_mass[tracer_name] = np.sum(coarse_emission_tensor, axis=(1, 2, 3))
        if tracer_name == smoke_tracer:
            coarse_emission_tensor_prob = get_emission_prob(coarse_emission_tensor)
    # mapping to final emission file
    CMAQ_emissions = {}
    for tracer_name in tracer_names:
        cur_specie_name = tracer_to_species(tracer_name)
        cur_mass = total_mass[tracer_name]
        cur_mapping = SFIRE_to_CMAQ(cur_specie_name)
        for key in cur_mapping:
            CMAQ_emissions[key] = coarse_emission_tensor_prob * cur_mass[:, np.newaxis, np.newaxis, np.newaxis] * cur_mapping[key][0]
    return CMAQ_emissions
