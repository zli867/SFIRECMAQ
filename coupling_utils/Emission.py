from MECHS.Mechanisms import SFIRE_to_F0AM, tracer_to_species
from utils import get_sfire_emissions


def calculate_mapped_emission(wrf_fire_data, tracer_names):
    emissions_dict = get_sfire_emissions(wrf_fire_data, tracer_names)
    emissions = {"time": emissions_dict["time"]}
    for tracer_name in tracer_names:
        specie_name = tracer_to_species(tracer_name)
        a_emissions = emissions_dict[tracer_name] * 1e6  # g to ug
        # map to F0AM species
        cur_mapping = SFIRE_to_F0AM(specie_name)
        for foam_species in cur_mapping.keys():
            # Currently, PNA in CMAQ emission file means sodium. However, it means peroxynitric acid in F0AM.
            if foam_species == "PNA":
                emissions[foam_species] = 0 * a_emissions
            else:
                cur_ratio = cur_mapping[foam_species][0]
                cur_emission = cur_ratio * a_emissions
                emissions[foam_species] = cur_emission
    return emissions
