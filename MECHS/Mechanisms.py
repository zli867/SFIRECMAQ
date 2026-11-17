def tracer_to_species(tracer_name):
    emission_mapping = {
        "tr17_1": "co",
        "tr17_2": "p25",
        "tr17_3": "no",
        "tr17_4": "no2",
        "tr17_5": "so2",
        "tr17_6": "nh3",
        "tr17_7": "ch4",
        "tr17_8": "NMOC"
    }
    return emission_mapping[tracer_name]


def SFIRE_to_F0AM(sfire_species):
    emission_mapping = {
        'CO': (1, 'co'),
        'SO2': (1, 'so2'),
        'NO2': (1, 'no2'),
        'NO': (1, 'no'),
        'CH4': (1, 'ch4'),
        'ALD2': (0.026318, 'NMOC'),
        'FORM': (0.115995, 'NMOC'),
        'ACET': (0.029243, 'NMOC'),
        'ALDX': (0.081322, 'NMOC'),
        'BENZ': (0.025761, 'NMOC'),
        'ETH': (0.059042, 'NMOC'),
        'ETHA': (0.03161, 'NMOC'),
        'ETHY': (0.014204, 'NMOC'),
        'ETOH': (0.000949, 'NMOC'),
        'IOLE': (0.012975, 'NMOC'),
        'ISOP': (0.00527, 'NMOC'),
        'KET': (0.007112, 'NMOC'),
        'MEOH': (0.105412, 'NMOC'),
        'OLE': (0.092184, 'NMOC'),
        'PAR': (0.200242, 'NMOC'),
        'PRPA': (0.013175, 'NMOC'),
        'TERP': (0.011595, 'NMOC'),
        'TOL': (0.03537, 'NMOC'),
        'XYLMN': (0.007905, 'NMOC'),
        'PAL': (4.6 * (10 ** (-4)), 'p25'),
        'PCA': (7.2 * (10 ** (-4)), 'p25'),
        'PCL': (0.00239, 'p25'),
        'PEC': (0.1093, 'p25'),
        'PFE': (4.45 * (10 ** (-4)), 'p25'),
        'PK': (0.00135, 'p25'),
        'PMN': (1.1 * (10 ** (-4)), 'p25'),
        'PMOTHR': (0.0125, 'p25'),
        'PNA': (0.00135, 'p25'),
        'PNCOM': (0.3513, 'p25'),
        'PNH4': (0.00341, 'p25'),
        'PNO3': (0.0107, 'p25'),
        'POC': (0.5019, 'p25'),
        'PSI': (1.0 * (10 ** (-4)), 'p25'),
        'PSO4': (0.0033, 'p25'),
        'PTI': (6.7 * (10 ** (-4)), 'p25'),
        # SERA
        'PMC': (0.2480483, 'p25')
    }
    res = {}
    for foam_species in emission_mapping.keys():
        if emission_mapping[foam_species][1] == sfire_species:
            res[foam_species] = emission_mapping[foam_species]
    return res


def CMAQ_to_F0AM(CMAQ_species):
    species_dict = {
        'O': 'O3P', 'BENZENE': 'BENZ', 'XYLMN': 'XYL',
        'O1D': 'O1D', 'OH': 'OH', 'HO2': 'HO2', 'C2O3': 'C2O3', 'XO2': 'XO2',
        'XO2N': 'XO2N', 'CXO3': 'CXO3', 'MEO2': 'MEO2', 'TO2': 'TO2', 'ROR': 'ROR',
        'HCO3': 'HCO3', 'CRO': 'CRO', 'BZO2': 'BZO2', 'EPX2': 'EPX2', 'ISO2': 'ISO2',
        'OPO3': 'OPO3', 'RO2': 'RO2', 'XLO2': 'XLO2', 'XO2H': 'XO2H', 'AACD': 'AACD',
        'ACET': 'ACET', 'ALD2': 'ALD2', 'ALDX': 'ALDX', 'CAT1': 'CAT1', 'CO': 'CO',
        'CRES': 'CRES', 'CRON': 'CRON', 'EPOX': 'EPOX', 'ETH': 'ETH', 'ETHA': 'ETHA',
        'ETHY': 'ETHY', 'ETOH': 'ETOH', 'FACD': 'FACD', 'FORM': 'FORM', 'GLY': 'GLY',
        'GLYD': 'GLYD', 'H2O2': 'H2O2', 'HNO3': 'HNO3', 'HONO': 'HONO', 'HPLD': 'HPLD',
        'INTR': 'INTR', 'IOLE': 'IOLE', 'ISOP': 'ISOP', 'ISPD': 'ISPD', 'ISPX': 'ISPX',
        'KET': 'KET', 'MEOH': 'MEOH', 'MEPX': 'MEPX', 'MGLY': 'MGLY', 'N2O5': 'N2O5',
        'NO': 'NO', 'NO2': 'NO2', 'NO3': 'NO3', 'NTR1': 'NTR1', 'NTR2': 'NTR2',
        'O3': 'O3', 'OLE': 'OLE', 'OPAN': 'OPAN', 'OPEN': 'OPEN', 'PACD': 'PACD',
        'PAN': 'PAN', 'PANX': 'PANX', 'PAR': 'PAR', 'PNA': 'PNA', 'PRPA': 'PRPA',
        'ROOH': 'ROOH', 'SO2': 'SO2', 'SULF': 'SULF', 'TERP': 'TERP', 'TOL': 'TOL',
        'XOPN': 'XOPN', 'ECH4': 'ECH4', 'H2': 'H2', 'CH4': 'CH4'
    }
    if CMAQ_species in species_dict.keys():
        return species_dict[CMAQ_species]
    else:
        return None


def F0AM_to_CMAQ(F0AM_species):
    # species_dict = {'O3P': 'O', 'BENZ': 'BENZ', 'XYL': 'XYLMN', 'O1D': 'O1D', 'OH': 'OH', 'HO2': 'HO2',
    #                 'C2O3': 'C2O3', 'XO2': 'XO2', 'XO2N': 'XO2N', 'CXO3': 'CXO3', 'MEO2': 'MEO2', 'TO2': 'TO2',
    #                 'ROR': 'ROR', 'HCO3': 'HCO3', 'CRO': 'CRO', 'BZO2': 'BZO2', 'EPX2': 'EPX2', 'ISO2': 'ISO2',
    #                 'OPO3': 'OPO3', 'RO2': 'RO2', 'XLO2': 'XLO2', 'XO2H': 'XO2H', 'AACD': 'AACD', 'ACET': 'ACET',
    #                 'ALD2': 'ALD2', 'ALDX': 'ALDX', 'CAT1': 'CAT1', 'CO': 'CO', 'CRES': 'CRES', 'CRON': 'CRON',
    #                 'EPOX': 'EPOX', 'ETH': 'ETH', 'ETHA': 'ETHA', 'ETHY': 'ETHY', 'ETOH': 'ETOH', 'FACD': 'FACD',
    #                 'FORM': 'FORM', 'GLY': 'GLY', 'GLYD': 'GLYD', 'H2O2': 'H2O2', 'HNO3': 'HNO3', 'HONO': 'HONO',
    #                 'HPLD': 'HPLD', 'INTR': 'INTR', 'IOLE': 'IOLE', 'ISOP': 'ISOP', 'ISPD': 'ISPD', 'ISPX': 'ISPX',
    #                 'KET': 'KET', 'MEOH': 'MEOH', 'MEPX': 'MEPX', 'MGLY': 'MGLY', 'N2O5': 'N2O5', 'NO': 'NO',
    #                 'NO2': 'NO2', 'NO3': 'NO3', 'NTR1': 'NTR1', 'NTR2': 'NTR2', 'O3': 'O3', 'OLE': 'OLE',
    #                 'OPAN': 'OPAN', 'OPEN': 'OPEN', 'PACD': 'PACD', 'PAN': 'PAN', 'PANX': 'PANX', 'PAR': 'PAR',
    #                 'PNA': 'PNA', 'PRPA': 'PRPA', 'ROOH': 'ROOH', 'SO2': 'SO2', 'SULF': 'SULF', 'TERP': 'TERP',
    #                 'TOL': 'TOL', 'XOPN': 'XOPN', 'ECH4': 'ECH4', 'H2': 'H2', 'CH4': 'CH4'}
    species_dict = {'O3P': 'O', 'BENZ': 'BENZ', 'XYL': 'XYLMN', 'O1D': 'O1D', 'OH': 'OH', 'HO2': 'HO2',
                    'C2O3': 'C2O3', 'XO2': 'XO2', 'XO2N': 'XO2N', 'CXO3': 'CXO3', 'MEO2': 'MEO2', 'TO2': 'TO2',
                    'ROR': 'ROR', 'HCO3': 'HCO3', 'CRO': 'CRO', 'BZO2': 'BZO2', 'EPX2': 'EPX2', 'ISO2': 'ISO2',
                    'OPO3': 'OPO3', 'RO2': 'RO2', 'XLO2': 'XLO2', 'XO2H': 'XO2H', 'AACD': 'AACD', 'ACET': 'ACET',
                    'ALD2': 'ALD2', 'ALDX': 'ALDX', 'CAT1': 'CAT1', 'CO': 'CO', 'CRES': 'CRES', 'CRON': 'CRON',
                    'EPOX': 'EPOX', 'ETH': 'ETH', 'ETHA': 'ETHA', 'ETHY': 'ETHY', 'ETOH': 'ETOH', 'FACD': 'FACD',
                    'FORM': 'FORM', 'GLY': 'GLY', 'GLYD': 'GLYD', 'H2O2': 'H2O2', 'HNO3': 'HNO3', 'HONO': 'HONO',
                    'HPLD': 'HPLD', 'INTR': 'INTR', 'IOLE': 'IOLE', 'ISOP': 'ISOP', 'ISPD': 'ISPD', 'ISPX': 'ISPX',
                    'KET': 'KET', 'MEOH': 'MEOH', 'MEPX': 'MEPX', 'MGLY': 'MGLY', 'N2O5': 'N2O5', 'NO': 'NO',
                    'NO2': 'NO2', 'NO3': 'NO3', 'NTR1': 'NTR1', 'NTR2': 'NTR2', 'O3': 'O3', 'OLE': 'OLE',
                    'OPAN': 'OPAN', 'OPEN': 'OPEN', 'PACD': 'PACD', 'PAN': 'PAN', 'PANX': 'PANX', 'PAR': 'PAR',
                    'PNA': 'HNO4', 'PRPA': 'PRPA', 'ROOH': 'ROOH', 'SO2': 'SO2', 'SULF': 'SULF', 'TERP': 'TERP',
                    'TOL': 'TOL', 'XOPN': 'XOPN', 'ECH4': 'ECH4', 'H2': 'H2', 'CH4': 'CH4'}
    if F0AM_species in species_dict.keys():
        return species_dict[F0AM_species]
    else:
        return None


def get_molecule_weight(specie_name):
    molecule_weights = {
        'CO': 28,
        'SO2': 64,
        'NH3': 17,
        'NO2': 46,
        'NO': 30,
        'CH4': 16.0425,
        'ALD2': 44.0526,
        'FORM': 30.026,
        'ACET': 58.0791,
        'ALDX': 40.7486,
        'BENZ': 78.1118,
        'ETH': 28.0532,
        'ETHA': 30.069,
        'ETHY': 26.0373,
        'ETOH': 46.0684,
        'IOLE': 55.6133,
        'ISOP': 68.117,
        'KET': 17.8646,
        'MEOH': 32.0419,
        'OLE': 31.8777,
        'PAR': 21.0484,
        'PRPA': 44.0956,
        'TERP': 136.234,
        'TOL': 93.0188,
        'XYLMN': 106.165,
        'PNA': 79.0
    }
    return molecule_weights[specie_name]


def convert_to_mole_s(molecule_weight):
    # convert to g/s -> moles/s
    return 1 / molecule_weight


def SFIRE_to_CMAQ(sfire_species):
    emission_mapping = {
        'CO': (1 * convert_to_mole_s(28), 'co'),
        # 'NH3': (1 * convert_to_mole_s(17), 'nh3'),
        'SO2': (1 * convert_to_mole_s(64), 'so2'),
        'NO2': (1 * convert_to_mole_s(46), 'no2'),
        'NO': (1 * convert_to_mole_s(30), 'no'),
        'CH4': (1 * convert_to_mole_s(16.0425), 'ch4'),
        'ALD2': (0.026318 * convert_to_mole_s(44.0526), 'NMOC'),
        'FORM': (0.115995 * convert_to_mole_s(30.026), 'NMOC'),
        'ACET': (0.029243 * convert_to_mole_s(58.0791), 'NMOC'),
        'ALDX': (0.081322 * convert_to_mole_s(40.7486), 'NMOC'),
        'BENZ': (0.025761 * convert_to_mole_s(78.1118), 'NMOC'),
        'ETH': (0.059042 * convert_to_mole_s(28.0532), 'NMOC'),
        'ETHA': (0.03161 * convert_to_mole_s(30.069), 'NMOC'),
        'ETHY': (0.014204 * convert_to_mole_s(26.0373), 'NMOC'),
        'ETOH': (0.000949 * convert_to_mole_s(46.0684), 'NMOC'),
        'IOLE': (0.012975 * convert_to_mole_s(55.6133), 'NMOC'),
        'ISOP': (0.00527 * convert_to_mole_s(68.117), 'NMOC'),
        'KET': (0.007112 * convert_to_mole_s(17.8646), 'NMOC'),
        'MEOH': (0.105412 * convert_to_mole_s(32.0419), 'NMOC'),
        'OLE': (0.092184 * convert_to_mole_s(31.8777), 'NMOC'),
        'PAR': (0.200242 * convert_to_mole_s(21.0484), 'NMOC'),
        'PRPA': (0.013175 * convert_to_mole_s(44.0956), 'NMOC'),
        'TERP': (0.011595 * convert_to_mole_s(136.234), 'NMOC'),
        'TOL': (0.03537 * convert_to_mole_s(93.0188), 'NMOC'),
        'XYLMN': (0.007905 * convert_to_mole_s(106.165), 'NMOC'),
        'PAL': (4.6 * (10 ** (-4)), 'p25'),
        'PCA': (7.2 * (10 ** (-4)), 'p25'),
        'PCL': (0.00239, 'p25'),
        'PEC': (0.1093, 'p25'),
        'PFE': (4.45 * (10 ** (-4)), 'p25'),
        'PK': (0.00135, 'p25'),
        'PMN': (1.1 * (10 ** (-4)), 'p25'),
        'PMOTHR': (0.0125, 'p25'),
        'PNA': (0.00135, 'p25'),
        'PNCOM': (0.3513, 'p25'),
        'PNH4': (0.00341, 'p25'),
        'PNO3': (0.0107, 'p25'),
        'POC': (0.5019, 'p25'),
        'PSI': (1.0 * (10 ** (-4)), 'p25'),
        'PSO4': (0.0033, 'p25'),
        'PTI': (6.7 * (10 ** (-4)), 'p25'),
        # SERA
        'PMC': (0.2480483, 'p25')
    }
    res = {}
    for foam_species in emission_mapping.keys():
        if emission_mapping[foam_species][1] == sfire_species:
            res[foam_species] = emission_mapping[foam_species]
    return res