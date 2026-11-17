import numpy as np
from utils import CMAQGridInfo
import netCDF4


def calculate_photo_rate(photo_diag_file, time):
    # Mapping CB6r2_AllRxns.m to https://github.com/USEPA/CMAQ/blob/main/CCTM/src/MECHS/cb6r3_ae7_aq/mech_cb6r3_ae7_aq.def
    # From integrateJ.m, it seems that F0AM uses J as s-1 unit.
    # J.JNO2        = Jmcm.J4;, NO2_IUPAC10
    # J.JNO3_NO     = Jmcm.J5; NO3NO_06
    # J.JNO3_NO2    = Jmcm.J6; NO3NO2_06
    # J.JHNO3       = Jmcm.J8; HNO3_IUPAC10
    # J.JO3P        = Jmcm.J2; O3_O3P_IUPAC10
    # J.JO1D        = Jmcm.J1; O3_O1D_IUPAC10
    # J.JHONO       = Jmcm.J7; HONO_IUPAC10
    # J.JH2O2       = Jmcm.J3; H2O2_IUPAC10
    # J.JNTR        = Jmcm.J54; NTR_IUPAC10
    # J.JCOOH       = Jmcm.J41; MEPX_IUPAC10
    # J.JCCHO_R     = Jmcm.J13; ALD2_R_IUPAC10
    # J.JC2CHO      = Jmcm.J14; ALDX_R_IUPAC10
    # J.JHCHO_M     = Jmcm.J12; FORM_M_IUPAC10
    # J.JHCHO_R     = Jmcm.J11; FORM_R_IUPAC10
    # J.JMGLY       = Jmcm.J34; MGLY_IUPAC10
    # J.JHPLD       = Jmcm.J20; HPALD
    # J.JGLY        = Jmcm.J31 + Jmcm.J32 + Jmcm.J33; GLY_R_IUPAC10
    # J.JACET       = Jmcm.J21; ACET_IUPAC10
    # J.JMEK        = Jmcm.J22; KET_IUPAC10
    #
    # % NO DIRECT MCM ANALOGUES
    # J.JN2O5       = Jmcm.Jn19 + Jmcm.Jn20; N2O5_IUPAC10
    # J.JHO2NO2     = Jmcm.Jn21 + Jmcm.Jn22; PNA_IUPAC10
    # J.JPAN        = Jmcm.Jn14 + Jmcm.Jn15; PAN_IUPAC10
    # J.JACRO       = Jmcm.Jn11; ISPD
    # J.JGLYD       = Jmcm.Jn9; GLYD_IUPAC10
    # J.JCRON       = Jmcm.Jn12 + Jmcm.Jn13; NTR_IUPAC10
    f0am_list = [
        "JNO2", "JNO3_NO", "JNO3_NO2", "JHNO3", "JO3P",
        "JO1D", "JHONO", "JH2O2", "JNTR", "JCOOH",
        "JCCHO_R", "JC2CHO", "JHCHO_M", "JHCHO_R", "JMGLY",
        "JHPLD", "JGLY", "JACET", "JMEK", "JN2O5",
        "JHO2NO2", "JPAN", "JACRO", "JGLYD", "JCRON"
    ]

    cmaq_list = [
        "NO2_IUPAC10", "NO3NO_06", "NO3NO2_06", "HNO3_IUPAC10", "O3_O3P_IUPAC10",
        "O3_O1D_IUPAC10", "HONO_IUPAC10", "H2O2_IUPAC10", "NTR_IUPAC10", "MEPX_IUPAC10",
        "ALD2_R_IUPAC10", "ALDX_R_IUPAC10", "FORM_M_IUPAC10", "FORM_R_IUPAC10", "MGLY_IUPAC10",
        "HPALD", "GLY_R_IUPAC10", "ACET_IUPAC10", "KET_IUPAC10", "N2O5_IUPAC10",
        "PNA_IUPAC10", "PAN_IUPAC10", "ISPD", "GLYD_IUPAC10", "NTR_IUPAC10"
    ]
    ds = netCDF4.Dataset(photo_diag_file)
    cmaq_info = CMAQGridInfo(ds)
    idxs = [cmaq_info["time"].index(t) for t in time]
    res = {}
    for i in range(0, len(f0am_list)):
        cur_cmaq_var = cmaq_list[i]
        cur_name = f0am_list[i]
        cur_val = np.mean(ds[cur_cmaq_var][:], axis=(1, 2, 3)) / 60.0
        res[cur_name] = cur_val[idxs]
    return res