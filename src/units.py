"""Konversi satuan lapangan (field) <-> SI. Semua rumus internal memakai SI."""

G = 9.80665                  # percepatan gravitasi (m/s2)

PPG_TO_KGM3 = 119.8264       # 1 lb/gal dalam kg/m3
PSI_TO_PA = 6894.757         # 1 psi dalam Pa
FT_TO_M = 0.3048             # 1 ft dalam m
IN_TO_M = 0.0254             # 1 inci dalam m
GPM_TO_M3S = 6.30902e-5      # 1 gal/menit dalam m3/s


def ppg_to_kgm3(ppg):
    return ppg * PPG_TO_KGM3


def kgm3_to_ppg(kgm3):
    return kgm3 / PPG_TO_KGM3


def psi_to_pa(psi):
    return psi * PSI_TO_PA


def pa_to_psi(pa):
    return pa / PSI_TO_PA


def ft_to_m(ft):
    return ft * FT_TO_M


def in_to_m(inch):
    return inch * IN_TO_M


def gpm_to_m3s(gpm):
    return gpm * GPM_TO_M3S