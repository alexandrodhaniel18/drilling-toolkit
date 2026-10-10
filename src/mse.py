import numpy as np


def _require_positive(x, name):
    """Ubah ke array float dan tolak nilai <= 0."""
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        raise ValueError(f'{name} harus > 0')
    return x


def _require_nonnegative(x, name):
    """Ubah ke array float dan tolak nilai < 0 (nol diperbolehkan)."""
    x = np.asarray(x, dtype=float)
    if np.any(x < 0):
        raise ValueError(f'{name} tidak boleh negatif')
    return x


def bit_area(d_bit):
    """Luas dasar lubang (m2) dari diameter bit d_bit (m)."""
    d_bit = _require_positive(d_bit, 'd_bit')
    return np.pi / 4.0 * d_bit ** 2


def mse_axial(wob, d_bit):
    """Suku aksial MSE (Pa) = WOB / A. wob dalam N."""
    return _require_nonnegative(wob, 'wob') / bit_area(d_bit)


def mse_rotary(torque, rps, rop, d_bit):
    """Suku rotari MSE (Pa) = 2*pi*N*T / (A*ROP).

    torque: N.m; rps: putaran/detik; rop: m/s (harus > 0); d_bit: m.
    """
    torque = _require_nonnegative(torque, 'torque')
    rps = _require_nonnegative(rps, 'rps')
    rop = _require_positive(rop, 'rop')
    return 2.0 * np.pi * torque * rps / (bit_area(d_bit) * rop)


def mse(wob, torque, rps, rop, d_bit):
    """MSE total (Pa), rumus Teale (1965). Semua masukan SI; array didukung."""
    return mse_axial(wob, d_bit) + mse_rotary(torque, rps, rop, d_bit)


def mse_breakdown(wob, torque, rps, rop, d_bit):
    """Rincian MSE: dict berisi axial, rotary, total (Pa), dan rotary_share."""
    axial = mse_axial(wob, d_bit)
    rotary = mse_rotary(torque, rps, rop, d_bit)
    total = axial + rotary
    return {'axial': axial, 'rotary': rotary, 'total': total,
            'rotary_share': rotary / total}


def mse_efficiency(mse_value, ccs):
    """Efisiensi mekanis = CCS / MSE (tanpa satuan). Keduanya dalam Pa."""
    mse_value = _require_positive(mse_value, 'mse')
    ccs = _require_positive(ccs, 'ccs')
    return ccs / mse_value