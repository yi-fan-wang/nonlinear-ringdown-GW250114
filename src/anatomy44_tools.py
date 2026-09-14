"""Helpers for the (4,4)-anatomy / referee-systematics notebook
(7-44anatomy_4Mf.ipynb).

Data flow: load waveform (SXS:BBH:3617 in geometric units, or NRSur7dq4 in
strain units) -> QNM table -> parent (2,2,n) WLS fits (direct at t0, or the
paper's 6-10 Mf multistart mean, or the production 2.0-3.67 ms grid mean) ->
QQNM amplitudes/series in (4,4) -> residual fits with 440/441(/442) ->
mismatch tables.  Everything is expressed per unit M_f so SXS and NRSur are
directly comparable.

Conventions
-----------
* times are seconds for NRSur (1 Mf = M_f * MTSUN_SI) and units of the
  total mass M for SXS (1 Mf = m_f * M); each Waveform carries `Mf` = the
  duration of one M_f in its own time unit;
* the QQNM enters the (4,4) multipole with sign -1 for NRSur7dq4 ('m/2 pi +
  pi' phase convention, as in tgr.gen_nrsur_remove_qqnm) and +1 for SXS
  (as in 20260730_NR_10qmodes.ipynb);
* `conv` converts a product of two parent amplitudes into the QQNM
  amplitude: A_q = R(chi_f) * A1 * A2 * conv, with conv = D/(G M_f/c^2) for
  strain-unit NRSur and 1/m_f for rh/M-unit SXS.
"""
import os
import sys
import numpy as np
import lal
from scipy.interpolate import interp1d
from pycbc.types import TimeSeries
from pycbc.conversions import get_lm_f0tau

import tgr.nrsurqnm as nr                                             # noqa: E402
from tgr.nrsurqnm import (QNMTable, least_square_qnmfitting,          # noqa: E402
                          weighted_least_square_qnmfitting,
                          load_interpolation_function, get_qnm_freqtau)

MODE22 = ['220', '221', '222', '223']
QUAD10 = ['220220', '220221', '221221', '220222', '221222', '222222',
          '220223', '221223', '222223', '223223']
QUAD6 = QUAD10[:6]
QUAD3 = ['220220', '220221', '221221']
QUAD1 = ['220220']
# 224 parents (Neev 2026-08 delivery, tgr branch add-224-quadratic-modes)
QUAD15 = QUAD10 + ['220224', '221224', '222224', '223224', '224224']
EXTRA_MODES = ['224', '225', '320', '440', '441', '442', '443']
DT = 1.0 / 4096
PAPER_STARTS_MF = [6., 7., 8., 9., 10.]     # Appendix C multistart of the paper


class Waveform:
    """h22, h44 (complex pycbc TimeSeries, t=0 at the peak), QNM table, units."""

    def __init__(self, label, h22, h44, tab, Mf, conv, sign, extra=None):
        self.label, self.h22, self.h44, self.tab = label, h22, h44, tab
        self.Mf, self.conv, self.sign = Mf, conv, sign
        self.extra = extra or {}

    def t(self, k_mf):
        return k_mf * self.Mf


def load_sxs3617(download=False):
    import sxs
    sim = sxs.load("SXS:BBH:3617", extrapolation='N2', ignore_deprecation=True,
                   download=download)
    h = sim.h
    h.t = h.t - h.max_norm_time()
    mf = sim.metadata['remnant_mass']
    af = sim.metadata['remnant_dimensionless_spin'][2]

    def mode_ts(mode):
        hm = h[h.index_closest_to(0.0):, h.index(*mode)]
        f = interp1d(hm.t, hm, kind='cubic')
        t = np.arange(hm.t[0], hm.t[-1], np.max(hm.t[1:] - hm.t[:-1]))
        return TimeSeries(f(t), delta_t=t[1] - t[0], epoch=t[0])
    freq, tau = {}, {}
    for mode in MODE22 + EXTRA_MODES + QUAD15:
        if len(mode) == 3:
            f0, t0 = get_lm_f0tau(mf, af, int(mode[0]), int(mode[1]), int(mode[2]))
        else:
            f1, t1 = get_lm_f0tau(mf, af, int(mode[0]), int(mode[1]), int(mode[2]))
            f2, t2 = get_lm_f0tau(mf, af, int(mode[3]), int(mode[4]), int(mode[5]))
            f0, t0 = f1 + f2, 1.0 / (1.0 / t1 + 1.0 / t2)
        freq[mode] = f0 * lal.MTSUN_SI      # 1/M
        tau[mode] = t0 / lal.MTSUN_SI       # M
    tab = QNMTable(mf, af, freq, tau)
    return Waveform('SXS:BBH:3617', mode_ts((2, 2)), mode_ts((4, 4)), tab,
                    Mf=mf, conv=1.0 / mf, sign=+1.0,
                    extra={'h32': mode_ts((3, 2)), 'h54': mode_ts((5, 4))})


def load_nrsur(label, **p):
    import pycbc.waveform
    hlm = pycbc.waveform.get_td_waveform_modes(
        approximant='NRSur7dq4', delta_t=DT, f_lower=20., f_ref=20.,
        mode_array=['22', '44', '32'], **p)
    h22 = hlm[(2, 2)][0] + 1j * hlm[(2, 2)][1]
    h44 = hlm[(4, 4)][0] + 1j * hlm[(4, 4)][1]
    tab = get_qnm_freqtau(MODE22 + EXTRA_MODES + QUAD15, f_ref=20., **p)
    Mf = tab.final_mass * lal.MTSUN_SI
    conv = p['distance'] * 1e6 * lal.PC_SI / (tab.final_mass * lal.MRSUN_SI)
    return Waveform(label, h22, h44, tab, Mf=Mf, conv=conv, sign=-1.0,
                    extra={'h32': hlm[(3, 2)][0] + 1j * hlm[(3, 2)][1],
                           'h4m4': hlm[(4, -4)][0] + 1j * hlm[(4, -4)][1],   # true (4,-4) (R2 symmetry test)
                           'params': p})


def load_bilby_maxl():
    """maxL sample parameters.  Primary source: zenodo posterior h5 (cluster;
    1.8G, not in the repo).  Fallback (laptop): attrs of the maxL injection
    file injection/GW250115-maxlpar-injection.hdf -- same maxL sample
    (ra/dec/psi/tc identical to the run config), but it stores NO spins, so
    the fallback returns zero spins: SXS-based results are unaffected, the
    NRSur 'maxL (precessing)' / 'aligned' rows are then nonspinning stand-ins
    and only the cluster run is authoritative for them.
    Also returns 'inclination' and 'coa_phase' for the (4,4) projection."""
    import h5py
    base = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base, 'injection', 'zenodo_download',
                        'posterior_samples_NRSur7dq4_no_calibration.h5')
    if os.path.exists(path):
        p = h5py.File(path, 'r')['bilby-NRSur7dq4_prod_nocal-reweighted']['posterior_samples']
        i = p['log_likelihood'][()].argmax()
        return dict(mass1=float(p['mass_1'][i]), mass2=float(p['mass_2'][i]),
                    spin1x=float(p['spin_1x'][i]), spin1y=float(p['spin_1y'][i]),
                    spin1z=float(p['spin_1z'][i]), spin2x=float(p['spin_2x'][i]),
                    spin2y=float(p['spin_2y'][i]), spin2z=float(p['spin_2z'][i]),
                    distance=float(p['luminosity_distance'][i]),
                    inclination=float(p['iota'][i]), coa_phase=float(p['phase'][i]))
    a = dict(h5py.File(os.path.join(base, 'injection',
                                    'GW250115-maxlpar-injection.hdf'), 'r').attrs)
    print('WARNING: zenodo posterior h5 not found -- maxL params from the '
          'injection-file attrs, SPINS SET TO ZERO (NRSur maxL/aligned rows '
          'are nonspinning stand-ins; SXS rows unaffected).')
    return dict(mass1=float(a['mass1']), mass2=float(a['mass2']),
                spin1x=0., spin1y=0., spin1z=0.,
                spin2x=0., spin2y=0., spin2z=0.,
                distance=float(a['distance']),
                inclination=float(a['inclination']), coa_phase=float(a['coa_phase']))


# ---------------------------------------------------------------- parent fits
def wls(w, t, fit=MODE22, omit=('224',), **kw):
    A, _, _, info = weighted_least_square_qnmfitting(list(fit), w.tab, t, w.h22,
                                                     list(omit), **kw)
    return A


def ols(w, t, fit=MODE22):
    A, _, _ = least_square_qnmfitting(list(fit), w.tab, t, w.h22)
    return A


def propagate(A, tab, t_from, t_to):
    return {m: A[m] * np.exp(-1j * tab.omega(m) * (t_to - t_from)) for m in A}


def parents_direct(w, k_mf, **kw):
    """single WLS fit at t0 = k_mf * Mf (parent_fit_draw = 'direct')."""
    return wls(w, w.t(k_mf), **kw)


def parents_direct5(w, k_mf, **kw):
    """5-parent WLS fit 220-224 (225 envelope) at t0 -- the parent set of the
    15-QQNM truncation order (224 couplings, Neev 2026-08)."""
    return wls(w, w.t(k_mf), fit=MODE22 + ['224'], omit=('225',), **kw)


def parents_multistart_mean(w, k_mf, starts_mf=PAPER_STARTS_MF, **kw):
    """WLS at each start (Mf), propagated to t0, complex mean (paper App. C recipe)."""
    t0 = w.t(k_mf)
    acc = {m: [] for m in MODE22}
    for s in starts_mf:
        A = propagate(wls(w, w.t(s), **kw), w.tab, w.t(s), t0)
        for m in MODE22:
            acc[m].append(A[m])
    return {m: complex(np.mean(acc[m])) for m in MODE22}


def parents_production_grid(w, k_mf, **kw):
    """production grid FIT_TSTART_GRID (2.0-3.67 ms, GW250114 Mf) -> t0, mean."""
    Mf_gw = 67.68 * lal.MTSUN_SI
    starts_mf = nr.FIT_TSTART_GRID / Mf_gw
    return parents_multistart_mean(w, k_mf, starts_mf=starts_mf, **kw)


# ---------------------------------------------------------------- QQNM sector
def qqnm_amps(w, A22, modes=QUAD10):
    return {q: A22[q[:3]] * A22[q[3:]] * complex(load_interpolation_function(q)(w.tab.final_spin)) * w.conv
            for q in modes}


def h44_from(w, k_mf):
    """slice of h44 from the last sample <= t0 to the end."""
    t0 = w.t(k_mf)
    idx = int(np.floor(float(t0 - w.h44.start_time) * w.h44.sample_rate))
    return w.h44.time_slice(w.h44.sample_times[idx], w.h44.sample_times[-1])


def qqnm_series(w, k_mf, Q, modes=None):
    """{mode: TimeSeries} of the QQNMs on the h44 grid from t0 (each mode
    referenced to exactly t0)."""
    sl = h44_from(w, k_mf)
    t = sl.sample_times.numpy()
    t0 = w.t(k_mf)
    out = {}
    for q in (modes or Q):
        out[q] = TimeSeries(w.sign * Q[q] * np.exp(-1j * w.tab.omega(q) * (t - t0)),
                            delta_t=sl.delta_t, epoch=sl.start_time)
    return out


def sum_series(series, modes):
    s = None
    for m in modes:
        s = series[m] if s is None else s + series[m]
    return s


def mismatch(a, b):
    """1 - <a,b>/sqrt(<a,a><b,b>) with the plain time-domain inner product
    over the common length (no PSD), as in 20260730_NR_10qmodes.ipynb."""
    a = np.asarray(a.numpy() if hasattr(a, 'numpy') else a)
    b = np.asarray(b.numpy() if hasattr(b, 'numpy') else b)
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    return 1 - np.vdot(a, b).real / np.sqrt(np.vdot(a, a).real * np.vdot(b, b).real)


def energy(a):
    a = np.asarray(a.numpy() if hasattr(a, 'numpy') else a)
    return np.vdot(a, a).real


def lin_table(w, modes):
    return QNMTable(w.tab.final_mass, w.tab.final_spin,
                    {m: w.tab.freq[m] for m in modes}, {m: w.tab.tau[m] for m in modes})


def fit_linear44(w, target, grid, lin_modes=('440', '441'), method='ols', omitted=('442',)):
    """Fit lin_modes to `target` (a slice starting at t0); reconstruct on `grid`
    (same epoch).  Returns (A dict, TimeSeries model on grid)."""
    t0 = float(target.start_time)
    ltab = lin_table(w, list(lin_modes) + list(omitted))
    if method == 'wls':
        A, _, _, _ = weighted_least_square_qnmfitting(list(lin_modes), ltab, t0, target, list(omitted))
    else:
        A, _, _ = least_square_qnmfitting(list(lin_modes), ltab, t0, target)
    t = grid.sample_times.numpy()
    lin = TimeSeries(sum(A[m] * np.exp(-1j * ltab.omega(m) * (t - t0)) for m in lin_modes),
                     delta_t=grid.delta_t, epoch=grid.start_time)
    return A, lin


def anatomy(w, k_mf, A22, lin_modes=('440', '441'), method='ols', omitted=('442',)):
    """Full (4,4) anatomy at t0: partial QQNM sums, linear-only fit, residual fits."""
    h44 = h44_from(w, k_mf)
    Q = qqnm_amps(w, A22, QUAD10)
    ser = qqnm_series(w, k_mf, Q)
    sums = {1: sum_series(ser, QUAD1), 3: sum_series(ser, QUAD3),
            6: sum_series(ser, QUAD6), 10: sum_series(ser, QUAD10)}
    # 15-QQNM truncation order: its own 5-parent (220-224, 225 envelope)
    # direct fit at t0, independent of the `A22` recipe of the 6/10 sums
    A22_5 = parents_direct5(w, k_mf)
    ser15 = qqnm_series(w, k_mf, qqnm_amps(w, A22_5, QUAD15))
    sums[15] = sum_series(ser15, QUAD15)
    A_lin_only, lin_only = fit_linear44(w, h44, h44, lin_modes, method, omitted)
    A_res, lin_res = fit_linear44(w, h44 - sums[10], h44, lin_modes, method, omitted)
    model = lin_res + sums[10]
    resid = h44 - model
    return dict(h44=h44, Q=Q, series=ser, sums=sums,
                A_lin_only=A_lin_only, lin_only=lin_only,
                A_res=A_res, lin_res=lin_res, model=model, resid=resid,
                mm_sums={n: mismatch(sums[n], h44) for n in sums},
                energy_frac_sums={n: energy(sums[n]) / energy(h44) for n in sums},
                mm_lin_only=mismatch(lin_only, h44),
                mm_model=mismatch(model, h44),
                energy_frac_resid=energy(resid) / energy(h44),
                energy_frac_lin=energy(lin_res) / energy(h44),
                A22_5=A22_5,
                qsum_over_h44=abs(sum(Q[q] for q in QUAD10)) / abs(np.interp(w.t(k_mf), w.h44.sample_times.numpy(),
                                                                              np.abs(w.h44.numpy()))))


def load_sxs3617_npz(path):
    """Load the released extrapolation-N2 SXS:BBH:3617 mode arrays."""
    data = np.load(path)
    remnant_mass = float(data['remnant_mass'])
    remnant_spin = float(data['remnant_spin_z'])
    raw_time = data['t']

    def mode_ts(key):
        values = data[key]
        interpolation = interp1d(raw_time, values, kind='cubic')
        time = np.arange(raw_time[0], raw_time[-1], np.max(np.diff(raw_time)))
        return TimeSeries(interpolation(time), delta_t=time[1] - time[0], epoch=time[0])

    frequency, damping_time = {}, {}
    for mode in MODE22 + EXTRA_MODES + QUAD15:
        if len(mode) == 3:
            f0, tau0 = get_lm_f0tau(
                remnant_mass, remnant_spin, int(mode[0]), int(mode[1]), int(mode[2]))
        else:
            f1, tau1 = get_lm_f0tau(
                remnant_mass, remnant_spin, int(mode[0]), int(mode[1]), int(mode[2]))
            f2, tau2 = get_lm_f0tau(
                remnant_mass, remnant_spin, int(mode[3]), int(mode[4]), int(mode[5]))
            f0, tau0 = f1 + f2, 1.0 / (1.0 / tau1 + 1.0 / tau2)
        frequency[mode] = f0 * lal.MTSUN_SI
        damping_time[mode] = tau0 / lal.MTSUN_SI
    table = QNMTable(remnant_mass, remnant_spin, frequency, damping_time)
    return Waveform('SXS:BBH:3617', mode_ts('h22'), mode_ts('h44'), table,
                    Mf=remnant_mass, conv=1.0 / remnant_mass, sign=+1.0)
