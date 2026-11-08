"""Protocol catalogue of the coexistence simulator: specs, frame layouts, access exchanges, period choice.

Purpose: turn parameters.json into one ProtocolSpec per (protocol, variant,
PHY): every stack at its full capability, the FULL "Simulation recipe" of each protocol in PARAMETERS.md section 3
(DESIGN.md 1), and the as-built models (variant 'thesis': the as-built BCC stack and the plain BLE model
ble4thesis); and give the engine what it needs per frame and per period (DESIGN.md 3 and 4).
- get_spec(pid, variant, phy=None, **flags) -> ProtocolSpec; variant 'full', or 'thesis' for bcc and ble4thesis.
- frame_layout(spec, L) -> blocks for phy.frame_outcome: dict(t0_us, t1_us, kind 'sync'|'header'|'payload',
  n_bits = hard decisions (bits, coded bits after repetition combining, or chips), info_bits (application payload
  bits, 8 L summed over the payload blocks), code None | ('bch', n, k, t) | ('walsh16',) | ('conv', gain_db),
  optional rep = SNR factor of the block when it differs from phy['rep'] (802.15.6 NB header at S=4 and preamble at
  S=1, HBC header and sync at their own FSC spreading); it also sets the block's processing gain against the C/I
  floor (phy.ci_gain)). ('bch', 64, 0, 7) is the 64-chip BCC/HBC sync correlator (<= 7 errors).
  airtime_us(spec, L) -> int, the layout end rounded up. Multi-tone BCC ('bcc', 'full', 'multitone') stripes a frame
  over its N tones: tone_octets(spec, L) -> the octets of each tone
  packet, tone_layouts(spec, L) -> their layouts; frame_layout and airtime_us give the longest tone packet.
- exchange(spec, L) -> [Step(actor, kind, dur_us, tx)], one access exchange; 'ifs' steps are idle gaps named
  after the node that transmits next. Backoff (CSMA) and retries are drawn by the engine, not listed here.
- reference_period(r_kbps, mode) -> the period P of a cell, the same for every protocol (DESIGN 3).
- select_period(r_kbps, mode, n_dev, spec) -> dict(period_us, payload_bits, frames_per_period,
  frame_payload_octets, r_actual_kbps, deadline_us, optional slot_us and reserved_us) or None (DESIGN 3; see its
  docstring for the scan).
- BCC full ladder (flag opt 0..5, DESIGN 4.3): pool_size(frames) -> exchanges of the L3 retry pool,
  poll_octets(spec) -> payload octets of a pool poll, beacon_layout(spec, n_dev, n) -> the beacon with its L2 ack
  field; select_period adds beacon_us, reserved_us and from L3 slot_us, pool and pool_us (one pool exchange).
- catalogue() -> the (pid, variant, phy) triples the sweep runs by default: every rate or PHY of each standard.
- acr_matrix(spec) -> ACR (dB) between every pair of channels, from their centre-frequency offset.
Spec conventions: channels = centre frequencies (MHz) indexed 0..n-1; acr_db[d] = adjacent-channel rejection in dB
at a centre-frequency offset of d grid spacings (the smallest channel spacing; the last table entry holds for
every larger offset, except 802.15.6 NB, whose table ends at offset 1); SINR of BCC and HBC is per chip (Ec/N0),
hence noise_bw_hz = chip rate for BCC; period_grid = allowed periods in us, ascending.
Inputs:  parameters.json (read once by load_params).
Outputs: none (library); tests in test_protocols.py.
"""
from __future__ import annotations
import bisect
import heapq
import json
import math
from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple, Optional

import channel_curves

HERE = Path(__file__).resolve().parent


@lru_cache(maxsize=None)
def load_params():
    """parameters.json as a dict, read once (shared: callers must not mutate it)."""
    return json.loads((HERE / "parameters.json").read_text(encoding="utf-8"))


P = load_params()
BLE, B154, CSMA, TSCH = P["ble5"], P["ieee802154_phy"], P["ieee802154_csma"], P["ieee802154_tsch"]
NB, HBC, SB, ANT, BCC, EN = (P["ieee802156_nb"], P["ieee802156_hbc"], P["smartban"], P["ant"], P["bcc"],
                             P["energy"])

# Values PARAMETERS.md gives that parameters.json lacks as numbers (several sit there only inside a formula string).
MISSING = dict(
    coded_gain_db={"coded_s2": 4.67, "coded_s8": 10.69},  # 3.1 BER model: soft-model gain over 1M at BER 1e-3
    coded_term_bits=3,              # 3.1 coding: TERM1 / TERM2 3-bit terminations (on air, not at risk)
    phr_octets=1,                   # 3.2 frame: PHR 1 octet (shr_phr_octets 6 = SHR 5 + PHR 1)
    beacon154_octets=(13, 1, 3),    # 3.2 frame: beacon PSDU 13 + 1 + 3 n_gts octets (formula string in the json)
    max_bo=14,                      # 3.2 beacon: BI = 15.36 ms x 2^BO, BO <= 14 (formula string in the json)
    nb_alloc_slot_us=500,           # 3.4 MAC: allocation slot 500 + L x 500 us (formula string in the json)
    smartban_sync_bits=(16, 32, 1),  # 3.5: preamble 16 + sync 32 bits; sync lost above 1 error ([50], BER model)
    smartban_phy_crc_bits=16,       # 3.5 PPDU: PHY CRC-16 after the PSDU
    smartban_mpdu_overhead_bits=64,  # 3.5 MAC frame: MPDU = 64 + 8 LF bits (formula string in the json)
    ant_frame_octets=(1, 5, 2, 2),  # 3.6 frame: preamble 1, address 5, control 2, CRC 2 (+ 8 data = 18 octets)
    ant_clock_hz=32768,             # 3.6 timing: channel period N / 32768 s
    ant_ack_rate_factor=1.25,       # 3.6 FULL step 2: acknowledged messages at f_ch = 1.25 f_data (A)
    hbc_sync_sf=8,                  # 3.7 PHY: preamble, SFD and pilots at FSC SF 8 (2048 = 256 x 8 chips)
    smartban_lost_beacons=10,       # OI-33: after 10 missed D-Beacons the node rescans ([51])
)
# Assumptions of this module, no value in PARAMETERS.md:
ARQ_QUEUE = CSMA["queue_depth"]     # queue of every ARQ variant: the 802.15.4 / TSCH FIFO depth 4 (OI-18)
SPARE_EXCHANGES = 1                 # retry room of a scheduled ARQ allocation per frame (DESIGN 3)
# DESIGN 3, 4.2 csma_max_load: the scheduled CSMA load holds the same retry room as the allocations and BLE events,
# one spare exchange per frame
CSMA_MAX_LOAD = 1 / (1 + SPARE_EXCHANGES)
# galvanic-optimisation ladder of BCC full: L3 retry pool of max(1, ceil(share x frames))
# exchanges per superframe, share = flag pool: 0.05 by default, 0.43 = 0.3 / 0.7 for the measured
# contact loss of about 30 % (both design points); L4 proposal: a 2-octet MAC header (frame type
# bit, then a data frame's 15-bit sequence number or a poll's recipient NID)
POOL_SHARE = 0.05
# touching bodies share the galvanic medium (author's demonstration: a 7-person hand-holding chain passed FM at
# maximum HackRF power; no loss figure exists, so the loss is a sweep); mean contact 5 s (A: a handshake to a
# caregiver's touch); the partner reaches a receiver at its own-link Ec/N0 less TOUCH_DB by default
TOUCH_DB, TOUCH_MEAN_US = 10.0, 5e6
# 15 bits: a sensor holds 4 periods (OI-18) of up to 2020 frames each (capacity grid), 8080 < 2^15;
# the beacon's base sequence number from L4 and the poll's body carry the same 15-bit number in 2 octets
LEAN_MAC_OCTETS = SEQ_OCTETS = 2
# block-ack reference (flag back): B-Ack payload octets, as claim4.back_cost (PARAMETERS 3.7b, 6.2.9.4);
# flag span: a block of this many frames or more never closes (the limit of an unlimited block, below every block the
# standard allows, 6.2.9.4)
BACK_OCTETS, BLOCK_UNLIMITED = 3, 10**9
# payload spreading of BCC full. 'walsh' the HBC Walsh-16 symbols (4 bits in 16 chips), 'none' unspread BFSK
# (1 bit per chip, 4 times the payload rate), 'adapt' the per-link switch: unspread, Walsh after ADAPT_DOWN
# consecutive frames the beacon reports lost, unspread again after ADAPT_UP consecutive delivered ones (A: hysteresis
# 2 / 16, so one bad contact dwell does not toggle it); the allocations are sized for Walsh frames (DESIGN 4.3)
SPREADS, ADAPT_DOWN, ADAPT_UP = ("walsh", "none", "adapt"), 2, 16
# rate rule: the weakest path of the link budget, RHW1-LHW1, needs A_95 = 15.7 mV at 500 kchip/s
# and NF 20 dB for Ec/N0 8 dB (the 95th-percentile participant, cleaned 50 Ohm means at 8 MHz)
A95_WEAKEST = dict(amp_v=15.7e-3, rc=500e3, nf_db=20.0, ec_n0_db=8.0)
# (results/bcc_rate_stability.md): the highest chip rate at which one BCC tone makes no chip error on any
# measured path at NF 10 to 30 in real passband; from 9.5 Mchip/s its lower tone f - R_c (1 MHz) fades with its image,
# which the engine's one-chip rule does not model, and 9.95 Mchip/s loses every packet
BCC_RC_MAX_CPS = 9e6
PERIOD_CAP_US = BLE["ci_max_us"]    # 4 s, the BLE CI maximum, for protocols without a period bound
P_LOSS_1M = P["common"]["legacy_ble_sim_residual_per"]["1M"]   # DESIGN 5.1: 1M-class residual loss

FLAGS = {"ble": set(), "ble4thesis": set(), "ieee802154": set(), "ieee802156nb": set(), "tsch": set(),
         "smartban": set(), "ant": set(), "hbc": {"ri_mode"},
         "bcc": {"rc_cps", "fcs", "contact", "coupling_db", "power", "opt", "pool", "fill", "touch", "touch_db", "fallback",
                 "spread", "lmax", "back", "burst", "span", "pilot128"}}
# 'thesis': the as-built models, the as-built BCC stack (HBCProtocol.py) and the plain BLE model (BLE_SIM.py)
VARIANTS = {"bcc": ("full", "thesis"), "ble4thesis": ("thesis",)}

# channel plans: the formula strings of parameters.json (ble5.channel_freq_MHz, ieee802154_phy.channels,
# ieee802156_nb.channels, ant.freq_pool_MHz)
BLE_CHANNELS = [2404 + 2 * i if i <= 10 else 2428 + 2 * (i - 11) for i in range(BLE["n_data_channels"])]
CH154 = [2405 + 5 * (k - 11) for k in range(11, 27)]
NB_CHANNELS = [2402 + nc for nc in range(79)]
ANT_POOL = [f for f in range(2403, 2481) if f not in (2450, ANT["freq_default_MHz"])]
# ANT channel periods N / 32768 s, N within the nRF24AP2 clamp; ponytail: rounded to whole us, the engine adds drift
ANT_PERIODS = [round(n * 1e6 / MISSING["ant_clock_hz"]) for n in range(ANT["n_clamp"]["nrf24ap2"][0],
                                                                        ANT["n_clamp"]["nrf24ap2"][1] + 1)]


@dataclass
class ProtocolSpec:
    """Everything the engine and phy.py need about one protocol, variant and PHY (see the module docstring)."""
    pid: str
    variant: str
    phy: dict               # name, mod, bit_rate_bps, noise_bw_hz, rep, fec, payload_crc, ci_co_lin, nf_db,
    #                         l_impl_db, tx_power_dbm (None: drive amplitude amp_v in V, BCC / HBC), sync_sinr_db,
    #                         sens_dbm (None: no sensitivity floor; 802.15.6 NB FULL: the rate's draft value)
    channels: list          # centre frequencies, MHz
    channel_bw_mhz: float
    acr_db: list            # by channel-index offset
    access: str             # 'tdma' | 'csma' | 'ant'
    hop: str                # 'perm' | 'csa2' (with AFH) | 'tsch' (with blacklisting) | 'fixed'
    period_grid: object     # allowed periods (us), ascending (range or list)
    max_payload_octets: int  # application payload per frame; larger L is fragmented (DESIGN 3)
    p_loss: float           # residual loss of a surviving frame in the collision model (DESIGN 5.1)
    inter_body: str         # 'rf' | 'isolated' | 'coupled'
    energy: dict            # tx_ma, rx_ma, ramp_tx_ma, ramp_rx_ma, ramp_us, wake_uc, sleep_ua, vdd_v
    flags: dict             # the get_spec flags (CSV label, BCC ladder level)
    arq_max_retx: Optional[int] = 0     # 0 none, None unlimited (BLE ACL)
    queue_len: Optional[int] = None     # frames a sensor may hold; None: unbounded (no ARQ, each frame sent once)
    ed_scan: bool = False               # least-used channel at BAN start (F5)
    channel_switch: bool = False        # SmartBAN DCH migration (F5)
    adapt: Optional[dict] = None        # window, threshold (failure ratio), min_channels of the F5 rule (DESIGN 4.1)
    csma: Optional[dict] = None         # 802.15.4 CSMA/CA constants, or SmartBAN MCA (mode 'mca')
    tsch: Optional[dict] = None         # TSCH slot template; hopping_sequence as indices into channels
    beacon: Optional[dict] = None       # airtime_us (+ airtime_per_dev_us), reserved_us before the first slot,
    #                                     max_lost: missed beacons after which a sensor loses its frames, energy:
    #                                     its radio currents when not spec.energy's (multi-tone BCC: one front end)
    ant: Optional[dict] = None          # ANT slot timing, burst and channel-rate constants
    sync_on_foreign: bool = False
    max_dev: Optional[int] = None       # sensors per hub the protocol allows (None: no cap)
    chip_rate_cps: Optional[float] = None
    tones: Optional[list] = None        # multi-tone BCC: tone centres (Hz) a frame is striped over (None: one tone)
    drift_ppm: float = 0.0              # clock error bound per device (OI-8)

    @property
    def period_min_us(self):
        return self.period_grid[0]

    @property
    def period_max_us(self):
        return self.period_grid[-1]

    @property
    def period_step_us(self):
        """Grid step; 0 when the grid is not uniform (802.15.4 2^BO, ANT N / 32768 s): use period_grid."""
        return getattr(self.period_grid, "step", 0)


class Step(NamedTuple):
    actor: str      # 'sensor' | 'hub'
    kind: str       # 'data' | 'burst' (an RI-mode data frame of a burst) | 'poll' | 'ack' | 'cca' | 'ifs'
    dur_us: float
    tx: bool


def _lin(db):
    return 10 ** (db / 10)


def _acr(table_mhz, spacing_mhz, channels=()):
    """ACR (dB) by offset 0..reach in units of spacing_mhz from a {MHz: dB} table; the largest key <= offset applies,
    so the last entry holds up to the span of `channels` (PARAMETERS 5.5 gives it as a floor, '>= 3 MHz'); without
    channels the table ends at its last key (802.15.6 NB: |dnc| >= 2 ignored, 3.4)."""
    keys = sorted((float(k), v) for k, v in table_mhz.items())
    reach = max(math.ceil(keys[-1][0] / spacing_mhz - 1e-9),
                round((max(channels) - min(channels)) / spacing_mhz) if channels else 0)
    return [0.0] + [float([v for k, v in keys if k <= d * spacing_mhz + 1e-9][-1]) for d in range(1, reach + 1)]


def acr_matrix(spec):
    """ACR (dB) between every pair of spec.channels, looked up by their centre-frequency offset in units of the
    smallest channel spacing (BLE data indices 10 and 11 sit 4 MHz apart), inf beyond the table (ignored)."""
    f = spec.channels
    if len(f) == 1:
        return [[0.0]]
    step = min(b - a for a, b in zip(sorted(f), sorted(f)[1:]))
    return [[spec.acr_db[d] if (d := round(abs(a - b) / step)) < len(spec.acr_db) else math.inf for b in f]
            for a in f]


def _energy(tx_ma, rx_ma, ramp_us, ramp_tx_ma=None, ramp_rx_ma=None):
    w = EN["wake_event"]
    return dict(tx_ma=tx_ma, rx_ma=rx_ma, ramp_tx_ma=tx_ma if ramp_tx_ma is None else ramp_tx_ma,
                ramp_rx_ma=rx_ma if ramp_rx_ma is None else ramp_rx_ma, ramp_us=ramp_us,
                wake_uc=w["us"] * w["mA"] / 1000, sleep_ua=EN["sleep_default_uA"], vdd_v=P["common"]["vdd_V"])


def _nrf(tx, rx, ramp):
    """nRF52840 compounded currents at 0 dBm; the ramp draws the start-up current of the target state (4.1)."""
    n = EN["nrf52840"]
    return _energy(n[tx], n[rx], n["ramp_us"][ramp], n["startup_tx"], n["startup_rx"])


def _rf_phy(name, mod, bit_rate, noise_bw, ci_co_db, nf_db, l_impl_db=0.0, rep=1, fec=None, sync_sinr_db=None,
            sens_dbm=None):
    return dict(name=name, mod=mod, bit_rate_bps=bit_rate, noise_bw_hz=noise_bw, rep=rep, fec=fec,
                payload_crc=True, ci_co_lin=_lin(ci_co_db), nf_db=nf_db, l_impl_db=l_impl_db,
                tx_power_dbm=BLE["tx_power_dBm"], sync_sinr_db=sync_sinr_db,   # 0 dBm for every RF stack (DESIGN 5.2)
                sens_dbm=sens_dbm)


def _coded_gain(name):
    """Coded-PHY gain of the stated per-bit model: soft-model gain minus the extra loss (3.1, OI-12)."""
    return MISSING["coded_gain_db"][name] - BLE["phy"][name]["extra_loss_dB"]


# --- builders: (variant, phy, flags) -> ProtocolSpec --------------------------------------------------------------

def _ble(variant, phy, flags):
    """ACL with poll, LL ACK and retransmission in the same event, CSA#2 with AFH; PHY 1M, 2M, coded_s2 or coded_s8."""
    name = phy or "1M"
    ph, coded = BLE["phy"][name], name.startswith("coded")
    rate = 2e6 if name == "2M" else 1e6     # symbol rate; the coded gain sits in the ('conv', gain) blocks
    legacy = P["common"]["legacy_ble_sim_residual_per"]
    return ProtocolSpec(
        pid="ble", variant=variant,
        phy=_rf_phy(name, "gfsk", rate, ph["noise_bw_Hz"], ph["ci_co_spec_dB"],
                    BLE["ber"]["nf_plus_limpl_dB"]["curve_a_nrf52840"],
                    fec=("conv", _coded_gain(name)) if coded else None),
        channels=BLE_CHANNELS, channel_bw_mhz=ph["noise_bw_Hz"] / 1e6,
        acr_db=_acr(ph["acr_spec_dB"], 2, BLE_CHANNELS),
        access="tdma", hop="csa2",
        adapt=dict(window=BLE["afh"]["window"], threshold=1 - BLE["afh"]["pdr_threshold"],
                   min_channels=BLE["afh"]["c_min"]),
        arq_max_retx=None, queue_len=ARQ_QUEUE,
        period_grid=range(BLE["ci_min_us"], BLE["ci_max_us"] + 1, BLE["ci_step_us"]),
        max_payload_octets=BLE["payload_max_octets"], p_loss=legacy.get(name, P_LOSS_1M), inter_body="rf",
        energy=_nrf("tx_0dBm", "rx_2M" if name == "2M" else "rx_1M", "ble"),
        drift_ppm=P["common"]["clock_ppm"]["ble"], flags=flags)


def _ble4thesis(variant, phy, flags):
    """BLE_SIM.py rules: payload only on air, one shuffled 37-channel permutation per body, no ARQ (3.1 legacy)."""
    name = phy or "1M"
    rate = {"1M": 1e6, "2M": 2e6}[name]
    s = _ble(variant, name, flags)
    s.pid, s.acr_db, s.hop, s.adapt, s.drift_ppm = "ble4thesis", [0.0], "perm", None, 0.0
    s.arq_max_retx, s.queue_len = 0, None
    s.phy.update(bit_rate_bps=rate, noise_bw_hz=rate, ci_co_lin=None)   # collision model only
    return s


def _phy154():
    return _rf_phy("oqpsk", "oqpsk_dsss", B154["bit_rate_bps"], B154["noise_bw_Hz"], B154["ci_co_dB"],
                   B154["nf_plus_limpl_dB"]["nrf52840"], sync_sinr_db=B154["preamble_drop_sinr_dB"])


def _ieee802154(variant, phy, flags):
    """Imm-Ack, 3 retries, ED-scan channel, on the beacon-interval grid 15.36 ms x 2^BO (OI-19). Phy 'oqpsk' (default):
    mode A, unslotted CSMA/CA. Phy 'gts' (the MAC mode in the phy column): beacon mode B, one GTS per sensor at the end
    of the superframe, retries inside the GTS (it holds one spare exchange per frame, as the 802.15.6 allocations),
    else in the next superframe's; a sensor that missed aMaxLostBeacons 4 beacons loses its frames (3.2 mode B)."""
    gts, bc = phy == "gts", CSMA["beacon"]
    o, (b0, b1, b2) = B154["octet_us"], MISSING["beacon154_octets"]
    return ProtocolSpec(
        pid="ieee802154", variant=variant, phy=dict(_phy154(), name="gts") if gts else _phy154(), channels=CH154,
        channel_bw_mhz=B154["noise_bw_Hz"] / 1e6, acr_db=_acr(B154["acr_chip_dB"], 5, CH154),
        access="tdma" if gts else "csma", hop="fixed", ed_scan=True,
        arq_max_retx=CSMA["mac_max_frame_retries"], queue_len=ARQ_QUEUE,
        csma=None if gts else dict({k: v for k, v in CSMA.items() if k not in ("beacon", "full", "fair")},
                                   ack_us=B154["imm_ack_us"], max_load=CSMA_MAX_LOAD),
        # GTS beacon: PSDU 13 + 1 + 3 n_gts octets; the GTSs follow the CAP (select_period sets reserved_us)
        beacon=dict(airtime_us=(B154["shr_phr_octets"] + b0 + b1) * o, airtime_per_dev_us=b2 * o,
                    reserved_us=bc["min_cap_symbols"] * B154["symbol_us"], max_lost=bc["max_lost_beacons"])
        if gts else None,
        max_dev=bc["max_gts"] if gts else None,
        period_grid=[bc["base_superframe_us"] << bo for bo in range(MISSING["max_bo"] + 1)],
        max_payload_octets=B154["payload_max_octets"], p_loss=P_LOSS_1M, inter_body="rf",
        energy=_nrf("tx_154", "rx_154", "154"), sync_on_foreign=True,
        drift_ppm=P["common"]["clock_ppm"]["ieee802154"], flags=flags)


def _tsch(variant, phy, flags):
    """Enh-Ack, 3 retries in dedicated and shared cells, 16-channel hopping with blacklisting (3.3 FULL+: a channel
    with PDR < 0.8 over its last 100 outcomes leaves the BAN's sequence, >= 4 kept), CCA before every data frame."""
    ts, b = TSCH["timeslot_us"], TSCH["full"]["blacklist"]
    t = {k: v for k, v in TSCH.items() if k not in ("full", "fair")}
    t["hopping_sequence"] = [c - 11 for c in TSCH["hopping_sequence"]]    # channel 11 = index 0
    return ProtocolSpec(
        pid="tsch", variant=variant, phy=_phy154(), channels=CH154, channel_bw_mhz=B154["noise_bw_Hz"] / 1e6,
        acr_db=_acr(B154["acr_chip_dB"], 5, CH154), access="tdma", hop="tsch",
        adapt=dict(window=b["window"], threshold=1 - b["pdr_threshold"], min_channels=b["min_channels"]),
        arq_max_retx=TSCH["max_retx"], queue_len=TSCH["queue_depth"], tsch=t,
        # slot 0 of the slotframe is the shared EB cell; an EB goes out every eb_period_s only
        beacon=dict(airtime_us=(B154["shr_phr_octets"] + TSCH["eb_psdu_octets"]) * B154["octet_us"],
                    reserved_us=ts, period_us=TSCH["eb_period_s"] * 1_000_000),
        period_grid=range(2 * ts, PERIOD_CAP_US + 1, ts), max_payload_octets=B154["payload_max_octets"],
        p_loss=P_LOSS_1M, inter_body="rf", energy=_nrf("tx_154", "rx_154", "154"),
        drift_ppm=P["common"]["clock_ppm"]["ieee802154"], flags=flags)


_NB_MS = {k: (v["M"], v["S"]) for k, v in NB["rates_kbps"].items()}


def _nb(variant, phy, flags):
    """BCH(63,51) PSDU + BCH(31,19) header at S=4, the rate's modulation and PSDU spreading (121.4 kb/s by default;
    242.9, 485.7, 971.4), I-Ack with 3 retransmissions, ED-scan channel, the rate's draft sensitivity."""
    name = phy or NB["full_default_rate_kbps"]
    m, s = _NB_MS[name]
    c = NB["bch_psdu"]
    spec = ProtocolSpec(
        pid="ieee802156nb", variant=variant,
        phy=_rf_phy(name, "dbpsk" if m == 2 else "dqpsk", NB["symbol_rate_Hz"] * math.log2(m), NB["symbol_rate_Hz"],
                    BLE["phy"]["1M"]["ci_co_spec_dB"], NB["nf_dB"], NB["impl_loss_dB"], rep=s,   # C/I: BLE analogue
                    fec=("bch", c["n"], c["k"], c["t"]),                                        # (OI-3)
                    # 3.4 FULL: 'frames below the rate's sensitivity (TI draft values) are lost'
                    sens_dbm=NB["sens_draft_dBm"][name]),
        channels=NB_CHANNELS, channel_bw_mhz=1.0, acr_db=_acr({"1": NB["acr_adjacent_dB"]}, 1),
        access="tdma", hop="fixed", ed_scan=True, arq_max_retx=NB["retry_limit"], queue_len=ARQ_QUEUE,
        period_grid=range(MISSING["nb_alloc_slot_us"], PERIOD_CAP_US + 1, MISSING["nb_alloc_slot_us"]),
        max_payload_octets=NB["body_max_octets"], p_loss=P_LOSS_1M, inter_body="rf",
        energy=_nrf("tx_0dBm", "rx_1M", "ble"), sync_on_foreign=True, flags=flags)
    spec.beacon = dict(airtime_us=airtime_us(spec, NB["beacon_body_octets"]))   # 17-octet body at the data PHY (A)
    spec.beacon["reserved_us"] = spec.beacon["airtime_us"] + NB["psifs_us"]
    return spec


# SmartBAN PHY schemes (3.5, TS 103 326 Table 3), assigned statically by the hub: PSDU BCH(127,113) or uncoded, the
# whole PPDU sent Nrep = 1, 2 or 4 times
SB_PHYS = {f"{c}{n}": (c == "bch", n) for c in ("bch", "unc") for n in (1, 2, 4)}


def _smartban(variant, phy, flags):
    """Header BCH, PSDU BCH ('bch1', default; 'bch2', 'bch4') or uncoded ('unc1', 'unc2', 'unc4') with the PPDU
    repeated Nrep times (ideal EGC in SINR, each copy decided alone in the collision model), ACK with 1 retransmission
    (in the C/M slots by the MCA backoff, else in the node's next slot), MCA sensing, scanned DCH with migration."""
    name = phy or "bch1"
    bch, nrep = SB_PHYS[name]
    c = SB["bch_psdu"]
    spec = ProtocolSpec(
        pid="smartban", variant=variant,
        phy=dict(_rf_phy(name, "gfsk", SB["symbol_rate_Hz"], SB["symbol_rate_Hz"], SB["ber"]["ci_co_dB_assumed"],
                         BLE["ber"]["nf_plus_limpl_dB"]["curve_a_nrf52840"],      # common GFSK curve, as BLE (3.5)
                         rep=nrep, fec=("bch", c["n"], c["k"], c["t"]) if bch else None), nrep=nrep),
        channels=BLE_CHANNELS, channel_bw_mhz=1.0,
        acr_db=_acr(BLE["phy"]["1M"]["acr_spec_dB"], 2, BLE_CHANNELS),
        access="tdma", hop="fixed", ed_scan=True, channel_switch=True,
        adapt=dict(window=SB["full"]["migration"]["window_ibi"], threshold=SB["full"]["migration"]["fail_ratio"],
                   min_channels=1),
        arq_max_retx=SB["retry_limit"], queue_len=ARQ_QUEUE,
        csma=dict(mode="mca", cca_us=SB["tmua_us"]["owner"], ed_threshold_dbm=SB["ed_threshold_dBm"]),
        # D-Beacon in slot 0 (reserved_us: the smallest slot; select_period sets the LSLOT slot of the cell)
        beacon=dict(airtime_us=SB["d_beacon_us"], reserved_us=SB["tmin_us"], max_lost=MISSING["smartban_lost_beacons"]),
        max_dev=SB["max_nodes"], period_grid=range(SB["tmin_us"], PERIOD_CAP_US + 1, SB["tmin_us"]),
        max_payload_octets=SB["lf_max_octets_nrep1"]["32"][bch] // nrep, p_loss=P_LOSS_1M, inter_body="rf",
        energy=_nrf("tx_0dBm", "rx_1M", "ble"), sync_on_foreign=True,
        drift_ppm=P["common"]["clock_ppm"]["smartban"], flags=flags)
    # LF,max of the table holds for TMUA 0 and Nrep 1 (its Nrep-th part is above the Nrep maximum); the owner's MCA
    # sensing (16 us) takes 2 octets off at LSLOT 32
    while _sb_slot(spec, spec.max_payload_octets) is None:
        spec.max_payload_octets -= 1
    return spec


def _ant(variant, phy, flags):
    """Acknowledged messages of 8 octets at 1.25 f_data, 6 attempts, listen before channel open (the least-occupied
    phase when none is free), advanced burst above 8 octets (packets of 8, 16 or 24 octets every 3.2 ms, 60 kb/s, 20
    retries; its optional frequency hopping not modelled), adaptive isochronous re-timing; a private network: one
    random frequency per BAN from the 76 of the pool (2450 and 2457 MHz excluded), with frequency agility over it and
    two more of the pool (a master moves on after 8 unacknowledged messages; DESIGN 4.1)."""
    return ProtocolSpec(
        pid="ant", variant=variant,
        phy=_rf_phy("gfsk", "gfsk", 1e6, 1e6, ANT["ci_co_dB"], ANT["nf_plus_limpl_dB"]["nrf52832"]),
        channels=ANT_POOL, channel_bw_mhz=1.0, acr_db=_acr(ANT["acr_dB"], 1, ANT_POOL), access="ant", hop="fixed",
        arq_max_retx=ANT["full"]["attempts"] - 1, queue_len=ARQ_QUEUE,
        ant=dict(master_tx_us=ANT["slot_offsets_us"]["master_tx"], slave_reply_us=ANT["slot_offsets_us"]["slave_reply"],
                 turnaround_us=ANT["turnaround_us"], rx_window_us=ANT["rx_window_us"],
                 burst_spacing_us=ANT["burst_spacing_us"], burst_retries=ANT["adv_burst_retries"],
                 ch_rate_factor=MISSING["ant_ack_rate_factor"], guard_us=ANT["full"]["guard_us"],
                 iso_shift_us=ANT["full"]["iso_shift_us"], agility_freqs=len(ANT["agility_freqs_MHz"]),
                 agility_unacked=ANT["full"]["agility_trigger_unacked"]),
        max_dev=ANT["max_channels"]["s212"],
        period_grid=ANT_PERIODS,
        # an advanced-burst packet of 24 octets (3.6 [53] 9.5.7.8); a message carries 8
        max_payload_octets=max(map(int, ANT["adv_burst_airtime_us"])), p_loss=P_LOSS_1M, inter_body="rf",
        energy=_nrf("tx_0dBm", "rx_1M", "ble"), sync_on_foreign=True, drift_ppm=P["common"]["clock_ppm"]["ant"],
        flags=flags)


def _hbc(variant, phy, flags, walsh_cps=None):
    """FSC SF 8 (phy 'sf8'..'sf64'), header at 164 kb/s (SF 64) unless ri_mode, I-Ack with 3 retransmissions (3.7).
    walsh_cps: the Walsh-chip rate of the SF 8 mode when it is not 5.25 Mcps (BCC full, _bcc)."""
    sf = int((phy or "sf8")[2:])
    sync_sf = MISSING["hbc_sync_sf"]
    hdr_sf = sf if flags.get("ri_mode") else HBC["sf_by_rate_kbps"]["164.0625"]
    rc, half_mw = HBC["chip_rate_cps"], HBC["power_total_mW"] / 2
    if walsh_cps is not None:
        rc = MISSING["hbc_sync_sf"] * walsh_cps
    spec = ProtocolSpec(
        pid="hbc", variant=variant,
        phy=dict(name=f"sf{sf}", mod="bfsk_chip_v2", bit_rate_bps=rc / sf, noise_bw_hz=HBC["bw_MHz"] * 1e6, rep=1,
                 fec=None, payload_crc=True, ci_co_lin=_lin(BLE["phy"]["1M"]["ci_co_spec_dB"]),  # OI-3 analogue
                 nf_db=BCC["nf_dB"][1], l_impl_db=0.0, tx_power_dbm=None, amp_v=BCC["drive_V"],
                 sync_sinr_db=None, sf=sf, sync_sf=sync_sf, header_sf=hdr_sf,
                 # processing gain against the C/I floor: the SF over SF 8, which is spectral shaping only (3.7:
                 # 'gamma_psdu = gamma + 10 log10(S) at the lower rates'); the blocks' rep carries header and sync
                 ci_pg=sf / sync_sf),
        channels=[HBC["center_MHz"]], channel_bw_mhz=HBC["bw_MHz"], acr_db=[0.0], access="tdma", hop="fixed",
        arq_max_retx=HBC["retry_limit"], queue_len=ARQ_QUEUE, chip_rate_cps=rc,
        period_grid=range(MISSING["nb_alloc_slot_us"], PERIOD_CAP_US + 1, MISSING["nb_alloc_slot_us"]),
        max_payload_octets=HBC["body_max_octets"], p_loss=P_LOSS_1M, inter_body="coupled",
        # 5.2 mW IC, TX = RX = 2.6 mW; settling time: the BCC typical AFE value (A, no HBC figure)
        energy=_energy(half_mw / P["common"]["vdd_V"], half_mw / P["common"]["vdd_V"], EN["bcc_wake_us"]["typical"]),
        sync_on_foreign=True, flags=flags)
    spec.beacon = dict(airtime_us=airtime_us(spec, NB["beacon_body_octets"]))   # NB MAC reused (3.7, A)
    spec.beacon["reserved_us"] = spec.beacon["airtime_us"] + HBC["psifs_us"]
    return spec


def _bcc(variant, phy, flags):
    """THESIS: the as-built stack as coded (3.8): Walsh-16 chips, header CRC-8, no payload CRC, no ACK, hub TDMA, one
    tone. FULL: the HBC standard stack (3.7 FULL, phy 'sf8'..'sf64') on the same galvanic channel :
    R_c plays the Walsh-chip rate of the SF 8 mode (3.7 last paragraph), so the payload keeps the as-built rate and the
    stack adds the SF 64 header, MAC header and FCS, I-Ack, 3 retransmissions and the NB MAC allocations.
    Flags: rc_cps (default 500 kchip/s), fcs (+2-octet FCS; the as-built stack only, refused for BCC full, which
    carries the 802.15.6 FCS), contact (measured-hardware case: the two-state contact chain at the measured 50 kchip/s,
    the only rate allowed with it), coupling_db (40 or 60: inter-body sensitivity case), power ('low' | 'typical' |
    'high' transceiver case). The one tone of both stacks sits at channel_curves.tone_centre(R_c) (the engine
    takes the path gains at its f -/+ R_c); above BCC_RC_MAX_CPS (9 Mchip/s): refused.
    MULTITONE (full, phy 'multitone'): BCC full on the N tones of the plan (bcc_fdma.plan: spacing 4 R_c, guard
    0, over the 3 dB span of the cohort 50 Ohm curve, 6.48 to 20 MHz with the 7.48 to 9.49 MHz dip; no tone above
    3.38 Mchip/s: refused), a frame striped over the tones (tone_octets), drive amplitude split by sqrt(N) (power
    limited). No contact chain or inter-body coupling per tone: those flags are refused.
    opt (full only, 0..5, default 0): the galvanic-optimisation ladder (DESIGN 4.3), cumulative: 1 RI mode and the
    '128 B' pilot period, 2 the beacon acknowledges the superframe's frames (no I-Ack), 3 the retry pool instead of the
    spare exchanges, 4 the 1-octet MAC header, 5 one preamble, SFD and header per stripe (multi-tone; one tone: as
    4). pool: the share of the superframe's frames the retry pool holds (default POOL_SHARE from opt 3); given below
    opt 3 it replaces the spare exchange per frame by the hub-polled pool, the reference stack of (the
    standard's improvised access, DESIGN 4.3). fill: the pool also takes every exchange that fits between the allocations' end and the next beacon;
    the schedule is admitted on the pool share alone, so capacity stays that of the share.
    touch: the share of time that the bodies of each pair of BANs (0-1, 2-3, ..)
    touch and share one galvanic medium; while they touch, the partner's frames reach a receiver at the partner node's
    own-link Ec/N0 less touch_db (default TOUCH_DB). fallback (with touch): the hub that hears a foreign beacon
    interleaves the pair's active periods (802.15.6 6.13.3) from one superframe after the contact starts, when both
    fit in one period (DESIGN 4.3). spread: the payload spreading, SPREADS; 'none'
    sends the PSDU unspread, 'adapt' (from opt 2) switches per link on the beacon's loss report (engine). lmax (capacity
    short frames, single-tone BCC full): the largest payload per frame, octets (default and upper limit 255); longer
    payloads are fragmented into frames of lmax. back (reference, BCC full at L0 or L1 with the pool): the
    standard's block acknowledgement instead of the I-Ack, data frames with L-Ack (data + pSIFS) and one B-Ack frame
    (BACK_OCTETS, the level's framing) pSIFS after the sensor's last frame of its allocation; the pool exchanges keep
    their I-Ack (claim4.back_cost). burst (with back at L0): the block in the HBC burst mode (6.2.9.4, 10.5.2,
    10.11.3): the sensor's first frame of its allocation in DRF mode, every later one pMIFS after the one before in RI
    mode (header at the payload rate, 608-chip SFD/RI), pSIFS after the last, then the B-Ack; spec.ri is the RI-mode
    frame (the default pilot period). span (with burst): the block of span frames spans
    allocations (6.2.9.4: a block may span more than one allocation interval); the B-Ack follows the allocation in
    which the sensor's frames since the last B-Ack reach span. The allocation keeps room for the B-Ack, as the block
    may close in any of them, unless span >= BLOCK_UNLIMITED: then no B-Ack is sent and only the pSIFS after the
    burst stays. pilot128 (at L0): the '128 B' pilot period of 10.7.3 on every frame,
    the lead reference; from L1 it is part of step 1."""
    contact, opt, pool = flags.get("contact", False), flags.get("opt", 0), flags.get("pool")
    if opt not in range(6) or opt and variant != "full":
        raise ValueError("opt is the ladder level 0..5 of BCC full (DESIGN 4.3), not of the as-built stack")
    if pool is not None and not 0 < pool <= 1:
        raise ValueError("pool is the retry-pool share (0, 1] of BCC full (DESIGN 4.3)")
    touch, fallback = flags.get("touch"), flags.get("fallback", False)
    if touch is not None and (variant != "full" or phy == "multitone" or not 0 < touch <= 1):
        raise ValueError("touch is the contact share (0, 1] of single-tone BCC full")
    if (fallback or flags.get("touch_db") is not None) and touch is None:
        raise ValueError("touch_db and fallback go with touch")
    spread = flags.get("spread", "walsh")
    if spread not in SPREADS or spread != "walsh" and (variant != "full" or phy == "multitone"):
        raise ValueError(f"spread is one of {SPREADS}, of single-tone BCC full")
    if spread == "adapt" and opt < 2:
        raise ValueError("spread 'adapt' switches on the beacon's loss report: BCC full from ladder level 2")
    lmax = flags.get("lmax")
    if lmax is not None and (variant != "full" or phy == "multitone" or not 1 <= lmax <= HBC["body_max_octets"]):
        raise ValueError(f"lmax is the frame payload 1..{HBC['body_max_octets']} octets of single-tone BCC full")
    if flags.get("back") and (variant != "full" or opt >= 2 or pool is None):
        raise ValueError("back is the block-ack reference: BCC full at ladder level 0 or 1 with the pool flag")
    if flags.get("burst") and (not flags.get("back") or opt):
        raise ValueError("burst is the block-ack reference in burst mode: with back at ladder level 0")
    if flags.get("span") is not None and (not flags.get("burst") or flags["span"] < 1
                                            or flags["span"] != int(flags["span"])):
        raise ValueError("span is the block size (frames >= 1) of the block-ack reference: with back and burst")
    if flags.get("pilot128") and (variant != "full" or opt):
        raise ValueError("pilot128 is the '128 B' pilot period of BCC full at ladder level 0 (from L1 it is the default)")
    if flags.get("fill") and opt < 3:
        raise ValueError("fill extends the retry pool of BCC full from ladder level 3")
    rc = flags.get("rc_cps", BCC["rc_sweep_cps"][0 if contact else 2])   # 50k with contact, else 500k (highest)
    if contact and rc != BCC["rc_sweep_cps"][0]:
        raise ValueError("the contact case is the 50 kchip/s hardware measurement (PARAMETERS 3.8, case b)")
    fcs, power, coupling = flags.get("fcs", False), flags.get("power", "typical"), flags.get("coupling_db")
    if coupling is not None and coupling not in BCC["inter_body_sensitivity_dB"]:
        raise ValueError(f"coupling_db must be one of {BCC['inter_body_sensitivity_dB']}")
    if rc > BCC_RC_MAX_CPS:
        raise ValueError(f"no stable BFSK tone at R_c {rc:g} chip/s: above {BCC_RC_MAX_CPS:g} its lower tone nears "
                         "0 Hz")
    tone = channel_curves.tone_centre(rc)
    sync_us = BCC["tdma"]["sync_chips"] * 1e6 / rc
    mw, v = EN["bcc_mW"][power], P["common"]["vdd_V"]
    inter = "coupled" if coupling is not None else "isolated"
    energy = _energy(mw["tx"] / v, mw["rx"] / v, EN["bcc_wake_us"][power])
    if variant == "full":   # SNR per Walsh chip (noise_bw = R_c) and the BCC 0 dB chip rule, as the as-built stack
        multi = phy == "multitone"
        if multi and (contact or coupling is not None):
            raise ValueError("multi-tone BCC has no contact chain or inter-body coupling per tone (new physics): run "
                             "--bcc-contact and --bcc-coupling on the single-tone BCC stacks")
        if fcs:
            raise ValueError("BCC full carries the 802.15.6 FCS already: the fcs flag is for the as-built stack")
        s = _hbc("full", "sf8" if multi else phy, {k: True for k, v in (("ri_mode", opt), ("pilot128",
                 flags.get("pilot128"))) if v}, walsh_cps=rc)
        s.pid, s.flags, s.chip_rate_cps, s.p_loss, s.inter_body, s.energy = "bcc", flags, rc, 0.0, inter, energy
        s.channels, s.channel_bw_mhz = [tone / 1e6], 4 * rc / 1e6
        # no C/I credit for the SF 64 header (or SF 16..64 payload): the energy detector integrates a Walsh chip of 8
        # R_c chips, and a same-tone interferer whose 8 chips all sit on the other tone still wins at SIR < 0 dB
        s.phy.update(noise_bw_hz=rc, ci_co_lin=1.0, ci_credit=False)
        s.max_payload_octets = lmax or s.max_payload_octets
        if multi:
            import bcc_fdma     # loads HBCProtocol and PHY2, so only for this entry
            s.tones = bcc_fdma.plan(bcc_fdma.Geometry("coexistence", "cohort", "50", rc, 0.0)).f.tolist()
            if not s.tones:
                raise ValueError(f"no multi-tone plan at R_c {rc:g} chip/s: 4 R_c exceeds the 3 dB span")
            n = len(s.tones)
            s.phy.update(name="multitone", amp_v=s.phy["amp_v"] / math.sqrt(n))
            s.max_payload_octets *= n       # one full-size BCC full packet per tone
            # ponytail: N parallel front ends, each at the single-tone transceiver power in every radio state
            # (conservative: a shared wideband front end may draw less); MCU wake-up and sleep current once. The
            # beacon stays a one-tone packet, sent and heard on one front end (DESIGN 6)
            s.energy = {x: v * n if x in ("tx_ma", "rx_ma", "ramp_tx_ma", "ramp_rx_ma") else v
                        for x, v in energy.items()}
            s.beacon["energy"] = energy
        if flags.get("burst"):
            s.ri = replace(s, phy=dict(s.phy, header_sf=s.phy["sf"]), flags=dict(flags, ri_mode=True))
        return s
    return ProtocolSpec(
        pid="bcc", variant=variant,
        phy=dict(name="bfsk", mod="bfsk_chip_v2", bit_rate_bps=rc, noise_bw_hz=rc, rep=1, fec=None, payload_crc=fcs,
                 ci_co_lin=1.0, nf_db=BCC["nf_dB"][1], l_impl_db=0.0, tx_power_dbm=None, amp_v=BCC["drive_V"],
                 sync_sinr_db=None),   # ci_co 0 dB: a chip under SIR < 0 dB is a guess (3.8 interference rule)
        channels=[tone / 1e6], channel_bw_mhz=4 * rc / 1e6, acr_db=[0.0], access="tdma", hop="fixed",
        beacon=dict(airtime_us=sync_us, reserved_us=sync_us), chip_rate_cps=rc,
        period_grid=range(BLE["ci_min_us"], PERIOD_CAP_US + 1, BLE["ci_step_us"]),   # the plain BLE CI grid
        max_payload_octets=BCC["payload_max_octets"],
        p_loss=0.0,     # noise model; the contact case runs the chain in both reception models (DESIGN 5.2)
        inter_body=inter, energy=energy, sync_on_foreign=True, flags=flags)


_BUILD = {"ble": _ble, "ble4thesis": _ble4thesis, "ieee802154": _ieee802154, "tsch": _tsch,
          "ieee802156nb": _nb, "smartban": _smartban, "ant": _ant, "hbc": _hbc, "bcc": _bcc}


def get_spec(pid, variant, phy=None, **flags):
    """Spec of protocol pid at its full capability (variant 'full'), or of an as-built model (variant 'thesis': bcc, the
    stack as coded, and ble4thesis, whose only variant it is); PHY name (None: default; see catalogue for the names).
    bcc 'full' is the HBC stack on the galvanic channel."""
    if pid not in _BUILD or variant not in VARIANTS.get(pid, ("full",)):
        raise ValueError(f"unknown protocol or variant: {pid!r}, {variant!r}")
    if set(flags) - FLAGS[pid]:
        raise ValueError(f"unknown flags for {pid}: {sorted(set(flags) - FLAGS[pid])}")
    return _BUILD[pid](variant, phy, dict(flags))


def catalogue():
    """(pid, variant, phy) run by default: every rate or PHY each standard offers (BLE 1M, 2M, Coded S=2 and S=8;
    802.15.6 NB 121.4, 242.9, 485.7 and 971.4 kb/s; HBC SF 8, 16, 32 and 64; 802.15.4 CSMA and GTS; SmartBAN BCH and
    uncoded at Nrep 1, 2 and 4; TSCH, ANT), BCC full
    on one tone and on the N tones of the plan, and the as-built BCC stack as coded ('thesis', the as-built
    reference). The report takes the best rate of a stack per cell as its link adaptation."""
    return ([("ble", "full", p) for p in BLE["phy"]]
            + [("ieee802154", "full", "oqpsk"), ("ieee802154", "full", "gts"), ("tsch", "full", "oqpsk")]
            + [("ieee802156nb", "full", p) for p in NB["rates_kbps"]]
            + [("smartban", "full", p) for p in SB_PHYS] + [("ant", "full", "gfsk")]
            + [("hbc", "full", f"sf{sf}") for sf in sorted(HBC["sf_by_rate_kbps"].values())]
            + [("bcc", "full", "sf8"), ("bcc", "full", "multitone"), ("bcc", "thesis", "bfsk")])


# --- frame layouts ------------------------------------------------------------------------------------------------

def _stack(parts):
    """(kind, dur_us, n_bits, info_bits, code[, rep]) in air order -> blocks with t0_us / t1_us."""
    out, t = [], 0.0
    for kind, dur, n, info, code, *rep in parts:
        b = dict(t0_us=t, t1_us=t + dur, kind=kind, n_bits=n, info_bits=info, code=code)
        if rep and rep[0] is not None:
            b["rep"] = rep[0]
        out.append(b)
        t += dur
    return out


def _shares(total, weights):
    """Split the integer total in proportion to weights, exactly."""
    out, acc, prev = [], 0, 0
    for w in weights:
        acc += w
        cur = total * acc // sum(weights)
        out.append(cur - prev)
        prev = cur
    return out


def _layout_ble(spec, L):
    name = spec.phy["name"]
    ph = BLE["phy"][name]
    aa, hdr, crc = BLE["access_address_bits"], BLE["header_bits"], BLE["crc_bits"]
    if not name.startswith("coded"):
        us = 1e6 / spec.phy["bit_rate_bps"]
        pre = ph["airtime_fixed_us"] - (aa + hdr + crc) * us      # 8 us: 1 octet (1M) or 2 octets (2M)
        return _stack([("sync", pre, 0, 0, None), ("sync", aa * us, aa, 0, None), ("header", hdr * us, hdr, 0, None),
                       ("payload", (8 * L + crc) * us, 8 * L + crc, 8 * L, None)])
    b1, term = BLE["coded_block1_bits"], MISSING["coded_term_bits"]     # block 1: AA + CI, always S=8
    us1, us = BLE["phy"]["coded_s8"]["airtime_per_octet_us"] / 8, ph["airtime_per_octet_us"] / 8
    g1, g2 = ("conv", _coded_gain("coded_s8")), ("conv", _coded_gain(name))
    pre = ph["airtime_fixed_us"] - (b1 + term) * us1 - (hdr + crc + term) * us   # 80 us uncoded preamble
    return _stack([("sync", pre, 0, 0, None), ("sync", (b1 + term) * us1, b1, 0, g1), ("header", hdr * us, hdr, 0, g2),
                   ("payload", (8 * L + crc + term) * us, 8 * L + crc, 8 * L, g2)])


def _layout_ble4thesis(spec, L):
    bits = round(8 * L)     # L may be fractional: BLE_SIM payloads are whole bits
    return _stack([("payload", max(1, math.ceil(bits * 1e6 / spec.phy["bit_rate_bps"])), bits, bits, None)])


def _layout_154(spec, L):
    o, phr = B154["octet_us"], MISSING["phr_octets"]
    shr, psdu = B154["shr_phr_octets"] - phr, B154["mhr_octets"] + L + B154["fcs_octets"]
    return _stack([("sync", shr * o, 8 * shr, 0, None), ("header", phr * o, 8 * phr, 0, None),
                   ("payload", psdu * o, 8 * psdu, 8 * L, None)])


def _layout_nb(spec, L):
    ph, us = spec.phy, 1e6 / NB["symbol_rate_Hz"]
    bps = ph["bit_rate_bps"] / NB["symbol_rate_Hz"]
    n_psdu = 8 * (NB["mac_header_octets"] + L + NB["fcs_octets"])
    n, k, t = NB["bch_psdu"]["n"], NB["bch_psdu"]["k"], NB["bch_psdu"]["t"]
    h, s_h = NB["bch_header"], NB["header_spread"]
    coded = n_psdu + (n - k) * math.ceil(n_psdu / k)
    # ponytail: the 971.4 kb/s header is pi/2-DBPSK on air but evaluated with the payload's DQPSK curve
    return _stack([("sync", NB["preamble_bits"] * us, 0, 0, None, 1),    # 90 bits at S = 1 (3.4 preamble)
                   ("header", h["n"] * s_h * us, h["n"], 0, ("bch", h["n"], h["k"], h["t"]), s_h),
                   ("payload", math.ceil(coded * ph["rep"] / bps) * us, coded, 8 * L, ("bch", n, k, t))])


def _layout_smartban(spec, L):
    """PPDU of 3.5: airtime Nrep (104 + L_PSDU) us. ponytail: every block is stretched Nrep times (its decisions spread
    over the copies' airtime, combined ideally by the phy's rep), so an overlap of one copy hits every block in
    proportion instead of that copy's blocks; the collision model decides copy by copy (engine)."""
    us = spec.phy["nrep"] * 1e6 / SB["symbol_rate_Hz"]
    pre, sync, tol = MISSING["smartban_sync_bits"]
    crc, mpdu = MISSING["smartban_phy_crc_bits"], MISSING["smartban_mpdu_overhead_bits"] + 8 * L
    hb, h, fec = SB["header_bits"], SB["bch_header"], spec.phy["fec"]
    # BCH: last codeword shortened to (mpdu mod 113) + 14
    coded = mpdu + (fec[1] - fec[2]) * math.ceil(mpdu / fec[2]) if fec else mpdu
    return _stack([("sync", pre * us, 0, 0, None), ("sync", sync * us, sync, 0, ("bch", sync, 0, tol)),
                   ("header", h["n"] * us, h["n"], 0, ("bch", h["n"], h["k"], h["t"])),
                   ("header", (hb - h["n"]) * us, hb - h["n"], 0, None),      # CRC-4, outside the BCH
                   ("payload", coded * us, coded, 8 * L, fec),
                   ("payload", crc * us, crc, 0, None)])


def _layout_ant(spec, L):
    """A message or burst packet carries 8 octets, an advanced-burst packet 16 or 24 (3.6 recipe: 8 (1 + 5 + L + 2 +
    2) us = 144 / 208 / 272 us); a short payload is padded to the next size."""
    pre, addr, ctrl, crc = MISSING["ant_frame_octets"]
    size = ANT["payload_octets"] * max(1, -(-L // ANT["payload_octets"]))
    pay = 8 * (ctrl + size + crc)    # at 1 Mb/s, 1 us per bit
    return _stack([("sync", 8 * pre, 0, 0, None), ("sync", 8 * addr, 8 * addr, 0, None),
                   ("payload", pay, pay, 8 * L, None)])


def _walsh_layout(counts, n_oct, info_bits, us, rep, plain=False):
    """BCC / HBC frame (3.7, 3.8) in Walsh chips: preamble of four 64-chip copies (the correlator needs the last copy
    and the SFD, <= 7 chip errors each), SFD, 8-symbol header, PSDU of n_oct octets, 2 symbols per octet, carrying
    info_bits application bits, with a 64-chip pilot after every full pilot period. counts = (preamble, SFD, header,
    pilot) in Walsh chips and the pilot period in symbols (preamble 0: the PSDU alone, an L5 stripe's later tones);
    us, rep = (sync, header, payload) chip duration in us and SNR factor (None: phy['rep']); plain: the PSDU
    unspread, 8 chips per octet, one bit per chip, no Walsh decoding (a pilot still after every period's octets)."""
    pre, sfd, hdr, pil, per = counts
    w = BCC["walsh"]["n"]
    q = ("bch", pil, 0, BCC["sync_max_errors_of_64"])
    parts = []
    if pre:
        parts = [("sync", (pre - pil) * us[0], 0, 0, None, rep[0]), ("sync", pil * us[0], pil, 0, q, rep[0]),
                 ("sync", sfd * us[0], sfd, 0, q, rep[0]), ("header", hdr * us[1], hdr, 0, ("walsh16",), rep[1])]
    cpo = BCC["chips_per_octet"]
    per = per * w // cpo    # the pilot period in octets
    sizes = [min(per, n_oct - i) for i in range(0, n_oct, per)]
    for m, info in zip(sizes, _shares(info_bits, sizes)):
        parts.append(("payload", 8 * m * us[2], 8 * m, info, None, rep[2]) if plain
                     else ("payload", m * cpo * us[2], m * cpo, info, ("walsh16",), rep[2]))
        if m == per:
            parts.append(("payload", pil * us[0], pil, 0, q, rep[0]))   # pilot: a miss corrupts the rest
    return _stack(parts)


def _opt(spec):
    """Ladder level of BCC full (DESIGN 4.3); 0 for every other stack."""
    return spec.flags.get("opt", 0)


def _pooled(spec):
    """BCC full with the hub-polled retry pool instead of a spare exchange per frame: from L3, and at L0 to L2 when
    the pool flag is given."""
    return _opt(spec) >= 3 or spec.flags.get("pool") is not None


def _hbc_framing(spec):
    """(counts, us, rep, MAC header octets) of the HBC frame for _walsh_layout (3.7); BCC full from L1 the RI mode with
    its 608-chip SFD/RI field (10.5.2, Eq. 117; 3.7a) and a pilot every 256 symbols (4096 Walsh chips, the HBC pilot
    option '128 B' of 3.7; at L0 with the flag pilot128, the lead reference), from L4 the
    lean MAC header."""
    ph, s, o = spec.phy, MISSING["hbc_sync_sf"], _opt(spec)
    ri = o or spec.flags.get("ri_mode")             # HBC's own RI flag; BCC full from L1 (its beacon is built with it)
    chip = 1e6 / (ph["bit_rate_bps"] * ph["sf"])     # the 42 Mcps chip; BCC full: 1 / (8 R_c)
    counts = (HBC["preamble_chips"] // s, HBC["sfd_ri_chips" if ri else "sfd_chips"] // s,
              HBC["header_chips"] // HBC["sf_by_rate_kbps"]["164.0625"],   # header chips at SF 64
              HBC["pilot_chips"] // s, HBC["pilot_period_symbols"] * (2 if o or spec.flags.get("pilot128") else 1))
    sf, ssf, hsf = ph["sf"], ph["sync_sf"], ph["header_sf"]
    return (counts, (ssf * chip, hsf * chip, sf * chip), (ssf / sf, hsf / sf, None),
            LEAN_MAC_OCTETS if o >= 4 else NB["mac_header_octets"])


def _layout_hbc(spec, L, plain=None):
    """plain None: the spec's own spreading; 'adapt' specs are sized as Walsh."""
    counts, us, rep, mac = _hbc_framing(spec)
    plain = spec.flags.get("spread") == "none" if plain is None else plain
    return _walsh_layout(counts, mac + L + NB["fcs_octets"], 8 * L, us, rep, plain)


@lru_cache(maxsize=None)
def _balance(n, octets, head_us, octet_us, pilot_us, pilot_octets):
    """PSDU octets per tone of an L5 stripe: each octet goes to the tone where it ends earliest (ties: the lower tone),
    tone 0 starting after head_us of preamble, SFD and header, a pilot after every pilot_octets. The per-tone end times
    grow with the octets, so this takes the octets smallest end times and minimises the stripe airtime."""
    end = lambda i, k: head_us * (i == 0) + k * octet_us + pilot_us * (k // pilot_octets)   # noqa: E731
    x, heap = [0] * n, [(end(i, 1), i) for i in range(n)]
    heapq.heapify(heap)
    for _ in range(octets):
        i = heapq.heappop(heap)[1]
        x[i] += 1
        heapq.heappush(heap, (end(i, x[i] + 1), i))
    return tuple(x)


def _stripe_layouts(spec, L):
    """L5 (proposal, DESIGN 4.3): one preamble, SFD and header per stripe, on tone 0; the stripe's PSDU (lean MAC
    header, L octets, FCS) split over the tones by _balance, the other tones' packets PSDU only from the stripe start;
    a tone left without octets stays idle (the last ones: equal tones fill in order). Application bits in proportion
    to the octets (ponytail: they only weight post_ber)."""
    counts, us, rep, mac = _hbc_framing(spec)
    pre, sfd, hdr, pil, per = counts
    octs = _balance(len(spec.tones), mac + L + NB["fcs_octets"], (pre + sfd) * us[0] + hdr * us[1],
                    BCC["chips_per_octet"] * us[2], pil * us[0], per * BCC["walsh"]["n"] // BCC["chips_per_octet"])
    return [_walsh_layout(counts if i == 0 else (0,) + counts[1:], x, info, us, rep)
            for i, (x, info) in enumerate(zip(octs, _shares(8 * L, octs))) if i == 0 or x]


def _layout_bcc(spec, L):
    c = 1e6 / spec.chip_rate_cps
    fcs = BCC["fcs_switch_chips"] // BCC["chips_per_octet"] if spec.phy["payload_crc"] else 0
    counts = (BCC["preamble_chips"], BCC["sfd_chips"], BCC["header_chips"], BCC["pilot_chips"],
              BCC["pilot_period_chips"] // BCC["walsh"]["n"])
    return _walsh_layout(counts, L + fcs, 8 * L, (c, c, c), (None, None, None))


_LAYOUT = {"ble": _layout_ble, "ble4thesis": _layout_ble4thesis, "ieee802154": _layout_154, "tsch": _layout_154,
           "ieee802156nb": _layout_nb, "smartban": _layout_smartban, "ant": _layout_ant, "hbc": _layout_hbc,
           "bcc": _layout_bcc}


def _recipe(spec):
    """Frame and exchange recipe of a spec: BCC full runs the HBC stack."""
    return "hbc" if spec.pid == "bcc" and spec.variant == "full" else spec.pid


def tone_octets(spec, payload_octets):
    """Application octets of each tone packet of a frame: multi-tone BCC stripes it over its N tones, ceil(L / N) on
    the first L mod N tones and floor(L / N) on the others (a packet of MAC header and FCS only when L < N, so every
    tone is active); every other stack sends one packet."""
    if not spec.tones:
        return [payload_octets]
    q, r = divmod(payload_octets, len(spec.tones))
    return [q + 1] * r + [q] * (len(spec.tones) - r)


def tone_layouts(spec, payload_octets, plain=None):
    """Layouts of the tone packets of one frame (tone_octets order: the longest first), sent in parallel; an L5
    stripe: tone 0 with the preamble, SFD and header (not always the longest), then the tones that carry PSDU octets
    (_stripe_layouts).
    plain True / False: the frame unspread or Walsh-spread, whatever the spec's flag."""
    if plain is not None:
        return [_layout_hbc(spec, payload_octets, plain)]
    if spec.tones and _opt(spec) >= 5:
        return _stripe_layouts(spec, payload_octets)
    return [_LAYOUT[_recipe(spec)](spec, o) for o in tone_octets(spec, payload_octets)]


def frame_layout(spec, payload_octets, plain=None):
    """Blocks of one frame carrying payload_octets application octets (phy.frame_outcome layout); multi-tone BCC: its
    longest tone packet, which sets the airtime (tone_layouts gives every tone's); plain as tone_layouts."""
    if plain is not None:
        return _layout_hbc(spec, payload_octets, plain)
    if spec.tones and _opt(spec) >= 5:
        return max(_stripe_layouts(spec, payload_octets), key=lambda lay: lay[-1]["t1_us"])
    return _LAYOUT[_recipe(spec)](spec, tone_octets(spec, payload_octets)[0])


def airtime_us(spec, payload_octets, plain=None):
    """Frame airtime in whole us (rounded up)."""
    return math.ceil(frame_layout(spec, payload_octets, plain)[-1]["t1_us"] - 1e-6)


# --- access exchange and period choice -----------------------------------------------------------------------------

def exchange(spec, payload_octets, pool=False, cont=False):
    """One access exchange of a sensor frame, as PARAMETERS prescribes (the data step is the frame airtime); pool: an
    exchange of the retry pool, which keeps the I-Ack under the flag back; cont (flag burst only): a frame after the
    first of the sensor's allocation, pMIFS after the frame before, in RI mode (step kind 'burst')."""
    if cont and spec.flags.get("burst"):
        return [Step("sensor", "ifs", HBC["pmifs_us"], False),
                Step("sensor", "burst", airtime_us(spec.ri, payload_octets), True)]
    data = Step("sensor", "data", airtime_us(spec, payload_octets), True)
    pid = _recipe(spec)
    if pid == "ble":                # ACL: hub poll (empty PDU), T_IFS, data, T_IFS; the ACK rides the next poll
        t = BLE["t_ifs_us"]
        return [Step("hub", "poll", airtime_us(spec, 0), True), Step("sensor", "ifs", t, False), data,
                Step("hub", "ifs", t, False)]
    if pid in ("ieee802154", "tsch"):
        if pid == "tsch":           # one timeslot: CCA at macTsCcaOffset, data at macTsTxOffset, Enh-Ack after it
            t = spec.tsch
            steps = [Step("sensor", "ifs", t["cca_offset_us"], False), Step("sensor", "cca", t["cca_us"], False),
                     Step("sensor", "ifs", t["tx_offset_us"] - t["cca_offset_us"] - t["cca_us"], False), data,
                     Step("hub", "ifs", t["tx_ack_delay_us"], False), Step("hub", "ack", t["enh_ack_us"], True)]
            return steps + [Step("sensor", "ifs", t["timeslot_us"] - sum(s.dur_us for s in steps), False)]
        mpdu = B154["mhr_octets"] + payload_octets + B154["fcs_octets"]
        ifs = CSMA["lifs_us"] if mpdu > CSMA["max_sifs_frame_octets"] else CSMA["sifs_us"]
        head = ([Step("sensor", "cca", CSMA["cca_us"], False), Step("sensor", "ifs", CSMA["turnaround_us"], False)]
                if spec.access == "csma" else [])     # a GTS: no CSMA in the CFP
        return head + [data, Step("hub", "ifs", CSMA["ack_delay_us"], False),
                       Step("hub", "ack", B154["imm_ack_us"], True), Step("sensor", "ifs", ifs, False)]
    if pid in ("ieee802156nb", "hbc"):  # data, pSIFS, I-Ack (the L = 0 frame at the data rate, A), pSIFS
        t = (NB if pid == "ieee802156nb" else HBC)["psifs_us"]
        # from L2 the next beacon acknowledges (DESIGN 4.3), under back a B-Ack per allocation: data and pSIFS
        if _opt(spec) < 2 and (pool or not spec.flags.get("back")):
            return [data, Step("hub", "ifs", t, False), Step("hub", "ack", airtime_us(spec, 0), True),
                    Step("sensor", "ifs", t, False)]
        if spec.flags.get("burst"):     # the burst's first frame; the pSIFS after the last is in back_us
            return [data]
        return [data, Step("sensor", "ifs", t, False)]
    if pid == "smartban":           # data, TIFS, ACK, TIFS
        t = SB["tifs_us"]
        return [data, Step("hub", "ifs", t, False), Step("hub", "ack", SB["ack_us"], True),
                Step("sensor", "ifs", t, False)]
    if pid == "ant":                # master TX 0..144, turnaround, slave ACK 274..418
        return [data, Step("hub", "ifs", ANT["turnaround_us"], False), Step("hub", "ack", ANT["ack_airtime_us"], True)]
    if pid == "bcc":                # MAC_layer.TDMA slot: frame + guard chips
        return [data, Step("sensor", "ifs", BCC["tdma"]["slot_guard_chips"] * 1e6 / spec.chip_rate_cps, False)]
    return [data]                   # ble4thesis


def _ex_us(spec, octets, pool=False, cont=False):
    return sum(st.dur_us for st in exchange(spec, octets, pool, cont))


def back_us(spec):
    """Allocation time of the block-ack reference's B-Ack (flag back): the frame and its pSIFS, under burst also the
    pSIFS before it; else 0. Under span >= BLOCK_UNLIMITED the pSIFS after the burst alone."""
    if not spec.flags.get("back"):
        return 0
    if spec.flags.get("span", 0) >= BLOCK_UNLIMITED:
        return HBC["psifs_us"]
    return airtime_us(spec, BACK_OCTETS) + HBC["psifs_us"] * (2 if spec.flags.get("burst") else 1)


def _frame_us(spec, octets, cont=False):
    """Schedule time one frame needs per period: its exchange, x (1 + SPARE_EXCHANGES) for the stacks with ARQ that
    retry in the same allocation or event: the 802.15.6 allocations (NB, HBC, BCC full; the author's rule of 2026-09-24
    replaces 3.4 'sized for (R + 1) full exchanges'), the 802.15.4 GTS (3.2 mode B 'retry inside the GTS if it fits',
    the same rule, A) and the BLE ACL connection event (a failed data PDU is resent in the same event, 3.1 FULL
    step 3); later retries go to the next allocation, GTS or event. x 1.25 for ANT
    acknowledged messages (3.6 FULL step 2). BCC full with the pool (_pooled): one exchange, the superframe's retry
    pool instead. cont: as exchange."""
    x = _ex_us(spec, octets, False, cont)
    if (_recipe(spec) in ("ieee802156nb", "hbc", "ble") or spec.phy["name"] == "gts") and not _pooled(spec):
        x *= 1 + SPARE_EXCHANGES
    return x * (spec.ant["ch_rate_factor"] if spec.ant else 1)


def pool_size(frames, share=POOL_SHARE):
    """L3 retry pool (DESIGN 4.3): max(1, ceil(share x the superframe's frames)) hub-polled exchanges."""
    return max(1, math.ceil(round(frames * share, 9)))   # round: 0.05 x 60 is 3.0000000000000004 in floats


def required_ec_n0_db(plain):
    """: Ec/N0 (dB) at which the largest PSDU (lean MAC header, 255 octets, FCS) is delivered as often as the
    Walsh-16 PSDU at the as-built error-free 8 dB (A95_WEAKEST): 8 dB for Walsh, for unspread BFSK the Ec/N0 at which
    all its 8 n bits are right with that probability (0.5 exp(-g / 2) per chip). Both sides stay on the unshifted curve,
    not the engine's 'bfsk_chip_v2': the rule compares two frames on one curve, and with the shift it gives
    13.09 dB for 13.02 and 62.8 Mchip/s for 63.9, far above the 9 Mchip/s cap that binds at the 1 V drive."""
    import phy
    e0 = A95_WEAKEST["ec_n0_db"]
    if not plain:
        return e0
    n = LEAN_MAC_OCTETS + BCC["payload_max_octets"] + NB["fcs_octets"]
    psr = phy.walsh16_symbol_ok_prob(phy.ber("bfsk_chip", 10 ** (e0 / 10))) ** (2 * n)
    return 10 * math.log10(-2 * math.log(2 * -math.expm1(math.log(psr) / (8 * n))))


def margin_rate(spread, nf_db, amp_v=BCC["drive_V"]):
    """ rate rule (DESIGN 4.3): the highest chip rate at which the 95th-percentile participant on the weakest
    path of the link budget meets required_ec_n0_db at drive amp_v and noise figure nf_db, capped at BCC_RC_MAX_CPS
    (the waveform limit); 'adapt' runs unspread while its link is clean. Ec/N0 = (A |H|)^2 / (2 R_c N0 F), so from
    the table row R_c = rc0 (A / A95)^2 10^((nf0 - NF) / 10) 10^((8 - required) / 10). Returns (rate, uncapped rate)."""
    t = A95_WEAKEST
    r = (t["rc"] * (amp_v / t["amp_v"]) ** 2 * 10 ** ((t["nf_db"] - nf_db) / 10)
         * 10 ** ((t["ec_n0_db"] - required_ec_n0_db(spread != "walsh")) / 10))
    return min(r, BCC_RC_MAX_CPS), r


def poll_octets(spec):
    """Payload octets of an L3 pool poll (DESIGN 4.3): none at L3, whose 7-octet header names the sensor and, in its
    sequence-number field, the polled frame (A); from L4 SEQ_OCTETS, the polled frame's 15-bit sequence number after
    the 2-octet header, which holds the frame type and the sensor's NID."""
    return SEQ_OCTETS if _opt(spec) >= 4 else 0


def _pool_us(spec, octets):
    """One exchange of the L3 retry pool: hub poll, pSIFS, the sensor's frame, pSIFS (802.15.6 improvised access: the
    polled allocation starts pSIFS after the poll)."""
    return airtime_us(spec, poll_octets(spec)) + HBC["psifs_us"] + _ex_us(spec, octets, True)


def beacon_layout(spec, n_dev, n):
    """One-tone layout of the BCC full beacon (the NB 17-octet body, A) at its ladder level: from L2 with the ack field
    of n_dev sensors sending n frames per superframe, per sensor a base sequence number (1 octet, from L4 the 15-bit
    number in SEQ_OCTETS: the hub holds every earlier frame of the sensor) and one bit for each of the n sequence
    numbers after it, base x n_dev + ceil(n_dev n / 8) octets; from L4 it keeps the 7-octet MAC header, which names the
    hub and the BAN (5 octets over the lean one).
    ponytail: n bits per sensor; the engine acknowledges every frame the hub holds, also one past them."""
    ack = (SEQ_OCTETS if _opt(spec) >= 4 else 1) * n_dev - (-n_dev * n // 8) if _opt(spec) >= 2 else 0
    return _layout_hbc(spec, NB["beacon_body_octets"] + ack
                       + (NB["mac_header_octets"] - LEAN_MAC_OCTETS if _opt(spec) >= 4 else 0))


def _sb_slot(spec, octets):
    """SmartBAN slot TS (us) of a frame: LSLOT = the smallest of 1, 2, 4, .., 32 with Nrep (104 + L_PSDU) <= TS - TMUA
    - 168 - 300 (3.5; TMUA = Tsch 16 us, the MCA sensing; the ACK is not repeated), None above 32."""
    room = spec.csma["cca_us"] + SB["ack_us"] + 2 * SB["tifs_us"] + airtime_us(spec, octets)
    return next((n * SB["tmin_us"] for n in SB["lslot"] if room <= n * SB["tmin_us"]), None)


def _fit(bits, t, n_dev, spec):
    """Frames of one period of L = bits and their schedule at period t (DESIGN 3), or None. Protocols with quantised
    slots add slot_us (one sensor's slot) and reserved_us (start of sensor 1's slot in the superframe): TSCH one cell
    per frame (3.3), SmartBAN one LSLOT slot per frame after the D-Beacon slot and before the 2 C/M slots, T a
    multiple of TS (3.5), 802.15.4 GTSs of whole slots (SD / 16) at the end of the SD, the smallest SO leaving the
    440-symbol CAP (3.2 mode B)."""
    lmax = spec.max_payload_octets
    if spec.pid == "ble4thesis":    # BLE_SIM: one PDU of the exact bit count, no fragmentation
        if bits > 8 * lmax:
            return None
        n, last = 1, bits / 8
    else:
        octets = -(-bits // 8)
        n = -(-octets // lmax)
        last = octets - lmax * (n - 1)
    frames = [lmax] * (n - 1) + [last]
    sel = dict(period_us=t, payload_bits=bits, frames_per_period=n, frame_payload_octets=frames,
               r_actual_kbps=bits * 1000.0 / t)
    if spec.access == "ant":        # the hub serves all masters: their slot time must fit one period
        if octets > ANT["payload_octets"]:      # advanced burst: one packet (+ ACK) per spacing, 60 kb/s at most
            sel["burst"] = True
            need = n * spec.ant["burst_spacing_us"]
        elif t / spec.ant["ch_rate_factor"] < spec.period_grid[0]:
            return None     # acknowledged messages at T / 1.25: the channel period stays >= N 171 (3.6, 12.3 / 1.25)
        else:
            need = _frame_us(spec, last)
        return sel if n_dev * need <= t else None
    need = sum(_frame_us(spec, o, i > 0) for i, o in enumerate(frames)) + back_us(spec)
    if spec.access == "csma":
        return sel if n_dev * need <= spec.csma["max_load"] * t else None
    if spec.pid == "smartban":
        ts = _sb_slot(spec, frames[0])
        if ts is None or t % ts or (1 + n_dev * n + SB["cm_slots"]) * ts > t:
            return None
        return dict(sel, slot_us=n * ts, reserved_us=ts)
    if spec.pid == "ieee802154":    # beacon mode B
        bc = CSMA["beacon"]
        for so in range(spec.period_grid.index(t) + 1):
            sd = bc["base_superframe_us"] << so
            gts = math.ceil(need * bc["slots"] / sd) * sd // bc["slots"]
            if n_dev * gts <= sd - bc["min_cap_symbols"] * B154["symbol_us"]:
                return dict(sel, slot_us=gts, reserved_us=sd - n_dev * gts)
        return None
    if _opt(spec) >= 2 or _pooled(spec):   # BCC full ladder (DESIGN 4.3): from L2 the beacon carries the ack field;
        # the retry pool (from L3, or the pool flag) follows the allocations of one exchange per frame (slots packed
        # from the beacon on)
        b_us = math.ceil(beacon_layout(spec, n_dev, n)[-1]["t1_us"] - 1e-6)
        pool = pool_size(n_dev * n, spec.flags.get("pool", POOL_SHARE)) if _pooled(spec) else 0
        p_us = _pool_us(spec, frames[0]) if pool else 0
        if n_dev * need + pool * p_us > t - b_us - HBC["psifs_us"]:
            return None
        sel.update(beacon_us=b_us, reserved_us=b_us + HBC["psifs_us"])
        if pool and spec.flags.get("fill"):     # the slack before the next beacon joins the pool
            pool = max(pool, int((t - b_us - HBC["psifs_us"] - n_dev * need) // p_us))
        return dict(sel, slot_us=need, pool=pool, pool_us=p_us) if pool else sel
    if n_dev * need > t - (spec.beacon["reserved_us"] if spec.beacon else 0):
        return None
    return dict(sel, slot_us=n * spec.tsch["timeslot_us"]) if spec.tsch else sel


def reference_period(r_kbps, mode):
    """Period P of a cell, the same for every protocol (PARAMETERS OI-4: payload and period from the BLE_SIM axis
    R = L / CI; OI-19 and 3.2 mode B: every protocol runs at P rounded to its grid), on the BLE CI grid: latency =
    the shortest CI carrying L >= 8 bits (BLE_SIM's 7.5 ms from 1.07 kb/s up), energy = the longest CI with L within
    one 251-octet PDU, at least 7.5 ms, at most 4 s."""
    g = range(BLE["ci_min_us"], BLE["ci_max_us"] + 1, BLE["ci_step_us"])
    if mode == "latency":
        return next((t for t in g if round(r_kbps * t / 1000) >= 8), g[-1])
    t = int(8000.0 * BLE["payload_max_octets"] / r_kbps) // g.step * g.step
    return min(max(t, g[0]), g[-1])


def select_period(r_kbps, mode, n_dev, spec):
    """Period and payload for per-sensor rate r_kbps (DESIGN 3), or None when infeasible.

    The shortest period of the protocol's grid at or above the cell's reference period P (reference_period) whose
    schedule holds the N_dev sensors (_fit), scanned upward; L = round(R T) bits (at least 8), carried in whole
    octets, fragmented into frames of max_payload_octets (ANT: advanced-burst packets, sel['burst']). deadline_us = P,
    the on-time deadline of every protocol in the cell. None also when N_dev exceeds the protocol's cap. ble4thesis
    keeps BLE_SIM.select_ci exactly: CI 7.5 ms (latency) or the nearest-quantised longest CI (energy), no scan, no
    fragmentation."""
    if mode not in ("latency", "energy"):
        raise ValueError(f"unknown mode {mode!r}")
    if r_kbps <= 0 or (spec.max_dev is not None and n_dev > spec.max_dev):
        return None
    g, lmax = spec.period_grid, spec.max_payload_octets
    if spec.pid == "ble4thesis":
        t = g[0] if mode == "latency" else min(max(round(int(8000.0 * lmax / r_kbps) / g.step) * g.step, g[0]), g[-1])
        bits = round(r_kbps * t / 1000)
        sel = _fit(bits, t, n_dev, spec) if bits >= 8 else None
        return sel and dict(sel, deadline_us=t)
    # every exchange costs at least its share of a full frame's, so N_dev sensors need at least
    # T x n_dev x R x t(L_max) / (8 L_max): above the medium's capacity no period fits (skips the scan); under burst
    # t(L_max) of a later frame of the burst, the cheapest
    share = n_dev * r_kbps / 1000 * _frame_us(spec, lmax, True) / (8 * lmax)
    if share > (spec.csma["max_load"] if spec.access == "csma" else 1.0) * (1 + 1e-9):
        return None
    p = reference_period(r_kbps, mode)
    for t in g[bisect.bisect_left(g, p):]:
        bits = round(r_kbps * t / 1000)
        sel = _fit(bits, t, n_dev, spec) if bits >= 8 else None
        if sel:
            return dict(sel, deadline_us=p)
    return None
