"""Discrete-event engine of the coexistence simulator: k body-area networks of one protocol and variant.

Purpose: DESIGN.md sections 3 to 7. One hub and n_dev sensors per body, periodic uplink traffic (select_period), three
access families: scheduled (TDMA superframes with hub beacon, BLE connection events, TSCH cells; FULL retries in the
TSCH shared cell and the SmartBAN C/M slots; the BCC full ladder, DESIGN 4.3: from L2 the beacon acknowledges and a
retry goes in the next allocation, from L3 in the superframe's retry pool), unslotted CSMA/CA (802.15.4 full) and ANT
periodic masters; channel rules perm / csa2 / tsch / fixed with AFH, TSCH blacklisting, ED scan, SmartBAN channel
switching and ANT adaptive isochronous re-timing; two reception models:
- 'collision': any overlap with another transmission on the same channel, or a receiver that transmits meanwhile,
  destroys the frame; a survivor is lost with spec.p_loss; BER: overlapped bits at 0.5, the rest at the BER that gives
  p_loss for the frame (expected values);
- 'sinr': per-segment SINR at the receiver from channel.py (room, on-body and body-to-body gains, per-packet fade
  redrawn after the coherence time, adjacent channels weighted by the ACR), phy.frame_outcome per frame (sync
  threshold, C/I floor, codewords), half duplex, sync lock on an earlier co-channel frame (sync_on_foreign), BCC
  galvanic isolation or coupling and the two-state contact chain, multi-tone BCC (every tone packet of a striped frame
  at its own Ec/N0, the frame needs them all), HBC coupling between bodies; ci_floor=False drops the C/I floor for
  every stack (interference as Gaussian noise, a sensitivity run).
Energy per node and state (TX, RX, ramps, wake-ups, sleep) and the DESIGN 7 metrics per BAN.
Mechanics: every link (a sensor and its hub side) and every hub beacon is a generator; it yields a Tx (start >= now,
registered at once in per-channel and per-node tracks sorted by start) and resumes at the Tx end with the outcome
decided there, when every overlapping transmission is known; or it yields ('cca' | 'sense', t0, t1, key, node) and
resumes at t1 with busy True/False, or ('wait', t, t) and resumes at t. One heap orders the resumptions. A frame counts
when its first access opportunity falls before duration_us; one still queued and undelivered at the end (waiting for a
retry or behind a backlog) counts as not delivered (drop_end).
Inputs:  protocols.ProtocolSpec, the grid point, reception model, duration, seed, venue, density, case, ci_floor.
Outputs: simulate_trial(...) -> one dict per BAN with the keys COLUMNS (DESIGN 7 plus venue, density, case).
argv:    none (library); sweep.py runs it, test_coexistence.py checks it.
"""
from __future__ import annotations
import heapq
import math
from bisect import bisect_left, bisect_right
from collections import defaultdict, deque
from itertools import count

import numpy as np

import channel as chn
import phy
import protocols as pr

COLUMNS = ("protocol variant phy reception strip_hopping ci_mode period_us r_target_kbps r_actual_kbps r_agg_kbps "
           "n_dev k payload_bits frames_per_period duration_us trial seed network_id feasible generated tx_attempts "
           "delivered_ok delivered_corrupt drop_collision drop_noise drop_sync drop_csma drop_retry drop_queue "
           "drop_end psr pdr per_first retx_per_frame lat_mean_ms lat_p50_ms lat_p95_ms lat_p99_ms lat_max_ms on_time "
           "raw_ber post_ber goodput_kbps airtime_frac e_sensor_j e_hub_j p_sensor_uw jpb_sensor jpb_ban lifetime_days "
           "collision_frac cca_busy_frac jain_fairness channel_occupancy venue density case").split()
COHERENCE_US = pr.P["common"]["fade_redraw_after_us"]   # per-packet fade redrawn after 20 ms (PARAMETERS 5.1)
CONTACT = pr.BCC["contact_model"]
BATT = pr.EN["battery"]
# buck converter: the cell current rises by 3.0 / V_cell late in life (PARAMETERS 4.4, 5-7 % at 2.8 V); ponytail:
# the 2.8 V end applied over the whole life (conservative), a discharge curve would give the mean factor
BUCK = pr.P["common"]["vdd_V"] / 2.8
ED_DBM = pr.CSMA["ed_threshold_dBm"]                    # -75 dBm, 802.15.4 and SmartBAN CCA
REV8 = [int(f"{i:08b}"[::-1], 2) for i in range(256)]
# SmartBAN C/M retries by the MCA backoff (TS 103 325 7.3.3.1.1; BEmin / BEmax not in the TS: 802.15.4 defaults, A)
BE_MIN, BE_MAX = pr.SB["be_min"], pr.SB["be_max"]


def spread_step(plain, run, ok):
    """ per-link spreading switch of BCC full 'adapt', one reported frame: run counts consecutive delivered (> 0)
    or lost (< 0) frames; unspread (plain) turns to Walsh at ADAPT_DOWN losses, Walsh back at ADAPT_UP deliveries.
    Returns (plain, run)."""
    run = max(run, 0) + 1 if ok else min(run, 0) - 1
    if plain and run <= -pr.ADAPT_DOWN or not plain and run >= pr.ADAPT_UP:
        return not plain, 0
    return plain, run


def csa2(counter, chan_id, used):
    """BLE CSA#2 data channel index (Core 5.4 Vol 6 Part B 4.5.8.3) of event `counter` over the ascending used map."""
    prn = counter ^ chan_id
    for _ in range(3):
        prn = (17 * (REV8[prn >> 8] << 8 | REV8[prn & 255]) + chan_id) & 0xFFFF
    prn ^= chan_id
    ch = prn % 37
    return ch if len(used) == 37 or ch in used else used[len(used) * prn >> 16]


class Tx:
    """One transmission; ok / reason / overlap / corrupt and the payload error counts are set at its end."""
    __slots__ = ("s", "e", "key", "node", "rx", "ban", "lay", "data", "ok", "reason", "overlap", "corrupt",
                 "raw", "post", "coded", "info", "w")

    def __init__(self, s, dur, key, node, rx, ban, lay, data):
        self.s, self.e, self.key, self.node, self.rx, self.ban, self.lay, self.data = (s, s + dur, key, node, rx,
                                                                                      ban, lay, data)
        self.ok, self.reason, self.overlap, self.corrupt = True, None, False, False
        self.raw = self.post = self.coded = self.info = 0.0
        self.w = s      # start of the receiver's listening for this frame (Sim._tx)


class Track:
    """Transmissions sorted by start; over(a, b) = those with start < b and end > a (maxd bounds the scan)."""
    __slots__ = ("s", "x")

    def __init__(self):
        self.s, self.x = [], []

    def add(self, tx):
        i = bisect_right(self.s, tx.s)
        self.s.insert(i, tx.s)
        self.x.insert(i, tx)

    def over(self, a, b, maxd):
        out, s, x = [], self.s, self.x
        i, lo = bisect_left(s, b), a - maxd
        while i:
            i -= 1
            if s[i] <= lo:
                break
            if x[i].e > a:
                out.append(x[i])
        return out


class Frame:
    __slots__ = ("gen", "octets", "bits", "tries", "delivered")

    def __init__(self, gen, octets, bits):
        # octets on air (airtime, BER); bits = its share of the L application bits (goodput, J per bit)
        self.gen, self.octets, self.bits, self.tries, self.delivered = gen, octets, bits, 0, False


class Adapt:
    """F5 rule (DESIGN 4.1): per channel the last `window` outcomes; a used channel with more than threshold x window
    failures among them is removed; a removal below min_channels resets the map to all channels (BLE AFH) or is
    refused (TSCH keeps >= min_channels). ponytail: a young window needs the full-window failure count."""

    def __init__(self, channels, a, reset):
        self.all, self.used, self.reset = list(channels), list(channels), reset
        # round: 1 - 0.8 is 0.19999999999999996 in floats, which removed a TSCH channel at PDR 0.8, not below it
        self.lim, self.min = round(a["threshold"] * a["window"], 9), a["min_channels"]
        self.hist = defaultdict(lambda: deque(maxlen=a["window"]))

    def record(self, ch, fail):
        h = self.hist[ch]
        h.append(fail)
        if sum(h) <= self.lim or ch not in self.used:
            return False
        if len(self.used) > self.min:
            self.used.remove(ch)
        elif self.reset:
            self.used = list(self.all)
            self.hist.clear()
        else:
            return False
        return True


class Sim:
    """One trial: k BANs of spec at the selection sel (protocols.select_period)."""

    def __init__(self, spec, sel, k, n_dev, reception, duration_us, seed, venue="ward", density="typical",
                 case="typical", ci_floor=True):
        sp = self.sp = spec
        self.ci_floor = ci_floor    # False: no C/I_co floor, interference as Gaussian noise (DESIGN 5.2)
        self.sel, self.k, self.n_dev, self.m, self.n = sel, k, n_dev, n_dev + 1, k * (n_dev + 1)
        self.T, self.D, self.case, self.sinr = sel["period_us"], duration_us, case, reception == "sinr"
        rng = self.rng = np.random.default_rng(seed)
        # the room, the link gains, the superframe and sample-clock phases and the clock errors come from their own
        # substreams in a fixed number of draws, so every protocol and variant of a cell sees the same ones for one seed
        self.room_rng, cell = np.random.default_rng([seed, 1]), np.random.default_rng([seed, 2])
        self.heap, self.seq, self.maxd = [], 0, 0.0
        self.nch, self.iso = len(sp.channels), sp.inter_body == "isolated"
        self.contact = sp.pid == "bcc" and sp.flags.get("contact", False)   # the chain decides in both models
        self.tracks, self.ntrack, self.atrack = {}, [Track() for _ in range(self.n)], Track()
        # linear ACR factor between channel indices (0: ignored), from the centre frequencies (PARAMETERS 5.5)
        self.acrlin = [[10 ** (-a / 10) for a in row] for row in pr.acr_matrix(sp)]
        self.c = [defaultdict(float) for _ in range(k)]
        self.lat = [[] for _ in range(k)]
        # t_tx, t_rx, TX ramps, RX ramps, wake-ups, then t_tx, t_rx, TX ramps, RX ramps of a beacon with its own
        # currents (spec.beacon['energy']: the one-tone beacon of multi-tone BCC on one front end, DESIGN 6)
        self.en = [[0.0] * 9 for _ in range(self.n)]
        self.queues = [[] for _ in range(k)]
        fo = sel["frame_payload_octets"]    # (octets, application bits) per frame: the last one holds the rest of L
        self.plan = list(zip(fo, [8 * o for o in fo[:-1]] + [sel["payload_bits"] - 8 * sum(fo[:-1])]))
        self.lays, self.clean, self.fades, self.cstate = {}, {}, {}, {}
        drift = sp.drift_ppm * 1e-6
        self.eps = drift * (2 * cell.random(self.n) - 1)
        # receive window: widening of two clocks over one period plus 2 x 16 us jitter (PARAMETERS 4.1)
        self.guard = 4 * sp.drift_ppm * 1e-6 * self.T + 2 * pr.BLE["jitter_us"]
        self.phase = np.floor(cell.random(k) * self.T)
        # sample clock of each sensor: period T, random phase (DESIGN 3); a scheduled sensor samples soff before its
        # slot starts (ponytail: locked to the hub's superframe clock, no slip between sampling and slots)
        self.soff = cell.random(self.n) * self.T
        # contact intervals of each BAN pair (0-1, 2-3, ..), alternating exponential on / off with mean contact
        # TOUCH_MEAN_US and the flag's share of time in contact; their own substream, so the other draws stay as they are
        self.touch = {}
        if sp.flags.get("touch"):
            trng, s, mc = np.random.default_rng([seed, 3]), sp.flags["touch"], pr.TOUCH_MEAN_US
            for b in range(0, k - 1, 2):
                iv, t, on = [], 0.0, trng.random() < s
                while t < duration_us:
                    d = trng.exponential(mc if on else mc * (1 - s) / s) if s < 1 else duration_us
                    if on:
                        iv.append((t, t + d))
                    t, on = t + d, not on
                self.touch[b] = self.touch[b + 1] = iv
        self.deadline = sel.get("deadline_us", self.T)   # on-time deadline D of the cell (DESIGN 3, OI-4)
        bc = sp.beacon
        self.reserved = math.ceil(sel.get("reserved_us", bc["reserved_us"] if bc else 0))
        self.bdur = sel.get("beacon_us", bc["airtime_us"] + bc.get("airtime_per_dev_us", 0) * n_dev) if bc else 0
        # BCC full ladder (DESIGN 4.3): from L2 the beacon's ack field acknowledges, so every sensor decodes the beacon
        # (its real one-tone layout; multi-tone: on the first tone); from L3 the hub polls failed frames in the pool
        self.late, self.pool_n, self.pend = sp.flags.get("opt", 0) >= 2, sel.get("pool", 0), [[] for _ in range(k)]
        if self.late:
            bl = pr.beacon_layout(sp, n_dev, sel["frames_per_period"])
            self.blay = ([bl] if sp.tones else bl, self.bdur, 0, 0, 0)
        self.lost_lim, self.missed = (bc or {}).get("max_lost"), [0] * self.n   # consecutive missed beacons
        self.blk = [0] * self.n     # flag span: frames each sensor sent since its last B-Ack
        # 'adapt': per sensor unspread or Walsh, its run of reported outcomes and the outcomes the next heard
        # beacon reports (the allocations stay sized for Walsh frames)
        self.sadapt = sp.flags.get("spread") == "adapt"
        self.plain, self.srun, self.sfb = [True] * self.n, [0] * self.n, [[] for _ in range(self.n)]
        self.bcol = 5 if bc and "energy" in bc else 0   # en column of the beacon's TX time (RX: + 1)
        if "slot_us" in sel:    # quantised slots: one exchange per TSCH cell or SmartBAN slot; a GTS holds its retries
            self.slot = sel["slot_us"]
            self.cap_ex = math.inf if sp.pid == "ieee802154" else sel["frames_per_period"]
        elif sp.pid == "ble4thesis":
            self.slot, self.cap_ex = self.T // n_dev, math.inf   # BLE_SIM: ci_us // num_devices
        else:
            self.slot, self.cap_ex = (self.T - self.reserved) / n_dev, math.inf
        if sp.tsch:
            self.S = self.T // sp.tsch["timeslot_us"]
        # channel state
        if sp.hop == "fixed":
            self.chist = [[(0, c)] for c in self._pick_channels()]
            if sp.channel_switch:   # per IBI [index, failed, attempted] uplink exchanges
                self.win = [deque(maxlen=sp.adapt["window"]) for _ in range(k)]
        elif sp.hop == "perm":
            self.perm = [rng.permutation(self.nch) for _ in range(k)]
        elif sp.hop == "csa2":
            aa = rng.integers(0, 2 ** 32, (k, self.m))
            self.cid = (aa >> 16) ^ (aa & 0xFFFF)
            self.afh = [[Adapt(range(37), sp.adapt, True) for _ in range(self.m)] for _ in range(k)]
            self.cmap = [[list(range(37)) for _ in range(self.m)] for _ in range(k)]
            self.pending = [[None] * self.m for _ in range(k)]
        else:   # tsch
            hs = sp.tsch["hopping_sequence"]
            self.bl = [Adapt(hs, sp.adapt, False) for _ in range(k)]
            self.hs = [a.used for a in self.bl]
            self.asn0 = rng.integers(0, 2 ** 16, k)     # own ASN origin per BAN (PARAMETERS 3.3)
        # retry room of the FULL schedules: TSCH retries in the slotframe's shared cell, SmartBAN
        # in the C/M slots of the IBI; bo: SmartBAN MCA backoff per node, [head frame, BE, BC]
        self.shared, self.cm, self.bo = bool(sp.tsch), sp.pid == "smartban", {}
        self.floor = -math.inf     # sensitivity in dB re the noise: frames below it are lost (NB FULL only)
        if self.sinr:
            self._setup_sinr(venue, density)
        if sp.access == "ant":
            self.ant_phase = self._ant_phases()
            # frequency agility (PARAMETERS 3.6 FULL step 4): per BAN its private frequency and two more of the pool
            self.agile = [[h[0][1]] + [int(x) for x in rng.choice([c for c in range(self.nch) if c != h[0][1]],
                                                                  sp.ant["agility_freqs"] - 1, replace=False)]
                          for h in self.chist]
        for b in range(k):
            if bc and sp.access == "tdma":
                self._send(self._beacon(b), None)
            if self.pool_n:
                self._send(self._pool(b), None)
            for j in range(1, n_dev + 1):
                self._send({"tdma": self._tdma, "csma": self._csma, "ant": self._ant}[sp.access](b, j), None)

    # --- event loop ---------------------------------------------------------------------------------------------

    def _send(self, gen, val):
        try:
            req = gen.send(val)
        except StopIteration:
            return
        self.seq += 1
        if req.__class__ is Tx:
            tr = self.tracks.get(req.key)
            if tr is None:
                tr = self.tracks[req.key] = Track()
            tr.add(req)
            self.ntrack[req.node].add(req)
            if self.sinr:
                self.atrack.add(req)
            self.maxd = max(self.maxd, req.e - req.s)
            self.c[req.ban]["airtime"] += req.e - req.s
            heapq.heappush(self.heap, (req.e, self.seq, gen, req))
        else:
            heapq.heappush(self.heap, (req[2], self.seq, gen, req))

    def run(self):
        heap = self.heap
        while heap:
            _, _, gen, arg = heapq.heappop(heap)
            if arg.__class__ is Tx:
                if arg.rx >= 0:
                    self._eval(arg)
                self._send(gen, arg)
            else:
                self._send(gen, getattr(self, "_" + arg[0])(*arg[1:]))   # _cca, _sense, _wait
        return self

    # --- helpers ------------------------------------------------------------------------------------------------

    def lay(self, octets, plain=None, ri=False):
        """(layout, airtime_us, bits before FEC, payload bits before FEC (phy.decisions), payload info bits) of a frame
        carrying `octets`; multi-tone BCC: the list of its tone packets' layouts, the counts summed over them.
        plain True / False: the frame unspread or Walsh, None: the spec's own. ri: the RI-mode frame
        of a burst (flag burst, spec.ri)."""
        key = ("ri", octets) if ri else octets if plain is None else (octets, plain)
        v = self.lays.get(key)
        if v is None:
            sp = self.sp.ri if ri else self.sp
            ls = pr.tone_layouts(sp, octets, plain)
            blocks = [b for L in ls for b in L]
            pay = [b for b in blocks if b["kind"] == "payload"]
            v = self.lays[key] = (ls if sp.tones else ls[0], pr.airtime_us(sp, octets, plain),
                                     sum(phy.decisions(b) for b in blocks), sum(phy.decisions(b) for b in pay),
                                     sum(b["info_bits"] for b in pay))
        return v

    def dlay(self, node, octets):
        """Layout of a data frame of sensor node: under 'adapt' in the sensor's current spreading."""
        return self.lay(octets, self.plain[node]) if self.sadapt else self.lay(octets)

    def _report(self, node, tx):
        if self.sadapt:
            self.sfb[node].append(tx.ok)

    def ctl(self, dur):
        """Layout of a control frame (poll, ACK) of airtime dur: the L = 0 frame when its airtime matches, else one
        uncoded block of dur x bit rate decisions (802.15.4 Imm-Ack, TSCH Enh-Ack, SmartBAN ACK)."""
        key = ("ctl", dur)
        if key not in self.lays:
            if pr.airtime_us(self.sp, 0) == dur:
                self.lays[key] = self.lay(0)
            else:
                n = round(dur * self.sp.phy["bit_rate_bps"] / 1e6)
                self.lays[key] = ([dict(t0_us=0.0, t1_us=float(dur), kind="header", n_bits=n, info_bits=0,
                                        code=None, rep=1)], dur, n, 0, 0)     # not repeated (SmartBAN Nrep)
        return self.lays[key]

    def ex(self, octets, pool=False, cont=False):
        key = ("ex", octets, pool, cont)
        if key not in self.lays:
            st = pr.exchange(self.sp, octets, pool, cont)
            self.lays[key] = (st, sum(x.dur_us for x in st))
        return self.lays[key]

    def _key(self, b, ch):
        return ch + b * self.nch if self.iso else ch

    def _tx(self, t, dur, ch, b, a, r, lay=None, data=False, lead=None, book=True, col=0):
        """Transmission node a -> node r (-1: nobody decodes it) on channel ch of BAN b; books TX time at a and RX time
        lead + dur at r, which listens from t - lead: lead None = the receive guard, booked once per event at its first
        reception (PARAMETERS 4.1); 0 for a response whose timing r knows (ACK, BLE T_IFS turn, the next exchange of a
        slot); book=False leaves r's RX time to the caller; col: the en column of the TX time (RX: col + 1)."""
        lead = self.guard if lead is None else lead
        self._on(a, col, dur)
        if r >= 0 and book:
            self._on(r, col + 1, dur + lead)
        tx = Tx(t, dur, self._key(b, ch), a, r, b, lay, data)
        tx.w = t - lead
        return tx

    def _on(self, node, col, dur):
        e = self.en[node]
        e[col] += dur
        e[col + 2] += 1

    def _wake(self, *nodes):
        for x in nodes:
            self.en[x][4] += 1

    def _enqueue(self, b, q, cap, f):
        self.c[b]["generated"] += 1
        if cap is not None and len(q) >= cap:
            self.c[b]["drop_queue"] += 1
        else:
            q.append(f)

    def _arrive(self, b, q, cap, ph, tg, kg, t):
        """Enqueue the frames generated at ph + i tg <= t (before the end); return the next generation index."""
        while True:
            g = ph + kg * tg
            if g > t or g >= self.D:
                return kg
            for o, nb in self.plan:
                self._enqueue(b, q, cap, Frame(g, o, nb))
            kg += 1

    def _queue(self, b):
        q = deque()
        self.queues[b].append(q)
        n, ql = len(self.sel["frame_payload_octets"]), self.sp.queue_len
        return q, None if ql is None else ql * n   # queue_len counts periods of frames (A)

    def _attempt(self, b, f, tx):
        c = self.c[b]
        c["tx_attempts"] += 1
        if f.tries == 0:
            c["first"] += 1
            c["first_fail"] += not tx.ok
        f.tries += 1
        c["coll"] += tx.overlap
        c["raw"] += tx.raw
        c["post"] += tx.post
        c["coded"] += tx.coded
        c["info"] += tx.info
        if tx.ok and not f.delivered:
            f.delivered = True
            c["delivered_ok"] += 1
            c["bits_ok"] += f.bits
            lat = tx.e - f.gen
            self.lat[b].append(lat)
            c["on_time"] += lat <= self.deadline   # D = the cell's reference period (OI-4), same for all protocols

    def _lost(self, b, tx):
        """Final failure of a frame sent once (no ARQ): corrupt delivery (no payload CRC) or its drop reason."""
        self.c[b]["delivered_corrupt" if tx.corrupt else "drop_" + tx.reason] += 1

    def _retry(self, b, f, q, done):
        """ARQ bookkeeping after an attempt of the head frame f: done (acked) or the retry limit ends it."""
        lim = self.sp.arq_max_retx
        if done or (lim is not None and f.tries > lim):
            q.popleft()
            if not f.delivered:
                self.c[b]["drop_retry"] += 1

    # --- channel rules ------------------------------------------------------------------------------------------

    def _pick_channels(self):
        """Fixed channel per BAN: uniform, or with ED scan (BANs formed in random order, ties random) the channel with
        the fewest BANs (802.15.4, SmartBAN: 'fewest active BANs', 'free DCH') or, 802.15.6 NB, the least energy (3.4
        'takes the minimum'): BANs weighted by the ACR of their channel offset, so a 0 dB neighbour counts in full
        (ponytail: BAN counts stand for the measured energy, no geometry)."""
        rng = self.rng
        if not self.sp.ed_scan:
            return [int(x) for x in rng.integers(0, self.nch, self.k)]
        w = np.array(self.acrlin) if self.sp.pid == "ieee802156nb" else np.eye(self.nch)
        load, out = np.zeros(self.nch), [0] * self.k
        for b in rng.permutation(self.k):
            e = w @ load
            out[b] = int(rng.choice(np.flatnonzero(e <= e.min() + 1e-9)))
            load[out[b]] += 1
        return out

    def _chan(self, b, j, n, cell=0):
        """Channel index of link j of BAN b in superframe (event) n; cell = timeslot of the TSCH slotframe."""
        hop = self.sp.hop
        if hop == "fixed":
            h = self.chist[b]
            i = len(h) - 1
            while h[i][0] > n:
                i -= 1
            return h[i][1]
        if hop == "perm":
            return int(self.perm[b][n % self.nch])
        if hop == "csa2":
            p = self.pending[b][j]
            if p and n >= p[0]:
                self.cmap[b][j], self.pending[b][j] = p[1], None
            return csa2(n, int(self.cid[b, j]), self.cmap[b][j])
        hs = self.hs[b]
        return hs[(self.asn0[b] + n * self.S + cell + j % 16) % len(hs)]   # HS[(ASN + chOff) mod len], chOff j mod 16

    def _adapt(self, b, j, n, ch, fail):
        sp = self.sp
        if sp.hop == "csa2":    # BLE AFH: per-link map, in force 6 events later (PARAMETERS 3.1 FULL step 6)
            a = self.afh[b][j]
            if a.record(ch, fail):
                self.pending[b][j] = (n + pr.BLE["afh"]["instant_offset_events"], list(a.used))
        elif sp.hop == "tsch":  # blacklisting (PARAMETERS 3.3 FULL+)
            self.bl[b].record(ch, fail)
        elif sp.channel_switch:   # SmartBAN migration (3.5, A): uplink failure ratio > 20 % over 10 IBIs
            w = self.win[b]
            if not w or w[-1][0] != n:      # first outcome of IBI n: judge the 10 completed IBIs before it
                if (len(w) == w.maxlen and sum(x[1] for x in w) > sp.adapt["threshold"] * sum(x[2] for x in w)
                        and self.chist[b][-1][0] <= n):     # to the least-loaded channel, 2 beacons later
                    load = np.bincount([h[-1][1] for h in self.chist], minlength=self.nch)
                    load[ch] = self.k + 1
                    new = int(self.rng.choice(np.flatnonzero(load == load.min())))
                    self.chist[b].append((n + pr.SB["full"]["migration"]["switch_after_beacons"], new))
                    w.clear()
                w.append([n, 0, 0])
            w[-1][1] += fail
            w[-1][2] += 1

    # --- access processes ---------------------------------------------------------------------------------------

    def _beacon(self, b):
        """Hub beacon at each superframe start (TSCH: an EB in cell 0 every eb period); sensors listen to it. With a
        beacon-loss rule (802.15.4 mode B: 4, SmartBAN: 10) the beacon is decided at every sensor and a sensor counts
        its consecutive misses (ponytail: evaluated as one uncoded block, no separate sync part)."""
        sp, hub = self.sp, b * self.m
        f1 = 1 + self.eps[hub]
        for n in count():
            t0 = self.phase[b] + n * self.T * f1
            if t0 >= self.D:
                return
            if sp.tsch:
                per, ts = self.sp.beacon["period_us"], sp.tsch
                if n and (n * self.T) // per == ((n - 1) * self.T) // per:
                    if self.shared:     # no EB: the RFC 8180 minimal cell (DESIGN 6)
                        end = t0 + ts["timeslot_us"]
                        yield ("wait", end, end)
                        self._cell_rx(b, t0, t0 + ts["tx_offset_us"] + 1, self._chan(b, 0, n),
                                      ts["tx_offset_us"] - ts["rx_offset_us"], ts["rx_wait_us"], hub, True)
                    continue
                if self.shared:     # EB in the minimal cell: the sensors that send nothing there receive it
                    self._wake(hub)
                    yield self._tx(t0 + ts["tx_offset_us"], self.bdur, self._chan(b, 0, n), b, hub, -1,
                                   self.ctl(self.bdur))
                    self._cell_rx(b, t0, t0 + ts["tx_offset_us"] + 1, self._chan(b, 0, n),
                                  ts["tx_offset_us"] - ts["rx_offset_us"], ts["rx_wait_us"], hub, True)
                    continue
                t0 += ts["tx_offset_us"]
            self._wake(hub)
            tx = yield self._tx(t0, self.bdur, self._chan(b, 0, n), b, hub, -1,
                                self.blay if self.late else self.ctl(self.bdur), col=self.bcol)
            for x in range(hub + 1, hub + self.m) if self.lost_lim or self.late else ():
                tx.rx, tx.ok, tx.reason = x, True, None
                self._eval(tx)
                self.missed[x] = 0 if tx.ok else self.missed[x] + 1
            if not sp.tsch:
                for x in range(hub + 1, hub + self.m):
                    self._on(x, self.bcol + 1, self.bdur + self.guard)
                    self._wake(x)
            if self.cm:     # C/M slots after the scheduled ones (_cm): the hub's receive guard and the longest frame
                ts, d = self.reserved, pr.SB["tmua_us"]["cm"]
                c0 = t0 + ts + self.n_dev * self.slot + d     # data after the MCA sensing
                w = c0 + (pr.SB["cm_slots"] - 1) * ts + 1   # after the last C/M data start: before the next beacon
                yield ("wait", w, w)
                for i in range(pr.SB["cm_slots"]):
                    self._cell_rx(b, c0 + i * ts, c0 + i * ts + 1, self._chan(b, 1, n), self.guard, self.guard, hub,
                                  False, (hub,))

    def _cell_rx(self, b, t0, t1, ch, lead, idle, hub, wake, nodes=None):
        """Listening of BAN b's nodes (default: hub and sensors) in a contention cell whose frames start in [t0, t1)
        (TSCH shared cell, SmartBAN C/M slot), booked once the frames are known: a node that sent one of them does not
        listen (its ACK wait is booked in _exchange); any other listens lead + the longest own-BAN frame it did not
        send (the EB, a sensor's frame), or idle when none came. One radio listens once however many frames collide
        there, so the frames of such a cell book no RX at their receiver (_tx book=False). ponytail: a frame of
        another BAN in the cell is not received."""
        tr = self.tracks.get(self._key(b, ch))
        got = [x for x in tr.over(t0, t1, self.maxd) if x.ban == b and t0 <= x.s < t1] if tr else []
        for r in range(hub, hub + self.m) if nodes is None else nodes:
            if any(x.node == r for x in got):
                continue
            if wake:
                self._wake(r)
            self._on(r, 1, lead + max(x.e - x.s for x in got) if got else idle)

    def _tdma(self, b, j):
        """Sensor j of BAN b in scheduled access: its slot in every superframe, exchanges back to back while the
        next one fits (DESIGN 4.1); BLE full runs connection events (_ble_event). TSCH FULL sends a frame left over from
        an earlier slotframe (failed, or queued behind a failed one) in the slotframe's shared cell (slot 0, the RFC
        8180 minimal cell, where the EB also goes), which comes before the sensor's dedicated cells (ponytail: the
        period's new sample joins the queue after the shared cell, as it samples on the superframe clock), with the
        shared-cell backoff of 802.15.4-2015 6.2.5.3 (PARAMETERS 3.3): no
        backoff before the first shared-cell attempt; after a failed one BE = macMinBE, then min(BE + 1, macMaxBE),
        and the frame skips U{0..2^BE - 1} shared cells; a success in a shared cell, or one in a dedicated cell that
        empties the queue, resets it; a failure in a dedicated cell leaves it (Contiki-NG tsch-queue.c). The sensors
        of a BAN contend in the shared cell (same cell, same channel offset 0). SmartBAN FULL retries in the C/M
        slots of the IBI (_cm)."""
        sp = self.sp
        hub, node = b * self.m, b * self.m + j
        f1 = 1 + self.eps[hub]
        q, cap = self._queue(b)
        mca = sp.csma is not None
        ts = sp.tsch["timeslot_us"] if sp.tsch else 0
        be, win = None, 0     # TSCH shared-cell backoff exponent (None: no shared-cell failure yet) and window
        for n in count():
            # the hub clock drifts the superframe start; the layout inside the superframe stays nominal (ponytail:
            # a slot filled to the superframe end overlaps the next superframe by T x |drift|, < 0.5 us)
            s0 = self.phase[b] + n * self.T * f1 + self.reserved + (j - 1) * self.slot
            if s0 >= self.D:
                return
            s1 = s0 + self.slot
            if self.lost_lim:   # after the superframe's beacon was decided: a sensor out of sync loses its frames
                yield ("wait", s0, s0)
                if self.missed[node] >= self.lost_lim:
                    self.c[b]["generated"] += len(self.sel["frame_payload_octets"])
                    self.c[b]["drop_sync"] += len(self.sel["frame_payload_octets"])
                    continue
            if self.late or self.pool_n:    # after the superframe's beacon was decided: its ack field (L2), or the
                yield ("wait", s0, s0)      # I-Acks of the pooled reference stack, acknowledge
                self._acked(b, node, q)
            if self.shared and win:
                win -= 1
            elif self.shared and q:   # a frame left over from an earlier slotframe: the shared cell comes first
                self._wake(node)
                if (yield from self._try(b, j, n, q, self.phase[b] + n * self.T * f1, self._chan(b, 0, n), True,
                                         False)):
                    be = None
                else:
                    be = sp.tsch["mac_min_be"] if be is None else min(be + 1, sp.tsch["mac_max_be"])
                    win = int(self.rng.integers(0, 2 ** be))
            for o, nb in self.plan:     # generated at the sample instant before the slot
                self._enqueue(b, q, cap, Frame(s0 - self.soff[node], o, nb))
            self._wake(node, hub)
            t = s0
            if mca:   # SmartBAN MCA: the owner senses Tsch at the slot start and skips a busy slot (PARAMETERS 3.5)
                d = sp.csma["cca_us"]
                busy = yield ("cca", t, t + d, self._key(b, self._chan(b, j, n)), node)
                self._on(node, 1, d)
                self.c[b]["cca"] += 1
                if busy:
                    self.c[b]["cca_busy"] += 1
                    yield from self._cm(b, j, n, q)
                    continue
                t += d
            if sp.pid == "ble":
                yield from self._ble_event(b, j, n, q, t, s1)
                continue
            if self.late or self.pool_n:
                yield from self._late(b, j, n, q, t, s1)
                continue
            done = 0
            while q and done < self.cap_ex and t + self.ex(q[0].octets)[1] <= s1 + 1e-6:
                done += 1
                ch = self._chan(b, j, n, (1 + (j - 1) * (self.slot // ts) + round((t - s0) / ts)) if ts else 0)
                dur = self.ex(q[0].octets)[1]
                if (yield from self._try(b, j, n, q, t, ch, done == 1)) and not q:
                    be, win = None, 0
                t += dur
            if self.cm:
                yield from self._cm(b, j, n, q)

    def _acked(self, b, node, q):
        """BCC full from L2 (DESIGN 4.3), at the sensor's allocation: the ack field of the beacon it just heard (base
        sequence number and bitmap) acknowledges every frame the hub holds, a duplicate sent after a missed beacon too
        (ponytail: also one past the field's n bits, which takes consecutive losses); a frame past the retry limit is
        given up; after a missed beacon the sensor keeps its frames (L2 resends them). The pooled reference stack
        learns from its I-Acks, so it cleans up at every allocation."""
        heard = self.missed[node] == 0 or not self.late
        if heard and self.sadapt:   # the beacon's report drives the spreading switch
            for ok in self.sfb[node]:
                self.plain[node], self.srun[node] = spread_step(self.plain[node], self.srun[node], ok)
            self.sfb[node].clear()
        for f in [f for f in q if f.tries]:
            if heard and f.delivered or f.tries > self.sp.arq_max_retx:
                q.remove(f)
                self.c[b]["drop_retry"] += not f.delivered

    def _late(self, b, j, n, q, t, s1):
        """BCC full allocation from L2 (DESIGN 4.3): every unacknowledged frame once, oldest first, while its exchange
        (data, pSIFS: no I-Ack) fits; from L3 the frames not sent yet first, then, in the room they leave (their
        arrivals were dropped at a full queue), the frames the heard beacon reported missing; a frame that fails for
        the first time goes to the hub's list for the superframe's retry pool (_pool). The pooled reference stack
        runs this path with its I-Ack exchanges and knows its failures without the beacon; under the flag
        back its exchanges carry no I-Ack and the hub sends one B-Ack pSIFS after the sensor's last frame of the
        allocation, which the sensor hears (ponytail: as with the I-Ack reference, the sensor learns the outcomes from
        the hub's record at its next allocation, so a lost B-Ack costs only its air time). Under the flag burst the
        frames after the allocation's first go pMIFS apart in RI mode and the B-Ack follows the last pSIFS later
        (ponytail: an RI frame is decided on its own, whatever happened to the frame before it). Under the flag span the
        B-Ack follows only the allocation in which the sensor's frames since the last B-Ack reach span."""
        ch, first, fs, sent = self._chan(b, j, n), True, list(q), 0
        bk = pr.back_us(self.sp)
        s1 -= bk
        if self.pool_n:
            known = not self.late or not self.missed[b * self.m + j]
            fs = [f for f in fs if not f.tries] + [f for f in fs if f.tries and known]
        for f in fs:
            steps, dur = self.ex(f.octets, False, not first)
            if t + dur > s1 + 1e-6:
                break
            tx = (yield from self._exchange(b, j, f, t, ch, steps, first))[0]
            if self.pool_n and not tx.ok and f.tries == 1:
                self.pend[b].append((f, j))
            first, t, sent = False, t + dur, sent + 1
        span, node = self.sp.flags.get("span"), b * self.m + j
        if span:    # the block closes at the end of the allocation in which it reaches span frames
            self.blk[node] += sent
            if self.blk[node] < span:
                sent = 0
            else:
                self.blk[node] = 0
        if bk and sent:
            t += pr.HBC["psifs_us"] if self.sp.flags.get("burst") else 0
            lay = self.lay(pr.BACK_OCTETS)
            yield self._tx(t, lay[1], ch, b, b * self.m, b * self.m + j, lay, lead=0)

    def _pool(self, b):
        """BCC full retry pool from L3 (DESIGN 4.3; 802.15.6 improvised access): after the allocations of every
        superframe the hub polls the frames that failed there, in the order they failed (sensor order, earlier
        superframes first), each at most once per pool, sel['pool'] exchanges of poll, pSIFS, the sensor's frame,
        pSIFS, back to back from the pool start; a frame not served, or failed again, waits for the next pool. Every
        sensor listens (A: before the next beacon it does not know whether its frame failed): from the receive guard
        at the pool start, poll after poll, the guard alone when none comes; the hub waits for an answer to a poll the
        sensor missed for the frame's airtime (as a sensor for a missing I-Ack); a missed poll costs no attempt."""
        sp, hub = self.sp, b * self.m
        sensors, f1, pend = range(hub + 1, hub + self.m), 1 + self.eps[hub], self.pend[b]
        poll, ifs = self.lay(pr.poll_octets(sp)), pr.HBC["psifs_us"]
        for n in count():
            t = self.phase[b] + n * self.T * f1 + self.reserved + self.n_dev * self.slot
            if t >= self.D:
                return
            yield ("wait", t, t)
            pend[:] = [x for x in pend if not x[0].delivered and x[0].tries <= sp.arq_max_retx]
            self._wake(*sensors)
            lead = self.guard
            for f, j in pend[:self.pool_n]:
                ch = self._chan(b, j, n)
                if lead:
                    self._wake(hub)
                p = yield self._tx(t, poll[1], ch, b, hub, hub + j, poll, lead=lead, book=False)
                for x in sensors:
                    self._on(x, 1, lead + poll[1])
                if p.ok:    # the sensor's exchange at its level: the data, and below L2 the I-Ack too
                    yield from self._exchange(b, j, f, p.e + ifs, ch, self.ex(f.octets, True)[0], False)
                else:
                    self._on(hub, 1, self.dlay(hub + j, f.octets)[1])
                lead, t = 0, t + self.sel["pool_us"]
            if lead:
                for x in sensors:
                    self._on(x, 1, lead)

    def _try(self, b, j, n, q, t, ch, first, book=True):
        """One exchange of the head frame of q from t on channel ch and its bookkeeping; True when the frame is done
        (acknowledged, or sent once without ARQ). book=False: a contention cell, whose listening _cell_rx books."""
        f = q[0]
        tx, ack_ok, busy = yield from self._exchange(b, j, f, t, ch, self.ex(f.octets)[0], first, book)
        self._adapt(b, j, n, ch, busy or not tx.ok)
        if self.sp.arq_max_retx == 0:
            q.popleft()
            if not tx.ok:
                self._lost(b, tx)
            return True
        if busy:   # TSCH CCA busy counts as a failed attempt (PARAMETERS 3.3, OI-16)
            c = self.c[b]
            c["tx_attempts"] += 1
            c["first"] += f.tries == 0
            c["first_fail"] += f.tries == 0
            f.tries += 1
        self._retry(b, f, q, ack_ok)
        return ack_ok

    def _cm(self, b, j, n, q):
        """SmartBAN FULL (PARAMETERS 3.5 FULL, TS 103 325 7.3.2.0, 7.3.3): a frame still queued after the node's own
        slot (it failed there, or its slot went to a retry or was busy: 7.3.2.0 'when it has not been allocated
        sufficient time slots') goes in the IBI's C/M slots (after the scheduled slots) by the MCA backoff, the C/M
        access of a BAN with MCA on (7.3.3.2.2, 7.3.3.1.1): per head frame BE = BEmin and a backoff counter BC =
        U{0..2^BE - 1} of C/M slots; at BC 0 the node senses TMUA = Tsch (16 us) at the slot start and, idle, sends its
        data after it, the ACK in the same slot; busy, or sent and not acknowledged, BE = min(BE + 1, BEmax) and a new
        BC. The other nodes of the BAN contend in the same slots (they sense at the same instant and neither hears the
        other's frame, which starts as the sensing ends, so they meet in both reception models). A frame
        not sent there waits for the node's next own slot, ahead of new data (the recipe's 'else'); RT = 1 ends it after
        one retransmission (the recipe's limit, not the MCA's RTmax). The hub's listening, once per slot, is booked
        with the beacon (_cell_rx)."""
        node, ts, d = b * self.m + j, self.reserved, pr.SB["tmua_us"]["cm"]   # the D-Beacon slot is one TS
        t0 = self.phase[b] + n * self.T * (1 + self.eps[b * self.m]) + ts + self.n_dev * self.slot
        ch = self._chan(b, j, n)
        for i in range(pr.SB["cm_slots"]):
            if not q:
                return
            st = self.bo.get(node)
            if st is None or st[0] is not q[0]:
                st = self.bo[node] = [q[0], BE_MIN, int(self.rng.integers(0, 2 ** BE_MIN))]
            if st[2]:
                st[2] -= 1
                continue
            t = t0 + i * ts
            self._wake(node)
            busy = yield ("cca", t, t + d, self._key(b, ch), node)
            self._on(node, 1, d)
            self.c[b]["cca"] += 1
            self.c[b]["cca_busy"] += busy
            if busy or not (yield from self._try(b, j, n, q, t + d, ch, False, False)):
                st[1] = min(st[1] + 1, BE_MAX)
                st[2] = int(self.rng.integers(0, 2 ** st[1]))

    def _exchange(self, b, j, f, t, ch, steps, first, book=True):
        """One exchange of protocols.exchange for frame f from time t; returns (data Tx, ACK ok, CCA busy). The hub
        listens with the guard for the first exchange of a slot; a TSCH cell follows the timeslot template (3.3): the
        hub listens from macTsRxOffset (macTsRxWait when nothing comes), the sensor for the Enh-Ack from
        macTsRxAckDelay (macTsAckWait when none comes). book=False: the hub's listening is booked by _cell_rx."""
        hub, node, ts = b * self.m, b * self.m + j, self.sp.tsch
        tx, ack_ok = None, False
        for st in steps:
            if st.kind == "cca":
                busy = yield ("cca", t, t + st.dur_us, self._key(b, ch), node)
                self._on(node, 1, st.dur_us)
                self.c[b]["cca"] += 1
                if busy:
                    self.c[b]["cca_busy"] += 1
                    if ts and book:     # the hub listens macTsRxWait for nothing (a shared cell: _cell_rx)
                        self._on(hub, 1, ts["rx_wait_us"])
                    return None, False, True
            elif st.kind in ("data", "burst"):
                lay = self.lay(f.octets, ri=True) if st.kind == "burst" else self.dlay(node, f.octets)
                lead = ts["tx_offset_us"] - ts["rx_offset_us"] if ts else None if first else 0
                tx = yield self._tx(t, lay[1], ch, b, node, hub, lay, True, lead, book)
                self._attempt(b, f, tx)
                self._report(node, tx)
            elif st.kind == "ack":
                if tx.ok:       # the hub acknowledges only a frame it received
                    lead = ts["tx_ack_delay_us"] - ts["rx_ack_delay_us"] if ts else 0
                    ack = yield self._tx(t, st.dur_us, ch, b, hub, node, self.ctl(st.dur_us), lead=lead)
                    ack_ok = ack.ok
                else:
                    self._on(node, 1, ts["ack_wait_us"] if ts else st.dur_us)
            t += st.dur_us
        return tx, ack_ok, False

    def _ble_event(self, b, j, n, q, t, end):
        """BLE ACL connection event (PARAMETERS 3.1 FULL, DESIGN 4.1): poll C, T_IFS, sensor PDU P, T_IFS while the
        next pair fits; the poll's NESN acknowledges the last P the hub received. SINR model: a C the sensor missed
        (sync) ends the event; a C it synced to but whose CRC failed still gets its P, the unacknowledged head again
        (no NESN seen), and counts one CRC error at the sensor; a P the hub missed ends the event, a P with a CRC
        error counts one at the hub. Collision model (sync and CRC failures look alike): a failed C ends the event,
        a failed P counts a CRC error. The event also ends after 2 consecutive CRC errors at either side or when the
        sensor has no more data (MD = 0). AFH scores every expected data P, a missed C included. No retry limit."""
        hub, node = b * self.m, b * self.m + j
        ch = self._chan(b, j, n)
        ifs, poll = pr.BLE["t_ifs_us"], self.ctl(pr.airtime_us(self.sp, 0))
        crc_hub = crc_sensor = 0
        first = True
        while True:
            nxt = q[1] if len(q) > 1 and q[0].delivered else q[0] if q and not q[0].delivered else None
            p_us = max(self.lay(x.octets if x else 0)[1] for x in (nxt, q[0] if q else None))   # either answer
            if t + poll[1] + p_us + 2 * ifs > end + 1e-6:
                return
            c = yield self._tx(t, poll[1], ch, b, hub, node, poll, lead=None if first else 0)
            first = False
            t = c.e + ifs
            if not c.ok and (not self.sinr or c.reason == "sync"):   # missed C: no P
                if nxt is not None:
                    self._adapt(b, j, n, ch, True)
                return
            if c.ok:
                crc_sensor = 0
                if q and q[0].delivered:
                    q.popleft()
                f = nxt
            else:
                crc_sensor += 1
                f = q[0] if q else None
            lay = self.lay(f.octets if f else 0)
            d = yield self._tx(t, lay[1], ch, b, node, hub, lay, f is not None, lead=0)
            t = d.e + ifs
            if f is not None:
                self._attempt(b, f, d)
                self._adapt(b, j, n, ch, not d.ok)
            if not d.ok and self.sinr and d.reason == "sync":     # missed P
                return
            crc_hub = 0 if d.ok else crc_hub + 1
            if 2 in (crc_hub, crc_sensor) or not q or (len(q) == 1 and q[0].delivered):
                return

    def _csma(self, b, j):
        """802.15.4 unslotted CSMA/CA (DESIGN 4.2, PARAMETERS 3.2 mode A FULL): own periodic generation with random
        phase, backoff, CCA, turnaround, data, Imm-Ack after aTurnaroundTime, retries from CSMA, FIFO drop newest; at
        the end of the run it finishes the frames of the sample in service, as a slot runs to its end."""
        rng, c = self.rng, self.sp.csma
        hub, node = b * self.m, b * self.m + j
        q, cap = self._queue(b)
        ph, tg = self.soff[node], self.T * (1 + self.eps[node])
        ch = self._chan(b, j, 0)
        key = self._key(b, ch)
        kg, t, cur = 0, ph, None
        while True:
            kg = self._arrive(b, q, cap, ph, tg, kg, t)
            if not q:
                t = ph + kg * tg
                if t >= self.D:
                    return
                continue
            if t >= self.D and q[0].gen != cur:   # past the end only the rest of the sample in service, as a slot
                return
            f = q[0]
            cur = f.gen
            self._wake(node)
            nb, be = 0, c["mac_min_be"]
            while True:
                t += int(rng.integers(0, 2 ** be)) * c["unit_backoff_us"]
                busy = yield ("cca", t, t + c["cca_us"], key, node)
                self._on(node, 1, c["cca_us"])
                t += c["cca_us"]
                self.c[b]["cca"] += 1
                if busy:
                    self.c[b]["cca_busy"] += 1
                    nb, be = nb + 1, min(be + 1, c["mac_max_be"])
                    if nb > c["mac_max_csma_backoffs"]:
                        kg = self._arrive(b, q, cap, ph, tg, kg, t)
                        q.popleft()
                        if not f.delivered:
                            self.c[b]["drop_csma"] += 1
                        break
                    continue
                lay = self.lay(f.octets)
                # the coordinator receives all the time (rows() books it): no listen window, no guard
                d = yield self._tx(t + c["turnaround_us"], lay[1], ch, b, node, hub, lay, True, math.inf, False)
                self._attempt(b, f, d)
                ack_ok = False
                if d.ok:
                    a = yield self._tx(d.e + c["ack_delay_us"], c["ack_us"], ch, b, hub, node, self.ctl(c["ack_us"]),
                                       lead=c["ack_delay_us"], book=False)
                    ack_ok = a.ok
                # the sensor listens until the ACK ends or macAckWaitDuration runs out
                self._on(node, 1, c["ack_delay_us"] + c["ack_us"] if ack_ok else c["ack_wait_us"])
                ifs = self.ex(f.octets)[0][-1].dur_us     # SIFS / LIFS after the ACK or the ACK timeout
                t = d.e + (c["ack_delay_us"] + c["ack_us"] if ack_ok else c["ack_wait_us"]) + ifs
                kg = self._arrive(b, q, cap, ph, tg, kg, t)
                self._retry(b, f, q, ack_ok)
                if not q or q[0] is not f:
                    break
                nb, be = 0, c["mac_min_be"]

    def _ant(self, b, j):
        """ANT master j of BAN b (PARAMETERS 3.6 FULL): channel period T / 1.25 (acknowledged messages), clock drift,
        phase from the listen at channel open (_ant_phases); acknowledged message with host retry in the next slot (6
        attempts), old payload as broadcast when the FIFO is empty; advanced burst (channel period T, a packet of up to
        24 octets every 3.2 ms, 1 + 20 attempts, abort) when L exceeds 8 octets; adaptive isochronous re-timing (FULL
        step 4): the master senses a 500 us guard window after its slot and moves 250 us per period away (earlier)
        while a co-channel frame starts in it (ponytail: bursts not re-timed); frequency agility (step 4): after 8
        consecutive unacknowledged messages or burst packets the master moves to the BAN's next agility frequency
        (ponytail: the hub's channel follows at once and the master keeps its phase, no new listen)."""
        a = self.sp.ant
        hub, node = b * self.m, b * self.m + j
        q, cap = self._queue(b)
        burst = self.sel.get("burst", False)
        tg = self.T * (1 + self.eps[node])
        tc = tg / (1 if burst else a["ch_rate_factor"])
        ph, ch = self.ant_phase[node], self._chan(b, j, 0)
        kg, sent, slot, shift, fi, miss = 0, False, 0, 0.0, 0, 0
        while True:
            t = ph + shift + slot * tc     # products, not running sums
            if t >= self.D:
                return
            kg = self._arrive(b, q, cap, self.soff[node], tg, kg, t)   # the sample clock has its own phase
            self._wake(node, hub)
            if burst:
                tb = t
                while q and q[0].gen <= t:
                    f = q[0]
                    lay = self.lay(f.octets)
                    d = yield self._tx(tb, lay[1], ch, b, node, hub, lay, True)
                    self._attempt(b, f, d)
                    ok = yield from self._ant_ack(b, j, d, ch)
                    fi, miss, ch = self._agile(b, fi, 0 if ok else miss + 1, ch)
                    tb += a["burst_spacing_us"]
                    if ok:
                        q.popleft()
                    elif f.tries > a["burst_retries"]:   # transfer aborted: the rest of the burst is lost
                        while q and q[0].gen <= t:
                            self._retry(b, q[0], q, True)
                slot += max(1, math.ceil((tb - t) / tc - 1e-9))
                continue
            lay = self.lay(q[0].octets if q else 0)
            msg = bool(q)
            if msg:
                f = q[0]
                d = yield self._tx(t, lay[1], ch, b, node, hub, lay, True)
                self._attempt(b, f, d)
                sent = True
                ok = yield from self._ant_ack(b, j, d, ch)
                self._retry(b, f, q, ok)
                fi, miss, ch = self._agile(b, fi, 0 if ok else miss + 1, ch)
            elif sent:              # empty FIFO: the old payload goes out as a broadcast (air and energy only)
                yield self._tx(t, lay[1], ch, b, node, -1)
                self._on(hub, 1, lay[1] + self.guard)
            if sent:                # adaptive isochronous: guard window after the slot [t, t + 418 us], every period
                e0 = t + self.ex(0)[1]
                if (yield ("sense", e0, e0 + a["guard_us"], self._key(b, ch), node)):
                    shift -= a["iso_shift_us"]
                if msg:     # the radio stays in RX from the ACK window W (t + 274 to 450 us) on: no new start-up
                    w_end = t + a["slave_reply_us"][0] + a["rx_window_us"]
                    self.en[node][1] += max(0.0, e0 + a["guard_us"] - w_end)
                else:       # after an old-payload broadcast (no W) the receiver starts anew
                    self._on(node, 1, a["guard_us"])
            slot += 1

    def _agile(self, b, fi, miss, ch):
        """Frequency agility: (index, unacknowledged count, channel) after an outcome; the trigger moves on."""
        a = self.sp.ant
        if miss < a["agility_unacked"]:
            return fi, miss, ch
        fi = (fi + 1) % a["agility_freqs"]
        return fi, 0, self.agile[b][fi]

    def _ant_ack(self, b, j, d, ch):
        """The slave's ACK after the turnaround; the master listens for the RX window W (3.6) either way."""
        hub, node, a = b * self.m, b * self.m + j, self.sp.ant
        dur = pr.ANT["ack_airtime_us"]
        self._on(node, 1, a["rx_window_us"])
        if d.ok:
            ack = yield self._tx(d.e + a["turnaround_us"], dur, ch, b, hub, node, self.ctl(dur), lead=0, book=False)
            return ack.ok
        return False

    def _ant_phases(self):
        """Master phases from the listen at channel open (PARAMETERS 3.6 FULL step 1; masters open one after another
        in random order): up to 100 random tries for a phase whose occupancy (one exchange, or the whole burst of 3.2
        ms per packet), widened by the guard, overlaps no master already on the channel; when none is free, the
        least-occupied phase: midway between the two heard masters farthest apart (DESIGN 4.1). The listen (one
        channel period of RX) is not in the energy (DESIGN 6)."""
        a, rng = self.sp.ant, self.rng
        n, burst = len(self.sel["frame_payload_octets"]), self.sel.get("burst", False)
        tc = self.T / (1 if burst else a["ch_rate_factor"])
        d = n * a["burst_spacing_us"] if burst else self.ex(self.sel["frame_payload_octets"][0])[1]
        g = a["guard_us"]
        ph, occ = np.zeros(self.n), defaultdict(list)
        for x in rng.permutation([b * self.m + j for b in range(self.k) for j in range(1, self.m)]):
            key = self._key(x // self.m, self._chan(x // self.m, x % self.m, 0))
            p = rng.uniform(0, tc)
            o = np.array(occ[key]) - g
            for _ in range(100):
                if not (((p - o) % tc < d + 2 * g) | ((o - p) % tc < d)).any():
                    break
                p = rng.uniform(0, tc)
            else:
                s = np.sort(occ[key])
                gap = np.diff(np.append(s, s[0] + tc))
                i = int(np.argmax(gap))
                p = (s[i] + gap[i] / 2) % tc
            occ[key].append(p)
            ph[x] = p
        return ph

    # --- reception ----------------------------------------------------------------------------------------------

    def _setup_sinr(self, venue, density):
        """Static SNR (dB, re the receiver noise) between every pair of nodes, per-packet fade sigmas (RF only)."""
        sp, rng, k, m, ph = self.sp, self.room_rng, self.k, self.m, self.sp.phy
        body = np.repeat(np.arange(k), m)
        self.fsig = None
        self.noise_dbm = 0.0
        if sp.inter_body == "rf":
            room = self.room = chn.make_room(rng, k, self.n_dev, venue, density)
            self.noise_dbm = chn.noise_dbm(ph["noise_bw_hz"], ph["nf_db"], ph["l_impl_db"])
            if ph["sens_dbm"] is not None:
                self.floor = ph["sens_dbm"] - self.noise_dbm
            snr = ph["tx_power_dbm"] + chn.gain_matrix_db(room, rng, self.case) - self.noise_dbm
            on, off = chn.CH["on_body"]["fade_sigma_dB"][self.case], chn.CH["inter_body"]["fade_sigma_dB"]
            self.fsig = np.where(body[:, None] == body[None, :], on, off).tolist()
        elif sp.pid == "hbc":   # SNR in the 5.25 MHz band: ec_n0_db with rc = bandwidth gives P / (N0 B)
            self.room = chn.make_room(rng, k, self.n_dev, venue, density)
            c = chn.hbc_coupling_db(self.room, self.case)
            snr = chn.ec_n0_db(ph["amp_v"], -chn.HBC_INTRA_DB - c[body][:, body], ph["noise_bw_hz"], ph["nf_db"])
        else:                   # BCC: per-link Ec/N0 of a measured path; isolated bodies do not hear each other
            hubs = np.repeat(np.arange(k) * m, m - 1)
            sens = hubs + np.tile(np.arange(1, m), k)
            rc = sp.chip_rate_cps
            e = chn.bcc_link(rng, k * self.n_dev, ph["amp_v"], rc, ph["nf_db"],
                             f_hz=[f + d for f in sp.tones or [sp.channels[0] * 1e6] for d in (-rc, rc)])
            tone = chn.bfsk_ec_n0_db(e.reshape(len(e), -1, 2))     # every BFSK tone: its gains at f -/+ R_c
            self.tsnr = dict(zip(sens.tolist(), tone.tolist()))    # sensor node -> its link's Ec/N0 on every tone
            link = tone.min(axis=1).reshape(k, -1)                  # the weakest tone (single-tone BCC: its tone)
            own = np.column_stack([10 * np.log10(np.mean(10 ** (link / 10), axis=1)), link]).ravel()
            snr = np.full((self.n, self.n), -np.inf)
            if sp.flags.get("coupling_db") is not None:   # interferer's own-body strength minus 40 / 60 dB (A)
                room = chn.make_room(rng, k, self.n_dev, venue, density)
                snr = own[:, None] - chn.hbc_coupling_db(room, f"bcc{sp.flags['coupling_db']}")[body][:, body]
            snr[hubs, sens] = snr[sens, hubs] = link.ravel()
            if sp.flags.get("touch"):   # the partner body's nodes at their own-link Ec/N0 less touch_db
                x = sp.flags.get("touch_db", pr.TOUCH_DB)
                for b in range(0, k - 1, 2):
                    for a, rb in ((b, b + 1), (b + 1, b)):
                        na, nr = np.arange(a * m, a * m + m), np.arange(rb * m, rb * m + m)
                        lvl = CONTACT["ec_n0_good_dB"] if self.contact else own[na][:, None]   # the chain's level
                        snr[np.ix_(na, nr)] = lvl - x
        self.snr = snr.tolist()
        # lock threshold (OI-6: the receiver locks on the first packet above sensitivity): the recipe's sync rule
        # (802.15.4), else the sensitivity (NB FULL), else the SINR of a 1e-3 bit (chip) error, the BLE sensitivity
        # criterion (GFSK stacks: 10.94 dB, the calibrated chip sensitivity; HBC and BCC chips (A))
        self.sync_thr = (ph["sync_sinr_db"] if ph["sync_sinr_db"] is not None else self.floor
                         if self.floor > -math.inf else phy.sinr_at_ber(ph))

    def _fade(self, a, r, t):
        key = a * self.n + r
        v = self.fades.get(key)
        if v is None or t - v[0] > COHERENCE_US:
            v = self.fades[key] = (t, self.rng.normal(0.0, self.fsig[a][r]))
        return v[1]

    def _touching(self, bp, bq, q):
        """: q of BAN bq reaches BAN bp: the two are a touching pair, q starts inside a contact interval, and
        bp's hub has not yet interleaved the pair's active periods (fallback: from one superframe after the contact
        starts, when both active periods fit in one period)."""
        iv = self.touch.get(bp)
        if iv is None or bq != bp ^ 1:
            return False
        c = next((a for a, e in iv if a <= q.s < e), None)
        if c is None:
            return False
        if self.sp.flags.get("fallback") and q.s >= c + self.T:
            sp, fr = self.sp, self.n_dev * self.sel["frames_per_period"]
            pool = pr.pool_size(fr, sp.flags.get("pool", pr.POOL_SHARE)) if pr._pooled(sp) else 0
            active = self.reserved + self.n_dev * self.slot + pool * self.sel.get("pool_us", 0)   # fill's slack is given up
            return 2 * active > self.T
        return True

    def _gain(self, a, r, t):
        return self.snr[a][r] + (self._fade(a, r, t) if self.fsig else 0.0)

    def _contact(self, link, a, b):
        """(t0, t1, Ec/N0 dB) pieces of the contact chain of a BCC link over [a, b): continuous in time per link,
        exponential dwell 35 000 / 13 000 chips, 14 dB connected, 12 dB less degraded (DESIGN 5.2, PARAMETERS 3.8)."""
        rng, us = self.rng, 1e6 / self.sp.chip_rate_cps
        mg, mb = CONTACT["good_dwell_chips"] * us, CONTACT["bad_dwell_chips"] * us
        st = self.cstate.get(link)
        x = 0.0 if st is None else math.exp(-(a - st[0]) * (1 / mg + 1 / mb))
        bad = rng.random() < mb / (mg + mb) * (1 - x) + (x if st and st[1] else 0.0)
        out, t = [], a
        while t < b:
            t1 = min(b, t + rng.exponential(mb if bad else mg))
            out.append((t, t1, CONTACT["ec_n0_good_dB"] - (CONTACT["depth_dB"] if bad else 0.0)))
            if t1 < b:
                bad = not bad
            t = t1
        self.cstate[link] = (b, bad)
        return out

    def _eval(self, p):
        """Outcome of transmission p at its receiver, decided at its end (DESIGN 5)."""
        same = [q for q in self.tracks[p.key].over(p.s, p.e, self.maxd) if q is not p]
        hd = self.ntrack[p.rx].over(p.s, p.e, self.maxd)
        p.overlap = bool(same)
        L = p.lay
        if self.sinr and not hd or self.contact and not same and not hd:
            return self._eval_sinr(p, L)
        # copies decided alone (SmartBAN Nrep, PARAMETERS 3.5 BER model: PER = the product of the copies' PER)
        n = self.sp.phy.get("nrep", 1) if p.data else 1
        d = (p.e - p.s) / n
        clean = sum(not any(q.s < p.s + (i + 1) * d and q.e > p.s + i * d for q in same + hd) for i in range(n))
        if not clean:
            p.ok, p.reason = False, "collision"
        elif self.sp.p_loss and self.rng.random() < self.sp.p_loss ** clean:
            p.ok, p.reason = False, "noise"
        if p.data:     # expected payload errors: overlapped part at 0.5, the rest at the BER that gives p_loss
            fov, cur = 0.0, p.s
            for a, e in sorted((max(q.s, p.s), min(q.e, p.e)) for q in same + hd):
                fov += max(0.0, e - max(a, cur))
                cur = max(cur, e)
            fov /= p.e - p.s
            res = 1 - (1 - self.sp.p_loss) ** (1 / max(L[2], 1))
            p.coded, p.info = L[3], L[4]
            p.raw = L[3] * (0.5 * fov + (1 - fov) * res)
            p.post = p.raw * L[4] / L[3] if L[3] else 0.0

    def _eval_sinr(self, p, L):
        sp, r, nch = self.sp, p.rx, self.nch
        acr, ints, lock = self.acrlin[p.key % nch], [], False
        for q in self.atrack.over(p.s, p.e, self.maxd):
            f = acr[q.key % nch]
            if q is p or not f or q.key // nch != p.key // nch and not self._touching(p.key // nch, q.key // nch, q):
                continue    # key // nch: the BAN of an isolated key; touching partners share the medium
            g = self._gain(q.node, r, q.s)
            if g == -math.inf:
                continue
            if q.key == p.key and sp.sync_on_foreign and q.s < p.s and not lock and self._hears(r, q.s, p):
                lock = g - 10 * math.log10(1 + self._power(q.key, q.s, r, q)) >= self.sync_thr
            ints.append((max(q.s, p.s), min(q.e, p.e), f * 10 ** (g / 10)))
        pieces = (self._contact(p.node if p.node % self.m else r, p.s, p.e) if self.contact
                  else [(p.s, p.e, self._gain(p.node, r, p.s))])
        if p.data:
            p.coded, p.info = L[3], L[4]
        if lock or pieces[0][2] < self.floor:   # locked onto an earlier frame, or p below the sensitivity (NB, 3.4)
            p.ok, p.reason, p.raw, p.post = False, "sync" if lock else "noise", 0.5 * L[3], 0.5 * L[4]
            return
        if sp.tones:
            return self._eval_tones(p, L, ints)
        if not ints and len(pieces) == 1 and self._is_clean(L, pieces[0][2]):
            return
        cuts = sorted({x for a, e, _ in ints for x in (a, e)} | {x for a, e, _ in pieces for x in (a, e)})
        segs, i = [], 0
        for a, e in zip(cuts, cuts[1:]):
            while pieces[i][1] <= a:
                i += 1
            s_lin = 10 ** (pieces[i][2] / 10)
            itf = sum(v for x, y, v in ints if x <= a and y >= e)
            segs.append((a, e, s_lin / (1 + itf), s_lin / itf if itf and self.ci_floor else math.inf))
        out = phy.frame_outcome(self.rng, segs, L[0], sp.phy)
        p.ok, p.corrupt = out["success"], out["delivered_corrupt"]
        if not p.ok:
            p.reason = "sync" if not out["sync_ok"] else "collision" if ints else "noise"
        if p.data:
            p.raw, p.post = out["raw_bit_errors"], out["post_fec_bit_errors"]

    def _eval_tones(self, p, L, ints):
        """Multi-tone BCC stripe (DESIGN 5.2): each tone packet decided at its tone's Ec/N0 on the link (isolated, no
        contact chain), the stripe (a frame or the I-Acks) delivered only when every tone packet is. An interferer
        cannot occur in the TDMA of an isolated BAN; ponytail: one would fail the whole stripe."""
        if ints:
            p.ok, p.reason, p.raw, p.post = False, "collision", 0.5 * L[3], 0.5 * L[4]
            return
        for g, lay in zip(self.tsnr[p.node if p.node % self.m else p.rx], L[0]):
            if self._is_clean((lay,), g):
                continue
            out = phy.frame_outcome(self.rng, [(p.s, p.e, 10 ** (g / 10), math.inf)], lay, self.sp.phy)
            if not out["success"]:
                p.ok, p.reason = False, p.reason or ("noise" if out["sync_ok"] else "sync")
            if p.data:
                p.raw, p.post = p.raw + out["raw_bit_errors"], p.post + out["post_fec_bit_errors"]

    def _hears(self, r, t, p):
        """r was receiving for p when another frame started at t: its listen window for p was open (p.w: the guard
        or the IFS before p; the CSMA coordinator always) and it was not transmitting."""
        return t >= p.w and not any(x.s <= t for x in self.ntrack[r].over(t, t + 1e-6, self.maxd))

    def _is_clean(self, L, snr_db):
        """True when a lone frame at snr_db (floored to 0.5 dB) fails with probability < 1e-7 and expects < 1e-7
        payload decision errors: it then succeeds without a sampled draw (exact to 1e-7)."""
        key = (id(L[0]), math.floor(snr_db * 2))
        v = self.clean.get(key)
        if v is None:
            o = phy.frame_outcome(None, [(0.0, L[0][-1]["t1_us"], 10 ** (key[1] / 20), math.inf)], L[0],
                                  self.sp.phy, sampled=False)
            v = self.clean[key] = o["success"] > 1 - 1e-7 and o["raw_bit_errors"] < 1e-7
        return v

    def _power(self, key, t, node, skip=None, strict=False):
        """Power at node (noise units, other channels weighted by the ACR) of the transmissions on air at time t,
        except node's own and skip; strict: only those that started before t (a CCA ending at t cannot hear a frame
        that starts at t, as in the collision model)."""
        mw, nch = 0.0, self.nch
        acr = self.acrlin[key % nch]
        for x in self.atrack.over(t, t + 1e-6, self.maxd):
            f = acr[x.key % nch]
            if (f and (x.s < t if strict else x.s <= t) and x.node != node and x is not skip
                    and x.key // nch == key // nch):
                mw += f * 10 ** (self._gain(x.node, node, x.s) / 10)
        return mw

    def _wait(self, t0, t1):
        """('wait', t, t): resume a process at t (after the events before t are decided)."""
        return None

    def _sense(self, t0, t1, key, node):
        """A co-channel frame of another node starts in [t0, t1] (ANT adaptive isochronous guard window)."""
        tr = self.tracks.get(key)
        return bool(tr) and any(x.s >= t0 and x.node != node for x in tr.over(t0, t1 + 1e-6, self.maxd))

    def _cca(self, t0, t1, key, node):
        """Channel busy at the end of the CCA window: collision model any co-channel frame on air; SINR model the
        received power (adjacent channels weighted by the ACR) above the ED threshold."""
        if not self.sinr:
            tr = self.tracks.get(key)
            return bool(tr) and any(q.s < t1 and q.node != node for q in tr.over(t1, t1 + 1e-6, self.maxd))
        return self._power(key, t1, node, strict=True) > 10 ** ((ED_DBM - self.noise_dbm) / 10)

    # --- metrics ------------------------------------------------------------------------------------------------

    def _joules(self, node):
        """V x (I_tx t_tx + I_rx t_rx + ramps + wake-ups + I_sleep t_sleep) of one node (DESIGN 6), the beacon columns
        at the beacon's own currents (empty unless spec.beacon has them)."""
        E, (t_tx, t_rx, n_rt, n_rr, n_w, b_tx, b_rx, b_rt, b_rr) = self.sp.energy, self.en[node]
        B = self.sp.beacon["energy"] if self.bcol else E
        t_sleep = max(self.D - t_tx - t_rx - E["ramp_us"] * (n_rt + n_rr)
                      - b_tx - b_rx - B["ramp_us"] * (b_rt + b_rr), 0.0)
        nc = (E["tx_ma"] * t_tx + E["rx_ma"] * t_rx + E["ramp_us"] * (E["ramp_tx_ma"] * n_rt + E["ramp_rx_ma"] * n_rr)
              + E["wake_uc"] * 1e3 * n_w + E["sleep_ua"] * 1e-3 * t_sleep    # mA us = nC
              + B["tx_ma"] * b_tx + B["rx_ma"] * b_rx
              + B["ramp_us"] * (B["ramp_tx_ma"] * b_rt + B["ramp_rx_ma"] * b_rr))
        return E["vdd_v"] * nc * 1e-9

    def _occupancy(self):
        """Mean over channels (per BAN when isolated) of the busy-time fraction."""
        keys = range(self.nch * (self.k if self.iso else 1))
        fr = []
        for key in keys:
            busy, cur = 0.0, 0.0
            for q in self.tracks[key].x if key in self.tracks else ():
                a, e = max(q.s, cur), min(q.e, self.D)
                busy += max(0.0, e - a)
                cur = max(cur, e)
            fr.append(busy / self.D)
        return float(np.mean(fr))

    def rows(self, base):
        sp, D, m, sel = self.sp, self.D, self.m, self.sel
        if sp.access == "csma":   # non-beacon PAN: the coordinator receives all the time (DESIGN 4.2)
            for b in range(self.k):
                e = self.en[b * m]
                e[1], e[3] = D - e[0], e[2]
        out = []
        for b in range(self.k):
            c, lat = self.c[b], np.array(self.lat[b]) / 1000
            gen = c["generated"]    # every frame counts; one still queued at the end is not delivered (drop_end)
            c["drop_end"] = sum(1 for q in self.queues[b] for f in q if not f.delivered)
            es = [self._joules(b * m + j) for j in range(1, m)]
            eh = self._joules(b * m)
            e_s = float(np.mean(es))
            i_ma = e_s / sp.energy["vdd_v"] / (D * 1e-6) * 1e3 * BUCK
            kc = BATT["k"]["typical"] * BATT["capacity_mAh"]
            bits = c["bits_ok"]
            # x / 0: inf when x > 0 (energy spent, nothing delivered: J per bit), else undefined
            div = lambda a, z: a / z if z else math.inf if a else math.nan   # noqa: E731
            row = dict(base, network_id=b, feasible=1, period_us=self.T, r_actual_kbps=sel["r_actual_kbps"],
                       r_agg_kbps=self.n_dev * sel["r_actual_kbps"], payload_bits=sel["payload_bits"],
                       frames_per_period=sel["frames_per_period"], generated=int(gen),
                       psr=div(c["delivered_ok"], gen), pdr=div(c["delivered_ok"] + c["delivered_corrupt"], gen),
                       per_first=div(c["first_fail"], c["first"]),
                       retx_per_frame=div(c["tx_attempts"] - c["first"], c["first"]),
                       lat_mean_ms=float(lat.mean()) if lat.size else math.nan, on_time=div(c["on_time"], gen),
                       raw_ber=div(c["raw"], c["coded"]), post_ber=div(c["post"], c["info"]),
                       goodput_kbps=bits / D * 1000, airtime_frac=c["airtime"] / D, e_sensor_j=e_s, e_hub_j=eh,
                       p_sensor_uw=e_s / (D * 1e-6) * 1e6, jpb_sensor=div(sum(es), bits),
                       jpb_ban=div(sum(es) + eh, bits),
                       lifetime_days=kc / (i_ma + BATT["self_discharge_per_year"] * kc / 8760) / 24,
                       collision_frac=div(c["coll"], c["tx_attempts"]), cca_busy_frac=div(c["cca_busy"], c["cca"]))
            for col in ("tx_attempts delivered_ok delivered_corrupt drop_collision drop_noise drop_sync drop_csma "
                        "drop_retry drop_queue drop_end").split():
                row[col] = int(c[col])
            for p in (50, 95, 99):
                row[f"lat_p{p}_ms"] = float(np.percentile(lat, p)) if lat.size else math.nan
            row["lat_max_ms"] = float(lat.max()) if lat.size else math.nan
            out.append(row)
        psr = np.array([r["psr"] for r in out], float)
        jain = float(psr.sum() ** 2 / (len(psr) * (psr ** 2).sum())) if np.nansum(psr ** 2) > 0 else math.nan
        occ = self._occupancy()
        for r in out:
            r["jain_fairness"], r["channel_occupancy"] = jain, occ
        return out


def label(spec, ci_floor=True, reception="sinr"):
    """Variant label of the CSV: the variant plus the non-default flags (bcc: '+contact', '+fcs', '+opt3' for ladder
    level 3, ...), '+nocifloor' without the C/I_co floor in the SINR model (the collision model has no floor: its rows
    keep the default label)."""
    floor = "" if ci_floor or reception != "sinr" else "+nocifloor"
    return spec.variant + "".join(f"+{k}" if v is True else f"+opt{v}" if k == "opt" else f"+{k}={v}"
                                  for k, v in sorted(spec.flags.items())
                                  if v not in (None, False) or k == "touch_db" and v is not None) + floor


def simulate_trial(spec, *, k, n_dev, r_kbps, ci_mode, reception, duration_us, seed, venue="ward",
                   density="typical", case="typical", trial=0, ci_floor=True):
    """One trial of k BANs: a list of k row dicts (COLUMNS); feasible = 0 rows without metrics when select_period
    finds no period for (r_kbps, ci_mode, n_dev). ci_floor=False: no C/I_co floor in the SINR model (DESIGN 5.2).
    The BCC coupling cases are an SINR concept (a 40 / 60 dB weaker interferer) and are refused in the collision
    model, where any overlap would kill."""
    if spec.flags.get("touch") and reception == "collision":
        raise ValueError("BCC touch needs reception 'sinr': the partner body arrives touch_db weaker")
    if spec.flags.get("coupling_db") is not None and reception == "collision":
        raise ValueError("BCC coupling_db needs reception 'sinr': in the collision model any overlap kills, although "
                         "the interferer is 40 / 60 dB weaker")
    # strip_hopping: always 0 since (the no-hopping runs are the v1 archive); kept, so the raw rows of v1, v2
    # and later runs share one schema and a resumed sweep appends to them column for column
    base = dict(protocol=spec.pid, variant=label(spec, ci_floor, reception), phy=spec.phy["name"],
                reception=reception, strip_hopping=0, ci_mode=ci_mode,
                r_target_kbps=r_kbps, n_dev=n_dev, k=k, duration_us=duration_us, trial=trial, seed=seed, venue=venue,
                density=density, case=case)
    sel = pr.select_period(r_kbps, ci_mode, n_dev, spec)
    if sel is None:
        return [dict(dict.fromkeys(COLUMNS, math.nan), **base, network_id=b, feasible=0) for b in range(k)]
    sim = Sim(spec, sel, k, n_dev, reception, duration_us, seed, venue, density, case, ci_floor).run()
    return [{c: r.get(c, math.nan) for c in COLUMNS} for r in sim.rows(base)]
