"""Write bodies drawn from the population channel model as Castalia/OMNeT++ channel files.

usage: python export_castalia.py OUTDIR --freq-mhz F [--bodies N] [--seed S] [--paths TX-RX ...]

Per body i one file body_<i>_pathLossMap.txt in the pathLossMap format of Castalia's WirelessChannel
(Castalia 3.2 user's manual, section 4.1, parameter pathLossMapFile; Simulations/Parameters/WirelessChannel/
BANmodels): one line per transmitting node, "txId>rxId:pathloss_dB,rxId:pathloss_dB,...", '#' starts a comment.
Path loss = -gain (dB at 50 Ohm, referred to the transmit voltage) at the one frequency F (a Castalia map holds one
number per pair). The map is reciprocal: every measured path is written in both directions. Node ids: NODES below,
written to nodes.txt. Pairs within the transmitter sites or within the receiver sites were not measured and are not
written; Castalia then uses its own path-loss model for them, so give such pairs no traffic. The model is static
within a session: leave WirelessChannel.temporalModelParametersFile empty (Castalia then builds no temporal model).
pathLossMap format checked against Castalia's WirelessChannel.cc and BANmodels/pathLossMap.txt (boulis/Castalia,
GitHub, 2026-10-05); one sample map was run in OMNeT++ 4.6 with Castalia 3.3.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for d in (HERE, HERE.parent):              # the release package keeps this file one folder below popmodel.py
    sys.path.insert(0, str(d))
import popmodel as pm  # noqa: E402

TX = ("CHEST", "LOIN", "RHE", "RHT", "RHW1")
RX = ("LHW1", "LLAJ", "LLB", "LLK", "LLQ")
NODES = dict(enumerate(TX + RX))          # Castalia node id -> site code
NODE_OF = {s: i for i, s in NODES.items()}

def losses(body, f_hz):
    """{(txId, rxId): path loss dB} of a drawn body at f_hz, both directions."""
    out = {}
    for path, curve in body.items():
        tx, rx = path.split("-")
        loss = -float(pm.gain_at(curve, f_hz))
        out[(NODE_OF[tx], NODE_OF[rx])] = loss
        out[(NODE_OF[rx], NODE_OF[tx])] = loss
    return out


def write_map(path, body, f_hz, header=""):
    L = losses(body, f_hz)
    lines = [f"# {header}".rstrip()] if header else []
    for tx in sorted({a for a, _ in L}):
        lines.append(f"{tx}>" + ",".join(f"{rx}:{L[(tx, rx)]:.3f}" for (a, rx) in sorted(L) if a == tx))
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")


def read_map(path):
    """Inverse of write_map: {(txId, rxId): loss dB}."""
    L = {}
    for line in Path(path).read_text(encoding="ascii").splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        tx, rest = line.split(">")
        for item in rest.split(","):
            rx, loss = item.split(":")
            L[(int(tx), int(rx))] = float(loss)
    return L


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("outdir")
    ap.add_argument("--freq-mhz", type=float, required=True, help="frequency of the map, 0.6 to 20 MHz")
    ap.add_argument("--bodies", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--paths", nargs="+", default=pm.PATHS, help="TX-RX paths to draw (default: all 25)")
    a = ap.parse_args(argv)
    if not 0.6 <= a.freq_mhz <= 20:
        ap.error("--freq-mhz must lie in the measured range 0.6 to 20 MHz")
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    for i in range(a.bodies):
        write_map(out / f"body_{i}_pathLossMap.txt", pm.draw_body(rng, a.paths), a.freq_mhz * 1e6,
                  f"body {i}, seed {a.seed}, {a.freq_mhz:g} MHz, 50 Ohm, population model of; "
                  "loss = -gain, reciprocal; node ids in nodes.txt")
    (out / "nodes.txt").write_text("".join(f"{i} {s}\n" for i, s in NODES.items()), encoding="ascii")
    print(f"wrote {a.bodies} bodies at {a.freq_mhz:g} MHz to {out}")


if __name__ == "__main__":
    main()
