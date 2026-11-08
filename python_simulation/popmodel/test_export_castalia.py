"""Plain-assert test of export_castalia.py: a written pathLossMap reproduces the drawn gains.

Draws three bodies with a fixed seed, exports them, reads the files back and checks every pair against
-gain_at(curve, f) in both directions, plus the node table. argv: none.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for d in (HERE, HERE.parent):              # the release package keeps this file one folder below popmodel.py
    sys.path.insert(0, str(d))
import popmodel as pm  # noqa: E402
import export_castalia as ex  # noqa: E402

F = 12.5e6
with tempfile.TemporaryDirectory() as d:
    ex.main([d, "--freq-mhz", "12.5", "--bodies", "3", "--seed", "7"])
    rng = np.random.default_rng(7)
    for i in range(3):
        body = pm.draw_body(rng)
        L = ex.read_map(Path(d) / f"body_{i}_pathLossMap.txt")
        assert len(L) == 2 * len(pm.PATHS), len(L)
        for path, curve in body.items():
            tx, rx = (ex.NODE_OF[s] for s in path.split("-"))
            want = -pm.gain_at(curve, F)
            assert abs(L[(tx, rx)] - want) < 1e-3 and abs(L[(rx, tx)] - want) < 1e-3, (i, path, L[(tx, rx)], want)
    nodes = dict(line.split() for line in (Path(d) / "nodes.txt").read_text().splitlines())
    assert {int(k): v for k, v in nodes.items()} == ex.NODES
print("ok: 3 bodies, 50 pairs each, map = -gain at 12.5 MHz in both directions; node table")
