"""Plain-assert test of popmodel.py: the mean layer follows Eq. 4.5 and drawn bodies have the covariances of the table.

Mean layer: at 1 MHz (log10 f = 0) the site slope is b0, so the both-wrist path RHW1-LHW1 (16.5 cm each) has
mean = M(1 MHz) + b0 (log10(16.5^2 / 272) - lc_ref), computed here by hand from params.json; a reference mean passed
as mean_db shifts the mean layer by exactly the difference of the two means. Covariances: two links of one body, the
covariance of one link's deficit (y - mean layer) must be Su + Se and the cross-covariance of the two links Su (the body
term), Su and Se built from the W and nugget entries of params.json. Tolerance: 5 standard errors of a sample
covariance, sqrt((S_ii S_jj + S_ij^2) / n). Runtime a few seconds.
argv: none (python python_simulation/popmodel/test_popmodel.py).
"""
import math
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import popmodel as pm   # noqa: E402

t, c = pm.PARAMS["site_term"], pm.PARAMS["site_circumference_cm"]
i1 = list(pm.FREQS).index(1e6)
want = pm.PARAMS["mean_db"][i1] + t["b0_db"] * (math.log10(c["RHW1"] * c["LHW1"] / 272.0) - t["lc_ref"])
assert abs(pm.path_mean("RHW1-LHW1")[i1] - want) < 1e-9, (pm.path_mean("RHW1-LHW1")[i1], want)
shift = np.linspace(-1, 1, 15)
moved = pm.path_mean("CHEST-LLK", np.asarray(pm.PARAMS["mean_db"]) + shift) - pm.path_mean("CHEST-LLK")
assert np.allclose(moved, shift)
print("ok mean layer: site term at 1 MHz by hand, reference mean plugged in")

N = 20000
rng = np.random.default_rng(1)
a, b = pm.PATHS[0], pm.PATHS[-1]
D = np.array([[curve - pm.path_mean(k) for k, curve in pm.draw_body(rng, (a, b)).items()] for _ in range(N)])
Su, Se = pm.covariances()
for got, want, var, name in ((np.cov(D[:, 0].T), Su + Se, Su + Se, "link"),
                             (np.cov(D[:, 0].T, D[:, 1].T)[:15, 15:], Su, Su + Se, "cross-link (body term)")):
    se = np.sqrt((np.outer(np.diag(var), np.diag(var)) + want ** 2) / N)
    assert np.all(np.abs(got - want) <= 5 * se), (name, np.max(np.abs(got - want) / se))
    print(f"ok {name}: max |sample - table| = {np.max(np.abs(got - want) / se):.1f} standard errors")
assert abs(pm.gain_at(np.arange(15.0), pm.FREQS[3]) - 3) < 1e-12
print("ok interpolation at a tabulated frequency")
