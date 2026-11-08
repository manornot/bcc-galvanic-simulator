# popmodel: population channel model of the galvanic link at 50 Ohm

The generator: whole bodies for a network simulator, every link with a gain curve over
0.6-20 MHz, the links of one body failing together as often as on the 30 measured persons. It is the simulation form
of the population channel model: the population mean, the deviation of a path from it set by
the two electrode sites, and a body term and a link term for the spread of persons and links around it.
Three forms: the TG6-form parameter table (params.json, this README), the Python generator (popmodel.py) and
Castalia/OMNeT++ channel files (export_castalia.py).

| file | what |
|---|---|
| params.json | mean table, site term, site circumferences and the generator table, 15 frequencies, dB at 50 Ohm; from the fit of the population model on the cleaned 50 Ohm curves of the open dataset |
| popmodel.py | path_mean(path, mean_db), draw_body(rng, paths, mean_db), gain_at(curve, f_hz), covariances() |
| test_popmodel.py | mean layer by hand; sample covariances of drawn bodies equal the table |
| export_castalia.py | writes per-body Castalia pathLossMap files and nodes.txt |
| test_export_castalia.py | a written map reproduces the drawn gains in both directions |

Run everything from the repository root with the thread caps of the review (numpy otherwise takes all cores):

    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python python_simulation/popmodel/test_popmodel.py
    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python python_simulation/popmodel/test_export_castalia.py

## The model in the TG6 form

The IEEE 802.15.6 TG6 channel-model document (IEEE P802.15-08-0780-09-0006, 8.2.3) writes the path loss of its
radio channels as a mean plus one normal deviation, drawn independently for every link; its human body
communication entry (8.2.2) is one transfer function with noise. This model keeps
the mean-plus-deviation form, gives the mean a site term and the deviation a body term shared by the links of one
body and a covariance over frequency. The gain y (dB; 50 Ohm receiver termination; referred to the transmit voltage)
of link l of body p between the sites a and b, a vector over the 15 tabulated frequencies, is

    y = M + (b0 + b1 log10 f) (lc - lc_ref) + u_p + e_l,   lc = log10(c_a c_b / 272),  f in MHz
    u_p = W_u z + s_u eps        once per body,  z ~ N(0, I_4), eps ~ N(0, I_15)
    e_l = W_e z' + s_e eps'      once per link,  z', eps' of the same form

Mean layer. M: population-mean table (params.json mean_db), the population mean, valid for two sites of 35.1 cm
(lc_ref = 0.6553). c_a, c_b: typical circumferences of the two sites in cm (site_circumference_cm). Site
term: b0 = 1.025 dB, b1 = 5.409 dB per decade of f, per decade of c_a c_b. Any pair of the twelve sites can
be drawn; the 25 selected paths are the default. A reference mean measured elsewhere for two sites of 35.1 cm
replaces M through the argument mean_db; the site term and the deviation layer stay.

| f (MHz) | 0.6 | 0.8 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 12.5 | 15 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| M (dB) | -60.30 | -57.79 | -55.70 | -50.56 | -48.36 | -48.03 | -45.83 | -46.01 | -45.22 | -48.25 | -46.19 | -43.58 | -42.33 | -42.83 | -44.70 |

Typical circumferences (cm): CHEST 85.2, LHW1 16.5, LLAJ 22.8, LLB 50.3, LLC 37.6, LLK 36.2, LLQ 50.2, LOIN 80.0, RHD 41.9, RHE 25.9, RHT 29.7, RHW1 16.5.

Deviation layer. W_u, W_e: four scaled eigencurves each (15 x 4, params.json body.W_db and link.W_db);
s_u^2 = 0.092 dB^2, s_e^2 = 0.084 dB^2 (the uncorrelated remainder, key nugget_db2). Four eigencurves are needed for the table form
to score as well as the full covariances on persons left out of the fit. Between the
tabulated frequencies interpolate linearly in log f. A draw is static for a body and a session, and reciprocal.

The generator table, dB:

| f (MHz) | W_u,1 | W_u,2 | W_u,3 | W_u,4 | W_e,1 | W_e,2 | W_e,3 | W_e,4 |
|---|---|---|---|---|---|---|---|---|
| 0.6 | 0.32 | 0.71 | 0.04 | -0.03 | 0.99 | 0.76 | 0.14 | 0.52 |
| 0.8 | 0.34 | 0.71 | 0.05 | -0.02 | 0.97 | 0.73 | 0.22 | 0.33 |
| 1 | 0.36 | 0.70 | 0.05 | -0.02 | 0.99 | 0.72 | 0.24 | 0.21 |
| 2 | 0.50 | 0.61 | 0.04 | 0.00 | 1.05 | 0.65 | 0.22 | -0.11 |
| 3 | 0.75 | 0.68 | 0.17 | 0.09 | 1.13 | 0.54 | 0.14 | -0.24 |
| 4 | 0.82 | 0.53 | 0.19 | 0.06 | 1.21 | 0.43 | 0.07 | -0.33 |
| 5 | 0.98 | 0.34 | 0.16 | 0.07 | 1.36 | 0.27 | -0.15 | -0.28 |
| 6 | 1.09 | 0.09 | 0.10 | -0.06 | 1.45 | 0.17 | -0.25 | -0.26 |
| 7 | 1.34 | -0.08 | 0.18 | -0.06 | 1.56 | 0.00 | -0.53 | -0.11 |
| 8 | 1.41 | -0.07 | 0.45 | -0.28 | 1.71 | -0.17 | -0.77 | 0.07 |
| 9 | 1.29 | 0.15 | -0.60 | 0.61 | 1.65 | -0.11 | -0.16 | 0.09 |
| 10 | 1.35 | -0.10 | -0.22 | 0.28 | 1.61 | -0.30 | -0.01 | 0.00 |
| 12.5 | 1.60 | -0.41 | 0.14 | 0.11 | 1.67 | -0.64 | 0.12 | 0.23 |
| 15 | 1.97 | -0.87 | 0.20 | 0.04 | 1.70 | -0.96 | 0.22 | 0.42 |
| 20 | 1.46 | 0.09 | -0.71 | -0.71 | 1.67 | -0.57 | 0.82 | -0.41 |

(params.json holds the unrounded values.)
A simulator that wants the plain TG6 form takes sigma_N(f) = sqrt(diag(W_u W_u^T + W_e W_e^T) + s_u^2 + s_e^2)
and loses the joint outage of the links of one body.

Site codes. Transmitter: CHEST (upper chest), LOIN (lower back), RHD (right shoulder), RHE (right elbow), RHT (right
upper arm, triceps), RHW1 (right wrist). Receiver: LHW1 (left wrist), LLAJ (left Achilles region), LLB (back of the
left thigh), LLC (back of the left calf), LLK (left knee), LLQ (left thigh, front). RHD and LLC were measured on few
participants and are outside the 25 selected paths.

## Python

    import numpy as np, popmodel as pm
    rng = np.random.default_rng(1)
    body = pm.draw_body(rng)                            # {path: 15 gains (dB) at pm.FREQS}, the 25 selected paths
    body = pm.draw_body(rng, ("CHEST-LLK", "RHT-LHW1"))  # only the paths a scenario uses
    g = pm.gain_at(body["CHEST-LLK"], 12.5e6)           # any frequency in 0.6-20 MHz
    m = pm.path_mean("RHD-LLC")                         # mean layer of any pair of the twelve sites

One call is one body in one session; several sensors on one body are several paths of the same call.

## Castalia/OMNeT++

    OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python \
      python_simulation/popmodel/export_castalia.py OUT --freq-mhz 12.5 --bodies 30 --seed 1

writes to OUT/:

- body_<i>_pathLossMap.txt, one per drawn body: the pathLossMap format of Castalia's WirelessChannel
  (Castalia 3.2 user's manual, section 4.1, parameter pathLossMapFile, as in Simulations/Parameters/
  WirelessChannel/BANmodels), one line per transmitting node, `txId>rxId:pathloss_dB,rxId:pathloss_dB,...`,
  `#` a comment. Path loss = -gain at the one frequency given (a map holds one number per pair), written in
  both directions (the measured channel is reciprocal). Pairs within the transmitter sites or within the
  receiver sites were not measured and are not in the file.
- nodes.txt: node id -> site code (0-4 CHEST, LOIN, RHE, RHT, RHW1; 5-9 LHW1, LLAJ, LLB, LLK, LLQ).

omnetpp.ini, one body per run (Castalia's path loss is in dB and positive):

    SN.wirelessChannel.pathLossMapFile = "body_0_pathLossMap.txt"
    SN.wirelessChannel.temporalModelParametersFile = ""
    SN.wirelessChannel.onlyStaticNodes = true

An empty temporalModelParametersFile builds no temporal model in Castalia, which is what the static-per-session
model wants. The pathLossMap format was checked against Castalia's WirelessChannel.cc parser and
BANmodels/pathLossMap.txt (boulis/Castalia, GitHub, 2026-10-05).
One sample map (body 0, every loss + 18 dB so that Castalia's 2.4 GHz BAN radio, 62 dB of loss between its
lowest transmit level and its sensitivity, separates the links) was run in OMNeT++ 4.6 with Castalia 3.3: the link
at 63.8 dB received no packet, the links at 57.0 to 60.9 dB 96 to 98 %. Castalia reads the map
and applies the loss per pair and direction. The radio is not the galvanic PHY; the run checks the file and its
use, not galvanic link performance. That run used the earlier table with 25 path means.

A pair not in the map gets Castalia's own path-loss model, so give such pairs no traffic.

## Limits

50 Ohm termination only; the body and link terms were fitted on the 25 selected paths, the site term on all 36
measured paths; no right leg, head or implant; static per body and session (no motion or contact term); 30
participants and 521 links, so the body term has wide intervals (first eigenvalue 11.2 to 29.8 dB^2);
the interpolation between the 15 frequencies is a convention, not a measurement.
The 8, 9 and 10 MHz rows of the mean table carry a depression of about 3.5 dB read as a resonance of the body in circuit with
the electrode leads and the switching unit, varying with the person; a leadless device may see it elsewhere or not at all. Version 2.2.0 repairs the 3 MHz row: in the sessions of 6 February to 19 May 2023 the 3 MHz reading carries a receiver-side step of 0.9 to 1.7 dB, replaced by log-f interpolation; the row moved by 0.8 dB in the mean and from 1.44 to 1.05 dB in the body term. The body covariance is estimated from the 29 persons with two or more links; with the
single-link person included (version 2.0.0) the body correlation at 0.6-2 MHz was 0.46-0.53 against a counted 0.18-0.20.
