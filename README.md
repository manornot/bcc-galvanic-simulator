# Galvanic body-coupled channel: population model and simulator

Population channel model and coexistence simulator for galvanic body-coupled communication. This repository holds scripts. It holds no measurements and no results: the measurements are the open dataset at https://doi.org/10.5281/zenodo.8214497, and the simulator and scripts write their results when they run.

| folder | what it is |
|---|---|
| `python_simulation/popmodel/` | The population channel model: parameter table, Python generator, Castalia and OMNeT++ export, tests. Start here to draw channels for a body area network simulation. |
| `python_simulation/coexistence/` | The coexistence simulator: code, tests, `DESIGN.md`, `PARAMETERS.md`, `run_sweeps.sh`, `claim4.py`, `report.py`. `formal/` holds the model checking of its medium access; `tlc.sh` writes the TLC configuration of each run from its arguments. |
| `python_simulation/` (top level) | Three earlier modules that the simulator imports or compares itself with. |

The only data files are the ones the code cannot start without: `python_simulation/coexistence/parameters.json` and `channel_curves.csv` (simulator) and `params.json` (model).

`PARAMETERS.md` gives, as the source of some values, files of earlier code that are not part of this repository (`MAC_layer.py`, `PHY_layer_BCC_simulation.py`, `viz_headline.py`, `fdma_ceiling.py`, `link_budget.py`).

## Requirements

Python 3.12 or later with numpy, scipy, pandas, matplotlib and statsmodels. TLC (`tla2tools.jar`, path in `TLA2TOOLS`) and Java for the TLA+ runs in `formal/`.

## The model

    python python_simulation/popmodel/test_popmodel.py
    python python_simulation/popmodel/test_export_castalia.py

`python_simulation/popmodel/README.md` describes the parameter table and the generator.

## The simulator

Run from `python_simulation/coexistence/`. Results go to `results/`, created on the first run.

    python test_phy.py; python test_channel.py; python test_protocols.py; python test_bcc_fdma.py
    python sweep.py --help
    bash run_sweeps.sh main 4
    python report.py
    python claim4.py --workers 4

`test_coexistence.py` also compares the engine with two tables, `python bcc_fdma.py rates` and `python bcc_fdma.py isi`, each more than half an hour to write; it skips those two checks until they exist. The sweeps take hours; each is resumable. `report.py` and `claim4.py` read the sweep tables, so run the sweeps first; `report.py` leaves out table T9, which compares with the BLE landscape of the BLE model. The closed-form comparison on the allocation slot needs no sweep: `python -c "import claim4; print(claim4.part_g6())"`, about 15 minutes.

## Licence

MIT (`LICENSE`).
