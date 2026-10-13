# Filler-free control run

Reference simulation for the WSSe/PEO-LiTFSI composite paper. Supplies the
one thing the manuscript currently lacks: a system identical in every
respect except the filler, so that observations about the composite can be
*attributed* to the filler rather than merely reported.

## What it is

`prod_end.data` (end of the 30 ns composite production) with the 432 WSSe
atoms deleted. 6592 atoms remain: 40 PEO 20-mers, 62 Li+, 62 TFSI-,
EO:Li = 16:1. WSSe carries zero partial charge, so the system stays exactly
charge neutral on deletion.

Identical force field, thermostat, timestep, temperature, cell and dump
cadence. One variable changes.

## Run it

    sbatch submit_control.sh

10 ns NVT equilibration (void collapse) then 30 ns NVT production at 400 K.
~23 ns/day on 48 cores, so ~1.8 days total -- one job on the `small`
partition, inside the 4-day limit.

## Then analyse

    python3 compare_control.py

Emits `control_comparison.txt` and `table2_control.tex` (drop-in Table 2).

The composite-side scripts also run here unchanged, since the output file
names and the 50 ps dump spacing match:

    cp ../npt_shrunk_restart/{traj_io,analyze_transport,analyze_interface}.py .
    python3 analyze_transport.py
    python3 analyze_interface.py

Face-resolved analyses do not apply -- there is no sheet.

## What to expect, and what each outcome means

| Observable | Expectation | If it holds |
|---|---|---|
| Li-O peak, n_coord | unchanged (2.13 A, 6.1) | filler is a structural, not chemical, perturbation -- currently asserted, would become established |
| D_TFSI/D_Li | ~2.6 in both | validates the charge-scaled force field against neat PEO/LiTFSI, free of filler effects |
| D_xy(Li+) | composite similar or slightly lower | quantifies in-plane obstruction by the filler |
| PEO P2 near sheet | composite ~-0.2, control flat ~-0.04 everywhere | templating is caused by the filler, not an artefact of the slab geometry |

The last row is the important one. In the composite, P2 dips near the sheet;
but the composite slab also has free surfaces, which can order chains by
themselves. The control has the same free surfaces and no sheet, so it
separates the two causes. If the control is flat in its interior, the dip is
the filler's doing.
