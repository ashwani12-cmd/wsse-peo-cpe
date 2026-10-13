# Production extension: 30 ns -> 150 ns

## Why

The face-selectivity result is the paper's real contribution and its weakest
statistic: 1.5x Se enrichment at **3.4 sigma**, from a near-interface occupancy
whose integrated autocorrelation time is tau ~ 1.33 ns. Over 30 ns that is only
~23 effective independent samples, and the block-averaged error does not
plateau with block length — i.e. the run is too short to converge this
observable.

Sampling scales as sqrt(t):

| total | n_eff | sigma becomes |
|---|---|---|
| 30 ns (now) | 23 | 3.4 |
| 60 ns | 45 | ~4.8 |
| 100 ns | 75 | ~6.2 |
| **150 ns** | **113** | **~7.6** |

Extending is the cleanest fix. Note that velocity-reseeded *replicas* would NOT
work here: Li+ barely exchanges across the sheet (~1 crossing in 30 ns), so
replicas launched from this configuration would share its initial ion
partitioning and would not be independent for this observable. Continuous
extension of the same trajectory is the right move.

## Run

    ./SUBMIT_CHAIN.sh

Two chained 60 ns stages (B waits on A via `--dependency=afterok`), each ~2.9
days at 21.2 ns/day — each fits inside the 4-day `small` limit.

**Submitted 2026-10-06 as 428882 (stage A) -> 428883 (stage B).**

The earlier rule here was "submit only after the control has started", to avoid
MaxJobsPerAccount delaying the control. That gating was wrong twice over: the
account has no `MaxJobs` limit set, and control job 427502 then died 3 s in on
a node I/O fault, which silently held the extension back for ~11 h. Control
(428878) and both extension stages now run concurrently.

## Then analyse

LAMMPS dump frames are self-contained and the cadence is identical (50 ps), so
the trajectories concatenate directly with no duplicate frames:

    cat ../npt_shrunk_restart/prod.lammpstrj \
        prod_ext_a.lammpstrj prod_ext_b.lammpstrj > prod.lammpstrj
    cp ../npt_shrunk_restart/{traj_io,analyze_faces_3d,analyze_transport,analyze_interface}.py .
    rm -f traj_cache.npz faces_3d.npz      # force a re-parse
    python3 analyze_faces_3d.py

Expect ~3000 frames and a few minutes of parsing.

## Honest caveat

Extending may show the effect **weakening**, not strengthening — 3.4 sigma from
non-converged errors is exactly the regime where a fluctuation can look like a
signal. If it weakens, that is the single most valuable thing to learn before a
referee learns it. Decide what the paper claims after the data is in, not
before.
