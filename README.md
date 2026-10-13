# WSSe / PEO-LiTFSI composite polymer electrolyte

All-atom molecular dynamics of a Janus WSSe monolayer embedded in a molten
PEO/LiTFSI electrolyte (EO:Li = 16:1, 400 K), paired with a filler-free
reference system, two independently re-packed replicas, and two purpose-built
force-field controls — built to *test*, not assume, the causal origin of any
interfacial signal. Manuscript in preparation for *Journal of Molecular
Liquids*; source in [`paper_fixed/`](paper_fixed/).

## Current result (see `paper_fixed/main.tex` for the full writeup)

- **Established**: the WSSe sheet obstructs in-plane Li⁺/TFSI⁻ transport
  (−24%/−18% vs. the filler-free control) and templates PEO chains within
  ~7 Å, both confirmed against the matched control in
  [`control_nofiller/`](control_nofiller/). Cross-plane transport is blocked
  outright, since the sheet spans the full periodic cross-section.
- **Tested and falsified**: a Janus face-selectivity hypothesis (Se face
  preferentially enriched in Li⁺), motivated by gas-phase DFT on other
  Janus TMDs. A single production trajectory showed up to 1.5× enrichment,
  but a null control, a causal chemistry-swap control, and two repacked
  replicas (`ff_controls/`, `replicas/`) show the signal does not survive —
  see `paper_fixed/analysis/fig8_control_battery.pdf`.

## Layout

| Directory | What it is |
|---|---|
| `paper_fixed/` | Manuscript source (`main.tex`), figures, analysis scripts, and the submission package (cover letter, highlights, graphical abstract). **Start here.** |
| `npt_shrunk_restart/` | The base 30 ns composite production run — source of the published figures/tables. |
| `control_nofiller/` | Filler-free reference run (WSSe deleted, everything else identical) — isolates the filler's effect from a free-surface artefact. |
| `ff_controls/null_se/` | Null control: both WSSe faces made chemically identical. Should read zero face preference if the effect is real. |
| `ff_controls/swap/` | Causal control: S/Se electrolyte-facing dispersion exchanged. Enrichment should follow the swapped chemistry if it's real. |
| `replicas/rep1/`, `replicas/rep2/` | Two independently 600 K-repacked replicas of the composite, for sign/magnitude consistency of the face-selectivity signal. |
| `finedump/` | Short high-frequency-dump rerun (2 ps cadence) to resolve Li⁺ residence kinetics at the two faces. |
| `npt_shrunk/` | Earlier NPT equilibration/shrink stage preceding `npt_shrunk_restart/`. |
| `deform_z/`, `deform_z2/` | Uniaxial z-deformation / stress-strain runs, exploratory, not part of the current manuscript. |
| `extend_prod/` | 30→150 ns trajectory extension, originally queued to improve the face-selectivity statistic's sampling (see its own `README.md`). Since that hypothesis is now reported as falsified rather than extended, **check with Ashwani before assuming this still feeds the paper** — it predates the control-battery result that reframed the manuscript. |
| `guard_test/` | Trivial LAMMPS input used to sanity-check the SLURM pre-flight/requeue guard on a new node, not physics. |

Each run directory with analysis code carries its own copy of `traj_io.py`
and friends (not symlinked) so it's self-contained and reproducible in place.

## Reproducing a run

```bash
cd npt_shrunk_restart          # or any run directory
sbatch submit.sh                # LAMMPS production run, Param Brahma `small`/`medium` partition
python3 analyze_transport.py    # after it finishes: diffusivity, t+
python3 analyze_interface.py    # PEO templating (P2 order parameter)
python3 analyze_faces_3d.py     # face-resolved Li+ coordination (RDF-based)
```

Force field: modified Stillinger-Weber for WSSe + all-atom OPLS-AA for the
electrolyte (Li⁺ charge scaled to 0.75). LAMMPS build path and cluster quirks
are documented per-script; see `NOTES.md` for the fuller working log.

## Data note

Raw trajectory dumps (`*.lammpstrj`) and cached-numpy trajectories (`*.npz`)
are **not** tracked in this repo — several individual files exceed GitHub's
100 MB limit (up to ~560 MB). They live on Param Brahma
(`/home/IITB/multiscale-mechanics/ashwani12/polymer/`) and are regenerable
from the LAMMPS inputs and restart files that *are* tracked here. If a
permanent, citable copy of the trajectories is needed for publication,
archive them separately (e.g. Zenodo) and link the DOI here.

## Status

See `NOTES.md` for the detailed working log and `paper_fixed/NOTES.md` for
manuscript-specific history. Short version: manuscript text and figures are
done and reframed around the obstruction/templating result as of
2026-10-12; `main.tex` has not yet been compiled (no LaTeX on the cluster —
compile via Overleaf) and the submission package needs a final check against
the live JMolLiq Guide for Authors before submitting.
