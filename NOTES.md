# WSSe / PEO-LiTFSI CPE paper — working notes

Last updated: 2026-10-06

## Status

**No fixed deadline.** Submission timing is Ashwani's call, and the decision
(2026-10-05) is to submit *when the science is ready* rather than rush an
October date. The paper is gated on two runs, not on a calendar.

Target journal: **Journal of Molecular Liquids** (Q1, IF 5.73, Elsevier).

### Progress as of 2026-10-06

**Where the paper stands: the manuscript is finished; everything open is either
compute that is queued, or the one Overleaf compile that cannot be done on the
cluster.**

#### Manuscript — done

| Item | State |
|---|---|
| Manuscript audit against raw trajectory | **done** — three fatal problems found (section 1) |
| Face analysis redone in 3D | **done** — `analyze_faces_3d.py`; 3.5x/6.2 sigma -> 1.5x/3.4 sigma |
| Bibliography rebuilt, all 14 entries verified vs Crossref | **done** — 4 pointed at the wrong paper |
| Broken `??` figure references | **done** — undefined `fig:results` label |
| New results from existing data (TFSI transport, P2 templating) | **done** — section 2 |
| `main.tex` rewritten end to end + title + preamble | **done** — passes all structural checks |
| Figures (7) | **done** — `paper_fixed/figures/` |
| Submission package (highlights, graphical abstract, declarations, cover letter) | **done** — `paper_fixed/submission/` |

#### Blocked on compute — all queued 2026-10-06 ~14:34

| Job | What | State | Data lands |
|---|---|---|---|
| **428878** | filler-free control (~1.8 d) | PENDING `(Priority)`, est. start Oct 7 ~11:46 | ~Oct 9 |
| **428882** | extension stage A, 60 ns (~2.83 d) | PENDING `(Priority)`, est. start Oct 7 ~11:46 | ~Oct 10 |
| **428883** | extension stage B, 60 ns | PENDING, `afterok:428882` | **~Oct 13-14** |

The first control attempt, **427502, failed** on 2026-10-06 03:41 — a node I/O
fault, not an input error (section 11). It took the extension down with it,
because the extension was gated on the control having started. ~11 h lost. Both
are now queued concurrently and all three scripts carry a pre-flight + requeue
guard.

#### Blocked on Ashwani — not compute

| Item | State |
|---|---|
| **Compile `main.tex` on Overleaf** | **NOT DONE** — no LaTeX on the cluster, so it has *never been compiled*. The single highest-risk open item: `[H]` float placement is unverified, with 7 figures total (3 of them new). |
| Verify JMolLiq Guide for Authors on the live page | **NOT DONE** — the page returns HTTP 403 to automated fetching, so the highlight/graphical-abstract specs used are standard-Elsevier assumptions |
| Paste `declarations.tex` into `main.tex` | **NOT DONE** — mechanical, do it at submission |
| Suggest 3-4 reviewers | **NOT DONE** — pool listed in `submission/README.md` |

#### Then, once the runs land

| Item | Depends on |
|---|---|
| `compare_control.py` -> Table 2 | 428878 |
| Re-run `analyze_faces_3d.py` on concatenated trajectory; update Table 1, abstract, conclusion with the new sigma | 428883 |
| Submit to JMolLiq | all of the above |

#### Infrastructure (does not gate the paper)

| Item | State |
|---|---|
| Pre-flight guard, logic paths | **verified 4/4** — passes clean, requeues on fault, stops at 3 restarts instead of looping |
| Pre-flight guard, on a real compute node | **VERIFIED** 2026-10-06 — job **428900**, node `rmcn042`, COMPLETED exit `0:0` in 5 s. Guard passed, did not requeue (`restart_count: 0`), 48-rank `mpirun` launched, LAMMPS `run 0` finished. |

### The one thing to decide

The extension may show the face-selectivity effect **weakening** rather than
strengthening (section 9). 3.4 sigma from non-converged error bars is exactly
the regime where a fluctuation masquerades as a signal. **Decide what the paper
claims after the Oct 13-14 data is in, not before** — and note the paper's
novelty rests on that face selectivity, because the cross-plane blocking result
is close to tautological on its own (section 8, "the cross-plane blocking
result is close to tautological").

Nothing is blocked on anything except compute, the Overleaf compile, and the
four manual submission steps above.

---

## 1. What was wrong with the original manuscript

Three problems were found by checking the manuscript's claims against the raw
trajectory. All three were fatal as written.

### 1.1 The headline Janus result was a geometry artefact

Original claim: *"Se-face accumulates ~3.5x more Li+ than S-face within 3 A
(6.2 sigma)."*

`make_paper_figure.py:153-156` assigned ions to a face using a **1D z-window**
about the mean plane of each chalcogen sublayer — not a distance to the face:

```python
ns  = np.sum((li_z > s_face)  & (li_z < s_face + rc))
nse = np.sum((li_z < se_face) & (li_z > se_face - rc))
```

The sheet ripples (1.9 A rms instantaneous, 1.5 A static), so the S and Se
sublayers' z-ranges **overlap over ~10 A** — wider than every cutoff used
except 8 A. Ions the rc=3 A criterion assigned to the Se-face sit at a true
median Li-Se distance of **5.8 A**.

Facts from the trajectory:
- **No Li+ ever comes within 3 A of either face.** Closest approach over 30 ns:
  **3.40 A (S), 3.48 A (Se)**. Li+ carries ~6 ether oxygens and cannot reach
  contact. The correct count at rc = 3 A is zero on both faces.
- The original Table 1 **changed sign with cutoff** (Se>S at 3,4 A; S>Se at
  5,6 A; Se>S at 8 A). A real face preference cannot oscillate with cutoff —
  that was the tell.

Corrected 3D result (`analyze_faces_3d.py`): Se-face enrichment is
**1.6x within 5 A at 3.4 sigma**, monotonic in sign at every cutoff.
Independently confirmed by the Li-S / Li-Se RDFs (`rdf_LiWSSe.dat`), which the
run had been writing all along and the paper never used: ratio **1.52, 3.4 sigma**.

### 1.2 D_z is not a diffusion coefficient

Original claim: *"D_xy/D_z = 2.4, establishing the impermeable WSSe layer as a
directional ion-transport feature"*, with the mechanism *"every cross-plane
trajectory requires Li+ to detour around the periodic sheet edges."*

That mechanism is impossible in this cell, for two independent reasons:
- The sheet is 40x40 A in a 40x40 A cell — it **spans the full periodic
  cross-section and has no edges**.
- The cell has **~12 A of vacuum at each z boundary** (density -> 0 below
  z~44 and above z~108). The electrolyte is a free-standing slab, so the two
  regions are not connected through the periodic image either.

Consequences:
- Cross-plane transport is **blocked, not slowed**: Li+ populations go 31/31
  -> 32/30 over 30 ns, changing in only 12 of 600 frames (~1 net crossing).
- MSD_z is **saturating toward confinement**, not diffusing: beta_z falls
  0.65 (3-10 ns) -> 0.38 (>12 ns); MSD_z reaches 57 A^2, 43% of the
  confinement limit W^2/6 ~ 131 A^2 (slab widths 22.5 / 32.7 A).
- So D_xy/D_z grows without bound with trajectory length. The original
  defence ("both components share the same convergence factor, which cancels
  in the ratio") is wrong — xy is sub-diffusive, z is saturating.

Replaced with the stronger, defensible claim: **complete cross-plane blocking**.

### 1.3 Bibliography

6 of 14 entries had `Author, A.` placeholders, 4 had `journal = {To be filled}`,
5 had `10.XXXX` DOIs. Worse, three entries pointed at the **wrong paper**:

| Key | Problem |
|---|---|
| `marchiori_electrochim_2017` | DOI `...2017.03.034` resolves to a **lead-acid battery** paper. The cited title is **Ebadi, Costa, Araujo, Brandell**, Electrochim. Acta 234, 43-51, DOI `...03.030` |
| `he_jpcb_2020_polarity` | Authors are **Gudla, Zhang, Brandell** — not "He, Xiang". Pages 8124-8131, not 7531-7544 |
| `wettstein_crosslinked_peo_2022` | Authors are **Fischer, Heuer, Diddens**; published in Macromolecules 55, 10229-10242, no longer a preprint |
| `jiang_sw_tmd_2019` | Cited a resonator paper. `tmd.sw.mod` names its own source in its header: **J.-W. Jiang, Acta Mech. Solida Sin. 32, 17-28 (2019)**, DOI `10.1007/s10338-018-0049-z` |

Two of those are **named in the body text**, so prose was corrected too:
"Marchiori et al." -> "Ebadi et al.", "Wettstein et al." -> "Fischer et al."

Also: all four results figures rendered as `??` — the text referenced a
`fig:results` label that was never defined.

---

## 2. New results added (all from the EXISTING trajectory, no new runs)

| Observable | Value | Why it matters |
|---|---|---|
| D_xy(Li+) | 4.40e-7 cm^2/s | in-plane only; 2D Einstein, MSD=4Dt |
| D_xy(TFSI-) | 1.17e-6 cm^2/s | was never computed before |
| **D_TFSI/D_Li** | **2.65** | canonical PEO/LiTFSI ordering — second independent force-field validation |
| ideal NE t+ | 0.274 | cation-limited, matches experiment |
| sigma (NE) | 6.0e-3 S/cm | upper bound only |
| beta | Li+ 0.68, TFSI- 0.92 | anion less coupled to PEO segmental motion |
| **PEO P2 near sheet** | **-0.12 to -0.23 within ~7 A** (bulk -0.038) | the "orientational template" claim, measured for the first time |
| Li-O coordination | flat at 6.1 to the interface | filler is a structural, not chemical, perturbation |

Residence kinetics were computed but **left out of the paper**: the 50 ps dump
spacing leaves the median visit unresolved, so the frequency-vs-duration split
is not established. `run_finedump.in` fixes this (5 ns at 2 ps, ~6 h).

---

## 3. Files

```
polymer/
├── PEO_paper_fixed.zip          <- upload to Overleaf
├── paper_fixed/
│   ├── main.tex                 corrected manuscript
│   ├── main.tex.ORIGINAL        original, for diffing
│   ├── refrences.bib            all 14 entries verified against Crossref
│   ├── REWRITTEN_PASSAGES.tex   the three rewritten passages, standalone
│   ├── figures/                 8 figures incl. 3 new
│   └── analysis/                all scripts + result summaries
└── control_nofiller/            SLURM job 428878 (427502 failed, resubmitted)
    ├── run_control.in
    ├── submit_control.sh
    ├── compare_control.py       run after the job finishes
    └── README.md
```

Analysis scripts (in `npt_shrunk_restart/`, copies in `paper_fixed/analysis/`):

| Script | Purpose |
|---|---|
| `traj_io.py` | shared loader, caches `prod.lammpstrj` to `traj_cache.npz` |
| `analyze_faces_3d.py` | 3D face statistics, regenerates Table 1 + fig4 |
| `analyze_transport.py` | Li+/TFSI- diffusivities, ratio, t+, sigma -> fig5 |
| `analyze_interface.py` | density, P2 chain order, coordination vs z -> fig6 |
| `analyze_residence.py` | residence kinetics -> fig7 (not in paper yet) |
| `make_fig_msd.py` | regenerates fig1 without the bogus D_z annotation |

---

## 4. Manuscript changes applied

- **Title** -> *Solvation-Mediated Janus Face Selectivity and Cross-Plane
  Transport Blocking at a WSSe Monolayer in PEO/LiTFSI Composite Polymer
  Electrolytes: An Atomistic Molecular Dynamics Study*
  (old title promised "Binding", "Anisotropic Ion Transport" and "Reinforced
  by" — none of which the paper supports)
- **Preamble**: removed `\renewcommand{\refname}{Notes and references}` (RSC
  style; Elsevier wants "References")
- **S2.1**: geometry now disclosed (full-cross-section sheet, vacuum slab,
  electrolyte = 58% of cell volume)
- **S2.3**: 1.20 g/cm^3 corrected to a cell average; condensed-phase plateau
  is 1.15-1.20 g/cm^3, which is the right comparison. Dropped the
  "stress-strain in progress" sentence
- **S2.4**: 4 -> 7 observables; states why D is in-plane only and why the 1D
  z-window criterion is rejected
- **S3.1**: new TFSI paragraph + fig:transport
- **S3.2**: retitled "In-plane transport and cross-plane blocking"; rewritten
- **S3.3**: NEW "Interfacial polymer structure" + fig:interface
- **S3.4**: Janus section rewritten (4 paragraphs)
- **Design implications / Limitations**: 3 -> 5 caveats, now including the
  non-converging error bars, the complete-blocking geometry, the absent control
- **Conclusion + Abstract**: rewritten
- n_coord made consistent at 6.13 throughout

Validated: refs/labels resolve through the `\fig`/`\Fig`/`\sect` macros, all 8
figures present, braces and environments balanced, 14/14 citations with no
orphans. **Not compiled** — no LaTeX on the cluster. Compile on Overleaf.

---

## 5. Control run (job 428878)

`prod_end.data` with the 432 WSSe atoms deleted -> 6592 atoms. WSSe carries
**zero partial charge**, so deletion is exactly charge neutral. 10 ns NVT
equilibration (void collapse) + 30 ns production. ~23 ns/day on 48 cores,
~1.8 days total.

Three bugs caught in pre-flight (`run 0` syntax test) before submitting:
1. `write_data prod_end.data` would have **overwritten its own input file** —
   outputs renamed `control_prod_end.*`
2. Styles must precede `read_data` — the data file carries Bond/Angle/Dihedral
   Coeffs sections
3. **PPPM accuracy mismatch**: left to itself PPPM chose a 12x12x20 mesh
   (rel. force accuracy 6e-5) vs the composite's forced 32x32x64 order 6
   (2e-7). Added `kspace_modify mesh 32 32 64 order 6` -> now 2.15e-7 vs the
   composite's 2.09e-7. A control must not differ in numerical accuracy.

Velocities are **inherited** from `prod_end.data` rather than re-created, to
avoid a thermal transient at t=0.

When it finishes:
```bash
cd ~/polymer/control_nofiller && python3 compare_control.py
```
-> `control_comparison.txt` and `table2_control.tex` (drop-in Table 2).

**The row that matters most is PEO P2.** The composite slab has free
polymer/vacuum surfaces, and free surfaces order chains by themselves. The
control has the same free surfaces and no sheet, so it separates the two
causes. If the control's interior is flat at -0.04, the -0.23 dip is genuinely
the filler's doing.

---

## 6. Journal plan

| Target | IF | Verdict |
|---|---|---|
| **J. Energy Storage** | ~8.9 | First choice **with** the control run. Direct genre precedent: refs Bounbaa + Rivera are both in this journal |
| **J. Molecular Liquids** | ~5.3 | Best match to what the paper is; accepts single-system MD; no control demanded as a matter of course |
| Computational Materials Science | ~3.1 | Safe fallback |
| PCCP / JPCC | ~3 | Would demand the control + a reformat out of `elsarticle` |
| Solid State Ionics | ~3.0 | Poor fit — wants a transport property we can't deliver |

`elsarticle` class already fits all the Elsevier options.

---

## 7. Remaining work

1. Compile on Overleaf; check `[H]` float placement with 3 new figures
2. Read the rewritten sections end to end — they differ substantially from the
   original and must be defensible as the author's own
3. When 428878 finishes: `compare_control.py`, add Table 2
4. Optional, strengthens further:
   - fine-dump rerun (`run_finedump.in`, ~6 h) -> residence kinetics
   - **finite platelet in a space-filling melt** (~18k atoms, ~8 ns/day, split
     into two chained jobs) -> a *real* tortuosity-limited anisotropy, which
     would rescue the original headline claim properly
   - PMF / umbrella sampling of Li+ approach to each face -> free-energy
     difference instead of 3.4 sigma counting statistics

### Honest caveat to hold onto
The control may show the composite's in-plane D_Li is **not** meaningfully
different from the filler-free system. That is a real result, not a failure —
but it would mean the contribution rests entirely on the blocking and
face-selectivity findings, and the framing needs another pass. Worth knowing
before writing the cover letter.

---

## 8. Journal decision and submission package (2026-10-05)

**Target: Journal of Molecular Liquids** — Q1, IF 5.73 (2025), Elsevier.
`elsarticle` class already matches. Chosen over J. Energy Storage (Q1, IF ~11)
because the scope fit is real rather than a stretch, and because JMolLiq does
not require the control run to be publishable.

Ruled out — **Computational Materials Science (IF 3.7) and Surfaces and
Interfaces are Q2**, so they fail the Q1 requirement. Do not re-suggest.
Electrochimica Acta / J. Power Sources are Q1 but want an electrochemical
observable this geometry cannot give.

Package in `paper_fixed/submission/`: `highlights.txt` (5 bullets, all under the
85-char Elsevier limit), `graphical_abstract.png|pdf` (1600x640, above the
1328x531 minimum), `declarations.tex` (CRediT + competing interest + data
availability), `cover_letter.txt`, `README.md` with the pre-submission checklist.

The JMolLiq Guide for Authors page blocks automated access (HTTP 403), so the
highlight limits and graphical-abstract dimensions above are **standard Elsevier
specs, not confirmed JMolLiq ones**. Verify on the live page before submitting.

### Honest assessment of publication prospects

**Genuinely novel, and defensible:**
- Janus face selectivity surviving in a solvating matrix, acting on *solvated
  cation partitioning* rather than bare-surface adsorption. No Li+ comes within
  3.4 A of either face, yet the Se face is still enriched 1.5x. This is a
  different mechanism from the gas-phase DFT literature and is the real
  contribution.
- Quantified interfacial templating (P2 = -0.23 within ~7 A vs -0.04 bulk).
  The filler literature asserts this qualitatively; measuring it is a modest
  but genuine addition.
- Two-observable force-field validation (Li-O structure *and*
  D_TFSI/D_Li = 2.65, t+ = 0.27). Better practice than most comparable papers.

**Weak points a referee will find:**
- **The cross-plane blocking result is close to tautological.** A sheet that
  spans the full periodic cross-section must block transport; showing that it
  does is not surprising. It is correct and worth stating, but it cannot carry
  the paper. Framing it as the "complete-blocking limit" helps but does not
  remove the objection.
- The face selectivity is 3.4 sigma from a single 30 ns trajectory whose
  block-averaged error does not converge. The manuscript admits this, which
  helps, but a careful referee may still ask for replicates.
- No filler-free control at the time of writing (job 428878 pending).
- The vacuum-bounded slab geometry is unusual for a bulk electrolyte study and
  will draw questions even though it is now disclosed.
- Single trajectory, no replicas, no experimental contact.

**Realistic expectation:** reasonable odds of acceptance after *major revision*.
The most predictable referee demands, in order: (1) the filler-free control,
(2) longer or replicate sampling for the face selectivity, (3) justification of
the vacuum-slab geometry. (1) is already queued and (3) is already answered in
S2.1 and the Limitations section.

**Therefore: wait for control run 428878 before submitting.** It is ~1.8 days
of compute, already in the queue, and it pre-empts the single most predictable
referee objection. Submitting without it to save a few days would be a bad
trade.

---

## 9. Production extension to 150 ns (submitted 2026-10-06 as 428882 -> 428883)

No hard deadline, so the highest-value remaining run is simply **more of the
same trajectory**. The face-selectivity result is the paper's real contribution
and its weakest statistic: tau ~ 1.33 ns for the near-interface occupancy means
30 ns gives only ~23 effective independent samples.

| total | n_eff | sigma |
|---|---|---|
| 30 ns (now) | 23 | 3.4 |
| 60 ns | 45 | ~4.8 |
| 100 ns | 75 | ~6.2 |
| **150 ns** | **113** | **~7.6** |

**Velocity-reseeded replicas would NOT work here.** Li+ barely exchanges across
the sheet (~1 crossing in 30 ns), so replicas launched from this configuration
share its initial ion partitioning and are not independent for this observable.
Continuous extension is the correct approach. (Genuinely independent statistics
would otherwise need fully rebuilt systems, which is far more expensive.)

Staged in `~/polymer/extend_prod/`: `run_extend_a.in` and `run_extend_b.in`
(60 ns each, generated from `run_restart.in` so the physics is identical),
`SUBMIT_CHAIN.sh` submits both with `--dependency=afterok`. ~2.9 days per stage.

Pre-flighted with `run 0`: resumes at step 40000113, 401 K, density 1.20,
MSD 627 A^2, PPPM relative accuracy 2.0866e-07 — matching the original
production exactly. Remaining warnings (sw unit conversion, manybody +
special_bonds, recenter ordering) are all pre-existing in the original run.

~~**Submit only after control job 427502 has started**~~ — superseded 2026-10-06.
The account has no `MaxJobs` limit set (`sacctmgr show assoc` reports none), so
the gating bought nothing and cost a day: 427502 died at 03:41 and the
extension, waiting on it, never went in. Control and extension are now queued
concurrently as 428878 / 428882 / 428883.

Analysis afterwards: dump frames are self-contained and the cadence is identical,
so `cat prod.lammpstrj prod_ext_a.lammpstrj prod_ext_b.lammpstrj > prod_all.lammpstrj`
works with no duplicate frames; delete `traj_cache.npz` and `faces_3d.npz` to
force a re-parse.

### The honest risk
Extending may show the effect **weakening** rather than strengthening. 3.4 sigma
from non-converged error bars is exactly the regime where a fluctuation can
masquerade as a signal. If it weakens, that is the most valuable thing to learn
before a referee learns it. Decide what the paper claims after the data is in.

---

## 10. Recommended order of work

1. **Control 428878** (queued) — removes the most predictable referee demand
2. **Extension to 150 ns** — queued concurrently as 428882 -> 428883; takes the
   face result to ~7.6 sigma if it is real
3. Re-run `analyze_faces_3d.py` on the concatenated trajectory, update Table 1,
   the abstract and the conclusion with the new sigma
4. `compare_control.py` -> Table 2
5. Compile on Overleaf, paste `submission/declarations.tex`, submit to JMolLiq

Optional extras, in descending value: fine-dump rerun (~6 h, residence
kinetics); finite WSSe platelet in a space-filling melt (a genuine
tortuosity-limited anisotropy — this would defuse the "blocking is tautological"
criticism and is arguably a second paper).

---

## 11. Node I/O failure and the pre-flight guard (2026-10-06)

Control job **427502 failed 3 seconds after launch** on node `rmcn211`, exit
code 127:

    mpirun: error while loading shared libraries:
    /home/apps/openmpi-4.1.4/lib/libopen-pal.so.40:
    cannot read file data: Input/output error

Not an input error — `run_control.in` was never read and `out.dat` was 0 bytes.
The library is intact (`md5sum` reads it fine from the login node), so the node
hit a transient fault reading `/home/apps` at launch. `rmcn211` was never
drained for it and is back in `allocated` state.

The real damage was the knock-on: the 150 ns extension was gated on "submit
after the control has started", so when the control died at 03:41 the extension
never went in either. **~11 hours of idle wall time on the paper's critical
path, from one bad node.**

### Fix

All three submit scripts (`submit_control.sh`, `submit_extend_a.sh`,
`submit_extend_b.sh`) now carry `#SBATCH --requeue` and an inlined pre-flight
block that runs before `mpirun`:

- probes `mpirun --version`, readability of the LAMMPS binary, and `ldd` on it
- retries 3 times, 10 s apart, so a brief NFS blip clears without a requeue
- on persistent failure, appends to `preflight_failures.log` and calls
  `scontrol requeue` so the job bounces to another node instead of burning
  the slot
- gives up after 3 requeues (`SLURM_RESTART_COUNT`) so a real breakage cannot
  requeue-loop forever

The guard is inlined in each script rather than sourced from a shared file on
purpose: a shared file would itself live on `/home` and be vulnerable to the
exact fault it guards against.

Note that a requeued job keeps its job ID, so `%J` in the `--error`/`--output`
filenames means a retry **overwrites** the previous attempt's logs. That is why
the guard appends to `preflight_failures.log`, which survives.

### Validation (2026-10-06) — the guard is tested, both paths

**Logic paths, 4/4**, run locally with `scontrol` stubbed so nothing was really
requeued:

| Case | Expected | Result |
|---|---|---|
| healthy binary | pass, no requeue | pass, exit 0 |
| missing binary | 3 retries -> requeue | requeue called |
| unreadable binary (`chmod 000`) | 3 retries -> requeue | requeue called |
| `SLURM_RESTART_COUNT=3` | fail, **no** requeue | gave up, no loop |

The last case is the one that mattered: a guard that requeues unconditionally
would be worse than no guard at all.

**On a real compute node:** job **428900**, `small`, node `rmcn042`, COMPLETED
exit `0:0` in 5 s. `pre-flight OK on rmcn042 at 2026-10-06T15:17:50`,
`restart_count: 0` (no false-positive requeue), 48-rank `mpirun` launched and
LAMMPS `run 0` completed. Test harness is `guard_test/`, and its guard block was
extracted *programmatically* from the live `submit_control.sh` and verified
byte-identical, so the real code was tested rather than a retyped copy.

Note on partitions: the first attempt (428895) was submitted to `debug`, on the
assumption that a 1-hour limit means fast turnover. **Wrong** — `debug` is the
lowest-priority partition on the shared node pool, so it waited behind every
other partition with `Reason=Nodes_required_for_job_are_DOWN,_DRAINED_or_
reserved_for_jobs_in_higher_priority_partitions` and never started. Resubmitted
to `small` with a 10-minute limit, it ran within minutes. Short wall limit does
not buy priority access; see section 12.

### Harmless warning seen in the test

All 48 ranks emit `mca_btl_openib: libosmcomp.so.3: cannot open shared object
file ... (ignored)`. OpenMPI marks it ignored and falls back. **Irrelevant for
these runs** — every job here is `-N 1` with 48 ranks, so intra-node transport
is shared memory, not InfiniBand verbs. It would matter if a future run went
multi-node, where losing the `openib` BTL could force a TCP fallback and cost
real performance. (Not confirmable as pre-existing from local logs: the only
other `.err` on hand is 427502's, which died before `mpirun` ever started.)

---

## 12. Why all three jobs are on `small`, and when they land (2026-10-06)

`small` is not a preference here — it is the **only** partition these runs fit in.

`debug`, `small`, `medium` and `large` are all **the same 482 nodes** (23136
CPUs): identical nodelists, identical 99 unusable (79 `down*` since 2026-09-28,
12 `drain`, 8 `down`). They differ *only* in wall-clock limit:

| partition | MaxTime | fits stage A (~2.83 d)? |
|---|---|---|
| `debug` | 1 h | no |
| `large` | 1 d | no |
| `medium` | 2 d | **no** |
| `small` | **4 d** | yes |

So moving a job to `medium` or `large` buys **zero** extra hardware and
guarantees a `TIMEOUT` instead of a faster start. `gpu` is also 4 d but is 30
separate, mostly-down nodes of the wrong architecture. Do not "spread the jobs
across queues" — there is only one queue's worth of hardware.

### Reading the backlog

22 jobs were pending on `small` and ours sat **last on priority** (13206), which
looks worse than it is: only **7** were genuinely competing `(Priority)`. The
other 12 were parked on limits that waiting cannot clear
(`QOSGrpBillingMinutes`, `QOSMaxJobsPerUserLimit`, `MaxJobsPerAccount`,
`JobArrayTaskLimit`) and do not overtake us when nodes free up. The 7 real
competitors are all single-node, three of them only 10 h.

Judge the queue by `REASON`, not job count, and trust
`squeue -u $USER --start` over counting rows.

### Expected timeline (Slurm's own estimate)

| Job | Start estimate | Data lands |
|---|---|---|
| 428878 control | 2026-10-07 ~11:46 | ~Oct 9 |
| 428882 ext A | 2026-10-07 ~11:46 | ~Oct 10 |
| 428883 ext B | after A succeeds | **~Oct 13-14** |

The extension is ~6 days of *serial* compute (B is chained `afterok` on A), so
that end date is queue-independent — it cannot be compressed by resubmitting
elsewhere. No idle whole nodes at the time of writing; the 69 `mix` nodes have
free cores but a 48-task `--exclusive` job needs a wholly free node.

### Open item: the pre-flight guard is not yet validated in a batch job

The section 11 guard was verified probe-by-probe on the **login node** only. If
it false-positives on a compute node it will requeue a healthy job three times
and then exit — a protection turning into a new way to lose a run. The cheap
check is a 1-hour `debug` job that loads the modules, runs the guard, and does
`run 0`. **Not yet submitted.**


---

## 13. Both stage-A jobs are running; the pre-flight guard is validated (2026-10-08)

Section 12's open item is **closed**. Both jobs started 2026-10-08 08:45 and both
`.err` files open with the guard's own success line:

```
pre-flight OK on rmcn353 at 2026-10-08T08:45:16+05:30   (428882 ext A)
pre-flight OK on rmcn239 at 2026-10-08T08:45:14+05:30   (428878 control)
```

So the guard does **not** false-positive on a compute node, and the separate
1-hour `debug` validation job is unnecessary — two real 48-task jobs have now
exercised it. The only other `.err` content is the known harmless
`mca_btl_openib` / `libosmcomp.so.3` warning (section 11), once per rank.

A was 21 h later than Slurm's section-12 estimate (Oct 7 ~11:46 → Oct 8 08:45);
shift the whole chain accordingly.

### Measured progress and revised landing times

Both at 48 cores on one node, 1 fs timestep, 7024 atoms (A) / 6592 (control).

| Job | step now | target | rate | wall used | **finishes** | limit margin |
|---|---|---|---|---|---|---|
| 428882 ext A | 49.24 M | 100.00 M | 21.1 ns/day | 10.5 h | **Oct 11 ~05:00** | 1.2 d spare of 4 d |
| 428878 control | 8.43 M | 40.00 M | 19.3 ns/day | 10.5 h | **Oct 10 ~10:30** | 1.9 d spare |
| 428883 ext B | — | 60 M | — | `afterok` on A | **~Oct 14 + queue wait** | — |

Both comfortably inside the 4-day `small` limit; the earlier ~23 ns/day cost
estimate was right to within 10%. Control is still in its 10 M-step NVT
equilibration (production starts ~Oct 8 20:50); its `out.dat` has no 7-column
MSD rows yet, which is expected, not a failure.

### MSD continuity across the restart is confirmed

`compute msd` stores its reference coordinates in an internal fix whose state
**is** carried in the restart file, so the extension continues the original time
origin rather than resetting it:

```
npt_shrunk_restart/msd_Li.dat : 40000000  625.356   (end of original prod)
extend_prod/msd_Li_a.dat      : 40000500  625.423   (first new frame)
```

Production began at step 10 M, so `t = (step - 1e7) fs`. The chain gives one
continuous 150 ns MSD (A → t=30-90 ns, B → t=90-150 ns). Nothing needs stitching.

### Do not read convergence off the live thermo MSD

`c_msd_Li[4]` in `out_a.dat` is a **single-origin** MSD over only **62 Li ions**.
At late lag it is pure noise — 2 ns blocks of it give exponents from -4.5 to
+3.2, and 3 ns sub-blocks of the new data give slopes of +4.4, **-19.8**, +34.3
A^2/ns (a negative diffusivity). The apparent "beta = 1.11 over 30-40 ns" that
falls out of it is meaningless and is **not** evidence of reaching the diffusive
regime.

This is not a criticism of the paper: `analyze_transport.py` already uses
`T.msd_sliding_xy`, a sliding-origin estimator, which is the correct one. The
point is narrower and operational: **the live `.dat` stream cannot be used to
judge whether the extension has converged.** That judgement has to wait for a
sliding-origin re-analysis of `prod_ext_a.lammpstrj` + `prod_ext_b`.

---

## 14. The paper's transport error bars are understated ~30x (2026-10-08)

Re-ran `analyze_transport.py`'s pipeline on the cached 600-frame / 30 ns
trajectory and **reproduced every central value exactly**, then replaced the
`linregress` fit error with a bootstrap over ions/molecules (2000 resamples,
same 3-10 ns window, same 2D `MSD_xy = 4Dt` fit).

| quantity | value | paper's error | **bootstrap error** | understated |
|---|---|---|---|---|
| D_Li in-plane | 4.405e-07 | +/- 1.2e-09 (0.27%) | **+/- 3.9e-08 (8.7%)** | 32x |
| D_TFSI in-plane | 1.167e-06 | +/- 2.6e-09 (0.22%) | **+/- 1.7e-07 (14.6%)** | 66x |
| D_TFSI / D_Li | 2.649 | none quoted | **+/- 0.45**, 95% [1.92, 3.70] | — |
| t+ (ideal NE) | 0.2740 | +/- 0.001 | **+/- 0.033**, 95% [0.213, 0.342] | 33x |
| sigma NE (S/cm) | 6.029e-03 | none quoted | **+/- 6.6e-04** | — |

### Why the quoted errors are wrong

`fit_D2d` returns `scipy.stats.linregress(...).stderr / 4`. That is the
regression's uncertainty *about its own line*, computed over ~140 lag points
that are (a) overlapping sliding windows and (b) all drawn from the **same 62
ions**. Neither correlation is accounted for, so the number shrinks with the
density of lag points rather than with the size of the ensemble. The physical
uncertainty is set by having **62 Li and 62 TFSI**, and 1/sqrt(62) ~ 13% is the
right order — which is what the bootstrap recovers.

Independent confirmation from splitting the trajectory in half:
D_Li = 4.91e-07 (0-15 ns) vs 4.27e-07 (15-30 ns), a 14% spread. Consistent with
8.7%, impossible with 0.27%.

### What survives and what must change in the manuscript

* **Survives.** `D_TFSI > D_Li` holds in **100.0%** of bootstrap resamples. The
  ordering, and therefore the support it lends the 0.75 charge-scaling choice,
  is robust. The qualitative story is intact.
* **Must change.** Every error bar in `transport_summary.txt` and in the text.
  `t+ = 0.274 +/- 0.001` is the dangerous one: at +/- 0.001 it reads as a
  converged prediction and a referee will treat it as one. It is
  `t+ = 0.27 +/- 0.03`, and quoting two decimals of a +/- 0.001 error on an
  ideal-NE value computed from 62 ions invites exactly the objection the
  summary's own "upper bound, not converged" caveat was trying to pre-empt.
* **Fix in code, not prose.** `fit_D2d` should return the bootstrap SD.
  `traj_io.block_sem` already exists for the autocorrelation-aware case; the
  per-ion bootstrap is the cheaper and more honest route here.

Note this is an *ensemble* limitation, so the 150 ns extension does not fix it
on its own: 5x the duration buys more independent time origins (helping roughly
as sqrt(5) ~ 2.2x at fixed lag), but 62 ions stays 62 ions. Budget for the
longer run to tighten D_Li to ~4%, not to 0.3%.

---

## 15. The model has NO Janus dipole — the face asymmetry is pure LJ (2026-10-08)

Checked the charges on the filler. **Every WSSe atom carries q = 0.0000**:

```
type 10 (W)   q=+0.0000  count=144
type 11 (S_W) q=+0.0000  count=144
type 12 (Se)  q=+0.0000  count=144
```

and no `set type ... charge` command exists in `run.in`,
`npt_shrunk/run_npt.in`, `npt_shrunk_restart/run_restart.in`,
`extend_prod/run_extend_a.in` or `control_nofiller/run_control.in`.
`run_control.in` even documents it ("The WSSe atoms carry zero partial
charge, so deletion leaves the system exactly charge neutral").

So the **only** thing distinguishing the two faces in this model is the
Lennard-Jones pair:

| face | type | epsilon (kcal/mol) | sigma (A) |
|---|---|---|---|
| S  | 11 | 0.2740 | 3.5940 |
| Se | 12 | 0.2910 | 3.7460 |

Se is 4.2% larger and 6.2% deeper. That is the entire Janus character
available to the simulation. There is **no out-of-plane dipole, no charge
asymmetry and no polarisability difference** beyond that epsilon.

### How much of a problem is this for the manuscript?

Less than it first looks, because the results section never claims a dipole
mechanism. `main.tex:251` attributes the preference to "the larger radius,
higher polarisability and softer electronic environment of Se relative to S",
and `main.tex:294` already concedes the block error does not plateau. The
dipole sentence is at `main.tex:132`, in the **introduction**, as cited
background about real Janus monolayers — not as a mechanism claim for the model.

But the gap is real: a reader who meets "permanent out-of-plane electric
dipole" in the intro and "Se-face enrichment" in the results will join them,
and a referee who greps the data file for charges will object. Two fixes, both
one sentence:

1. In the methods, state that the filler is modelled as neutral and that the
   face asymmetry enters **only** through the S/Se LJ parameters.
2. Where the preference is interpreted, say the mechanism available to the
   model is steric/dispersive — "larger radius" is captured by sigma,
   "higher polarisability" only crudely by epsilon, and the intrinsic dipole
   not at all.

This *narrows* the claim but does not weaken the novelty: the conceptual
result is that a face preference acts on **solvated** cations, and that
survives whether its origin is dispersive or electrostatic.

---

## 16. Five jobs queued: two replicas, two force-field controls, one fine dump (2026-10-08)

`small` QOS limits confirmed by `sacctmgr`: **MaxJobsPU = 5 running,
MaxSubmitPU = 10 submitted**. Five slots were free, so five jobs went in.
`small` is now at **10/10 submitted, 4/5 running** — at the cap. Note a
`scontrol requeue` from the pre-flight guard keeps the job ID and does **not**
consume a submit slot, so the guard still works at the cap.

| job | ID | dir | what it answers | ~wall |
|---|---|---|---|---|
| rep1 | **430713** | `replicas/rep1` | packing systematic | 1.9 d |
| rep2 | **430714** | `replicas/rep2` | packing systematic | 1.9 d |
| swap | **430715** | `ff_controls/swap` | is the effect causally the LJ asymmetry? | 1.6 d |
| null_se | **430716** | `ff_controls/null_se` | does the estimator read zero when it should? | 1.6 d |
| finedump | **430717** | `finedump` | residence kinetics at 2 ps resolution | 8 h |

### Why these four and not "more nanoseconds"

Section 13 showed the extension cannot fix the face result's systematic:
**60 of 62 Li never change sides in 30 ns** (19 crossing events, all from just
2 ions — flickering at the dividing surface, not transport). Each face is
sampled by a disjoint set of ~31 ions, in one realization of the packing.
Extending to 150 ns multiplies time origins, not realizations. So:

* **rep1 / rep2** — branch from the same equilibrated state
  (`npt_shrunk_restart/equil.data`, timestep 10000113) and re-melt at 600 K for
  2 ns with new velocity seeds (90210, 31337), cool 600->400 K over 2 ns,
  re-equilibrate 5 ns, then 30 ns production. This is `run.in`'s own stage-2/3
  protocol, so it is describable as "the published preparation, reseeded".
  **This is not velocity reseeding** — the earlier note that reseeding "would
  not help" was about reseeding at 400 K, which preserves the configuration.
  A 600 K melt destroys it: at 600 K Li+ moves ~15 A rms in 2 ns, comparable to
  the reservoir half-thickness, so the interfacial PEO layer and the near-face
  ion set are genuinely new. Reservoir populations stay 31/31 (geometry forbids
  otherwise) but Li+ are indistinguishable, so a fresh configuration is exactly
  what independence means here. Three packings then give a
  **between-packing error bar**, which is the one a referee accepts.

* **swap** — the causal test. The LJ cross terms of types 11 and 12 with types
  1-9 are **exchanged**; the diagonals (11-11, 12-12) and every sheet-internal
  term are left alone, so the sheet's structure and internal energetics are
  bit-identical and only the *electrolyte-facing* dispersion swaps. Starts from
  the same `equil.data` with its stored velocities **not** re-drawn, so packing
  and momenta are identical to the published run. Prediction: if the enrichment
  is caused by Se dispersion it must **follow the swap** onto the +z
  (S-geometry) face. If it stays on the Se-geometry side, it was never
  chemistry — it was packing or sheet corrugation. Confirmed near energy-neutral
  at the start, as a relabeling must be: **+5.4 kcal/mol out of 17689 (0.031%)**.
  Note this system is a numerical control, not a physically realizable
  material — do not present it as one.

* **null_se** — the estimator-bias test. The type-11 cross terms are replaced
  by type 12's, so **both** faces present Se-like dispersion and the two faces
  are chemically identical to the electrolyte. The measured enrichment must
  then be **zero**; whatever it returns instead is the noise floor under the
  published 1.6x. Starting PE shifts -16.2 kcal/mol (0.092%), correctly more
  negative since Se is both larger and deeper. This is the same class of check
  as section 14: do not trust an error bar you have not seen return zero on a
  null.

* **finedump** — `run_finedump.in` already existed and is unchanged; it only
  needed the section-11 guard wrapped around it. Clears the RESOLUTION WARNING
  in `residence_summary.txt` blocking the frequency-vs-duration decomposition.

### All five were `run 0` pre-flighted before submission

Per the standing rule. Patched copies (`run` -> `run 0`, `processors` commented
so 1 rank works) in the scratchpad, all five exit 0 with no ERROR:

```
rep1/rep2  density 1.2000023, T=400 (fresh velocities), PE -17688.985
swap       density 1.2000023, T=403.71, KE 8451.46  (velocities preserved)
null_se    density 1.2000023, T=403.71, KE 8451.46  (velocities preserved)
finedump   step 40000113, 401.2 K, PE -17615.411    (matches section 10)
```

The identical KE in swap and null_se confirms the velocities were carried, not
re-drawn — which is what makes them controls rather than replicas. rep1 and
rep2 reproduce the reference PE exactly, confirming their force field is
untouched.

Inputs were generated from `npt_shrunk_restart/run_restart.in` by script, not
retyped: 39 `pair_coeff` and 25 bonded-coeff lines are harvested and
re-emitted, so the replica force field diffs **clean** against the published
one and the controls differ in exactly 18 and 9 `pair_coeff` lines
respectively. Generator kept at `paper_fixed/analysis/`-adjacent scratch; the
per-job inputs are self-contained.

---

## 17. Script inventory and what to run when each job lands (2026-10-08)

Everything used to reach sections 13-16 is in the repo, not in a scratch dir.
Four scripts, all verified to reproduce the numbers quoted above:

| script | reproduces | run it on |
|---|---|---|
| `paper_fixed/analysis/bootstrap_transport.py` | section 14 error bars | the parent run (`npt_shrunk_restart/`) |
| `paper_fixed/analysis/analyze_reservoirs.py` | section 13 exchange + symmetry | **every** new trajectory as it lands |
| `paper_fixed/analysis/check_thermo_msd.py` | section 13 single-origin noise | any `msd_*.dat` while a run is in flight |
| `build_replicas_and_controls.py` | section 16 job generation | regenerate inputs if a protocol changes |

`analyze_reservoirs.py` and `check_thermo_msd.py` take paths as arguments, so
they work on the new jobs without editing. Example:

```
python3 paper_fixed/analysis/analyze_reservoirs.py replicas/rep1 replicas/rep2 \
        ff_controls/swap ff_controls/null_se
```

### Checklist as each job finishes

1. **finedump 430717** (~Oct 9 04:00) — re-run `analyze_residence.py` against
   `prod_fine.lammpstrj`. This either confirms or kills the
   frequency-vs-duration split; the RESOLUTION WARNING in
   `residence_summary.txt` says the 50 ps data cannot decide it. If the split
   does not survive 2 ps resolution, **drop the claim** rather than soften it.
2. **control 428878** (~Oct 10 10:30) — run `control_nofiller/compare_control.py`
   to generate Table 2. Hold the section-10 caveat in mind: if in-plane D_Li is
   not meaningfully different from the filler-free system, that is a real
   result and the framing needs another pass.
3. **rep1 430713 / rep2 430714** (~Oct 12-14) — first
   `analyze_reservoirs.py` on each, to confirm the replica reservoirs are
   symmetric to the same ~4% as the parent (section 13). A badly asymmetric
   replica must be reported, not averaged in. Then `analyze_faces_3d.py` on
   all three packings and quote the **between-packing** spread as the error on
   the enrichment. That number replaces the 3.4 sigma, and it is the one a
   referee will accept.
4. **swap 430715** (~Oct 14) — does the enrichment follow the swapped
   dispersion onto the +z face? If yes, the effect is causally the LJ
   asymmetry. If it stays on the Se-geometry side, it is packing or sheet
   corrugation and the face-selectivity claim does not survive.
5. **null_se 430716** (~Oct 14) — the measured enrichment must be zero.
   Whatever it returns is the noise floor under the published 1.6x, and it
   bounds every other face number in the paper.
6. **ext_b 428883** (~Oct 14) — concatenate with `prod_ext_a.lammpstrj`
   (identical 50 ps cadence, continuous MSD origin, section 13) and re-run
   `analyze_transport.py` **with the bootstrap error from section 14**, not
   `linregress.stderr`.

### Read 4 and 5 before 3

`swap` and `null_se` are cheaper and land at the same time, and they can
invalidate the face result outright. If `null_se` returns a non-zero
enrichment comparable to 1.6x, the estimator is biased and the replicas'
between-packing spread is measuring the wrong thing. Check the null first.

---

## 18. Status check (2026-10-09 ~19:50) and decision to skip the optional strengtheners

Live `squeue`/`sacct` check, not just notes: `small` QOS is at its 5-running
cap. Running: `control` 428878 (1d11h in), `ext_a` 428882 (1d11h in), `rep1`
430713 (22h in), `rep2` 430714 (14h in), `swap` 430715 (12h in). Queued:
`null_se` 430716 (`QOSMaxJobsPerUserLimit` — will auto-start when a slot
frees, likely when `control` lands ~Oct 10 10:30) and `ext_b` 428883
(`Dependency` — waits on `ext_a`, ETA ~Oct 11 05:00 for eligibility). Nothing
stuck or failed.

**`finedump` 430717 completed 2026-10-09 01:20 but `analyze_residence.py` has
NOT been re-run against `prod_fine.lammpstrj` yet** — `residence_summary.txt`
in both `npt_shrunk_restart/` and `paper_fixed/analysis/` still reflects the
old 50 ps data. This is the one checklist item (section 17, #1) that's ready
to run right now, not blocked on anything else landing.

**Decision: skip the two optional strengtheners (finite WSSe platelet in a
space-filling melt for a genuine tortuosity anisotropy; PMF of Li+ approach to
each face) for this submission.** Both are new-system builds, not quick
add-ons, and there's no deadline forcing them. Plan is to submit with the
current five jobs + the 150 ns extension, and only build the platelet system
if a referee specifically asks for the anisotropy number.

**Submission-ready estimate (clean path, i.e. `null_se` reads ~zero and `swap`
follows the dispersion): ~Oct 17-18.** Breakdown: all compute lands ~Oct 14
(`ext_b` is the long pole) -> 2-3 days to re-run
`bootstrap_transport.py`/`analyze_reservoirs.py`/`analyze_faces_3d.py` and
fold real numbers into Tables 1-2 -> 1 day buffer for the first-ever Overleaf
compile of `main.tex` (untested, no LaTeX on the cluster) -> half a day final
read-through of cover letter and submission package. If `null_se` or `swap`
come back bad, this is not a schedule slip — it's a reframe of the
face-selectivity claim, timeline unknown until that's seen.

---

## 19. Checklist executed, controls came back against the Janus claim, manuscript rewritten (2026-10-12)

Ran every item in section 17's checklist against the real completed jobs
(control 428878, rep1 430713, rep2 430714, swap 430715, null_se 430716,
finedump 430717 — all COMPLETED; `ext_b` 428883 still RUNNING, ~1d6h in as
of this check). Full commands and raw output are reproducible from each
run directory; summary here.

### Results, in the "read 4 and 5 before 3" order

**null_se (estimator-bias test) — FAILED to read zero.** Both faces
chemically identical by construction; RDF cross-check still gives Se:S
ratio 1.06 / 1.08 / 1.16 at rc = 5 / 6 / 8 A, the last at 4.8 sigma and
close to the published 1.18 at the same cutoff. The noise floor overlaps
the claimed effect.

**swap (causal test) — enrichment did NOT follow the swapped chemistry.**
Electrolyte-facing S/Se dispersion exchanged, same packing and momenta.
Prediction was: if real, effect follows the swap onto the geometrically-S
face. It stayed on the geometrically-Se face instead, statistically
indistinguishable from the unswapped composite at every cutoff
(1.55/1.36/1.27 vs published 1.52/1.35/1.18). Conclusion per section 16's
own stated criterion: "it was never chemistry — it was packing or sheet
corrugation."

**rep1 / rep2 (independent 600 K-repack replicas) — disagree on sign.**
rep1 matches the published direction at rc=6,8 (up to 9.9 sigma) but
reverses at rc=4,5. rep2 reverses at rc=4,5,6 (S-face enriched instead,
up to 4.5 sigma) and is flat at rc=8. rep2 also failed its own
reservoir-symmetry check (14.4% vs ~4% for everyone else) via
`analyze_reservoirs.py` — individually a weaker data point, but it is the
one pointing the wrong way, which can't be waved off.

**finedump (2 ps residence rerun) — frequency/duration split still doesn't
resolve.** Copied `traj_io.py`+`analyze_residence.py` into `finedump/`,
patched `TRAJ='prod_fine.lammpstrj'` and `DT_FRAME_NS=0.002`. Median visit
duration is still 1 frame (2 ps) at rc=5,6 — the RESOLUTION WARNING still
fires — and the frequency-vs-duration conclusion flips between cutoffs
(frequency-dominant at rc=5/6, duration-dominant at rc=7). Per section 17
item 1: dropped, not softened.

**control 428878 (filler-free) — this one is GOOD news.** Ran
`compare_control.py`. D_xy(Li+) 24.0% slower and D_xy(TFSI-) 17.8% slower
in composite vs filler-free; Li-O peak and n_coord unchanged (structural,
not chemical, obstruction confirmed); PEO P2 near-sheet -0.117 (composite)
vs control's own near-range -0.003, matching its own bulk -0.080
(templating confirmed as real, not a free-surface artefact). This
resolves the old section-10 caveat about composite D_Li possibly not
differing from filler-free — it does, by a lot.

### Net verdict

The Janus face-selectivity claim — the paper's intended headline — does
not survive its own control battery: null reads non-zero at comparable
magnitude, causal swap doesn't move the effect, independent packings
disagree on sign. This is not an error-bar problem, it's a "the effect
may not be real" problem. The filler-free control, by contrast, is clean
and strong, and supports a real, attributable obstruction + templating
result. See [[wsse-face-selectivity-controls-failed]] (memory) for the
full writeup with all numbers.

### Manuscript rewritten same day (2026-10-12) in `paper_fixed/main.tex`

Given the above, the user asked for a full reframe rather than a quiet
fix. Changes made:

- `\journal{}` corrected from a stale "ACS Applied Energy Materials" to
  "Journal of Molecular Liquids" (the actual 2026-10-05 decision — this
  was a leftover from an earlier draft, never caught before).
- **Title** changed from "Solvation-Mediated Janus Face Selectivity and
  Cross-Plane Transport Blocking..." to "Obstruction, Chain Templating,
  and the Limits of Janus Face Selectivity...: An Atomistic Molecular
  Dynamics Study with Causal and Null Controls" — leads with what
  survived, not what didn't.
- **Abstract** rewritten to lead with the filler-free control numbers,
  then states plainly that the Janus hypothesis was tested and not
  confirmed, with the null/causal/replica numbers given explicitly rather
  than hedged.
- **New `\subsection{Causal attribution via a filler-free control}`**
  (`sec:res_control`) inserted after the interfacial-structure section,
  with a new Table 2 (`tab:control`, the `compare_control.py` output) and
  four-point discussion.
- **`sec:res_janus` retitled and rewritten**: "Testing a Janus
  face-selectivity hypothesis: null, causal, and replicate controls".
  Keeps the original single-trajectory Table 1 / Fig 4 (so a referee can
  see the uncontrolled result), adds new Table 3 (`tab:controls`, the
  5-system x 4-cutoff ratio/sigma grid) and a Verdict paragraph.
- **New Figure 8** (`figures/fig8_control_battery.pdf`, generated by new
  script `paper_fixed/analysis/make_fig_control_battery.py`): 4-panel dot
  plot of Se:S ratio vs system (composite/null/swap/rep1/rep2) per
  cutoff, null line at 1, sigma annotated per point, Okabe-Ito
  colorblind-safe palette. Regenerate with
  `cd paper_fixed/analysis && python3 make_fig_control_battery.py`, numbers
  are hardcoded from Table 3 (not re-parsed from trajectories).
- Design implications, Limitations, and Conclusion all rewritten to
  match: obstruction+templating promoted to the paper's principal result,
  Janus selectivity demoted to an explicitly tested-and-not-confirmed
  hypothesis, reported as a methodological worked example (why a
  single-trajectory face count on a corrugated periodic-spanning sheet
  needs null+causal controls) rather than deleted outright.
- Added `\newcommand{\Tab}` macro (didn't exist before; `\Fig`, `\figs`,
  `\sect` did).
- Sanity-checked: brace balance, figure/table env counts, and label
  uniqueness all verified by script after every edit — no duplicate
  labels, no unmatched environments. Still **never compiled** (no LaTeX
  on the cluster — Overleaf remains the only way to check `main.tex`
  actually builds).

### What's still open (superseded by 2026-10-12 follow-up below — kept for history)

- `ext_b` 428883 still running (150 ns extension, long pole). Doesn't
  change any of the above — it's an ensemble problem, not a sampling-time
  one — but still needed for the final transport/t+ numbers via
  `bootstrap_transport.py`.
- `main.tex` still never compiled on Overleaf — now more important than
  ever to check before anything else, since this session's edits were
  large.
- Cover letter / highlights / graphical abstract in
  `paper_fixed/submission/` were NOT touched and still describe the old
  Janus-selectivity framing — need a pass once `main.tex` compiles clean.

### Follow-up session, same day (2026-10-12): submission package reconciled

- **`ext_b` 428883 progress checked** via `extend_prod/log.lammps`: run
  started at step 100,000,000, target `run 60000000` → ends at
  160,000,000. Currently at step 132,150,000 → **32.15/60 ns done
  (53.6%)**. Job had run 1d 12h49m in squeue at check time; roughly
  another ~1 day 3h expected at the same rate. Still the long pole for
  final transport/t+ numbers via `bootstrap_transport.py`.
- **`paper_fixed/submission/cover_letter.txt` rewritten** to match the
  reframed manuscript: leads with the filler-free obstruction (24%/18%
  slower Li+/TFSI- mobility, unchanged Li-O coordination) and templating
  (P2 -0.117 near sheet vs control's own -0.003) results, then states the
  Janus face-selectivity hypothesis was tested and not confirmed, giving
  the null/causal-swap/replica numbers explicitly. D_TFSI/D_Li=2.65 and
  t+=0.27 force-field validation kept (framing-independent).
- **`paper_fixed/submission/highlights.txt` rewritten** (5 bullets, all
  ≤85 chars per Elsevier limit — checked): obstruction, mobility
  reduction, templating, "hypothesis tested and not confirmed", controls
  falsify the face effect. No longer leads with "Se face holds 1.5x more
  Li+" as an established result.
- **`paper_fixed/submission/declarations.tex`**: read, left unchanged —
  CRediT/competing-interest/data-availability statements are
  framing-independent.
- **`paper_fixed/submission/README.md` item 5 updated**: no longer
  describes the filler-free control (428878) as pending — it finished
  and is already Table 2 in `main.tex`. Added a note that
  `graphical_abstract.png`/`.pdf` still need a redesign.
- **`graphical_abstract.png` inspected directly** (not assumed): its
  third panel is literally titled "Se face enriched in Li+" showing the
  1.5x/3.4σ single-trajectory result as if established — exactly the
  claim the control battery falsified. This image was NOT regenerated
  (needs a real layout decision, e.g. swapping that panel for the
  filler-free D_Li/D_TFSI comparison or a null-control inset) — flagged
  in README instead of silently left stale.

### What's still open (current)

- `ext_b` 428883: 53.6% done, ~1 day 3h remaining (see above).
- `main.tex` still never compiled on Overleaf — highest-priority item,
  needed before anything else since both the manuscript and the
  submission package changed substantially.
- `graphical_abstract.png`/`.pdf` still depict the old Janus-selectivity
  framing (see above) — needs a redesign, not yet done.
