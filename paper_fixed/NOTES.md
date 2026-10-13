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

