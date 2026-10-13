# Submission package — Journal of Molecular Liquids

| File | Use |
|---|---|
| `highlights.txt` | 5 bullets, each <= 85 characters (Elsevier limit). Paste into the online system or upload as a separate "Highlights" file. |
| `graphical_abstract.png` / `.pdf` | 1600 x 640 px, above the 1328 x 531 minimum. Upload as "Graphical Abstract". |
| `declarations.tex` | CRediT, competing interest, data availability. Paste into `main.tex` before `\section*{Acknowledgments}`. |
| `cover_letter.txt` | Paste into the submission system. |

## Before submitting

1. Compile `main.tex` on Overleaf — it has never been compiled (no LaTeX on
   the cluster). Check float placement with the three new figures.
2. Paste `declarations.tex` into `main.tex`.
3. **Verify the current Guide for Authors.** The live page blocks automated
   access, so the requirements above are the standard Elsevier ones. Confirm
   the highlight count/length, graphical-abstract dimensions and whether
   JMolLiq wants a structured or unstructured abstract:
   https://www.sciencedirect.com/journal/journal-of-molecular-liquids/publish/guide-for-authors
4. Suggest 3-4 reviewers if asked. Sensible pool from the reference list:
   Brandell / Gudla (charge-scaled PEO-LiTFSI MD), Diddens / Heuer (PEO
   transport modelling), C.V. Nguyen (Janus WSSe DFT). Do not suggest anyone
   you have collaborated with.
5. The filler-free control (428878) finished and is already incorporated
   as Table 2 (`sec:res_control`) in `main.tex` — the paper was reframed
   around it on 2026-10-12 after the face-selectivity controls failed.
   `graphical_abstract.png`/`.pdf` still depict the old Janus-selectivity
   framing and need to be redesigned around the obstruction + templating
   result before submission.
