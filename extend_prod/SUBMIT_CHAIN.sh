#!/bin/bash
# Submit both extension stages, B dependent on A finishing cleanly.
set -e
A=$(sbatch --parsable submit_extend_a.sh)
echo "stage A submitted: $A  (60 ns, ~2.9 days)"
B=$(sbatch --parsable --dependency=afterok:$A submit_extend_b.sh)
echo "stage B submitted: $B  (60 ns, starts only after A succeeds)"
echo
echo "Total production after both: 30 + 60 + 60 = 150 ns"
