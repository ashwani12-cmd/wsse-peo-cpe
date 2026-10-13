#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time=4-00:00:00
#SBATCH --job-name=rep2
#SBATCH --partition=small
#SBATCH --error=rep2.%J.err
#SBATCH --output=rep2.%J.out
#SBATCH --mem=0
#SBATCH --requeue

module load openmpi/4.1.4
module load gcc/11.2.0
export OMP_NUM_THREADS=1
cd $SLURM_SUBMIT_DIR

LMP=/home/IITB/multiscale-mechanics/ashwani12/software/builds/manybody-LAMMPS/lammps-29Aug2024/src/lmp_mpi

# ---- pre-flight: MPI stack and LAMMPS binary must actually load ----
# Job 427502 died 3 s in with exit 127 because a node could not read
# libopen-pal.so.40 off /home. A transient fault must bounce the job,
# not silently burn a 4-day slot. Inlined, not sourced: a shared file
# would live on the filesystem this guard exists to distrust.
preflight() {
  timeout 60 mpirun --version >/dev/null 2>&1 || return 1
  timeout 60 head -c1 "$LMP"  >/dev/null 2>&1 || return 1
  timeout 60 ldd "$LMP" 2>&1 | grep -qiE "not found|error" && return 1
  return 0
}

ok=0
for attempt in 1 2 3; do
  if preflight; then ok=1; break; fi
  echo "pre-flight attempt $attempt failed on $(hostname), retrying in 10 s" >&2
  sleep 10
done

if [ "$ok" -ne 1 ]; then
  n=${SLURM_RESTART_COUNT:-0}
  echo "$(date -Is) job=$SLURM_JOB_ID node=$(hostname) pre-flight FAILED restart_count=$n" \
    >> preflight_failures.log
  if [ "$n" -lt 3 ]; then
    echo "pre-flight failed; requeueing job $SLURM_JOB_ID onto another node" >&2
    scontrol requeue "$SLURM_JOB_ID"
    sleep 60
    exit 1
  fi
  echo "pre-flight failed after $n requeues -- giving up so this does not loop" >&2
  exit 1
fi
echo "pre-flight OK on $(hostname) at $(date -Is)" >&2

time mpirun -np $SLURM_NTASKS "$LMP" -in run_rep2.in > out_rep2.dat
