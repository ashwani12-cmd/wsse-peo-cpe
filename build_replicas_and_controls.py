import os, re, shutil
ROOT='/home/IITB/multiscale-mechanics/ashwani12/polymer'
SRC=f'{ROOT}/npt_shrunk_restart/run_restart.in'
LMP='/home/IITB/multiscale-mechanics/ashwani12/software/builds/manybody-LAMMPS/lammps-29Aug2024/src/lmp_mpi'
lines=open(SRC).read().split('\n')

# ---- harvest the force-field block from the validated input -------------
pair, order = {}, []
for ln in lines:
    m=re.match(r'\s*pair_coeff\s+(\d+)\s+(\d+)\s+lj/cut/coul/long\s+([\d.]+)\s+([\d.]+)',ln)
    if m:
        i,j=int(m.group(1)),int(m.group(2)); pair[(i,j)]=(m.group(3),m.group(4)); order.append((i,j))
bonded=[l for l in lines if re.match(r'\s*(bond_coeff|angle_coeff|dihedral_coeff)\s',l)]
assert len(pair)==12+9*3 and len(bonded)==7+12+6, (len(pair),len(bonded))
print(f"harvested {len(pair)} pair_coeff, {len(bonded)} bonded coeff lines")

def pairblock(mode):
    p=dict(pair)
    if mode=='swap':      # 11<->12 ONLY for electrolyte cross terms; sheet-internal LJ untouched
        for i in range(1,10): p[(i,11)],p[(i,12)]=pair[(i,12)],pair[(i,11)]
    elif mode=='null_se': # both faces present Se-like dispersion to the electrolyte
        for i in range(1,10): p[(i,11)]=pair[(i,12)]
    return '\n'.join(f"pair_coeff {i:2d} {j:2d}  lj/cut/coul/long  {p[(i,j)][0]}  {p[(i,j)][1]}" for i,j in order)

HEAD="""units           real
atom_style      full
boundary        p p p
processors      4 4 3
timestep        1.0

# ---- ALL styles BEFORE read_data ----
pair_style      hybrid/overlay sw/mod lj/cut/coul/long 12.0 12.0
bond_style      harmonic
angle_style     harmonic
dihedral_style  opls
improper_style  none
comm_modify     cutoff 20.0
kspace_style    pppm 1.0e-4

read_data       equil.data
special_bonds   lj/coul 0.0 0.0 0.5

pair_coeff * * sw/mod tmd.sw.mod NULL NULL NULL NULL NULL NULL NULL NULL NULL W S Se
{PAIRS}
pair_modify mix arithmetic
kspace_modify mesh 32 32 64 order 6

{BONDED}
special_bonds lj/coul 0.0 0.0 0.5

group  WSSe     type 10 11 12
group  PEO      type 1 2 3
group  Li       type 4
group  TFSI     type 5 6 7 8 9
group  polymer  type 1 2 3 4 5 6 7 8 9

neighbor        2.0  bin
neigh_modify    every 1 delay 0 check yes

fix  recenter_wsse  WSSe  recenter  INIT  INIT  75.9  units box shift all
"""

PROD="""
print ">>> NVT Production 400 K (30 ns)"
fix  nvt_prod  all  nvt  temp  400.0  400.0  100.0
compute  msd_Li  Li  msd  com yes
fix  out_msd  Li  ave/time  500  1  500  c_msd_Li[4]  file  msd_Li.dat  title1 "# Step MSD"
compute  rdf_LiO  all  rdf  200  4 3
fix  out_rdf  all  ave/time  10000  100  1000000  c_rdf_LiO[*]  file  rdf_LiO.dat  mode vector
compute  rdf_LiWSSe  all  rdf  200  4 11  4 12
fix  out_rdfW  all  ave/time  10000  100  1000000  c_rdf_LiWSSe[*]  file  rdf_LiWSSe.dat  mode vector
compute  cZ  all  chunk/atom  bin/1d  z  0.0  1.0  units box
fix  out_rhoZ  all  ave/chunk  5000  200  1000000  cZ  density/mass  file  density_z.dat
dump  d_prod  all  custom  50000  prod.lammpstrj  id mol type x y z
dump_modify  d_prod  sort id
thermo  10000
thermo_style  custom  step  temp  press  density  pe  ke  c_msd_Li[4]
run  30000000
write_data    prod_end.data
write_restart prod_end.restart
print ">>> ALL DONE"
"""

def replica(seed,tag):
    return f"""
# =====================================================================
# {tag}: independent packing via a 600 K re-melt (seed {seed}).
# Same preparation protocol as the published run (run.in stages 2-3),
# branched from the same equilibrated state with a different seed, so
# the interfacial PEO layer and the near-face ion set are a NEW
# realization. Reservoir populations stay 31/31 -- the sheet spans the
# full cross-section, so no protocol can exchange them -- but Li+ are
# indistinguishable, so a fresh configuration is what independence means
# for this observable.
# =====================================================================
print ">>> STAGE A: re-melt 400 -> 600 K (2 ns), seed {seed}"
velocity  all  create  400.0  {seed}  dist gaussian  rot yes  mom yes
fix  melt  polymer  nvt  temp  400.0  600.0  100.0
fix  spr   WSSe     spring/self  10.0
thermo  2000
thermo_style  custom  step  temp  press  density  pe  ke
run  2000000
write_data  melt.data
unfix  melt
unfix  spr

print ">>> STAGE B: cool 600 -> 400 K (2 ns)"
fix  cool  polymer  nvt  temp  600.0  400.0  100.0
fix  spr   WSSe     spring/self  10.0
run  2000000
write_data  cool.data
unfix  cool
unfix  spr

print ">>> STAGE C: re-equilibrate 400 K (5 ns)"
fix  eq  all  nvt  temp  400.0  400.0  100.0
variable  rho  equal  density
fix  rho_log  all  ave/time  5000  1  5000  v_rho  file  density_equil.dat
run  5000000
write_data    equil_rep.data
write_restart equil_rep.restart
unfix  eq
unfix  rho_log
"""

CTL_HDR="""
# =====================================================================
# {tag}
# Starts from the IDENTICAL state the published 30 ns production started
# from (npt_shrunk_restart/equil.data, timestep 10000113, velocities
# included and NOT re-drawn). Packing, momenta, topology and the
# sheet-internal LJ are bit-identical to that run. {what}
# Only the electrolyte-facing dispersion changed, so this isolates
# {isolates}
# =====================================================================
print ">>> STAGE A: re-equilibrate 400 K under the modified dispersion (2 ns)"
fix  eq  all  nvt  temp  400.0  400.0  100.0
thermo  2000
thermo_style  custom  step  temp  press  density  pe  ke
run  2000000
write_data  equil_ctl.data
unfix  eq
"""

def submit(name,infile,outfile,wall='4-00:00:00'):
    return f"""#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time={wall}
#SBATCH --job-name={name}
#SBATCH --partition=small
#SBATCH --error={name}.%J.err
#SBATCH --output={name}.%J.out
#SBATCH --mem=0
#SBATCH --requeue

module load openmpi/4.1.4
module load gcc/11.2.0
export OMP_NUM_THREADS=1
cd $SLURM_SUBMIT_DIR

LMP={LMP}

# ---- pre-flight: MPI stack and LAMMPS binary must actually load ----
# Job 427502 died 3 s in with exit 127 because a node could not read
# libopen-pal.so.40 off /home. A transient fault must bounce the job,
# not silently burn a 4-day slot. Inlined, not sourced: a shared file
# would live on the filesystem this guard exists to distrust.
preflight() {{
  timeout 60 mpirun --version >/dev/null 2>&1 || return 1
  timeout 60 head -c1 "$LMP"  >/dev/null 2>&1 || return 1
  timeout 60 ldd "$LMP" 2>&1 | grep -qiE "not found|error" && return 1
  return 0
}}

ok=0
for attempt in 1 2 3; do
  if preflight; then ok=1; break; fi
  echo "pre-flight attempt $attempt failed on $(hostname), retrying in 10 s" >&2
  sleep 10
done

if [ "$ok" -ne 1 ]; then
  n=${{SLURM_RESTART_COUNT:-0}}
  echo "$(date -Is) job=$SLURM_JOB_ID node=$(hostname) pre-flight FAILED restart_count=$n" \\
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

time mpirun -np $SLURM_NTASKS "$LMP" -in {infile} > {outfile}
"""

JOBS=[
 dict(d='replicas/rep1', name='rep1', mode=None, seed=90210,
      tag='REPLICA 1'),
 dict(d='replicas/rep2', name='rep2', mode=None, seed=31337,
      tag='REPLICA 2'),
 dict(d='ff_controls/swap', name='swap', mode='swap',
      tag='CAUSAL TEST -- S/Se electrolyte-facing dispersion SWAPPED',
      what='The LJ cross terms of types 11 and 12 with atom types 1-9 are exchanged.',
      isolates='the LJ asymmetry from the packing and the sheet corrugation:\n# if the enrichment is caused by Se dispersion it must FOLLOW the swap and\n# reverse onto the +z (S-geometry) face; if it stays put it was never\n# chemistry.'),
 dict(d='ff_controls/null_se', name='null_se', mode='null_se',
      tag='NULL CONTROL -- both faces present Se-like dispersion',
      what='The LJ cross terms of type 11 with types 1-9 are replaced by those of type 12.',
      isolates='estimator bias: the two faces are now chemically identical to the\n# electrolyte, so the measured enrichment MUST be zero. Whatever it\n# returns instead is the floor on the published 1.6x.'),
]

for j in JOBS:
    d=os.path.join(ROOT,j['d']); os.makedirs(d,exist_ok=True)
    body=HEAD.replace('{PAIRS}',pairblock(j['mode'])).replace('{BONDED}','\n'.join(bonded))
    if j['mode'] is None:
        body += replica(j['seed'], j['tag'])
    else:
        body += CTL_HDR.format(tag=j['tag'], what=j['what'], isolates=j['isolates'])
    body += PROD
    inf=f"run_{j['name']}.in"
    open(os.path.join(d,inf),'w').write(body)
    open(os.path.join(d,f"submit_{j['name']}.sh"),'w').write(
        submit(j['name'],inf,f"out_{j['name']}.dat"))
    os.chmod(os.path.join(d,f"submit_{j['name']}.sh"),0o755)
    shutil.copy(f'{ROOT}/tmd.sw.mod', d)
    shutil.copy(f'{ROOT}/npt_shrunk_restart/equil.data', d)
    print(f"  built {j['d']:24s} {inf}")

# ---- finedump: reuse the already-written input, add the guard --------
d=os.path.join(ROOT,'finedump'); os.makedirs(d,exist_ok=True)
shutil.copy(f'{ROOT}/paper_fixed/analysis/run_finedump.in', d)
shutil.copy(f'{ROOT}/tmd.sw.mod', d)
shutil.copy(f'{ROOT}/npt_shrunk_restart/prod_end.restart', d)
open(os.path.join(d,'submit_finedump.sh'),'w').write(
    submit('finedump','run_finedump.in','out_finedump.dat',wall='12:00:00'))
os.chmod(os.path.join(d,'submit_finedump.sh'),0o755)
print("  built finedump/            run_finedump.in (guard added)")
