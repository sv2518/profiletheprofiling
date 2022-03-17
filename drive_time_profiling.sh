#!/bin/bash

MYDIR=$HOME/firedrakeexamples/profiletheprofiling/
FDDIR=$HOME/firedrakeinstalls/fd060322/firedrake/

source "$FDDIR"bin/activate
cd "$FDDIR"src/PyOP2
git fetch sv/log-LACallables
git checkout sv/log-LACallables

# Take care, rerunning will remove previous results
rm "$MYDIR"time_*.txt
rm "$MYDIR"results.png
rm "$MYDIR"flame*.txt

# Run timings on commit which does not allocate ids from python
git checkout 7a6893414d9311150cddcc7c65ffda44b1fc5416
for l in $(seq 20)
do
    firedrake-clean
    python3 $MYDIR/time_profiling.py "old" $MYDIR "compile"
    python3 $MYDIR/time_profiling.py "old" $MYDIR "compute"

    firedrake-clean
    python3 $MYDIR/time_profiling.py "old" $MYDIR "compile" -log_view :$MYDIR/flame_setidfromC_compile.txt:ascii_flamegraph
    python3 $MYDIR/time_profiling.py "old" $MYDIR "compute" -log_view :$MYDIR/flame_setidfromC_compute.txt:ascii_flamegraph
done

# Run timings on commit which allocate ids from python
git checkout sv/log-LACallables
for l in $(seq 20)
do
    firedrake-clean
    python3 $MYDIR/time_profiling.py "new" $MYDIR "compile"
    python3 $MYDIR/time_profiling.py "new" $MYDIR "compute"

    firedrake-clean
    python3 $MYDIR/time_profiling.py "new" $MYDIR "compile" -log_view :$MYDIR/flame_setidfrompython_compile.txt:ascii_flamegraph
    python3 $MYDIR/time_profiling.py "new" $MYDIR "compute" -log_view :$MYDIR/flame_setidfrompython_compute.txt:ascii_flamegraph
done

# gather data
python3 $MYDIR/vis_results_time_profiling.py $MYDIR