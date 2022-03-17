#!/bin/bash

MYDIR=$HOME/firedrakeexamples/profiletheprofiling/
FDDIR = /Users/sv2518/firedrakeinstalls/fd060322/firedrake/

cd $FDDIRsrc/PyOP2
git fetch sv/log-LACallables
git checkout sv/log-LACallables

# Take care, rerunning will remove previous results
rm $MYDIRtime_*.txt
rm $MYDIRresults.png
rm $MYDIRflame*.txt

# Run timings on commit which does not allocate ids from python
git checkout 7a6893414d9311150cddcc7c65ffda44b1fc5416
for l in $(seq 20)
do
    python3 $MYDIR/time_profiling.py "old" $MYDIR
    python3 $MYDIR/time_profiling.py "old" $MYDIR -log_view :$MYDIR/flame_setidfromC.txt:ascii_flamegraph
done

# Run timings on commit which allocate ids from python
git checkout sv/log-LACallables
for l in $(seq 20)
do
    python3 $MYDIR/time_profiling.py "new" $MYDIR
    python3 $MYDIR/time_profiling.py "new" $MYDIR -log_view :$MYDIR/flame_setidfrompython.txt:ascii_flamegraph
done

# gather data
python3 $MYDIR/vis_results_time_profiling.py $MYDIR