#!/bin/bash

# choose which problem to profile
PROBLEM="problem1"

# set working repo automatically, but you need to set the path to your firedrake install
MYDIR=$HOME/firedrakeexamples/profiletheprofiling/"$PROBLEM"/
mkdir $MYDIR
FDDIR=$HOME/firedrakeinstalls/fd060322/firedrake/
source "$FDDIR"bin/activate

# setups for the different problems
if [ "$PROBLEM" = "problem1" ]; then
    # problem1 profile only the profiling of the solve and inverse callables in PyOP2
    BASE_COMMIT=7a6893414d9311150cddcc7c65ffda44b1fc5416
    PERFORM_COMMIT=sv/log-LACallables
    cd "$FDDIR"src/PyOP2
    git fetch $PERFORM_COMMIT
    git checkout $PERFORM_COMMIT
    N=20
else
    # problem2 profiles also the profiling of all local kernels (e.g. forms for TSFC)
    # so the corresponding commits in firedrake and tsfc also have to be updated
    # if code is force pushed to the branches the base commits have to be updated
    BASE_COMMITS_PTF=(3e5a09c5edbbd9739edf29a5bdbb006c8d8a5537  9a9ab94f255b277c1b6b53e6f41ea02dd9df3926 9a9ab94f255b277c1b6b53e6f41ea02dd9df3926)
    PERFORM_COMMIT=sv/profiling-local-kernels-trackevents
    cd "$FDDIR"src/PyOP2
    git fetch origin $PERFORM_COMMIT
    git checkout $PERFORM_COMMIT
    cd "$FDDIR"src/firedrake
    git fetch origin $PERFORM_COMMIT
    git checkout $PERFORM_COMMIT 
    cd "$FDDIR"src/tsfc
    git fetch origin $PERFORM_COMMIT
    git checkout $PERFORM_COMMIT
    N=10
fi

# Take care, rerunning will remove previous results
rm "$MYDIR"time_*.txt
rm "$MYDIR"results.png
rm "$MYDIR"flame*.txt

# Run timings on commit which does not allocate ids from python
if [ "$PROBLEM" = "problem1" ]; then
    cd "$FDDIR"src/PyOP2
    git checkout $BASE_COMMIT
    cd "$FDDIR"src/tsfc
    git checkout master
    cd "$FDDIR"src/firedrake
    git checkout master
else
    cd "$FDDIR"src/PyOP2
    git checkout "${BASE_COMMITS_PTF[0]}"
    cd "$FDDIR"src/tsfc
    git checkout "${BASE_COMMITS_PTF[1]}"
    cd "$FDDIR"src/firedrake
    git checkout "${BASE_COMMITS_PTF[2]}"
fi

for l in $(seq $N)
do
    firedrake-clean
    python3 $MYDIR../time_profiling.py "old" $MYDIR "compile" $PROBLEM
    python3 $MYDIR../time_profiling.py "old" $MYDIR "compute" $PROBLEM

    firedrake-clean
    python3 $MYDIR../time_profiling.py "old" $MYDIR "compile" $PROBLEM -log_view :$MYDIR/flame_setidfromC_compile.txt:ascii_flamegraph
    python3 $MYDIR../time_profiling.py "old" $MYDIR "compute" $PROBLEM -log_view :$MYDIR/flame_setidfromC_compute.txt:ascii_flamegraph
done

# Run timings on commit which allocate ids from python
if [ "$PROBLEM" = "problem1" ]; then
    cd "$FDDIR"src/PyOP2
    git checkout $PERFORM_COMMIT
    cd "$FDDIR"src/tsfc
    git checkout master
    cd "$FDDIR"src/firedrake
    git checkout master
else
    cd "$FDDIR"src/PyOP2
    git checkout $PERFORM_COMMIT
    cd "$FDDIR"src/tsfc
    git checkout $PERFORM_COMMIT
    cd "$FDDIR"src/firedrake
    git checkout $PERFORM_COMMIT
fi

for l in $(seq $N)
do
    firedrake-clean
    python3 $MYDIR../time_profiling.py "new" $MYDIR "compile" $PROBLEM
    python3 $MYDIR../time_profiling.py "new" $MYDIR "compute" $PROBLEM

    firedrake-clean
    python3 $MYDIR../time_profiling.py "new" $MYDIR "compile" $PROBLEM -log_view :$MYDIR/flame_setidfrompython_compile.txt:ascii_flamegraph
    python3 $MYDIR../time_profiling.py "new" $MYDIR "compute" $PROBLEM -log_view :$MYDIR/flame_setidfrompython_compute.txt:ascii_flamegraph
done

# gather data
python3 $MYDIR../vis_results_time_profiling.py $MYDIR