#!/bin/bash
# Checks the "touch SAVE_SOLUTION" on-demand solution save.
#
# The trigger file is created BEFORE the run, so the solver finds it at its
# first check (iteration 50) and must then:
#   1) write RESULTS/SaveSolutionTouch_0000000050.hsol (autosave is disabled in
#      the control file, so this file can only come from the trigger), and
#   2) delete the trigger file.
#
# Usage: check_save_solution.sh <solver launcher and arguments...>
#   e.g. check_save_solution.sh mpiexec --oversubscribe -n 8 ./horses3d.ns SaveSolutionTouch.control
set -e

TRIGGER=SAVE_SOLUTION
EXPECTED=RESULTS/SaveSolutionTouch_0000000050.hsol

rm -f "$TRIGGER" "$EXPECTED"
mkdir -p RESULTS
touch "$TRIGGER"

"$@"

if [ -e "$TRIGGER" ]; then
   echo "FAIL: trigger file $TRIGGER was not removed by the solver"
   exit 1
fi
if [ ! -s "$EXPECTED" ]; then
   echo "FAIL: $EXPECTED was not written"
   ls -l RESULTS
   exit 1
fi
echo "OK: $EXPECTED written and $TRIGGER consumed"
