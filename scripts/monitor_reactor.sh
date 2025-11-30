#!/usr/bin/env bash
# Monitor script for reactor run
# Watches PID in /tmp/reactor_persistent.pid and collects results after exit.

LOG=/tmp/reactor_persistent.log
PIDFILE=/tmp/reactor_persistent.pid
FLAG=/tmp/reactor_finished.flag
RESULT_DIR=/workspaces/MEng172-OpenMC/results
SUMMARY=$RESULT_DIR/reactor_final_summary.txt

mkdir -p "$RESULT_DIR"

if [ ! -f "$PIDFILE" ]; then
  echo "No PID file $PIDFILE found. Exiting." >&2
  exit 1
fi

PID=$(cat "$PIDFILE")
if [ -z "$PID" ]; then
  echo "PID file empty. Exiting." >&2
  exit 1
fi

echo "Monitor starting for PID $PID (log: $LOG)"

# Poll until process exits
while kill -0 "$PID" 2>/dev/null; do
  # optionally, we could inspect the log for progress messages here
  sleep 30
done

echo "Process $PID has exited. Collecting outputs..."

# Collect final log tail
if [ -f "$LOG" ]; then
  tail -n 500 "$LOG" > "$SUMMARY"
else
  echo "Log $LOG not found" > "$SUMMARY"
fi

# Copy key output files if present
cp -n depletion_results.h5 "$RESULT_DIR/" 2>/dev/null || true
cp -n summary.h5 "$RESULT_DIR/" 2>/dev/null || true
# copy any statepoint files produced
mkdir -p "$RESULT_DIR/statepoints"
cp -n statepoint.*.h5 "$RESULT_DIR/statepoints/" 2>/dev/null || true

# mark finished
touch "$FLAG"

# final status
if [ -f "$RESULT_DIR/depletion_results.h5" ] || ls "$RESULT_DIR/statepoints" 2>/dev/null | grep -q .; then
  echo "Reactor run finished and outputs collected in $RESULT_DIR" >> "$SUMMARY"
  echo "DONE" > "$FLAG"
else
  echo "Reactor run finished but no outputs found. See $SUMMARY for final log tail." >> "$SUMMARY"
  echo "DONE_NO_OUTPUT" > "$FLAG"
fi

# Keep a copy of the full log (optional)
cp -n "$LOG" "$RESULT_DIR/reactor_full.log" 2>/dev/null || true

# Run post-processing `results.py` if outputs exist
REPO_ROOT="/workspaces/MEng172-OpenMC"
RESULTS_SCRIPT="models/FusionFissionReactor/Iteration1/results.py"
if [ -f "$RESULT_DIR/depletion_results.h5" ] || ls "$RESULT_DIR/statepoints" 2>/dev/null | grep -q .; then
  echo "Running post-processing script $RESULTS_SCRIPT" >> "$SUMMARY"
  cd "$REPO_ROOT" || true
  # Run inside the `openmc` micromamba environment
  if command -v micromamba >/dev/null 2>&1; then
    micromamba run -n openmc python "$RESULTS_SCRIPT" > "$RESULT_DIR/results_run.log" 2>&1 || \
      echo "results.py failed — see $RESULT_DIR/results_run.log" >> "$SUMMARY"
  else
    echo "micromamba not found; skipping results.py run" >> "$SUMMARY"
  fi
  touch /tmp/reactor_results_done.flag
fi

# Exit
exit 0
