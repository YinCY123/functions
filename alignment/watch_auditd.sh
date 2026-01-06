#!/usr/bin/env bash
set -euo pipefail

# Usage: sudo ./watch_auditd.sh [WATCH_DIR] [RULE_KEY] [INTERVAL_SECONDS]
# Example: sudo ./watch_auditd.sh /home/yincy/BioHome dir_monitor 10

WATCH_DIR="${1:-/home/yincy/BioHome}"
RULE_KEY="${2:-dir_monitor}"
INTERVAL="${3:-10}"

# Logs written next to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_LOG="$SCRIPT_DIR/dir_monitor.audit.log"
LAST_TS_FILE="$SCRIPT_DIR/.dir_monitor.last_ts"

trap 'echo "Received stop signal, exiting"; exit 0' INT TERM

require_tool() {
  command -v "$1" >/dev/null 2>&1 || { echo "Required tool '$1' not found. Install it and rerun." >&2; exit 2; }
}

require_tool auditctl
require_tool ausearch

if [ "$EUID" -ne 0 ]; then
  echo "This script must be run as root (sudo)." >&2
  exit 1
fi

add_rule_if_missing() {
  if ! auditctl -l | grep -F -- "$WATCH_DIR" >/dev/null 2>&1; then
    auditctl -w "$WATCH_DIR" -p wa -k "$RULE_KEY"
    echo "Added audit rule for $WATCH_DIR (key=$RULE_KEY)"
  fi
}

init_last_ts() {
  if [ ! -f "$LAST_TS_FILE" ]; then
    date -u '+%F %T' > "$LAST_TS_FILE"
  fi
}

echo "Starting auditd watcher for: $WATCH_DIR (key=$RULE_KEY)"
mkdir -p "$(dirname "$OUT_LOG")"
touch "$OUT_LOG"
init_last_ts

while true; do
  set +e
  add_rule_if_missing
  set -e

  LAST_TS="$(cat "$LAST_TS_FILE")"
  NOW="$(date -u '+%F %T')"

  # Query audit logs between LAST_TS (exclusive) and NOW (inclusive)
  EVENTS=$(ausearch -k "$RULE_KEY" -ts "$LAST_TS" -te "$NOW" -i 2>/dev/null || true)

  if [ -n "$EVENTS" ]; then
    printf "\n===== %s - events from %s to %s =====\n" "$(date '+%F %T')" "$LAST_TS" "$NOW" >> "$OUT_LOG"
    printf "%s\n" "$EVENTS" >> "$OUT_LOG"
  fi

  # update last timestamp to NOW to avoid re-processing
  echo "$NOW" > "$LAST_TS_FILE"

  sleep "$INTERVAL"
done
