#!/bin/bash
# watch.sh - keep watching a directory 24/7 using auditd
WATCH_DIR="/home/yincy/BioHome/"
AUDIT_KEY="dir_monitor"
OUTPUT_LOG="/home/yincy/BioHome/github/functions/alignment/dir_monitor.log"
AUDIT_LOG="/home/yincy/BioHome/github/functions/alignment/dir_monitor.log"

# checks
if ! command -v auditctl &>/dev/null; then
    echo "auditd not installed. Run: sudo apt install auditd audispd-plugins"
    exit 1
fi

# add audit rule (idempotent if re-run)
echo "Adding audit rule for $WATCH_DIR ..."
sudo auditctl -w "$WATCH_DIR" -p rwxa -k "$AUDIT_KEY"

# ensure output dir exists
mkdir -p "$(dirname "$OUTPUT_LOG")"
touch "$OUTPUT_LOG"

echo "Watching $WATCH_DIR (key=$AUDIT_KEY) - appending matches from $AUDIT_LOG to $OUTPUT_LOG"
# run pipeline as root so it can read /var/log/audit/audit.log continuously
sudo bash -c "tail -n0 -F '$AUDIT_LOG' | grep --line-buffered -E 'key=$AUDIT_KEY' >> '$OUTPUT_LOG'"
