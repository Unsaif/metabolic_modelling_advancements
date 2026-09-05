#!/bin/bash
# collect_browser_download.sh <dest-relative-to-data> [expected-bytes]
# Copies the newest not-yet-collected download made by the Claude desktop browser pane (which lands as a hidden
# temp file ~/Downloads/.<id>.com.anthropic.claudefordesktop.<6 chars>) into data/<dest> of the project folder and
# records the temp id in data/fitness_browser/used_temp_ids.txt so it is never collected twice. Runs inside the
# Cowork device shell where the project folder is mounted at $HOME/mnt/metabolic_modelling_advancements.
# With an expected size the copy is refused on mismatch; without one it waits until the file size is stable for 6 s.
R="$HOME/mnt/metabolic_modelling_advancements"; D="$HOME/mnt/Downloads"; USED="$R/data/fitness_browser/used_temp_ids.txt"; touch "$USED"
dest="$R/data/$1"; mkdir -p "$(dirname "$dest")"
for try in $(seq 1 26); do
  f=$(ls -t "$D"/.*claudefordesktop* 2>/dev/null | while read x; do id="${x##*.}"; grep -q "$id" "$USED" || { echo "$x"; break; }; done)
  if [ -n "$f" ]; then
    s1=$(stat -c %s "$f")
    if [ -n "$2" ]; then [ "$s1" = "$2" ] || { echo "SIZE MISMATCH $s1 != $2 ($f)"; exit 1; }; else sleep 6; s2=$(stat -c %s "$f"); [ "$s1" = "$s2" ] || continue; fi
    [ "$s1" -gt 0 ] && cp "$f" "$dest" && echo "${f##*.}" >> "$USED" && echo "OK $1 <- ${f##*.} ($s1 bytes)" && exit 0
  else sleep 6; fi
done; echo "FAILED $1"; exit 1
