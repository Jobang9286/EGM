#!/bin/bash
# Shared setup: load config, resolve paths, provide small helpers.
# Sourced by every step script; not meant to be run directly.
set -euo pipefail

HM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

hm_load_config() {
  local preset="$1"
  source "$HM_ROOT/config/common.conf"
  if [ -f "$preset" ]; then
    source "$preset"
  elif [ -f "$HM_ROOT/config/$preset.conf" ]; then
    source "$HM_ROOT/config/$preset.conf"
  else
    echo "no such preset: $preset (try: atac, h3k27me3, h3k56ac, or a file path)" >&2
    return 1
  fi
}

hm_need() {
  # Fail early and by name rather than deep inside a pipeline.
  local missing=0
  for t in "$@"; do
    command -v "$t" >/dev/null 2>&1 || { echo "missing tool: $t" >&2; missing=1; }
  done
  [ "$missing" -eq 0 ]
}

hm_peak_suffix() {
  [ "$PEAK_MODE" = "broad" ] && echo "broadPeak" || echo "narrowPeak"
}

hm_log() { echo "[$(date +%H:%M:%S)] $*"; }
