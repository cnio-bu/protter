#!/usr/bin/env bash

# bash strict mode - snakemake.readthedocs.io/en/stable/project_info/faq.html
set -euo pipefail
IFS=$'\n\t'

log () {
  echo "[$(date --utc +%Y-%m-%dT%H:%M:%SZ)] $1"
}

usage_error () {
  echo "usage: $0 file_url file_path" >&2
  exit 1
}

if [[ $# -ne 2 ]]
then usage_error
fi

in_file_url=$1
out_file_path=$2

if [[ -z "$in_file_url" || -z "$out_file_path" ]]
then usage_error
fi

out_file_name=$(basename "$out_file_path")

log "download process starting"

# Download file to temp directory to help
# reduce occurrence of partial downloads.
TWD=$(mktemp -d)
log "preparing temp directory: '$TWD'"
trap "{ rm -rf $TWD; }" EXIT

log "downloading file: '$in_file_url'"
tmp_file_path="$TWD/$out_file_name"
wget --quiet --limit-rate=10m --random-wait --wait=97 \
  "$in_file_url" -O "$tmp_file_path"

log "moving downloaded file to: '$out_file_path'"
out_dir=$(dirname "$out_file_path")
mkdir -p "$out_dir"
mv "$tmp_file_path" "$out_file_path"

log "download process done"
