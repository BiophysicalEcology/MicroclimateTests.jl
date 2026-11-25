#!/bin/sh
# Download the latest L3 NetCDF of every TERN OzFlux site into data/ as <Site>_L3.nc.
# Skips sites already on disk (any <Site>[n]_L3.nc), so it resumes and keeps existing files.
# Usage: ./download_ozflux.sh
set -u
BASE=https://dap.tern.org.au/thredds
CAT=$BASE/catalog/ecosystem_process/ozflux
DATA=$(dirname "$0")/data
mkdir -p "$DATA"
titles() { curl -sf "$1/catalog.xml" | grep -oE 'title="[^"]+"' | cut -d'"' -f2; }
for site in $(titles "$CAT"); do
    if ls "$DATA" | grep -qE "^${site}[0-9]*_L3\.nc$"; then
        echo "$site: on disk"; continue
    fi
    version=$(titles "$CAT/$site" | grep -E '^[0-9]{4}' | sort | tail -1)
    path=$(curl -sf "$CAT/$site/$version/L3/default/catalog.xml" | grep -oE 'urlPath="[^"]+_L3\.nc"' | head -1 | cut -d'"' -f2)
    if [ -z "$path" ]; then
        echo "$site: no L3 file in $version"; continue
    fi
    echo "$site: $version"
    curl -sf -o "$DATA/${site}_L3.nc.part" "$BASE/fileServer/$path" && mv "$DATA/${site}_L3.nc.part" "$DATA/${site}_L3.nc" ||
        echo "$site: download failed"
done
