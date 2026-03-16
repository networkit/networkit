#!/bin/bash

set -euo pipefail

: "${DOCS_CHANNEL:?DOCS_CHANNEL is required}"
: "${DOCS_LABEL:?DOCS_LABEL is required}"
: "${DOCS_BUILD_DIR:?DOCS_BUILD_DIR is required}"
: "${DOCS_PUBLISH_DIR:?DOCS_PUBLISH_DIR is required}"
: "${DOCS_BASE_PATH:?DOCS_BASE_PATH is required}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
METADATA_PATH="${DOCS_PUBLISH_DIR}/versions.json"

if [[ "${DOCS_CHANNEL}" == "stable" ]]; then
    find "${DOCS_PUBLISH_DIR}" \
        -mindepth 1 \
        -maxdepth 1 \
        ! -name '.git' \
        ! -name 'nightly' \
        ! -name 'CNAME' \
        -exec rm -rf {} +

    rsync -a "${DOCS_BUILD_DIR}/" "${DOCS_PUBLISH_DIR}/"
elif [[ "${DOCS_CHANNEL}" == "nightly" ]]; then
    mkdir -p "${DOCS_PUBLISH_DIR}/nightly"
    rsync -a --delete "${DOCS_BUILD_DIR}/" "${DOCS_PUBLISH_DIR}/nightly/"
else
    echo "Unsupported DOCS_CHANNEL: ${DOCS_CHANNEL}" >&2
    exit 1
fi

python3 "${SCRIPT_DIR}/update_docs_versions.py" \
    --file "${METADATA_PATH}" \
    --channel "${DOCS_CHANNEL}" \
    --label "${DOCS_LABEL}" \
    --base-path "${DOCS_BASE_PATH}" \
    --nightly-date "${DOCS_NIGHTLY_DATE:-}"
