#!/bin/bash
# KPP upstream sync script
# Usage: ./sync.sh /path/to/kpp/repo
set -euo pipefail

KPP_SRC="$1"
if [ ! -d "$KPP_SRC/src" ]; then
  echo "Error: $KPP_SRC does not contain KPP source tree"
  exit 1
fi

cd "$(dirname "$0")"

# Record version
TAG=$(cd "$KPP_SRC" && git describe --tags --always 2>/dev/null || echo "unknown")
COMMIT=$(cd "$KPP_SRC" && git rev-parse HEAD 2>/dev/null || echo "unknown")

cat > VERSION << EOF
# KPP upstream version anchor
KPP_UPSTREAM_TAG=$TAG
KPP_UPSTREAM_COMMIT=$COMMIT
SYNC_DATE=$(date +%Y-%m-%d)
LOCAL_CHANGES=no
EOF

echo "KPP upstream synced: $TAG ($COMMIT)"
echo "Note: This is a version anchor only. KPP source is NOT copied into the repository."
echo "To regenerate mechanism code, run KPP's bin/kpp on your .kpp file."
