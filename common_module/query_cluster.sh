#!/usr/bin/env bash
# ============================================================
#  query_cluster.sh — Probe a SLURM cluster and write cluster.conf
# ============================================================
#  Usage:
#      bash query_cluster.sh [partition_name] [output_path]
#
#  * partition_name — optional; force a specific partition as the
#    GPU/default target.  If the name does not exist on the cluster,
#    the script warns and exits without writing any config file.
#  * output_path    — optional; defaults to cluster.conf next to this
#                     script.
#
#  Must be run on (or SSH'd into) a SLURM-enabled host.
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORCE_PARTITION="${1:-}"
OUTPUT="${2:-${SCRIPT_DIR}/cluster.conf}"

die() { echo "[ERROR] $*" >&2; exit 1; }

# ---- Helpers -------------------------------------------------------------
run_cmd() {
    local out
    out=$(eval "$1" 2>/dev/null) || true
    echo "$out"
}

# Extract GPU model from a GRES string.
# Handles formats:
#   "a100:8"              → a100
#   "gpu:a100:4"          → a100
#   "gpu:v100:4(S:0)"     → v100
#   "gpu"                 → (empty — generic, not a model)
extract_gpu_model() {
    local gres="$1"
    IFS=':' read -ra parts <<< "$gres"
    for part in "${parts[@]}"; do
        # Strip metadata like (S:0-1) first
        clean=$(echo "$part" | sed 's/([^)]*)//g; s/^[[:space:]]*//; s/[[:space:]]*$//')
        local lower="${clean,,}"
        # Skip generic/non-model tokens
        if [[ -z "$lower" ]] || [[ "$lower" == "gpu" ]] || [[ "$lower" == "(null)" ]] || [[ "$lower" == "none" ]]; then
            continue
        fi
        # Strip trailing digits and check against known models
        stripped=$(echo "$clean" | sed 's/[0-9]*$//')
        local lower_stripped="${stripped,,}"
        case "$lower_stripped" in
            *a100*|*h100*|*v100*|*a6000*|*a40*|*a1000*|*l40*|*l4*|*t4*)
                echo "$stripped"; return 0 ;;
        esac
        # Also try without stripping (for names like a100_1g)
        case "$lower" in
            *a100*|*h100*|*v100*|*a6000*|*a40*|*a1000*|*l40*|*l4*|*t4*)
                echo "$clean"; return 0 ;;
        esac
    done
    return 1
}

# ---- SLURM version -------------------------------------------------------
SLURM_VER=$(run_cmd "sinfo --version" | head -1)
if [ -z "$SLURM_VER" ]; then
    die "Cannot run 'sinfo --version' — is SLURM installed?"
fi

echo "[*] SLURM: ${SLURM_VER}"

# ---- Partitions ----------------------------------------------------------
declare -a PART_NAMES=()
declare -A PART_TIME=()
declare -A PART_NODES=()
declare -A PART_STATE=()
declare -A PART_GRES=()
declare -A PART_ALLOWED_QOS=()
DEFAULT_PARTITION=""

while IFS= read -r line; do
    line="${line%\\}"
    line="${line//\\ / }"

    if [[ "$line" =~ ^PartitionName=([A-Za-z0-9_-]+) ]]; then
        pname="${BASH_REMATCH[1]}"
        PART_NAMES+=("$pname")
        PART_NODES["$pname"]=""
        PART_TIME["$pname"]="UNLIMITED"
        PART_STATE["$pname"]="UNKNOWN"
        PART_GRES["$pname"]=""
        PART_ALLOWED_QOS["$pname"]=""
    fi
    if [[ "$line" =~ TotalNodes=([0-9]+) ]]; then
        PART_NODES["$pname"]="${BASH_REMATCH[1]}"
    fi
    if [[ "$line" =~ ^Time=([A-Za-z0-9:-]+) ]]; then
        PART_TIME["$pname"]="${BASH_REMATCH[1]}"
    fi
    if [[ "$line" =~ State=([A-Z]+) ]]; then
        PART_STATE["$pname"]="${BASH_REMATCH[1]}"
    fi
    if [[ "$line" =~ GresTypes=(.+)$ ]]; then
        gres_val="${BASH_REMATCH[1]}"
        gres_val=$(echo "$gres_val" | sed 's/(null)//g; s/NONE//g; s/^[[:space:]]*//; s/[[:space:]]*$//')
        if [ -n "$gres_val" ]; then
            PART_GRES["$pname"]="$gres_val"
        fi
    fi
    if [[ "$line" =~ AllowedQos=(.+)$ ]]; then
        qos_val="${BASH_REMATCH[1]}"
        qos_val=$(echo "$qos_val" | sed 's/(null)//g; s/NONE//g; s/^[[:space:]]*//; s/[[:space:]]*$//')
        if [ -n "$qos_val" ]; then
            PART_ALLOWED_QOS["$pname"]="$qos_val"
        fi
    fi
    # SLURM 25.x may use AllowQos (without 'ed') or split across lines
    if [[ "$line" =~ AllowQos=(.+)$ ]]; then
        qos_val="${BASH_REMATCH[1]}"
        qos_val=$(echo "$qos_val" | sed 's/(null)//g; s/NONE//g; s/^[[:space:]]*//; s/[[:space:]]*$//')
        if [ -n "$qos_val" ]; then
            PART_ALLOWED_QOS["$pname"]="$qos_val"
        fi
    fi
done < <(run_cmd "scontrol show partition")

if [ ${#PART_NAMES[@]} -eq 0 ]; then
    die "No partitions found."
fi

echo "[*] Partitions found: ${PART_NAMES[*]}"

# ---- Validate forced partition (if given) --------------------------------
if [ -n "$FORCE_PARTITION" ]; then
    found=0
    for p in "${PART_NAMES[@]}"; do
        if [ "$p" = "$FORCE_PARTITION" ]; then
            found=1
            break
        fi
    done
    if [ "$found" -eq 0 ]; then
        echo "[ERROR] Partition '$FORCE_PARTITION' does not exist."
        echo "         Available partitions: ${PART_NAMES[*]}"
        exit 1
    fi
    echo "[*] Forced partition: $FORCE_PARTITION"
fi

# ---- GPU partition detection ---------------------------------------------
GPU_PARTITION=""
GPU_TYPE=""

# Pass 1: scontrol show partition GresTypes
for p in "${PART_NAMES[@]}"; do
    gres="${PART_GRES[$p]:-}"
    state="${PART_STATE[$p]:-UNKNOWN}"
    if [[ "$state" == "DOWN" ]] || [[ "$state" == "DRAIN" ]]; then continue; fi
    if [ -n "$gres" ]; then
        # Try each comma-separated GRES entry
        IFS=',' read -ra entries <<< "$gres"
        for entry in "${entries[@]}"; do
            entry=$(echo "$entry" | tr -d '[:space:]')
            model=$(extract_gpu_model "$entry") || continue
            GPU_PARTITION="$p"
            GPU_TYPE="$model"
            break 2
        done
    fi
done

# Pass 2: sinfo --noheader -o '%P|%G' — parse per-partition GRES
if [ -z "$GPU_PARTITION" ]; then
    while IFS='|' read -r pname gres; do
        pname=$(echo "$pname" | tr -d '[:space:]')
        gres=$(echo "$gres" | tr -d '[:space:]')
        if [ -z "$gres" ] || [[ "$gres" == "(null)" ]] || [[ "$gres" == "NONE" ]]; then
            continue
        fi
        # GRES can be comma-separated; try each entry
        IFS=',' read -ra entries <<< "$gres"
        for entry in "${entries[@]}"; do
            entry=$(echo "$entry" | tr -d '[:space:]')
            model=$(extract_gpu_model "$entry") || continue
            GPU_PARTITION="$pname"
            GPU_TYPE="$model"
            break 2
        done
    done < <(run_cmd "sinfo --noheader -o '%P|%G' 2>/dev/null") || true
fi

# Pass 3: scontrol show node — get GRES from individual nodes
if [ -z "$GPU_PARTITION" ]; then
    NODE_GRES_RAW=$(run_cmd "scontrol show node" | grep -oP 'GRES=\K\S+' | sort -u) || true
    if [ -n "$NODE_GRES_RAW" ]; then
        while IFS= read -r gres; do
            gres=$(echo "$gres" | tr -d '[:space:]')
            [ -z "$gres" ] && continue
            [[ "$gres" == "(null)" ]] && continue
            [[ "$gres" == "NONE" ]] && continue

            # Try to extract model from this GRES string
            model=$(extract_gpu_model "$gres") || continue

            GPU_TYPE="$model"
            echo "[*] Found GPU model on nodes: ${GPU_TYPE}"

            # Find which partition has nodes with this GPU
            for p in "${PART_NAMES[@]}"; do
                state="${PART_STATE[$p]:-UNKNOWN}"
                if [[ "$state" == "DOWN" ]] || [[ "$state" == "DRAIN" ]]; then continue; fi
                match=$(run_cmd "scontrol show node -p $p 2>/dev/null" | grep -i "${model}" | head -1) || true
                if [ -n "$match" ]; then
                    GPU_PARTITION="$p"
                    break
                fi
            done

            # If no specific partition found, pick the most specific gpu-* partition
            if [ -z "$GPU_PARTITION" ]; then
                for p in "${PART_NAMES[@]}"; do
                    lower="${p,,}"
                    # Prefer partitions that contain the model name (e.g. gpu-a100)
                    if [[ "$lower" == *"${model,,}"* ]]; then
                        GPU_PARTITION="$p"
                        break
                    fi
                done
            fi

            # Still nothing? Pick any gpu-* partition
            if [ -z "$GPU_PARTITION" ]; then
                for p in "${PART_NAMES[@]}"; do
                    if [[ "${p,,}" == *gpu* ]]; then
                        GPU_PARTITION="$p"
                        break
                    fi
                done
            fi

            break
        done <<< "$NODE_GRES_RAW"
    fi
fi

# Pass 4: Last resort — pick any partition with "gpu" in its name
if [ -z "$GPU_PARTITION" ]; then
    for p in "${PART_NAMES[@]}"; do
        if [[ "${p,,}" == *gpu* ]]; then
            GPU_PARTITION="$p"
            break
        fi
    done
fi

if [ -n "$GPU_PARTITION" ]; then
    echo "[*] GPU partition: ${GPU_PARTITION} (GPU type: ${GPU_TYPE:-auto-detect failed, will use generic})"
else
    echo "[*] No GPU partition detected"
fi

# ---- CPU-only partition detection ----------------------------------------
CPU_PARTITION=""
for p in "${PART_NAMES[@]}"; do
    state="${PART_STATE[$p]:-UNKNOWN}"
    if [[ "$state" == "DOWN" ]] || [[ "$state" == "DRAIN" ]]; then continue; fi
    gres="${PART_GRES[$p]:-}"
    if [ -z "$gres" ]; then
        CPU_PARTITION="$p"
        break
    fi
done
# If no GPU-less partition found, fall back to the first UP partition
if [ -z "$CPU_PARTITION" ]; then
    for p in "${PART_NAMES[@]}"; do
        state="${PART_STATE[$p]:-UNKNOWN}"
        if [[ "$state" != "DOWN" ]] && [[ "$state" != "DRAIN" ]]; then
            CPU_PARTITION="$p"
            break
        fi
    done
fi
if [ -n "$CPU_PARTITION" ]; then
    echo "[*] CPU-only partition: ${CPU_PARTITION}"
else
    echo "[*] No CPU-only partition detected"
fi

# ---- Override with forced partition (if given) ---------------------------
if [ -n "$FORCE_PARTITION" ]; then
    GPU_PARTITION="$FORCE_PARTITION"
    # Try to detect GPU type for the forced partition
    gpu_gres="${PART_GRES[$FORCE_PARTITION]:-}"
    if [ -n "$gpu_gres" ]; then
        IFS=',' read -ra entries <<< "$gpu_gres"
        for entry in "${entries[@]}"; do
            entry=$(echo "$entry" | tr -d '[:space:]')
            model=$(extract_gpu_model "$entry") || continue
            GPU_TYPE="$model"
            break
        done
    fi
    echo "[*] Overriding GPU partition to: $FORCE_PARTITION (GPU type: ${GPU_TYPE:-auto-detect failed, will use generic})"
fi

# ---- Default partition ---------------------------------------------------
CURRENT_USER=$(whoami 2>/dev/null || echo "${USER:-unknown}")
DEFAULT_PARTITION=$(run_cmd "scontrol show user ${CURRENT_USER}" | grep -oP 'DefaultPartition=\K\S+' | head -1) || true
if [ -z "$DEFAULT_PARTITION" ]; then
    DEFAULT_PARTITION="${PART_NAMES[0]}"
fi

# Override default partition if forced
if [ -n "$FORCE_PARTITION" ]; then
    DEFAULT_PARTITION="$FORCE_PARTITION"
fi

echo "[*] Default partition: ${DEFAULT_PARTITION}"

# ---- Accounts with balance check -----------------------------------------
declare -a ACCOUNTS=()
VALID_ACCOUNT=""

# Method 1: sacctmgr (pipe-delimited output)
# Output looks like: Account|State|
#                    my-account|ACTIVE|
while IFS='|' read -r acct state rest; do
    acct=$(echo "$acct" | tr -d '[:space:]')
    # Skip header line and empty entries
    if [ -z "$acct" ] || [[ "$acct" == "Account" ]]; then
        continue
    fi
    ACCOUNTS+=("$acct")
done < <(run_cmd "sacctmgr -p show assoc user=${CURRENT_USER} format=account%20,state%10")

# Method 2: scontrol show user fallback
if [ ${#ACCOUNTS[@]} -eq 0 ]; then
    while IFS= read -r acct; do
        acct=$(echo "$acct" | tr -d '[:space:]')
        if [ -n "$acct" ]; then
            ACCOUNTS+=("$acct")
        fi
    done < <(run_cmd "scontrol show user ${CURRENT_USER}" | grep -oP 'Account=\K\S+') || true
fi

# Deduplicate
if [ ${#ACCOUNTS[@]} -gt 0 ]; then
    readarray -t ACCOUNTS < <(printf '%s\n' "${ACCOUNTS[@]}" | sort -u)
fi

echo "[*] Accounts: ${ACCOUNTS[*]:-none}"

# ---- Check account balances using sbalance -------------------------------
if [ ${#ACCOUNTS[@]} -gt 0 ]; then
    BALANCE_OUTPUT=$(run_cmd "sbalance --all") || true

    if [ -n "$BALANCE_OUTPUT" ]; then
        for acct in "${ACCOUNTS[@]}"; do
            BALANCE_LINE=$(echo "$BALANCE_OUTPUT" | grep -i "$acct" | head -1) || true
            if [ -n "$BALANCE_LINE" ]; then
                BALANCE_VAL=$(echo "$BALANCE_LINE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+(\.[0-9]+)?$/) {print $i; exit}}') || true
                if [ -n "$BALANCE_VAL" ]; then
                    HAS_BALANCE=$(awk "BEGIN {print ($BALANCE_VAL > 0) ? 1 : 0}")
                    if [ "$HAS_BALANCE" -eq 1 ]; then
                        VALID_ACCOUNT="$acct"
                        echo "[*] Account with balance: ${VALID_ACCOUNT} (balance: ${BALANCE_VAL})"
                        break
                    fi
                fi
            fi
        done
    fi

    # Fallback: use first account if sbalance didn't work
    if [ -z "$VALID_ACCOUNT" ] && [ ${#ACCOUNTS[@]} -gt 0 ]; then
        VALID_ACCOUNT="${ACCOUNTS[0]}"
        echo "[*] sbalance unavailable; using first account: ${VALID_ACCOUNT}"
    fi
fi

echo "[*] Selected account: ${VALID_ACCOUNT:-none}"

# ---- QoS -----------------------------------------------------------------
declare -a QOS_LIST=()

# Try multiple commands for compatibility across SLURM versions:
#   SLURM <23:  scontrol show qos          (Qos=name)
#   SLURM >=23: sacctmgr show qos format=Name
qos_output=""
if command -v sacctmgr &>/dev/null; then
    qos_output=$(run_cmd "sacctmgr show qos format=Name") || true
fi
if [ -z "$qos_output" ] && command -v sacctmgr &>/dev/null; then
    qos_output=$(run_cmd "sacctmgr -p list qos format=Name") || true
fi
if [ -z "$qos_output" ]; then
    qos_output=$(run_cmd "scontrol show qos") || true
fi

while IFS= read -r q; do
    q=$(echo "$q" | tr -s ' ' | xargs)  # trim whitespace
    [ -z "$q" ] && continue
    # Skip headers, empty lines, and non-name fields
    [[ "$q" =~ ^(Name|\#|$) ]] && continue
    [[ "$q" =~ ^[0-9] ]] && continue
    # Try Qos= or Name= key=value format
    kv=$(echo "$q" | grep -oP '(?:Qos=|Name=)\K\S+') || true
    if [ -n "$kv" ]; then
        QOS_LIST+=("$kv")
        continue
    fi
    # Bare name (single word, looks like a QoS)
    if [[ "$q" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]]; then
        QOS_LIST+=("$q")
    fi
done <<<"$qos_output"

# Deduplicate
if [ ${#QOS_LIST[@]} -gt 0 ]; then
    readarray -t QOS_LIST < <(printf '%s\n' "${QOS_LIST[@]}" | sort -u)
fi

echo "[*] QoS: ${QOS_LIST[*]:-none}"

# ---- GRES (GPUs, etc.) ---------------------------------------------------
declare -a GRES_TYPES=()

while IFS= read -r g; do
    g=$(echo "$g" | tr -s ' ' | xargs)
    if [ -n "$g" ] && [[ "$g" != "(null)" ]] && [[ "$g" != "NONE" ]]; then
        GRES_TYPES+=("$g")
    fi
done < <(run_cmd "sinfo --noheader -o %G 2>/dev/null") || true

if [ ${#GRES_TYPES[@]} -eq 0 ]; then
    while IFS= read -r g; do
        for item in $(echo "$g" | tr ',' ' '); do
            if [ -n "$item" ] && [[ "$item" != "(null)" ]]; then
                GRES_TYPES+=("$item")
            fi
        done
    done < <(run_cmd "scontrol show partition" | grep -oP 'GresTypes=\K\S+') || true
fi

if [ ${#GRES_TYPES[@]} -gt 0 ]; then
    readarray -t GRES_TYPES < <(printf '%s\n' "${GRES_TYPES[@]}" | sort -u)
fi

echo "[*] GRES types: ${GRES_TYPES[*]:-none}"

# ---- Node features (constraints) -----------------------------------------
declare -a FEATURES=()
while IFS= read -r f; do
    for item in $(echo "$f" | tr ',' ' '); do
        if [ -n "$item" ] && [[ "$item" != "(null)" ]]; then
            FEATURES+=("$item")
        fi
    done
done < <(run_cmd "scontrol show node" | grep -oP 'Features=\K\S+') || true

if [ ${#FEATURES[@]} -gt 0 ]; then
    readarray -t FEATURES < <(printf '%s\n' "${FEATURES[@]}" | sort -u)
fi

echo "[*] Node features: ${FEATURES[*]:-none}"

# ---- Write cluster.conf --------------------------------------------------
{
    echo "# ============================================================"
    echo "#  cluster.conf — Auto-generated by query_cluster.sh"
    echo "#  Generated: $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
    echo "#  Hostname:  $(hostname 2>/dev/null || echo unknown)"
    echo "#  SLURM:     ${SLURM_VER}"
    echo "# ============================================================"
    echo ""
    echo "[general]"
    echo "default_partition = ${DEFAULT_PARTITION}"
    if [ -n "$VALID_ACCOUNT" ]; then
        echo "default_account = ${VALID_ACCOUNT}"
    else
        echo "# default_account = <not detected or no balance>"
    fi
    if [ -n "$GPU_PARTITION" ]; then
        echo "gpu_partition = ${GPU_PARTITION}"
    else
        echo "# gpu_partition = <none detected>"
    fi
    if [ -n "$GPU_TYPE" ]; then
        echo "gpu_type = ${GPU_TYPE}"
    else
        echo "# gpu_type = <none detected>"
    fi
    if [ -n "$CPU_PARTITION" ]; then
        echo "cpu_partition = ${CPU_PARTITION}"
    else
        echo "# cpu_partition = <none detected>"
    fi
    echo ""

    echo "[partitions]"
    for p in "${PART_NAMES[@]}"; do
        echo "  ${p}:"
        echo "    nodes = ${PART_NODES[$p]:-unknown}"
        echo "    state = ${PART_STATE[$p]:-unknown}"
        echo "    time_limit = ${PART_TIME[$p]:-UNLIMITED}"
        gres="${PART_GRES[$p]:-}"
        if [ -n "$gres" ]; then
            echo "    gres = ${gres}"
            echo "    has_gpu = true"
        else
            echo "    has_gpu = false"
        fi
        # Allowed QoS for this partition
        allowed_qos="${PART_ALLOWED_QOS[$p]:-}"
        if [ -n "$allowed_qos" ]; then
            echo "    allowed_qos = ${allowed_qos}"
        else
            echo "    allowed_qos = (all)"
        fi
        echo ""
    done

    echo "[accounts]"
    if [ ${#ACCOUNTS[@]} -gt 0 ]; then
        for a in "${ACCOUNTS[@]}"; do
            marker=""
            if [ "$a" = "$VALID_ACCOUNT" ]; then
                marker=" (selected — has balance)"
            fi
            echo "  - ${a}${marker}"
        done
    else
        echo "  # No accounts detected"
    fi
    echo ""

    echo "[qos]"
    if [ ${#QOS_LIST[@]} -gt 0 ]; then
        for q in "${QOS_LIST[@]}"; do
            echo "  - ${q}"
        done
    else
        echo "  # No QoS detected"
    fi
    echo ""

    echo "[gres]"
    if [ ${#GRES_TYPES[@]} -gt 0 ]; then
        for g in "${GRES_TYPES[@]}"; do
            echo "  - ${g}"
        done
    else
        echo "  # No GRES resources detected"
    fi
    echo ""

    echo "[node_features]"
    if [ ${#FEATURES[@]} -gt 0 ]; then
        for f in "${FEATURES[@]}"; do
            echo "  - ${f}"
        done
    else
        echo "  # No features detected"
    fi
    echo ""

} > "$OUTPUT"

echo ""
echo "[+] Written to ${OUTPUT}"
