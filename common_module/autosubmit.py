#!/usr/bin/env python3
"""
autosubmit.py — Submit any bash script to SLURM using cluster.conf + template.

Usage:
    # Minimal — uses all defaults from cluster.conf
    python3 autosubmit.py my_script.sh

    # Override resources
    python3 autosubmit.py my_script.sh --ntasks 4 --cpus-per-task 2 --mem 16G --time 04:00:00

    # Manually set number of CPUs/threads (overrides default of 4)
    python3 autosubmit.py my_script.sh --num-cpus 32

    # Override partition / account
    python3 autosubmit.py my_script.sh --partition gpu --account my_account

    # No GPU — use a CPU-only partition and skip the --gpus line
    python3 autosubmit.py my_script.sh --nogpu

    # Add extra #SBATCH lines (e.g. --gres, --constraint, --qos)
    python3 autosubmit.py my_script.sh --extra "#SBATCH --gres=gpu:1" "#SBATCH --constraint=v100"

    # Load modules / activate conda
    python3 autosubmit.py my_script.sh --modules "CUDA/12.1 Python/3.11" --conda-env myenv

    # Set working directory and output directory
    python3 autosubmit.py my_script.sh --work-dir /scratch/runs/exp1 --output-dir /scratch/logs

    # Dry-run — print the generated .slurm file without submitting
    python3 autosubmit.py my_script.sh --dry-run

    # Custom job name
    python3 autosubmit.py my_script.sh --job-name "my_experiment"

    # Use a manual SLURM header file (bypasses cluster.conf entirely)
    python3 autosubmit.py my_script.sh --manual-header my_header.slurm

Paths:
    cluster.conf       : auto-detected next to this script, or via --config
    slurm_template.slurm: auto-detected next to this script, or via --template
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
#  Defaults (used when cluster.conf is missing a key)
# ---------------------------------------------------------------------------
DEFAULTS = {
    "JOB_NAME": "dockflow_job",
    "NTASKS_PER_NODE": "4",
    "MEM": "8000",
    "TIME": "24:00:00",
    "ACCOUNT": "",
    "MODULES": "",
    "CONDA_ENV": "",
    "EXTRA_SBATCH": "",
}


# ---------------------------------------------------------------------------
#  Config loader
# ---------------------------------------------------------------------------
def run_cmd(cmd):
    """Run a shell command and return stripped stdout."""
    r = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, timeout=15
    )
    return r.stdout.strip()


def load_cluster_conf(path):
    """
    Parse the simple INI-like cluster.conf into a dict.

    The file is *not* strict INI (it uses 'key = value' under sections with
    indented sub-keys), so we do lightweight manual parsing.
    """
    cfg = {}
    p = Path(path)
    if not p.exists():
        print(f"[!] cluster.conf not found at {p}", file=sys.stderr)
        return cfg

    current_section = None
    for line in p.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        # Section header
        m = re.match(r"^\[(\w+)\]$", stripped)
        if m:
            current_section = m.group(1)
            continue

        # Key = value (top-level in [general])
        m = re.match(r"^(\w+)\s*=\s*(.+)$", stripped)
        if m and current_section == "general":
            cfg[m.group(1)] = m.group(2).strip()
            continue

        # Sub-key under a partition: "  key = value"
        m = re.match(r"^\s+(\w+)\s*=\s*(.+)$", stripped)
        if m and current_section == "partitions":
            cfg[current_section + "_" + m.group(1)] = m.group(2).strip()

    # Also parse top-level [qos] list for fallback
    current_qos_list = []
    in_qos = False
    for line in p.read_text().splitlines():
        stripped = line.strip()
        if re.match(r"^\[qos\]$", stripped):
            in_qos = True
            continue
        if stripped.startswith("[") and in_qos:
            in_qos = False
            continue
        if in_qos:
            m = re.match(r"^-\s+(.+)$", stripped)
            if m:
                current_qos_list.append(m.group(1).strip())
    cfg["qos_list"] = current_qos_list

    return cfg


# ---------------------------------------------------------------------------
#  Template substitution
# ---------------------------------------------------------------------------
def fill_template(template_text, values):
    """Replace <<KEY>> placeholders with values."""
    result = template_text
    for key, val in values.items():
        result = result.replace(f"<<{key}>>", str(val))
    # If any placeholders remain unfilled, warn but don't fail
    remaining = re.findall(r"<<(\w+)>>", result)
    if remaining:
        print(f"[!] Unfilled placeholders: {', '.join(remaining)}", file=sys.stderr)
    return result


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------
def main():
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="Submit a bash script to SLURM via cluster.conf + template.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("user_script", help="Path to the bash script to submit")
    parser.add_argument("--config", default=str(script_dir / "cluster.conf"), help="Path to cluster.conf (auto-detected by default)")
    parser.add_argument("--template", default=str(script_dir / "slurm_template.slurm"), help="Path to SLURM template (auto-detected by default)")
    parser.add_argument("--job-name", default=None, help="Job name (overrides default)")
    parser.add_argument("--partition", default=None, help="Partition (overrides cluster.conf default)")
    parser.add_argument("--account", default=None, help="Account (overrides cluster.conf default)")
    parser.add_argument("--ntasks", type=int, default=None, help="Number of tasks per node")
    parser.add_argument("--cpus-per-task", type=int, default=None, help="CPUs per task (if template supports it)")
    parser.add_argument("--mem", default=None, help="Memory in MB (e.g. 4096, 8000, 16384)")
    parser.add_argument("--time", default=None, help="Time limit (HH:MM:SS or D-HH:MM)")
    parser.add_argument("--modules", default=None, help="Space-separated module names to load")
    parser.add_argument("--conda-env", default=None, help="Conda environment to activate")
    parser.add_argument("--work-dir", default=None, help="Working directory for the job")
    parser.add_argument("--output-dir", default=None, help="Directory for SLURM output/error files")
    parser.add_argument("--extra", action="append", default=[], help="Extra #SBATCH lines (may be repeated)")
    parser.add_argument("--dry-run", action="store_true", help="Print the generated script without submitting")
    parser.add_argument("--no-submit", action="store_true", help="Write the .slurm file but do not sbatch")
    parser.add_argument("--output-file", default=None, help="Write generated .slurm to this path (default: <job_name>.slurm)")
    parser.add_argument("--nogpu", action="store_true", help="Do not request a GPU; use a CPU-only partition if available")
    parser.add_argument("--num-cpus", type=int, default=None, help="Number of CPUs/threads (overrides default of 4)")
    parser.add_argument("--manual-header", default=None, help="Path to a manual SLURM header file; bypasses cluster.conf and template entirely")

    args = parser.parse_args()

    # ---- Resolve user script early (needed by manual-header mode too) ----
    user_script = Path(args.user_script).resolve()
    if not user_script.exists():
        print(f"[!] User script not found: {user_script}", file=sys.stderr)
        sys.exit(1)

    # ---- Manual header mode: bypass cluster.conf + template -------------
    if args.manual_header:
        header_path = Path(args.manual_header)
        if not header_path.exists():
            print(f"[!] Manual header file not found: {header_path}", file=sys.stderr)
            sys.exit(1)

        header_text = header_path.read_text()
        work_dir = args.work_dir or str(user_script.parent)

        # Build the combined script: header + user script execution
        combined = header_text.rstrip() + "\n\n"
        combined += f"# ---- Run the user script -------------------------------------------------\n"
        combined += f'echo "=== Job started at $(date -u \'+%Y-%m-%d %H:%M:%S UTC\') ==="\n'
        combined += f'echo "=== SLURM_JOB_ID  = $SLURM_JOB_ID ==="\n'
        combined += f'echo "=== Working dir   = {work_dir} ==="\n\n'
        combined += f'cd "{work_dir}"\n'
        combined += f'bash "{user_script}"\n'
        combined += f"EXIT_CODE=$?\n\n"
        combined += f'echo "=== Job finished at $(date -u \'+%Y-%m-%d %H:%M:%S UTC\') ==="\n'
        combined += f'echo "=== Exit code     = $EXIT_CODE ==="\n\n'
        combined += f"exit $EXIT_CODE\n"

        if args.dry_run:
            print("=" * 60)
            print("DRY RUN — generated SLURM script (manual header mode):")
            print("=" * 60)
            print(combined)
            print("=" * 60)
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        job_name = args.job_name or DEFAULTS["JOB_NAME"]
        safe_name = re.sub(r"[^\w\-]", "_", job_name)
        default_out = f"{safe_name}_{timestamp}.slurm"
        out_path = Path(args.output_file) if args.output_file else Path(default_out)

        out_path.write_text(combined)
        print(f"[+] Generated: {out_path.resolve()}")

        if args.no_submit:
            print("[i] --no-submit set; not calling sbatch.")
            return

        # ---- Submit ----------------------------------------------------------
        has_slurm = shutil.which("sbatch") is not None
        if not has_slurm:
            print("[!] SLURM (sbatch) not found — cannot submit.", file=sys.stderr)
            sys.exit(1)

        try:
            result = subprocess.run(
                ["sbatch", str(out_path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                timeout=30,
            )
            if result.returncode == 0:
                m = re.search(r"Submitted batch job (\d+)", result.stdout)
                job_id = m.group(1) if m else "unknown"
                print(f"[+] Submitted! Job ID: {job_id}")
            else:
                print(f"[!] sbatch failed (rc={result.returncode})", file=sys.stderr)
                print(result.stderr, file=sys.stderr)
                sys.exit(result.returncode)
        except subprocess.TimeoutExpired:
            print("[!] sbatch timed out after 30s", file=sys.stderr)
            sys.exit(1)

        return

    # ---- Detect whether SLURM is available --------------------------------
    has_slurm = shutil.which("sbatch") is not None

    if not has_slurm:
        print("[i] SLURM (sbatch) not found — running script directly.")
        work_dir = args.work_dir or str(user_script.parent)
        if args.dry_run:
            print(f"  Would run: bash {user_script}")
            return
        os.chdir(work_dir)
        ret = os.system(f"bash {user_script}")
        if ret != 0:
            print(f"[!] Script exited with code {ret >> 8}", file=sys.stderr)
            sys.exit(ret >> 8)
        print("[+] Script finished successfully.")
        return

    # ---- Ensure cluster.conf exists (auto-generate if missing) -----------
    config_path = Path(args.config)
    query_script = script_dir / "query_cluster.sh"
    if not config_path.exists():
        if query_script.exists():
            print(f"[i] cluster.conf not found at {config_path} — running query_cluster.sh to generate it...")
            try:
                result = subprocess.run(
                    ["bash", str(query_script), str(config_path)],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    timeout=60,
                )
                if result.returncode != 0:
                    print(f"[!] query_cluster.sh failed (rc={result.returncode})", file=sys.stderr)
                    if result.stderr:
                        print(result.stderr, file=sys.stderr)
                    sys.exit(1)
                print(f"[+] Generated {config_path} from cluster probe.")
            except FileNotFoundError:
                print("[!] 'bash' not found — cannot run query_cluster.sh", file=sys.stderr)
                sys.exit(1)
            except subprocess.TimeoutExpired:
                print("[!] query_cluster.sh timed out after 60s", file=sys.stderr)
                sys.exit(1)
        else:
            print(f"[!] cluster.conf not found at {config_path}", file=sys.stderr)
            print(f"    Run 'bash {query_script}' first, or provide --config <path>", file=sys.stderr)
            sys.exit(1)

    # ---- Load cluster.conf -----------------------------------------------
    cfg = load_cluster_conf(args.config)

    # ---- Resolve values (CLI > cluster.conf > built-in defaults) ----------
    job_name = args.job_name or DEFAULTS["JOB_NAME"]

    # Partition: prefer GPU partition if detected, otherwise use default
    gpu_partition = cfg.get("gpu_partition", "")
    gpu_type = cfg.get("gpu_type", "")
    default_partition = cfg.get("default_partition", "")

    # ---- Partition selection -------------------------------------------
    want_gpu = not args.nogpu

    if args.partition:
        partition = args.partition
    elif want_gpu and gpu_partition:
        partition = gpu_partition
    elif not want_gpu and cfg.get("cpu_partition", ""):
        # --nogpu: prefer a CPU-only partition from cluster.conf
        partition = cfg.get("cpu_partition", "")
    elif not want_gpu:
        # --nogpu but no cpu_partition in config — scan for any has_gpu=false partition
        raw_conf = Path(args.config).read_text()
        in_part_block = False
        found_cpu = ""
        for line in raw_conf.splitlines():
            if re.match(r"^\s{2}\S+:\s*$", line):
                block_name = line.strip().rstrip(":")
                in_part_block = True
                continue
            if in_part_block:
                if line.strip().startswith("["):
                    in_part_block = False
                    continue
                if re.search(r"has_gpu\s*=\s*false", line):
                    found_cpu = block_name
                    break
        partition = found_cpu or default_partition
    elif want_gpu:
        # GPU requested but no GPU partition detected; fall back to default
        partition = default_partition
    else:
        partition = default_partition

    # ---- Resolve QoS for the chosen partition ---------------------------
    # Strategy: live-query SLURM for AllowedQos of the chosen partition.
    # This is more reliable than trusting a stale cluster.conf.
    part_allowed_qos = []
    if has_slurm:
        qos_raw = run_cmd(f"scontrol show partition {partition}")
        # Match both AllowedQos= (older SLURM) and AllowQos= (SLURM 25.x)
        for m in re.finditer(r"(?:AllowedQos|AllowQos)=(\S+)", qos_raw):
            val = m.group(1)
            if val not in ("(null)", ""):
                part_allowed_qos = [q.strip() for q in val.split(",")]
                break
        if part_allowed_qos:
            print(f"[i] Partition '{partition}' allows QoS: {', '.join(part_allowed_qos)}")

    # Fallback: scan cluster.conf for the allowed_qos line inside the chosen partition block
    if not part_allowed_qos:
        raw_conf = Path(args.config).read_text()
        in_part_block = False
        for line in raw_conf.splitlines():
            s = line.strip()
            if re.match(r"^\s{2}\S+:\s*$", line):
                block_name = s.rstrip(":")
                if block_name == partition:
                    in_part_block = True
                else:
                    in_part_block = False
                continue
            if in_part_block:
                if s.startswith("["):
                    in_part_block = False
                    continue
                m = re.match(r"^\s+allowed_qos\s*=\s*(.+)$", line)
                if m:
                    val = m.group(1).strip()
                    if val != "(all)" and val:
                        part_allowed_qos = [q.strip() for q in val.split(",")]

    # Pick best QoS: prefer the first allowed one, fallback to global list
    qos_to_use = ""
    all_qos = cfg.get("qos_list", [])

    if part_allowed_qos:
        for q in part_allowed_qos:
            if q in all_qos or not all_qos:
                qos_to_use = q
                break
        if not qos_to_use and part_allowed_qos:
            qos_to_use = part_allowed_qos[0]
    elif all_qos:
        for preferred in ["normal", "default"]:
            if preferred in all_qos:
                qos_to_use = preferred
                break
        if not qos_to_use:
            qos_to_use = all_qos[0]

    # GPU specification: only include if GPUs are wanted
    if want_gpu and gpu_type and gpu_type.lower() != "gpu":
        gpus_line = f"#SBATCH --gpus={gpu_type}:1"
    elif want_gpu:
        # Fallback: request 1 GPU without specifying type
        gpus_line = "#SBATCH --gpus=1"
    else:
        gpus_line = ""  # --nogpu: no GPU line

    # Account: only include if a valid account was detected
    if args.account:
        account = args.account
    else:
        account = cfg.get("default_account", "").strip()

    # Validate account: reject empty strings or obviously bogus values
    if not account or account.startswith("#") or "|" in account:
        account = ""

    ntasks_per_node = str(args.num_cpus) if args.num_cpus else (str(args.ntasks) if args.ntasks else DEFAULTS["NTASKS_PER_NODE"])
    mem = args.mem or DEFAULTS["MEM"]
    time_limit = args.time or DEFAULTS["TIME"]
    modules = args.modules or DEFAULTS["MODULES"]
    conda_env = args.conda_env or DEFAULTS["CONDA_ENV"]

    work_dir = args.work_dir or str(user_script.parent)
    output_dir = args.output_dir or work_dir

    extra_sbatch = "\n".join(args.extra) if args.extra else ""

    # ---- Read template ---------------------------------------------------
    template_path = Path(args.template)
    if not template_path.exists():
        print(f"[!] Template not found: {template_path}", file=sys.stderr)
        sys.exit(1)

    template_text = template_path.read_text()

    # ---- Account line: include only if account is non-empty ---------------
    if account:
        account_line = f"#SBATCH --account={account}"
    else:
        account_line = ""  # Placeholder removed — no account line

    # ---- QoS line: include only if a valid QoS was resolved ---------------
    if qos_to_use:
        qos_line = f"#SBATCH --qos={qos_to_use}"
    else:
        qos_line = ""

    # ---- Fill placeholders -----------------------------------------------
    values = {
        "JOB_NAME": job_name,
        "OUTPUT_DIR": output_dir,
        "PARTITION": partition,
        "GPUS_LINE": gpus_line,
        "NTASKS_PER_NODE": ntasks_per_node,
        "MEM": mem,
        "TIME": time_limit,
        "ACCOUNT_LINE": account_line,
        "QOS_LINE": qos_line,
        "MODULES": modules,
        "CONDA_ENV": conda_env,
        "WORK_DIR": work_dir,
        "USER_SCRIPT": str(user_script),
        "EXTRA_SBATCH": extra_sbatch,
    }

    filled = fill_template(template_text, values)

    # ---- Write / print ---------------------------------------------------
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_name = re.sub(r"[^\w\-]", "_", job_name)
    default_out = f"{safe_name}_{timestamp}.slurm"

    out_path = Path(args.output_file) if args.output_file else Path(default_out)

    if args.dry_run:
        print("=" * 60)
        print("DRY RUN — generated SLURM script:")
        print("=" * 60)
        print(filled)
        print("=" * 60)
        return

    out_path.write_text(filled)
    print(f"[+] Generated: {out_path.resolve()}")

    if args.no_submit:
        print("[i] --no-submit set; not calling sbatch.")
        return

    # ---- Submit ----------------------------------------------------------
    try:
        result = subprocess.run(
            ["sbatch", str(out_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=30,
        )
        if result.returncode == 0:
            m = re.search(r"Submitted batch job (\d+)", result.stdout)
            job_id = m.group(1) if m else "unknown"
            print(f"[+] Submitted! Job ID: {job_id}")
        else:
            print(f"[!] sbatch failed (rc={result.returncode})", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            sys.exit(result.returncode)

    except subprocess.TimeoutExpired:
        print("[!] sbatch timed out after 30s", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
