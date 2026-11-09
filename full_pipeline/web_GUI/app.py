#!/usr/bin/env python3
"""
ANTIOCH-VS Web GUI - Main Application

A web-based interface for preparing docking workflows.
Designed to be easily extensible with new modules.

Steps:
  1. Input: protein (PDB ID or upload) + chain selection + ligand file
  2. Visualize: 3D viewer to pick interface residues
  3. Output: downloadable dockflow_prepped.zip

Usage:
    source /scratch/software/micromamba_root/sourceme.sh
    micromamba activate chem_env
    python app.py

Architecture:
  - Each "step" is a separate route + template for easy extension.
  - Session state carries cleaned PDB data and ligand data between steps.
"""

import os
import sys
import io
import zipfile
import tempfile
from pathlib import Path

from flask import (
    Flask, render_template, request, send_file, flash, redirect, url_for, session
)
from werkzeug.utils import secure_filename

# ─── Add workspace to path for clean_pdb module ──────────────────────────
WORKSPACE = Path(os.environ.get("OPENCLAW_WORKSPACE", "~/.openclaw/workspace")).expanduser()
sys.path.insert(0, str(WORKSPACE / "scripts"))

from clean_pdb_for_tleap import clean_pdb, download_pdb, is_pdb_id

# ─── Flask App Setup ─────────────────────────────────────────────────────
app = Flask(__name__)
app.secret_key = "dockflow-secret-key-change-in-production"
app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024  # 50 MB max upload

# Allowed file extensions
ALLOWED_PROTEIN_EXTS = {".pdb", ".ent"}
ALLOWED_LIGAND_EXTS = {".smi", ".umf"}


def allowed_file(filename: str, allowed_exts: set) -> bool:
    """Check if file extension is allowed."""
    return Path(filename).suffix.lower() in allowed_exts


# ─── Step 1: Input Form ─────────────────────────────────────────────────

@app.route("/")
def index():
    """Main input form — protein + ligand."""
    return render_template("index.html")


@app.route("/prepare", methods=["POST"])
def prepare():
    """Process protein input, clean it, store in session → go to viewer."""

    protein_source = request.form.get("protein_source", "pdb_id").strip()
    chain_selection = request.form.get("chains", "").strip()
    chain_filter: set[str] | None = set(chain_selection.upper().split()) if chain_selection else None

    cleaned_pdb_data: bytes | None = None
    cleaned_pdb_filename: str = "protein_cleaned.pdb"

    if protein_source == "pdb_id":
        pdb_id = request.form.get("pdb_id", "").strip().upper()
        if not pdb_id or not is_pdb_id(pdb_id):
            flash("Invalid PDB ID. Please enter a 4-character alphanumeric ID.", "error")
            return redirect(url_for("index"))

        with tempfile.TemporaryDirectory() as tmpdir:
            raw_path = os.path.join(tmpdir, f"{pdb_id.lower()}.pdb")
            cleaned_path = os.path.join(tmpdir, f"cleaned_{pdb_id.lower()}.pdb")

            try:
                download_pdb(pdb_id, raw_path)
            except Exception as e:
                flash(f"Failed to download PDB {pdb_id}: {e}", "error")
                return redirect(url_for("index"))

            clean_pdb(raw_path, cleaned_path, chain_filter=chain_filter)
            with open(cleaned_path, "rb") as f:
                cleaned_pdb_data = f.read()
            cleaned_pdb_filename = f"cleaned_{pdb_id.lower()}.pdb"

    elif protein_source == "upload":
        if "protein_file" not in request.files:
            flash("No protein file uploaded.", "error")
            return redirect(url_for("index"))

        protein_file = request.files["protein_file"]
        if protein_file.filename == "" or not allowed_file(protein_file.filename, ALLOWED_PROTEIN_EXTS):
            flash("Invalid protein file. Please upload a .pdb or .ent file.", "error")
            return redirect(url_for("index"))

        with tempfile.TemporaryDirectory() as tmpdir:
            raw_path = os.path.join(tmpdir, secure_filename(protein_file.filename))
            protein_file.save(raw_path)
            cleaned_path = os.path.join(tmpdir, "protein_cleaned.pdb")

            clean_pdb(raw_path, cleaned_path, chain_filter=chain_filter)
            with open(cleaned_path, "rb") as f:
                cleaned_pdb_data = f.read()
    else:
        flash("Invalid protein source selected.", "error")
        return redirect(url_for("index"))

    if cleaned_pdb_data is None:
        flash("Protein cleaning failed.", "error")
        return redirect(url_for("index"))

    # ── Ligand Input ────────────────────────────────────────────────
    ligand_filename: str | None = None
    ligand_data: bytes | None = None

    if "ligand_file" in request.files:
        ligand_file = request.files["ligand_file"]
        if ligand_file.filename != "":
            if not allowed_file(ligand_file.filename, ALLOWED_LIGAND_EXTS):
                flash("Invalid ligand file. Please upload a .smi or .umf file.", "error")
                return redirect(url_for("index"))
            ligand_filename = secure_filename(ligand_file.filename)
            ligand_data = ligand_file.read()

    if ligand_filename is None:
        flash("No ligand file provided. Please upload a .smi or .umf file.", "error")
        return redirect(url_for("index"))

    # ── Store in session for next steps ─────────────────────────────
    # Flask sessions have size limits (~4KB signed cookie). For large PDBs
    # we store data in a server-side temp file and track the path in session.
    import uuid
    job_id = str(uuid.uuid4())[:8]

    # Create a persistent job directory
    job_dir = Path(app.instance_path) / "jobs" / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    pdb_out = job_dir / cleaned_pdb_filename
    pdb_out.write_bytes(cleaned_pdb_data)

    ligand_out = job_dir / ligand_filename
    ligand_out.write_bytes(ligand_data)

    session["dockflow_job"] = {
        "job_id": job_id,
        "pdb_filename": str(cleaned_pdb_filename),
        "ligand_filename": ligand_filename,
    }

    return redirect(url_for("visualize"))


# ─── Step 2: 3D Visualization & Residue Selection ───────────────────────

@app.route("/api_test")
def api_test():
    """3Dmol.js API exploration page."""
    return render_template("api_test.html")


@app.route("/test")
def test_viewer():
    """Minimal 3Dmol.js test page."""
    return render_template("test_viewer.html")


@app.route("/visualize")
def visualize():
    """3D protein viewer — user clicks residues to select interface."""
    job = session.get("dockflow_job")
    if not job:
        flash("No active job. Please start from the input form.", "error")
        return redirect(url_for("index"))

    pdb_path = Path(app.instance_path) / "jobs" / job["job_id"] / job["pdb_filename"]
    if not pdb_path.exists():
        flash("Job data expired. Please restart.", "error")
        session.pop("dockflow_job", None)
        return redirect(url_for("index"))

    return render_template("visualize.html", job=job)


import logging
logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO)


@app.route("/api/pdb")
def api_pdb():
    """Return the cleaned PDB file for the current job."""
    from flask import send_file as _send_file
    job = session.get("dockflow_job")
    logger.info(f"[API] /api/pdb requested, job in session: {job is not None}")
    if not job:
        logger.warning("[API] No active job in session!")
        return "No active job", 400
    pdb_path = Path(app.instance_path) / "jobs" / job["job_id"] / job["pdb_filename"]
    logger.info(f"[API] Looking for PDB at: {pdb_path}")
    if not pdb_path.exists():
        logger.warning(f"[API] PDB file not found: {pdb_path}")
        return "PDB file not found", 404
    size = pdb_path.stat().st_size
    logger.info(f"[API] Serving PDB ({size} bytes)")
    return _send_file(str(pdb_path), mimetype="chemical/x-pdb")


@app.route("/finalize", methods=["POST"])
def finalize():
    """Collect selected residues, build ZIP, offer download."""
    job = session.get("dockflow_job")
    if not job:
        flash("No active job. Please start from the input form.", "error")
        return redirect(url_for("index"))

    job_dir = Path(app.instance_path) / "jobs" / job["job_id"]

    # Get selected residues from form (JSON array of {chain, resseq, resname})
    import json
    raw_selections = request.form.get("selected_residues", "[]")
    try:
        selections = json.loads(raw_selections)
    except json.JSONDecodeError:
        selections = []

    # Build interface_residues.csv
    csv_lines = ["Chain,ResidueID,ResidueName"]
    for res in selections:
        chain = str(res.get("chain", "")).strip().upper()
        resseq = str(res.get("resseq", "")).strip()
        resname = str(res.get("resname", "")).strip().upper()
        if chain and resseq and resname:
            csv_lines.append(f"{chain},{resseq},{resname}")
    csv_content = "\n".join(csv_lines) + "\n"

    # Build ZIP
    pdb_path = job_dir / job["pdb_filename"]
    ligand_path = job_dir / job["ligand_filename"]

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(pdb_path, arcname=job["pdb_filename"])
        if ligand_path.exists():
            zf.write(ligand_path, arcname=job["ligand_filename"])
        zf.writestr("interface_residues.csv", csv_content)

    zip_buffer.seek(0)

    # Clean up job directory
    import shutil
    shutil.rmtree(job_dir, ignore_errors=True)
    session.pop("dockflow_job", None)

    return send_file(
        zip_buffer,
        mimetype="application/zip",
        as_attachment=True,
        download_name="dockflow_prepped.zip",
    )


# ─── Main ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Ensure instance/jobs directory exists
    Path(app.instance_path).mkdir(parents=True, exist_ok=True)
    (Path(app.instance_path) / "jobs").mkdir(exist_ok=True)

    print("Starting ANTIOCH-VS Web GUI...")
    print("Open http://localhost:5000 in your browser.")
    app.run(host="0.0.0.0", port=5000, debug=True)
