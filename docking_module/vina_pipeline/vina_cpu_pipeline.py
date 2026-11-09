import argparse
import json
import logging
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import yaml
from Bio.PDB import PDBParser


class MaxLogLevelFilter(logging.Filter):
    """Pass records at or below one logging level."""

    def __init__(self, maximum_level):
        super().__init__()
        self.maximum_level = maximum_level

    def filter(self, record):
        return record.levelno <= self.maximum_level


def configure_logging(log_file):
    formatter = logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S')
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.addFilter(MaxLogLevelFilter(logging.WARNING))
    stdout_handler.setFormatter(formatter)
    root_logger.addHandler(stdout_handler)

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.ERROR)
    stderr_handler.setFormatter(formatter)
    root_logger.addHandler(stderr_handler)

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)


def write_config_snapshot(log_file, config_data):
    with open(log_file, "a") as handle:
        handle.write("\n" + "=" * 60 + "\n")
        handle.write("--- CONFIG SNAPSHOT ---\n")
        handle.write("=" * 60 + "\n")
        handle.write(yaml.safe_dump(config_data, sort_keys=False))
        handle.write("\n")


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def write_json_summary(summary, json_out):
    json_path = Path(json_out)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with json_path.open("w") as handle:
        json.dump(_json_safe(summary), handle, indent=2, sort_keys=False)
        handle.write("\n")


class VinaCPUPipeline:
    def __init__(self, receptor_pdb=None, receptor_pdbqt=None, ligand_dir=None, output_dir=None, **kwargs):
        # Ensure paths are absolute
        self.receptor_pdb = os.path.abspath(receptor_pdb) if receptor_pdb else None
        self.ligand_dir = os.path.abspath(ligand_dir) if ligand_dir else None
        self.output_dir = os.path.abspath(output_dir) if output_dir else None

        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)

        # Handle explicit PDBQT bypass logic
        if receptor_pdbqt:
            self.receptor_pdbqt = os.path.abspath(receptor_pdbqt)
            self.skip_prep = True
            self.prefix = os.path.splitext(os.path.basename(self.receptor_pdbqt))[0]
        elif self.receptor_pdb:
            self.prefix = os.path.splitext(os.path.basename(self.receptor_pdb))[0]
            self.receptor_pdbqt = os.path.join(self.output_dir, f"{self.prefix}.pdbqt")
            self.skip_prep = False
        else:
            self.receptor_pdbqt = None
            self.skip_prep = False
            self.prefix = "docking_batch"

        # Handle directional or uniform padding (Ported from GPU version)
        raw_padding = kwargs.get("box_padding", 10.0)
        if isinstance(raw_padding, (int, float)):
            parsed_padding = [float(raw_padding)] * 3
        elif isinstance(raw_padding, (list, tuple)) and len(raw_padding) == 3:
            parsed_padding = [float(x) for x in raw_padding]
        else:
            raise ValueError("box_padding must be a single number or a list/tuple of three numbers (x, y, z).")

        self.config = {
            "box_center": kwargs.get("box_center", "auto"),
            "box_size": kwargs.get("box_size", "auto"),
            "box_padding": parsed_padding,
            "min_box_size": kwargs.get("min_box_size", 20.0),
            "mamba_env": kwargs.get("mamba_env", "analysis"),
            "mamba_exe": kwargs.get("mamba_exe", "micromamba"),
            "vina_exe": kwargs.get("vina_exe", "vina"),
            "vina_cmd": kwargs.get("vina_cmd", "vina"),
            "singularity_binds": kwargs.get("singularity_binds", [])
        }

        self.extra_receptor_args = kwargs.get("extra_receptor_args", {})
        self.extra_docking_args = kwargs.get("extra_docking_args", {})
        self.target_chain = kwargs.get("target_chain", None)
        self.target_residues = []

        target_csv_path = kwargs.get("target_residues_file", None)

        if target_csv_path and os.path.exists(target_csv_path):
            with open(target_csv_path) as f:
                for line in f.readlines()[1:]:
                    if not line.strip(): continue
                    parts = line.strip().split(',')
                    if len(parts) >= 2:
                        chain = parts[0].strip()
                        if self.target_chain and chain != self.target_chain: continue
                        self.target_residues.append((chain, int(parts[1].strip())))

        self.vina_config = os.path.join(self.output_dir, f"{self.prefix}_vina_config.txt")
        self.log_file = os.path.join(self.output_dir, "batch_docking.log")
        self.run_config = kwargs.get("run_config", {})

        open(self.log_file, "w").close()
        configure_logging(self.log_file)
        logging.info("Starting Vina CPU Pipeline")
        write_config_snapshot(self.log_file, self.run_config)

    def _run_command(self, cmd, step_name):
        logging.info(f"STARTING: {step_name}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable='/bin/bash')
        if result.returncode != 0:
            logging.error(f"FAILED: {step_name}\n{result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, cmd)
        logging.info(f"FINISHED: {step_name}")
        return result.stdout

    def _calculate_box_parameters(self):
        if not self.receptor_pdb:
            raise ValueError(
                "Auto-calculating the box requires the raw 'receptor_pdb' file. "
                "Please provide it in the pipeline initialization, or explicitly set hardcoded 'box_center' and 'box_size' values."
            )

        logging.info("Calculating binding site dimensions from raw PDB...")
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", self.receptor_pdb)
        coords = [atom.coord for model in structure for chain in model for res in chain
                  if (chain.id, res.id[1]) in self.target_residues for atom in res]

        if not coords:
            raise ValueError("No target residues found in the PDB!")

        coords = np.array(coords)
        centroid = coords.mean(axis=0)
        dimensions = coords.max(axis=0) - coords.min(axis=0)

        # Unpack directional padding
        pad_x, pad_y, pad_z = self.config['box_padding']
        min_s = float(self.config['min_box_size'])

        # Calculate size based on directional padding
        size_x = int(np.ceil(max(dimensions[0] + pad_x, min_s)))
        size_y = int(np.ceil(max(dimensions[1] + pad_y, min_s)))
        size_z = int(np.ceil(max(dimensions[2] + pad_z, min_s)))

        return f"{centroid[0]:.3f} {centroid[1]:.3f} {centroid[2]:.3f}", f"{size_x} {size_y} {size_z}"

    def step_1_prepare_receptor(self):
        if self.skip_prep:
            logging.info(f"Bypassing receptor preparation. Using provided PDBQT: {self.receptor_pdbqt}")
            return

        print("\n--- STEP 1: PREPARING RECEPTOR ---")
        extras = " ".join([f"--{k} {v}" if v else f"--{k}" for k, v in self.extra_receptor_args.items()])
        cmd = f"{self.config['mamba_exe']} run -p {self.config['mamba_env']} mk_prepare_receptor.py -i {self.receptor_pdb} --write_pdbqt {self.receptor_pdbqt} {extras}"
        self._run_command(cmd, "Receptor Preparation")

    def step_2_write_config(self):
        print("\n--- STEP 2: WRITING VINA CONFIGURATION ---")

        if self.config['box_center'] == "auto" or self.config['box_size'] == "auto":
            c, s = self._calculate_box_parameters()
            if self.config['box_center'] == "auto": self.config['box_center'] = c
            if self.config['box_size'] == "auto": self.config['box_size'] = s

        cx, cy, cz = self.config['box_center'].split()
        sx, sy, sz = self.config['box_size'].split()

        with open(self.vina_config, 'w') as f:
            f.write("# Input/Output\n")
            f.write(f"receptor = {self.receptor_pdbqt}\n")
            f.write(f"batch = {self.ligand_dir}\n")
            f.write(f"dir = {self.output_dir}\n\n")

            f.write("# Search Space\n")
            f.write(f"center_x = {cx}\ncenter_y = {cy}\ncenter_z = {cz}\n\n")
            f.write(f"size_x = {sx}\nsize_y = {sy}\nsize_z = {sz}\n\n")

            f.write("# Performance Settings\n")
            for k, v in self.extra_docking_args.items():
                f.write(f"{k} = {v}\n" if v is not None else f"{k}\n")

        logging.info(f"Configuration saved to {self.vina_config}")

    def step_3_run_docking(self):
        print("\n--- STEP 3: RUNNING BATCH DOCKING ---")

        ligand_count = len(list(Path(self.ligand_dir).glob("*.pdbqt")))
        logging.info(f"Passing directory with {ligand_count} .pdbqt files to Vina processor...")

        vina_path = self.config['vina_exe']
        is_sif = vina_path.endswith('.sif')

        if is_sif:
            internal_cmd = self.config['vina_cmd']
            logging.info(f"Detected Singularity container (.sif). Using internal command: {internal_cmd}")

            # Use absolute path if the user just provided the .sif name and it's in the current dir
            if not os.path.exists(vina_path):
                raise FileNotFoundError(f"Singularity image not found at: {vina_path}")

            # Dynamically build the bind string
            bind_list = self.config['singularity_binds']
            bind_flag = f"-B {','.join(bind_list)} " if bind_list else ""

            # Removed --nv since this is the CPU pipeline
            cmd = (
                f"singularity exec -u --fakeroot {bind_flag}"
                f"{vina_path} {internal_cmd} --config {self.vina_config} "
                f">> {self.log_file} 2>&1"
            )
        else:
            logging.info("Detected native execution environment (Mamba/Conda).")
            cmd = (
                f"{self.config['mamba_exe']} run -p {self.config['mamba_env']} "
                f"{vina_path} --config {self.vina_config} "
                f">> {self.log_file} 2>&1"
            )

        result = subprocess.run(cmd, shell=True, executable='/bin/bash')

        if result.returncode != 0:
            logging.error(f"Vina failed. Check the log file at {self.log_file}")
            raise subprocess.CalledProcessError(result.returncode, cmd)
        else:
            logging.info("Batch docking complete! All models saved to output directory.")

    def build_summary(self, status="success", error=None):
        ligand_count = len(list(Path(self.ligand_dir).glob("*.pdbqt"))) if self.ligand_dir else 0
        output_files = sorted(str(path) for path in Path(self.output_dir).glob("*.pdbqt")) if self.output_dir else []
        warnings = []
        if ligand_count == 0:
            warnings.append("No .pdbqt ligand files were found in ligand_dir.")

        return {
            "stage": "vina_cpu",
            "status": status,
            "timestamp": datetime.now().isoformat(),
            "inputs": {
                "receptor_pdb": self.receptor_pdb,
                "receptor_pdbqt": self.receptor_pdbqt,
                "ligand_dir": self.ligand_dir,
                "target_residues": self.target_residues,
            },
            "configuration": {
                "box_center": self.config.get("box_center"),
                "box_size": self.config.get("box_size"),
                "box_padding": self.config.get("box_padding"),
                "extra_docking_args": self.extra_docking_args,
                "skip_receptor_prep": self.skip_prep,
            },
            "outputs": {
                "output_dir": self.output_dir,
                "vina_config": self.vina_config,
                "log_file": self.log_file,
                "output_pdbqt_files": output_files,
            },
            "metrics": {
                "ligand_count": ligand_count,
                "output_pdbqt_count": len(output_files),
            },
            "warnings": warnings,
            "error": error,
        }

    def run_all(self):
        if not self.receptor_pdbqt or not self.ligand_dir or not self.output_dir:
            raise ValueError("Missing essential paths! You must provide at least a receptor (PDB or PDBQT), a ligand directory, and an output directory.")

        self.step_1_prepare_receptor()
        self.step_2_write_config()
        self.step_3_run_docking()
        return self.build_summary()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vina CPU Batch Docking Pipeline")
    parser.add_argument("--config", required=True, help="Path to the YAML configuration file")
    parser.add_argument("--json-out", help="Optional path for a machine-readable JSON run summary")
    args = parser.parse_args()

    try:
        with open(args.config) as file:
            config_data = yaml.safe_load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find the configuration file: {args.config}")
    except yaml.YAMLError as e:
        raise ValueError(f"Error reading the YAML file. Please check your formatting:\n{e}")

    paths = config_data.get('paths', {})
    box = config_data.get('box', {})
    pocket_definition = paths.get('pocket_definition')
    if pocket_definition:
        with open(pocket_definition) as handle:
            pocket_data = yaml.safe_load(handle) or {}
        pocket_box = pocket_data.get("box", {})
        center = pocket_box.get("center")
        size = pocket_box.get("size")
        if not isinstance(center, list) or len(center) != 3 or not isinstance(size, list) or len(size) != 3:
            raise ValueError("pocket_definition must contain box.center and box.size lists with three values each.")
        box = dict(box)
        box["center"] = " ".join(f"{float(value):.3f}" for value in center)
        box["size"] = " ".join(f"{float(value):.3f}" for value in size)
    docking = config_data.get('docking', {})
    env = config_data.get('environment', {})

    pipeline = None
    json_out = args.json_out or config_data.get("reporting", {}).get("json_out")
    try:
        pipeline = VinaCPUPipeline(
            receptor_pdbqt=paths.get('receptor_pdbqt'),
            receptor_pdb=paths.get('receptor_pdb'),
            ligand_dir=paths.get('ligand_dir'),
            output_dir=paths.get('output_dir'),
            target_residues_file=paths.get('target_residues_file'),

            box_center=box.get('center', 'auto'),
            box_size=box.get('size', 'auto'),
            box_padding=box.get('padding', [10.0, 10.0, 10.0]),
            target_chain=box.get('target_chain'),

            mamba_env=env.get('mamba_env', 'analysis'),
            mamba_exe=env.get('mamba_exe', 'micromamba'),
            vina_exe=env.get('vina_exe', 'vina'),
            vina_cmd=env.get('vina_cmd', 'vina'),
            singularity_binds=env.get('singularity_binds', []),
            run_config=config_data,

            extra_docking_args={
                "cpu": docking.get('cpu', 8),
                "exhaustiveness": docking.get('exhaustiveness', 32),
                "num_modes": docking.get('num_modes', 9),
                "min_rmsd": docking.get('min_rmsd', 3),
                "energy_range": docking.get('energy_range', 3),
                "scoring": docking.get('scoring', 'vina'),
            }
        )
        summary = pipeline.run_all()
        write_json_summary(summary, json_out or Path(pipeline.output_dir) / "dockflow_vina_cpu_summary.json")
    except Exception as exc:
        if pipeline:
            failure = pipeline.build_summary(status="failed", error=str(exc))
            write_json_summary(failure, json_out or Path(pipeline.output_dir) / "dockflow_vina_cpu_summary.json")
        raise
