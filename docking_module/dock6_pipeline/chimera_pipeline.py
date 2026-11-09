import argparse
import json
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path

import yaml


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


class ChimeraPipeline:
    def __init__(self, input_pdb, out_dir=None, density=5, 
                 execution_mode="local", chimera_cmd="chimera", 
                 sif_image=None, sif_binds=None, mamba_env=None, run_config=None):
                 
        self.input_pdb = os.path.abspath(input_pdb)
        self.density = density
        self.execution_mode = execution_mode.lower()
        self.chimera_cmd = chimera_cmd
        self.sif_image = sif_image
        self.sif_binds = sif_binds if sif_binds else []
        self.mamba_env = mamba_env
        self.run_config = run_config or {}
        
        # Output Directory Handling
        self.out_dir = os.path.abspath(out_dir) if out_dir else os.getcwd()
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
            
        # Dynamic Prefix Naming
        self.prefix = os.path.splitext(os.path.basename(self.input_pdb))[0]
        self.log_file = os.path.join(self.out_dir, f"{self.prefix}_chimera_pipeline.log")
        
        self.h_added_pdb = os.path.join(self.out_dir, f"{self.prefix}_wH.pdb")
        self.charged_mol2 = os.path.join(self.out_dir, f"{self.prefix}_wH_wcH.mol2")
        self.dms_out = os.path.join(self.out_dir, f"{self.prefix}_surf_d{self.density}.dms")
        
        self.script_dir = os.path.dirname(os.path.abspath(__file__))

    def _build_command(self, com_file):
        """Constructs the exact execution command based on the environment mode."""
        if self.execution_mode == "singularity":
            if not self.sif_image or not os.path.exists(self.sif_image):
                raise FileNotFoundError(f"Singularity image not found: {self.sif_image}")
            
            cmd = ["singularity", "exec"]
            
            # Combine user-provided binds with the script directory (so it finds companion scripts)
            all_binds = set(self.sif_binds)
            all_binds.add(self.script_dir)
            
            for b in all_binds:
                cmd.extend(["-B", b])
            
            # Equivalent to: singularity exec -B <binds> <sif> chimera --nogui <com_file>
            cmd.extend([self.sif_image, "chimera", "--nogui", com_file])
            return cmd
            
        elif self.execution_mode == "micromamba":
            if not self.mamba_env:
                raise ValueError("mamba_env must be provided in micromamba mode.")
            return ["micromamba", "run", "-n", self.mamba_env, self.chimera_cmd, "--nogui", com_file]
            
        else: # "local"
            return [self.chimera_cmd, "--nogui", com_file]

    def has_hydrogens(self):
        print(f"[*] Checking {os.path.basename(self.input_pdb)} for existing protein hydrogens...")
        with open(self.input_pdb) as f:
            for line in f:
                if line.startswith("ATOM") and (line[76:78].strip() == 'H' or line[12:16].strip().startswith('H')):
                    print("[+] Protein hydrogens detected. Skipping AddH step.")
                    return True
        print("[-] No standard protein hydrogens found. AddH step is required.")
        return False

    def generate_com_file(self, step_name, input_file, script_filename):
        com_filename = os.path.join(self.out_dir, f"{self.prefix}_{step_name}.com")
        absolute_script_path = os.path.join(self.script_dir, script_filename)
        
        if not os.path.exists(absolute_script_path):
            raise FileNotFoundError(f"Missing required script '{script_filename}' in {self.script_dir}")

        with open(com_filename, 'w') as f:
            f.write(f"open {input_file}\n")
            f.write(f"runscript {absolute_script_path}\n")
        return com_filename

    def execute_chimera(self, com_file, step_name):
        print(f"[*] Executing {step_name} via [{self.execution_mode.upper()}]...")
        
        env = os.environ.copy()
        env["SURF_DENSITY"] = str(self.density)
        env["CHIMERA_OUT_PDB"] = self.h_added_pdb
        env["CHIMERA_OUT_MOL2"] = self.charged_mol2
        env["CHIMERA_OUT_DMS"] = self.dms_out
        
        cmd = self._build_command(com_file)
        result = subprocess.run(cmd, env=env, capture_output=True, text=True)
        
        with open(self.log_file, 'a') as f:
            f.write(f"\n{'='*50}\n--- OUTPUT FOR {step_name.upper()} ---\n{'='*50}\n")
            f.write(f"Executed Command: {' '.join(cmd)}\n\n")
            f.write(result.stdout)
            if result.stderr: f.write("\n[STDERR]:\n" + result.stderr)
        
        if os.path.exists(com_file): os.remove(com_file)
        if result.returncode != 0:
            raise RuntimeError(f"Chimera crashed during {step_name}. Check log: {self.log_file}")
        return result.stdout

    def check_charge_integrity(self, stdout_text):
        match = re.search(r"Total charge for #[0-9]+:\s*(-?\d+\.\d+)", stdout_text)
        if match:
            charge_val = float(match.group(1))
            if charge_val.is_integer():
                print(f"[+] Integrity check passed: Total charge is {match.group(1)} (Perfect Integer).")
            else:
                print(f"\n[!] ALARM: Non-integral total charge detected ({match.group(1)})!")
        else:
            print("[-] WARNING: Could not parse total charge from Chimera output.")

    def run(self):
        with open(self.log_file, 'w') as f: f.write(f"Starting Pipeline for {self.input_pdb}\n")
        write_config_snapshot(self.log_file, self.run_config)

        current_input = self.input_pdb
        if not self.has_hydrogens():
            com_file = self.generate_com_file("step1_addH", current_input, "chimera_add_hydrogen.py")
            self.execute_chimera(com_file, "Add Hydrogens")
            current_input = self.h_added_pdb
        
        com_file = self.generate_com_file("step2_charge", current_input, "chimera_add_charge.py")
        self.check_charge_integrity(self.execute_chimera(com_file, "Add Charges"))
        current_input = self.charged_mol2

        com_file = self.generate_com_file("step3_dms", current_input, "chimera_dms.py")
        self.execute_chimera(com_file, "Surface Generation")
        print(f"\n[+] Pipeline complete! Final surface: {self.dms_out}")
        return self.build_summary()

    def build_summary(self, status="success", error=None):
        outputs = {
            "output_dir": self.out_dir,
            "log_file": self.log_file,
            "hydrogenated_pdb": self.h_added_pdb,
            "charged_mol2": self.charged_mol2,
            "surface_dms": self.dms_out,
        }
        return {
            "stage": "chimera_receptor_prep",
            "status": status,
            "timestamp": datetime.now().isoformat(),
            "inputs": {
                "input_pdb": self.input_pdb,
            },
            "configuration": {
                "density": self.density,
                "execution_mode": self.execution_mode,
                "chimera_cmd": self.chimera_cmd,
                "sif_image": self.sif_image,
                "sif_binds": self.sif_binds,
                "mamba_env": self.mamba_env,
            },
            "outputs": outputs,
            "metrics": {
                "hydrogenated_pdb_exists": os.path.exists(self.h_added_pdb),
                "charged_mol2_exists": os.path.exists(self.charged_mol2),
                "surface_dms_exists": os.path.exists(self.dms_out),
            },
            "warnings": [],
            "error": error,
        }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", default=None)
    parser.add_argument("-o", "--out_dir", default=None)
    parser.add_argument("-d", "--density", type=int, default=None)
    parser.add_argument("--config", default=None, help="Path to YAML config")
    
    # Execution Mode Overrides
    parser.add_argument("--mode", choices=["local", "singularity", "micromamba"], default=None)
    parser.add_argument("--sif", default=None, help="Path to .sif image")
    parser.add_argument("-b", "--bind", nargs="*", default=[], help="Directories to bind mount for Singularity")
    parser.add_argument("--env", default=None, help="Micromamba environment name")
    parser.add_argument("--json-out", help="Optional path for a machine-readable JSON run summary")
    
    args = parser.parse_args()
    kwargs = {"execution_mode": "local", "chimera_cmd": "chimera", "density": 5}
    
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f)
            kwargs.update({k: v for k, v in cfg.items() if k in 
                          ['input_pdb', 'out_dir', 'density', 'execution_mode', 
                           'chimera_cmd', 'sif_image', 'sif_binds', 'mamba_env']})

    if args.input: kwargs["input_pdb"] = args.input
    if args.out_dir: kwargs["out_dir"] = args.out_dir
    if args.density is not None: kwargs["density"] = args.density
    if args.mode: kwargs["execution_mode"] = args.mode
    if args.sif: kwargs["sif_image"] = args.sif
    if args.env: kwargs["mamba_env"] = args.env
    
    # Merge command-line binds with any YAML binds
    if args.bind:
        existing_binds = kwargs.get("sif_binds", [])
        kwargs["sif_binds"] = list(set(existing_binds + args.bind))

    if "input_pdb" not in kwargs or not kwargs["input_pdb"]:
        parser.error("Input PDB required via -i or inside YAML.")

    json_out = args.json_out or kwargs.get("reporting", {}).get("json_out")
    kwargs["run_config"] = kwargs.copy()
    pipeline = None
    try:
        pipeline = ChimeraPipeline(**kwargs)
        summary = pipeline.run()
        write_json_summary(summary, json_out or Path(pipeline.out_dir) / "dockflow_chimera_summary.json")
    except Exception as exc:
        if pipeline:
            failure = pipeline.build_summary(status="failed", error=str(exc))
            write_json_summary(failure, json_out or Path(pipeline.out_dir) / "dockflow_chimera_summary.json")
        raise
