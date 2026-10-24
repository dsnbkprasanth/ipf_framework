"""Configuration loading for the IPF framework."""
from __future__ import annotations
import os
from pathlib import Path
import yaml

PKG_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PKG_DIR.parent
CONFIG_DIR = PROJECT_DIR / "config"


def _load(name: str) -> dict:
    with open(CONFIG_DIR / name) as fh:
        return yaml.safe_load(fh)


def load_datasets() -> dict:
    return _load("datasets.yaml")


def load_modules() -> dict:
    """Return {module_id: {'label': str, 'genes': [str, ...]}}."""
    raw = _load("modules.yaml")["modules"]
    return {k: {"label": v.get("label", k), "genes": list(v["genes"])}
            for k, v in raw.items()}


def load_ligand_receptor() -> list:
    return _load("ligand_receptor.yaml")["axes"]


def load_drug_safety() -> dict:
    return _load("drug_safety.yaml")


def default_output_dir() -> Path:
    out = os.environ.get("IPFFX_OUTPUT", PROJECT_DIR / "results")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    (out / "tables").mkdir(exist_ok=True)
    return out
