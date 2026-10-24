"""Explicit run directory keeps controlled outputs outside the code package."""
import os
from pathlib import Path

def required_env(name):
    value = os.environ.get(name)
    if not value:
        raise SystemExit(f"Set {name}; see WORKFLOW.md. No default patient input is bundled.")
    return value

def project_root():
    root = Path(required_env("MVA_PROJECT_ROOT")).expanduser().resolve()
    package = Path(__file__).resolve().parents[1]
    if root == package or package in root.parents:
        raise SystemExit("Use a controlled run directory outside the code package.")
    if not root.is_dir():
        raise SystemExit("MVA_PROJECT_ROOT must be an existing controlled run directory.")
    return root
