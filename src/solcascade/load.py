"""YAML loading. Paths in a config file are relative to that file."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path} did not contain a mapping")
    return data


def flatten_material(raw: dict) -> dict:
    material = dict(raw)
    one_group = material.pop("one_group", None)
    if one_group:
        material["one_group_macroscopic_total_per_m"] = one_group["macroscopic_total_per_m"]
        material["one_group_kerma_fraction"] = one_group["kerma_fraction"]
        for source, target in (
            ("fission_fraction", "one_group_fission_fraction"),
            ("capture_fraction", "one_group_capture_fraction"),
            ("tritium_fraction", "one_group_tritium_fraction"),
        ):
            if source in one_group and one_group[source] is not None:
                material[target] = one_group[source]
    return material


def load_problem(path: Path) -> tuple[dict, dict[str, dict]]:
    config_path = Path(path)
    config = load_yaml(config_path)
    materials_name = config.get("materials_file", "materials.yaml")
    materials_path = (config_path.parent / materials_name).resolve()
    raw_materials = {
        name: flatten_material(body)
        for name, body in load_yaml(materials_path).items()
    }
    config["_config_path"] = str(config_path)
    config["_materials_path"] = str(materials_path)
    return config, raw_materials
