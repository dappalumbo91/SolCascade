"""End-to-end pipeline, scalar injection, and generated backend files."""

from pathlib import Path

import pytest

from solcascade.interfaces import MappingProvider, collect_slots, reconcile_energy_density, apply_provider
from solcascade.load import load_problem
from solcascade.pipeline import run_pipeline
from solcascade.transport.openmc_backend import generated_script_text_mentions

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs" / "minimal_pipeline.yaml"


def test_minimal_pipeline_writes_a_closed_ledger(tmp_path: Path):
    summary = run_pipeline(CONFIG, workdir=tmp_path)
    conversion = summary["conversion"]
    assert conversion["alpha_J"] + conversion["neutron_J"] == pytest.approx(conversion["fusion_J"])
    assert (
        conversion["alpha_escaped_J"] + conversion["alpha_wall_J"] + conversion["alpha_retained_J"]
        == pytest.approx(conversion["alpha_J"])
    )
    deposited = conversion["neutron_deposited_J"]
    missing = conversion["neutron_not_deposited_J"]
    backend = summary["transport"]["backend"]
    assert backend in {"analytical", "openmc"}
    # Analytical kerma cannot exceed the neutron kinetic energy. OpenMC
    # heating can, because the tally includes reaction Q.
    if backend == "analytical":
        assert deposited + missing == pytest.approx(conversion["neutron_J"])
    else:
        assert deposited + missing >= conversion["neutron_J"]
        assert summary["transport"]["openmc"]["executed"] is True
    assert conversion["electric_gross_J"] - conversion["supply_drawn_J"] == pytest.approx(
        conversion["net_electrical_J"]
    )
    assert conversion["electric_shell_J"] == 0.0
    assert (tmp_path / "radial_map.csv").is_file()
    assert (tmp_path / "plasma_history.csv").is_file()
    assert (tmp_path / "field_grid.npz").is_file()
    assert (tmp_path / "summary.json").is_file()
    heat = (tmp_path / "heat_conduction.i").read_text(encoding="utf-8")
    assert "coord_type = RSPHERICAL" in heat
    assert "HeatConduction" in heat
    script = (tmp_path / "openmc_model.py").read_text(encoding="utf-8")
    assert "fixed source" in script
    assert "Discrete" in script
    assert ".temperature = 600.0" in script
    assert generated_script_text_mentions(
        __import__("solcascade.assemble", fromlist=["build_shell_from_slots"]).build_shell_from_slots(
            * _resolved_shell_args()
        )
    ) in script
    assert (tmp_path / "cardinal_skeleton.i").is_file()
    assert (tmp_path / "tokamaker_equilibrium.py").is_file()
    assert (tmp_path / "griffin_interface.json").is_file()
    assert (tmp_path / "bison_interface.json").is_file()
    assert summary["plasma"]["peak_triple_product_ratio"] < 1.0
    liquid = next(layer for layer in summary["stress"]["layers"] if layer["layer"] == "absorber")
    assert liquid["structural"] is False


def test_injected_energy_density_rescales_temperature():
    config, materials = load_problem(CONFIG)
    slots = collect_slots(config, materials)
    original = slots["plasma.thermal_energy_density_J_m3"]
    provider = MappingProvider({"plasma.thermal_energy_density_J_m3": original * 4.0})
    resolved, provenance = apply_provider(slots, provider)
    resolved, provenance = reconcile_energy_density(resolved, provenance)
    assert resolved["plasma.ion_temperature_keV"] == pytest.approx(
        slots["plasma.ion_temperature_keV"] * 4.0
    )
    assert provenance["plasma.ion_temperature_keV"].endswith("->temperature")


def test_provider_json_skips_comments_and_nulls():
    import json

    from solcascade.cli import provider_scalars

    payload = json.loads((ROOT / "examples" / "fsot_scalars.example.json").read_text(encoding="utf-8"))
    values = provider_scalars(payload)
    assert values["plasma.ion_density_m3"] == 1.0e22
    assert "plasma.thermal_energy_density_J_m3" not in values
    assert "_comment" not in values
    config, materials = load_problem(CONFIG)
    apply_provider(collect_slots(config, materials), MappingProvider(values))
    with pytest.raises(ValueError):
        provider_scalars({"plasma.ion_density_m3": "high"})


def test_unknown_provider_key_is_rejected():
    config, materials = load_problem(CONFIG)
    slots = collect_slots(config, materials)
    provider = MappingProvider({"plasma.not_a_slot": 1.0})
    with pytest.raises(KeyError):
        apply_provider(slots, provider)


def _resolved_shell_args():
    from solcascade.interfaces import NullProvider

    config, materials = load_problem(CONFIG)
    slots, _provenance = apply_provider(collect_slots(config, materials), NullProvider())
    return config, slots
