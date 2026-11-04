"""Fertile blanket, fission ledger, depletion, schema, and yield-scale injection."""

import math
from pathlib import Path

import pytest

from solcascade.geometry.shell import build_shell
from solcascade.hybrid.blanket import blanket_from_config
from solcascade.hybrid.depletion import deplete_pulses
from solcascade.hybrid.rprocess import score_fission_energy_fraction
from solcascade.interfaces import MappingProvider, collect_slots
from solcascade.load import load_problem
from solcascade.materials.library import Material
from solcascade.pipeline import run_pipeline
from solcascade.transport.analytical import first_collision
from solcascade.transport.openmc_backend import write_model_script

ROOT = Path(__file__).resolve().parents[1]
MINIMAL = ROOT / "configs" / "minimal_pipeline.yaml"
HYBRID_U238 = ROOT / "configs" / "hybrid_u238.yaml"
FUSION_J = 0.4977033387900156


def _metal(name: str, nuclide: str) -> Material:
    return Material(
        name=name,
        density_kg_m3=19050.0,
        thermal_conductivity_W_mK=27.0,
        specific_heat_J_kgK=116.0,
        youngs_modulus_Pa=2.08e11,
        poissons_ratio=0.23,
        thermal_expansion_per_K=1.39e-5,
        yield_strength_Pa=2.0e8,
        resistivity_ohm_m=None,
        resistivity_temperature_coeff_per_K=None,
        resistivity_reference_K=None,
        one_group_macroscopic_total_per_m=10.0,
        one_group_kerma_fraction=0.0,
        composition=((nuclide, 1.0),),
        note="",
        one_group_fission_fraction=0.25,
        one_group_capture_fraction=0.05,
    )


def test_fertile_fuel_tally_matches_first_collision(tmp_path: Path):
    shell = build_shell(
        0.2,
        [{"name": "fertile", "material": "uranium_238", "thickness_m": 0.1}],
        radial_cells=1,
    )
    material = _metal("uranium_238", "U238")
    deposition = first_collision(shell, {"uranium_238": material})
    collision = 1.0 - math.exp(-1.0)
    assert float(deposition.fission_per_source.sum()) == pytest.approx(collision * 0.25)
    assert float(deposition.capture_per_source.sum()) == pytest.approx(collision * 0.05)
    script = tmp_path / "openmc_model.py"
    write_model_script(
        script,
        shell,
        {"uranium_238": material},
        neutron_energy_ev=1.41e7,
        particles=10,
        batches=2,
        photon_transport=False,
        seed=1,
        source_rate_per_s=1.0,
        temperature_k=600.0,
    )
    text = script.read_text(encoding="utf-8")
    assert "add_nuclide('U238'" in text
    assert "fission-rate" in text
    assert "H3-production" in text


def test_fission_ledger_is_zero_on_the_pure_fusion_baseline(tmp_path: Path):
    summary = run_pipeline(MINIMAL, workdir=tmp_path)
    conversion = summary["conversion"]
    hybrid = summary["hybrid"]
    assert hybrid["fuel_type"] == "none"
    assert hybrid["fission_thermal_J"] == 0.0
    assert conversion["fission_thermal_J"] == 0.0
    assert conversion["electric_fission_J"] == 0.0
    assert conversion["fusion_J"] == pytest.approx(FUSION_J, rel=1.0e-12)
    assert conversion["electric_shell_J"] == 0.0
    assert conversion["electric_gross_J"] - conversion["supply_drawn_J"] == pytest.approx(
        conversion["net_electrical_J"]
    )


def test_depletion_step_conserves_heavy_atoms():
    result = deplete_pulses("u238", 1.0e6, 100.0, 1000.0, 1, 86400.0)
    assert result["inventory"]["Pu239"] > 0.0
    heavy = sum(result["inventory"].values())
    assert heavy + result["fissioned"] == pytest.approx(result["conserved_atoms"])
    assert result["conserved_atoms"] == pytest.approx(result["initial_parent_atoms"])
    assert result["history"][-1]["heavy_plus_fissioned"] == pytest.approx(1.0e6)


def test_yaml_schema_accepts_hybrid_slots_and_rejects_unknown_fuel():
    config, materials = load_problem(HYBRID_U238)
    slots = collect_slots(config, materials)
    for key in (
        "blanket.thickness_m",
        "blanket.density_g_cm3",
        "blanket.enrichment",
        "blanket.neutron_yield_scale",
        "conversion.thermal_efficiency",
    ):
        assert key in slots
    assert config["blanket"]["fuel_type"] == "u238"
    assert config["fission"]["backend"] == "serpent"
    assert slots["conversion.thermal_efficiency"] == pytest.approx(0.33)
    bad = dict(config)
    bad["blanket"] = dict(config["blanket"])
    bad["blanket"]["fuel_type"] = "pu239"
    with pytest.raises(ValueError, match="fuel_type"):
        blanket_from_config(bad, slots)


def test_provider_neutron_yield_scale_leaves_multiplication_fixed(tmp_path: Path):
    base = run_pipeline(HYBRID_U238, workdir=tmp_path / "base")
    scaled = run_pipeline(
        HYBRID_U238,
        workdir=tmp_path / "scaled",
        provider=MappingProvider({"blanket.neutron_yield_scale": 2.0}),
    )
    assert base["hybrid"]["fission_thermal_J"] > 0.0
    assert scaled["hybrid"]["energy_multiplication_M"] == pytest.approx(
        base["hybrid"]["energy_multiplication_M"]
    )
    assert scaled["conversion"]["fission_thermal_J"] == pytest.approx(
        2.0 * base["conversion"]["fission_thermal_J"]
    )
    assert scaled["conversion"]["electric_fission_J"] == pytest.approx(
        0.33 * scaled["conversion"]["fission_thermal_J"]
    )
    assert scaled["conversion"]["fusion_J"] == pytest.approx(base["conversion"]["fusion_J"])
    assert scaled["conversion"]["fusion_J"] == pytest.approx(FUSION_J, rel=1.0e-12)
    assert scaled["conversion"]["supply_drawn_J"] == pytest.approx(
        base["conversion"]["supply_drawn_J"]
    )
    assert scaled["overrides"]["blanket.neutron_yield_scale"] == "mapping"
    assert scaled["provider"] == "mapping"
    assert scaled["conversion"]["electric_gross_J"] - scaled["conversion"]["supply_drawn_J"] == pytest.approx(
        scaled["conversion"]["net_electrical_J"]
    )
    assert base["hybrid"]["depletion"]["pulses"] == 100
    assert base["hybrid"]["fission_backend"]["executed"] is False
    assert (tmp_path / "base" / "isotopic_inventory.json").is_file()
    assert (tmp_path / "base" / "blanket_serpent.txt").is_file()
    assert score_fission_energy_fraction(0.10)["residual"] == 0.0
