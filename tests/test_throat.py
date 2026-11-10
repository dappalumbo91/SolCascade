"""Bore throat stronger than the frozen internal field."""

import math

import numpy as np
import pytest

from solcascade.geometry.coils import Coil, CoilSet
from solcascade.production.nature import isotropic_loss_fraction
from solcascade.production.refine import load_checked_in, project_root
from solcascade.production.throat import (
    build_throat,
    default_output_path,
    main,
    mirror_ratio_for_loss,
    render,
)


def test_half_trapped_ratio_inverts_the_loss_cone():
    assert mirror_ratio_for_loss(1.0) == 1.0
    assert mirror_ratio_for_loss(0.5) == pytest.approx(4.0 / 3.0)
    assert mirror_ratio_for_loss(0.1) == pytest.approx(1.0 / 0.19)
    assert isotropic_loss_fraction(mirror_ratio_for_loss(0.5)) == pytest.approx(0.5)
    assert isotropic_loss_fraction(mirror_ratio_for_loss(0.1)) == pytest.approx(0.1)
    with pytest.raises(ValueError):
        mirror_ratio_for_loss(0.0)
    with pytest.raises(ValueError):
        mirror_ratio_for_loss(1.5)


@pytest.fixture(scope="module")
def report():
    return build_throat()


def test_bore_pair_is_stronger_than_the_internal_field(report):
    assert report.radius_m == math.nextafter(report.wire_m, math.inf)
    assert report.design.mirror_ratio == pytest.approx(4.0 / 3.0)
    assert report.design.loss_fraction == pytest.approx(0.5)
    assert report.design.throat_field_t > report.design.internal_field_t
    assert 2.0 < report.design.throat_field_t < 2.5
    assert 0.2 < report.design.stored_j < 0.4
    assert report.design.turns == 1.0
    assert report.design.current_a == pytest.approx(report.design.ampere_turns)
    assert report.design.current_a < report.conductor_a
    assert report.design.persistent_draw_j == 0.0
    assert report.design.pulsed_supply_j == pytest.approx(
        (1.0 - report.recovery) * report.design.stored_j
    )
    assert report.design.trapped_electric_j == pytest.approx(0.25 * report.walk_electric_j)
    assert report.design.trapped_electric_j < report.walk_electric_j
    assert report.design.trapped_electric_j < report.design.pulsed_supply_j
    assert report.threshold.mirror_ratio == 1.0
    assert report.threshold.loss_fraction == 1.0
    assert report.threshold.trapped_electric_j == 0.0
    assert report.threshold.stored_j < report.design.stored_j < report.tenth.stored_j
    assert report.tenth.loss_fraction == pytest.approx(0.1)
    assert report.tenth.mirror_ratio == pytest.approx(mirror_ratio_for_loss(0.1))
    assert report.tenth.current_a == report.conductor_a
    assert report.tenth.turns * report.tenth.current_a == pytest.approx(report.tenth.ampere_turns)
    assert report.bare_internal_t == pytest.approx(1.6696763518751803)
    assert report.bare_ratio < 1.0
    assert report.vacuum_rise > report.design.mirror_ratio * report.compression
    assert report.external.steady_reaches is False
    assert report.external.asymptotic_ratio < 1.0
    assert report.external.vacuum_rise == pytest.approx(6.68753568754885)
    assert report.external.stored_j > 1000.0 * report.design.stored_j
    assert report.checked_fusion_j == 0.4977033387900156
    problem = load_checked_in()
    assert [coil.name for coil in problem.coils.coils] == [
        "mirror_inner_minus",
        "mirror_inner_plus",
        "mirror_outer_minus",
        "mirror_outer_plus",
    ]
    winding = CoilSet(
        (
            Coil(
                "throat_minus",
                report.radius_m,
                -report.z_m,
                report.design.turns,
                report.design.current_a,
                report.wire_m,
                "copper",
            ),
            Coil(
                "throat_plus",
                report.radius_m,
                report.z_m,
                report.design.turns,
                report.design.current_a,
                report.wire_m,
                "copper",
            ),
        )
    )
    currents = np.array([report.design.current_a, report.design.current_a], dtype=float)
    assert winding.magnetic_energy_j(currents) == pytest.approx(report.design.stored_j)
    base = problem.coils.base_currents_a()
    scales = problem.schedule.current_scales(base.size)
    peak_index = int(np.argmin(problem.schedule.plasma_radius_m))
    checked = float(problem.coils.on_axis_field(report.z_m, base * scales[peak_index]))
    throat = float(winding.on_axis_field(report.z_m, currents))
    assert checked + throat == pytest.approx(report.design.throat_field_t)


def test_report_names_the_bore_pair_and_writes_under_runs(report):
    text = render()
    assert text == render()
    assert "The gyro-walk remains the cross-field clock." in text
    assert "Half of an isotropic distribution is trapped" in text
    assert "The bore pair is stronger than the internal field." in text
    assert "The checked-in coil list is unchanged." in text
    assert "A path whose strongest field is the occupied field has no trapped set." in text
    assert "A persistent throat current draws no magnetic energy from the supply." in text
    assert "Trapped particles stay on the gyro-walk." in text
    assert "Trapped electricity versus the gyro-walk electricity: below." in text
    assert "Trapped electricity versus the pulsed supply: below." in text
    assert "A steady pair outside the shell stays below the internal field." in text
    assert "Bore stored energy versus that pulsed store: below." in text
    assert "The checked-in schedule is not rewritten." in text
    assert f"{report.design.stored_j:.6g}" in text
    assert f"{report.checked_fusion_j:.6g}" in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    nature = readme.index("## Nature")
    throat = readme.index("## Throat")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    tests = readme.index("## Tests")
    assert nature < throat < plants < observables < tests
    assert "py -3.14 -m solcascade.production.throat" in readme
    source = (project_root() / "src" / "solcascade" / "production" / "throat.py").read_text(
        encoding="utf-8"
    )
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "throat" not in package
