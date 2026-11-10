"""Pinned-engine nuclear readings against the stored targets and published bars."""

import pytest

from solcascade.constants import FISSION_RECOVERABLE_ENERGY_MEV
from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_observables import (
    GREEN_GATE_PERCENT,
    NUCLEAR_NAMES,
    PUBLISHED_BARS,
    SUMMARY_GATE_PERCENT,
    build_observables,
    default_output_path,
    main,
    render,
)
from solcascade.production.refine import project_root


@pytest.fixture(scope="module")
def report():
    return build_observables()


def _by_name(rows, name):
    return next(row for row in rows if row.name == name)


def test_wave_readings_stay_on_the_engine_targets(report):
    assert tuple(row.name for row in report.readings) == NUCLEAR_NAMES
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX
    assert report.quantum_mev == FISSION_RECOVERABLE_ENERGY_MEV
    assert report.quantum_mev == 200.0
    assert report.quantum_emitted is False
    assert report.bindings_inside_one is True
    assert report.dialed is False
    assert report.nuclear_scalar == pytest.approx(report.thermodynamics_scalar)
    assert report.nuclear_scalar < 1.0
    assert report.valve_beta < 1.0e-15
    for row in report.readings:
        assert row.error_percent == pytest.approx(
            abs(row.computed - row.wave_target) / abs(row.wave_target) * 100.0
        )
        assert row.inside_green is (row.error_percent <= GREEN_GATE_PERCENT)
        assert row.inside_summary is (row.error_percent < SUMMARY_GATE_PERCENT)
        assert row.inside_green is True
        assert row.inside_summary is True
        if row.percent_tag is None:
            assert row.inside_tag is None
        else:
            assert row.inside_tag is (row.error_percent <= row.percent_tag)
            assert row.inside_tag is True
    assert _by_name(report.readings, "Deuteron_mu_muN").percent_tag is None
    assert "fission" not in " ".join(row.name for row in report.readings).lower()


def test_published_bars_use_the_wave_value(report):
    assert tuple((row.name, row.central, row.uncertainty) for row in report.published) == PUBLISHED_BARS
    for row in report.published:
        reading = _by_name(report.readings, row.name)
        bare = (reading.computed - row.central) / row.uncertainty
        assert row.bare_residual == pytest.approx(bare)
        assert row.bare_inside is (abs(bare) <= 1.0)
        if row.leaf_value is None:
            assert row.leaf_residual is None
            assert row.leaf_inside is None
            assert row.inside_one is row.bare_inside
        else:
            leaf = (row.leaf_value - row.central) / row.uncertainty
            assert row.leaf_residual == pytest.approx(leaf)
            assert row.leaf_inside is True
            assert row.inside_one is True
            assert row.leaf_value != reading.computed
    helium = _by_name(report.published, "He4_binding_MeV")
    triton = _by_name(report.published, "Triton_binding_MeV")
    moment = _by_name(report.published, "mu_p_muN")
    deuteron = _by_name(report.published, "Deuteron_binding_MeV")
    assert helium.bare_inside is False
    assert triton.bare_inside is False
    assert moment.bare_inside is False
    assert deuteron.leaf_formula == "bare*(1+alpha3*(sqrt(e)+phi^-4))"
    assert deuteron.leaf_value == pytest.approx(report.deuteron_rounded.value)
    assert deuteron.inside_one is True
    assert tuple((hit.factor, hit.seed) for hit in report.deuteron_hits) == (
        ("alpha2", "pi^{-4}"),
        ("alpha3", "phi"),
        ("alpha3", "sqrt(e)"),
    )
    assert all(hit.same_side for hit in report.deuteron_hits)
    assert all(abs(hit.residual_uncertainties) <= 1.0 for hit in report.deuteron_hits)
    assert report.deuteron_rounded.factor == "alpha3"
    assert report.deuteron_rounded.piece == "sqrt(e)+phi^-4"
    assert report.deuteron_rounded.same_side is True
    assert report.deuteron_rounded.residual_uncertainties == pytest.approx(-0.00043684852049348424)
    assert report.deuteron_excess.factor == "alpha3"
    assert report.deuteron_excess.piece == "e^3/phi^5"
    assert report.deuteron_excess.same_side is True
    assert report.deuteron_excess.residual_uncertainties == pytest.approx(-0.03386969152956236)
    assert report.helium_mismatch == pytest.approx(4.234313898927058e-06)
    assert report.deuteron_rounded.relative_mismatch > report.helium_mismatch
    assert (report.deuteron_rounded.piece, report.deuteron_rounded.factor) != (
        report.deuteron_excess.piece,
        report.deuteron_excess.factor,
    )
    assert _by_name(report.published, "m_n-m_p_MeV").inside_one is True
    assert _by_name(report.published, "r_p_fm").inside_one is True
    assert _by_name(report.published, "r_p_fm").leaf_value is None
    deuteron_moment = _by_name(report.published, "Deuteron_mu_muN")
    assert deuteron_moment.leaf_formula == "bare*(1+alpha2*sqrt(K)/psi)"
    assert deuteron_moment.leaf_value == pytest.approx(report.deuteron_moment.value)
    assert deuteron_moment.inside_one is True
    assert deuteron_moment.leaf_inside is True
    assert deuteron_moment.bare_residual < -1000.0
    assert report.deuteron_moment.piece == "sqrt(K)/psi"
    assert report.deuteron_moment.factor == "alpha2"
    assert report.deuteron_moment.same_side is False
    assert report.deuteron_moment.residual_uncertainties == pytest.approx(0.11111486972275693)
    assert abs(report.deuteron_moment.residual_uncertainties) < 1.0
    assert report.deuteron_moment.relative_mismatch == pytest.approx(5.2216334752203434e-06)
    assert report.deuteron_moment.relative_mismatch > report.helium_mismatch
    assert abs(helium.leaf_residual) < 1.0
    assert abs(triton.leaf_residual) < 1.0
    assert abs(moment.leaf_residual) < 1.0


def test_report_leaves_the_scalar_out_of_the_quantum(report):
    text = render()
    assert text == render()
    assert "The pinned engine emits the reading." in text
    assert "Every wave reading is inside the green gate." in text
    assert "Every tagged reading is inside its printed percent tag." in text
    assert "A leaf is installed where the accuracy note names one formula." in text
    assert "A derived dressing is installed when that dressing sits inside the published bar." in text
    assert "The helium-4 leaf multiplies the wave formula by one minus alpha squared times pi plus P_base." in text
    assert "The triton leaf multiplies the wave formula by one plus (Poof times Suction) squared times gamma times psi_con squared." in text
    assert "The proton moment leaf is the g-factor leaf divided by two." in text
    assert "The alpha-squared window holds pi to the minus four." in text
    assert "The alpha-cubed window holds phi and the square root of e." in text
    assert "The suction window holds no named seed." in text
    assert "Those are different dressings." in text
    assert "The scored deuteron binding is the rounded same-side nearest." in text
    assert "That dressing sits inside the mass-excess bar." in text
    assert (
        "The helium quotient matches pi plus P_base at "
        f"`{report.helium_mismatch:.6g}`."
    ) in text
    assert (
        "The rounded same-side nearest is "
        f"`{report.deuteron_rounded.piece}` on {report.deuteron_rounded.factor_label} "
        f"at `{report.deuteron_rounded.residual_uncertainties:.6g}`."
    ) in text
    assert (
        "The mass-excess nearest is "
        f"`{report.deuteron_excess.piece}` on {report.deuteron_excess.factor_label} "
        f"at `{report.deuteron_excess.residual_uncertainties:.6g}`."
    ) in text
    assert "The two centrals return different products." in text
    assert "The rounded match is wider than the helium standard." in text
    assert "sqrt(e)+phi^-4" in text
    assert "e^3/phi^5" in text
    for hit in report.deuteron_hits:
        assert f"`{hit.seed}` on {hit.factor_label} reads `{hit.residual_uncertainties:.6g}`." in text
    assert "The helium-4 leaf sits inside one published uncertainty." in text
    assert "The triton leaf sits inside one published uncertainty." in text
    assert "The proton-moment leaf sits inside one published uncertainty." in text
    assert "The deuteron binding sits inside one published uncertainty." in text
    assert "`sqrt(K)/psi` on alpha squared reads " in text
    assert f"The moment quotient match is `{report.deuteron_moment.relative_mismatch:.6g}`." in text
    assert "That dressing crosses the central." in text
    assert "Its quotient match is wider than the helium standard." in text
    assert "The scored deuteron moment is `sqrt(K)/psi` on alpha squared." in text
    assert "The deuteron moment sits inside one published uncertainty." in text
    assert "A domain scalar is not a megaelectronvolt." in text
    assert "The nuclear scalar stays out of the fission joules." in text
    assert "gamma_FSOT does not multiply the plant." in text
    assert "The engine valve stays out of the fission joules." in text
    assert "The 200 MeV ledger quantum is the plant energy." in text
    assert "No wave result carries that quantum." in text
    assert "The channel search stays closed." in text
    assert "The checked-in schedule is not rewritten." in text
    scaled = report.nuclear_scalar * report.quantum_mev
    assert f"{scaled:.6g}" not in text
    assert f"{scaled:.6e}" not in text
    reciprocal = 1.0 / report.valve_beta
    assert f"{reciprocal:.6e}" not in text
    assert f"{reciprocal:.6g}" not in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    assert "137.035999177" not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    plants = readme.index("## Working plants")
    observables = readme.index("## Fission observables")
    tests = readme.index("## Tests")
    assert plants < observables < tests
    assert "py -3.14 -m solcascade.production.fission_observables" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_observables" not in package
