"""Neutron, gamma, and delayed balance on the ranked U-235 kilogram."""

import pytest

from solcascade.energy_baseline.phase2.map_fsot import PIN_PREFIX
from solcascade.production.fission_balance import (
    BALANCE_SHARE,
    build_balance,
    default_output_path,
    delayed_proportions,
    main,
    render,
)
from solcascade.production.fission_paths import ADDED_MEV, LEDGER_MEV
from solcascade.production.refine import project_root


def test_delayed_proportions_count_neutrons():
    rows = delayed_proportions()
    assert [row.nucleus for row in rows] == ["U-235", "U-233", "Pu-239"]
    assert [row.percent for row in rows] == [0.66, 0.27, 0.23]
    assert rows[0].proportion == pytest.approx(0.0066)
    assert BALANCE_SHARE == ADDED_MEV / LEDGER_MEV
    assert BALANCE_SHARE == pytest.approx(0.15)


@pytest.fixture(scope="module")
def report():
    return build_balance()


def test_balance_stays_inside_the_rankine_product(report):
    assert report.balance_share + 0.85 == pytest.approx(1.0)
    assert report.balance_inventory_j == pytest.approx(report.balance_share * report.density_j)
    assert report.fragment_inventory_j + report.balance_inventory_j == pytest.approx(report.density_j)
    assert report.balance_electricity_j == pytest.approx(report.efficiency * report.balance_inventory_j)
    assert report.balance_electricity_j < report.delivered_j
    assert report.ceiling_headroom_j == pytest.approx(
        (report.ceiling - report.efficiency) * report.density_j
    )
    assert report.upper_bound_increment_j == pytest.approx(
        (1.0 - report.efficiency) * report.balance_inventory_j
    )
    assert report.ceiling_headroom_j > report.upper_bound_increment_j
    assert report.balance_nearest.name == "gamma/pi"
    assert report.balance_nearest.value == pytest.approx(0.18373345, rel=1.0e-5)
    assert report.balance_nearest.relative_gap == pytest.approx(0.2248897, rel=1.0e-4)
    assert report.complement.name == "1-G^2"
    assert report.complement.value == pytest.approx(0.16100703028357413)
    assert report.complement.relative_gap < report.balance_nearest.relative_gap
    assert report.complement.relative_gap > 1.0e-3
    assert report.delayed_nearest.name == "gamma/pi"
    assert report.delayed_nearest.relative_gap > 10.0
    assert report.fine_structure == pytest.approx(0.007297352564920234)
    assert report.fine_structure_gap == pytest.approx(0.10565947953336877)
    assert report.fine_structure_gap > report.complement.relative_gap
    assert report.nuclear_scalar == pytest.approx(report.thermodynamics_scalar)
    assert report.nuclear_scalar != pytest.approx(report.balance_share)
    assert report.pin_prefix == PIN_PREFIX
    assert report.digest_prefix == PIN_PREFIX
    assert "physics-of-nuclear-energy" in report.governing_source


def test_report_keeps_the_balance_out_of_a_second_plant(report):
    text = render()
    assert text == render()
    assert "That page assigns the balance, after the fragment kinetic energy, to gamma rays and to the kinetic energy of the neutrons." in text
    assert "The balance is already inside the ranked density." in text
    assert "The Rankine product already includes it." in text
    assert "Scoring that electricity again would count the same joules twice." in text
    assert "The fragment inventory and the balance inventory add to the ranked density." in text
    assert "The catalog files no shield-recovery efficiency." in text
    assert "Neutrons and gamma rays that stop outside the coolant are a shield load." in text
    assert "The shield stays a shield." in text
    assert "Those proportions count neutrons." in text
    assert "They are not an energy share of the ledger." in text
    assert "Delayed neutrons enable control of the chain." in text
    assert "The control margin stays in the core." in text
    assert "That bound is not a filed plant." in text
    assert "The file-ceiling headroom is larger than that bound." in text
    assert "The refinement inside the operating plant is the file efficiency ceiling." in text
    assert "A domain scalar is not the balance share." in text
    assert "The nuclear scalar stays out of the fission joules." in text
    assert "gamma_FSOT does not multiply the plant." in text
    assert "The engine valve stays out of the fission joules." in text
    assert "The balance grammar's nearest seed is `gamma/pi` at " in text
    assert "The complement of the fragment nearest is `1-G^2` at " in text
    assert "The complement sits closer than the balance grammar." in text
    assert "The complement is not installed as the balance." in text
    assert "The fine-structure leaf is not the delayed-neutron proportion." in text
    assert "No seed is installed as the balance or as the delayed proportion." in text
    assert "The checked-in schedule is not rewritten." in text
    counted_twice = report.delivered_j + report.balance_electricity_j
    assert f"{counted_twice:.6g}" not in text
    assert f"{counted_twice:.6e}" not in text
    for scaled in (
        report.nuclear_scalar * report.balance_inventory_j,
        report.nuclear_scalar * report.delivered_j,
        report.gamma_fsot * report.balance_inventory_j,
    ):
        assert f"{scaled:.6g}" not in text
        assert f"{scaled:.6e}" not in text
    reciprocal = 1.0 / report.valve_beta
    assert f"{reciprocal:.6e}" not in text
    assert f"{reciprocal:.6g}" not in text
    alpha_reciprocal = 1.0 / report.fine_structure
    assert f"{alpha_reciprocal:.6g}" not in text
    assert f"{alpha_reciprocal:.12f}" not in text
    assert "3.328" not in text
    assert "3.703" not in text
    assert "warp_shell_energy_j" not in text
    assert "137.035999177" not in text
    output = default_output_path()
    assert output.parent.parent == project_root() / "runs"
    assert main() == 0
    assert output.read_text(encoding="utf-8") == text
    readme = (project_root() / "README.md").read_text(encoding="utf-8")
    paths = readme.index("## Fission paths")
    balance = readme.index("## Fission balance")
    tests = readme.index("## Tests")
    assert paths < balance < tests
    assert "py -3.14 -m solcascade.production.fission_balance" in readme
    source = (
        project_root() / "src" / "solcascade" / "production" / "fission_balance.py"
    ).read_text(encoding="utf-8")
    assert "warp_shell_energy_j" not in source
    for needle in ("137.035999177", "1.602176634", "8071.31713", "1.17302e-9"):
        assert needle not in source
    observables_source = (
        project_root() / "src" / "solcascade" / "production" / "fission_observables.py"
    ).read_text(encoding="utf-8")
    assert "The channel search stays closed." in observables_source
    package = (project_root() / "src" / "solcascade" / "production" / "__init__.py").read_text(
        encoding="utf-8"
    )
    assert "fission_balance" not in package
