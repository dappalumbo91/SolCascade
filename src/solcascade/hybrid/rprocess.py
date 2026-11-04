"""Comparison library for r-process fission cycling. Not a nucleosynthesis run.

Stored numbers are ones read from the papers named on each record. A SolCascade
blanket multiplication is not one of these yields. ``score_fission_energy_fraction``
returns the residual of a caller-supplied fraction against the Barnes benchmark.
It does not invent an FSOT neutron yield.
"""

from __future__ import annotations

# Barnes, Kasen, Wu, and Martinez-Pinedo, ApJ 829, 110 (2016),
# DOI 10.3847/0004-637X/829/2/110, arXiv:1605.07218.
# Read from the paper text. Table 1 fit coefficients are not copied: the
# extracted table lost its row labels.
BARNES_2016 = {
    "citation": "Barnes, Kasen, Wu, and Martinez-Pinedo, ApJ 829, 110 (2016)",
    "fiducial_ejecta_mass_solar": 5.0e-3,
    "fiducial_velocity_c": 0.2,
    "representative_initial_electron_fraction": 0.04,
    "fission_energy_fraction_before_1_day": 0.10,
    "fission_fraction_statement": (
        "On their representative merger trajectory, fission contributes about "
        "10 percent of the radioactive energy at times up to about 1 day."
    ),
    "beta_energy_to_beta_particles": 0.20,
    "beta_energy_to_gamma": 0.45,
    "beta_energy_to_neutrinos": 0.35,
    "heating_rate_scale_erg_per_g_s": 1.0e11,
    "heating_rate_index_range": (1.1, 1.4),
    "heating_rate_statement": (
        "Barnes et al. quote a heating rate near 1e11 erg/s/g times t_day to "
        "the power minus 1.1 to 1.4, citing Metzger et al. 2010, Roberts et al. "
        "2011, Goriely et al. 2011, and Korobkin et al. 2012."
    ),
    "thermalization_fit_prefactor": 0.36,
    "thermalization_fit": (
        "f_tot(t) = 0.36 * [exp(-a*t) + ln(1 + 2*b*t^d) / (2*b*t^d)]. "
        "a, b, and d depend on the ejecta model and are Table 1 of that paper."
    ),
    "peak_luminosity_reduction_factor": 2.0,
    "late_time_days": 15.0,
    "late_time_statement": (
        "Thermalization reduces kilonova luminosity by a factor of about 2 at "
        "peak and by an order of magnitude or more at 15 days or later."
    ),
}

# Korobkin, Rosswog, Arcones, and Winteler, MNRAS 426, 1940 (2012).
# The abstract and the opening of that paper state the robustness result below.
# Barnes 2016 cites this paper for the heating-rate scale.
KOROBKIN_2012 = {
    "citation": (
        "Korobkin, Rosswog, Arcones, and Winteler, "
        "MNRAS 426, 1940 (2012)"
    ),
    "electron_fraction": 0.04,
    "binary_systems": 23,
    "abundance_statement": (
        "The 23 investigated binary systems yield practically identical "
        "abundance patterns, from the second r-process peak through the third, "
        "because the ejecta are extremely neutron rich."
    ),
    "ejecta_mass_variation_factor": 5.0,
}

# The pivot request cites Korobkin et al., MNRAS 448, 2544 (2015). That page
# was not opened here, so this record stores the citation and no yield number.
KOROBKIN_2015_CITATION = {
    "citation": "Korobkin et al., MNRAS 448, 2544 (2015)",
    "yields_stored": False,
    "note": (
        "No numerical fission-cycling yield from MNRAS 448, 2544 (2015) is "
        "stored. The comparison numbers in this module are the Barnes 2016 "
        "and Korobkin 2012 records above."
    ),
}

FREIBURGHAUS_1999 = {
    "citation": "Freiburghaus, Rosswog, and Thielemann, ApJ 525, L121 (1999)",
    "role": "r-process nucleosynthesis in neutron star merger ejecta",
}

GW170817 = {
    "citation": "Abbott et al., Physical Review Letters 119, 161101 (2017)",
    "event": "GW170817 / GRB 170817A",
    "kilonova": "AT2017gfo",
    "role": "Observational anchor that a neutron-star merger produced an r-process kilonova.",
}


def reference_records() -> dict:
    return {
        "barnes_2016": BARNES_2016,
        "korobkin_2012": KOROBKIN_2012,
        "korobkin_2015_citation": KOROBKIN_2015_CITATION,
        "freiburghaus_1999": FREIBURGHAUS_1999,
        "gw170817": GW170817,
    }


def score_fission_energy_fraction(predicted_fraction: float) -> dict:
    """Residual against the Barnes 2016 fission share of radioactive energy.

    ``predicted_fraction`` is supplied by the caller. This function does not
    read the SolCascade ledger and does not turn a domain scalar into a yield.
    """
    if not 0.0 <= predicted_fraction <= 1.0:
        raise ValueError("fission energy fraction must be in [0, 1]")
    benchmark = float(BARNES_2016["fission_energy_fraction_before_1_day"])
    return {
        "benchmark_fraction": benchmark,
        "predicted_fraction": predicted_fraction,
        "residual": predicted_fraction - benchmark,
        "citation": BARNES_2016["citation"],
        "statement": BARNES_2016["fission_fraction_statement"],
    }
