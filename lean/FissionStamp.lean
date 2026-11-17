/-
Copyright 2026 Damian Arthur Palumbo
SPDX-License-Identifier: Apache-2.0
-/

import Mathlib.Tactic.NormNum

/-
Fission stamp. Not part of `lean_lib FSOT`.

The rational plant figures are proved here. Each installed closing, each held
pair, and the four-hour year share are Octave obligations. The Octave log
recomputes the expressions from the pinned seeds. A page formula is not this
stamp. The tails percent and the tails fraction are installed closings. The
steam-share quotients and the electricity gains are rationals of the named
steps. The operating feed quotients are installed closings. The comparison
column's feed multiple is a rational of the printed masses and stays unrun.
The file ledger quantum is an installed closing. The operating-column
remainder and the cooling-page extra heat are rationals of the named scores.
The 55 GWd/t step and the 70 GWd/t step stay unrun.
The catalog ceiling, the 6.5 percent page, and the 60 GWd/t stack stay
unrun. Domain scalars stay out of the joules.
-/

namespace FissionStamp

def demonstrated : ℚ := 33 / 100

def builtEff : ℚ := 9 / 25

def carnot : ℚ := 1 - 310 / 560

def fabricationCap : ℚ := 5 / 100

def catalogCeiling : ℚ := 45 / 100

def fastHeadroom : ℚ := (1 / 2) / (2 / 5)

theorem carnot_bound : carnot = 25 / 56 := by
  unfold carnot
  norm_num

theorem fabrication_cap : fabricationCap = 1 / 20 := by
  unfold fabricationCap
  norm_num

theorem catalog_ceiling : catalogCeiling = 9 / 20 := by
  unfold catalogCeiling
  norm_num

theorem fast_headroom : fastHeadroom = 5 / 4 := by
  unfold fastHeadroom
  norm_num

theorem shared_headroom : catalogCeiling / builtEff = fastHeadroom := by
  unfold catalogCeiling builtEff fastHeadroom
  norm_num

theorem built_step : builtEff / demonstrated = 12 / 11 := by
  unfold builtEff demonstrated
  norm_num

theorem efficiency_step : ((2 : ℚ) / 5) / demonstrated = 40 / 33 := by
  unfold demonstrated
  norm_num

theorem steam_order : demonstrated < builtEff ∧ builtEff < carnot ∧ carnot < catalogCeiling := by
  unfold demonstrated builtEff carnot catalogCeiling
  norm_num

theorem built_difference :
    builtEff - demonstrated = demonstrated * (builtEff / demonstrated - 1) := by
  unfold builtEff demonstrated
  norm_num

theorem carnot_difference : carnot - builtEff = builtEff * (carnot / builtEff - 1) := by
  unfold carnot builtEff
  norm_num

theorem rejected_fraction : (310 : ℚ) / 560 = 31 / 56 := by norm_num

theorem kept_built_share : builtEff / carnot = 504 / 625 := by
  unfold builtEff carnot
  norm_num

theorem kept_demonstrated_share : demonstrated / carnot = 462 / 625 := by
  unfold demonstrated carnot
  norm_num

theorem closed_share :
    (builtEff - demonstrated) / (carnot - demonstrated) = 42 / 163 := by
  unfold builtEff demonstrated carnot
  norm_num

theorem remaining_share :
    1 - (builtEff - demonstrated) / (carnot - demonstrated) = 121 / 163 := by
  unfold builtEff demonstrated carnot
  norm_num

theorem unused_built_share : 1 - builtEff / carnot = 121 / 625 := by
  unfold builtEff carnot
  norm_num

theorem unused_demonstrated_share : 1 - demonstrated / carnot = 163 / 625 := by
  unfold demonstrated carnot
  norm_num

theorem published_gain : builtEff / demonstrated - 1 = 1 / 11 := by
  unfold builtEff demonstrated
  norm_num

theorem remaining_gain : carnot / builtEff - 1 = 121 / 504 := by
  unfold carnot builtEff
  norm_num

theorem gain_ratio :
    (carnot / builtEff - 1) / (builtEff / demonstrated - 1) = 1331 / 504 := by
  unfold carnot builtEff demonstrated
  norm_num

def burnupStep : ℚ := 60 / 45

theorem burnup_step : burnupStep = 4 / 3 := by
  unfold burnupStep
  norm_num

theorem burnup_gain : burnupStep - 1 = 1 / 3 := by
  unfold burnupStep
  norm_num

theorem burnup_ratio : (burnupStep - 1) / (builtEff / demonstrated - 1) = 11 / 3 := by
  unfold burnupStep builtEff demonstrated
  norm_num

theorem gain_order :
    builtEff / demonstrated - 1 < carnot / builtEff - 1
      ∧ carnot / builtEff - 1 < burnupStep - 1 := by
  unfold builtEff demonstrated carnot burnupStep
  norm_num

theorem feed_multiple : (211 : ℚ) * 10 / 243 = 2110 / 243 := by norm_num

theorem feed_per_twh : (211 : ℚ) * 1000000 / (1000 * 8760) = 5275 / 219 := by norm_num

theorem feed_concentrate : (211 : ℚ) / 249 = 211 / 249 := by norm_num

theorem feed_hex : (211 : ℚ) / 312 = 211 / 312 := by norm_num

theorem feed_enriched_hex : (243 : ℚ) / 359 = 243 / 359 := by norm_num

theorem feed_oxide : (243 : ℚ) / 276 = 81 / 92 := by norm_num

theorem feed_core : (243 : ℚ) / 730 = 243 / 730 := by norm_num

theorem feed_transuranic : (280 : ℚ) / 27600 = 7 / 690 := by norm_num

theorem feed_compare : (163 : ℚ) * 20 / 337 = 3260 / 337 := by norm_num

theorem feed_compare_order : (3260 : ℚ) / 337 > 2110 / 243 := by norm_num

theorem column_swu_year : (182000 : ℚ) * 1000 / 8760 = 4550000 / 219 := by norm_num

theorem column_product : (243 : ℚ) / 2110 * (2110 / 243) = 1 := by norm_num

theorem column_tails : (1 : ℚ) - 243 / 2110 = 1867 / 2110 := by norm_num

theorem column_swu_natural : (1820 : ℚ) / 243 * 1000 / (2110 / 243) = 182000 / 211 := by norm_num

theorem column_enriched : (5275 : ℚ) / 219 / (2110 / 243) = 405 / 146 := by norm_num

theorem column_fissile :
    ((45 : ℚ) / 1000 - (22 : ℚ) / 10000) / (2110 / 243) = 26001 / 5275000 := by
  norm_num

theorem column_cap_step : (55 : ℚ) / 45 = 11 / 9 := by norm_num

theorem column_cap_keeps_built : (60 : ℚ) / 55 = 12 / 11 := by norm_num

theorem column_compare_product : (337 : ℚ) / 3260 * (3260 / 337) = 1 := by norm_num

theorem column_sighted_step : (70 : ℚ) / 45 = 14 / 9 := by norm_num

theorem reject_ratio :
    (12 : ℚ) / 11 * (((1 : ℚ) - 33 / 100) / ((1 : ℚ) - 9 / 25)) = 201 / 176 := by
  norm_num

theorem reject_extra : (201 : ℚ) / 176 - 1 = 25 / 176 := by norm_num

theorem reject_over_gain : (25 : ℚ) / 176 / (1 / 11) = 25 / 16 := by norm_num

/-- Octave discharges one closing string per installed unique match. -/
def installedClosings : List String := [
  "(phi^-6/sqrt(G))*(1-alpha*e^-4)",
  "(1-(gamma^6+sqrt(G)))*(1-alpha*pi^-5)",
  "(1-(gamma^3/sqrt(G)))*(1+alpha*pi^-5)",
  "(((pi^-5+phi^-6)*(1-alpha*phi^-3))*(1-alpha*pi^-4))*(1-alpha*pi^-6)",
  "((pi^-4*G^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)",
  "(((e^-3*phi^-4)*(1+alpha*phi^-2))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
  "((pi^-5*phi)*(1-alpha*e^-4))*(1-alpha*pi^-6)",
  "((gamma^3*G^6)*(1+alpha*gamma^5))*(1+alpha*e^-5)",
  "((pi-gamma^4)*(1-alpha*pi^-4))*(1-alpha*e^-6)",
  "((pi^-3*G^-5)*(1-alpha*phi^-6))*(1-alpha*e^-6)",
  "(((pi^-2*gamma^6)*(1-alpha*gamma^2))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
  "(((gamma^-5-gamma^6)*(1+alpha*pi^-3))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
  "((((G^3+G^5)*(1+alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
  "(((sqrt(e)*e^-5)*(1-alpha*sqrt(gamma)))*(1+alpha*e^-4))*(1-alpha*e^-6)",
  "((pi^-4-pi^-5)*(1-alpha*pi^-4))*(1-alpha*pi^-6)",
  "(((phi^-2*gamma^6)*(1+alpha*G^-2))*(1-alpha*gamma^6))*(1+alpha*e^-6)",
  "(((e^-3-e^-6)*(1+alpha*phi^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
  "(sqrt(gamma)*G^-3)*(1+alpha*e^-6)",
  "((pi^-2+e^2)*(1-alpha*pi^-4))*(1-alpha*e^-6)",
  "((phi^5+gamma^-5)*(1+alpha*pi^-3))*(1+alpha*pi^-6)",
  "(((pi^6-gamma^-4)*(1-alpha*pi^-2))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
  "(((pi^-3*phi^5)*(1+alpha*G))*(1-alpha*e^-4))*(1-alpha*e^-5)",
  "(((phi^-2+gamma^5)*(1+alpha*gamma^4))*(1+alpha*e^-5))*(1+alpha*pi^-6)",
  "(((e^-3+G^-3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*e^-6)",
  "(((e^-3*gamma^-2)*(1+alpha*gamma))*(1-alpha*phi^-6))*(1+alpha*pi^-6)",
  "((phi^-5+sqrt(gamma))*(1+alpha*pi^-4))*(1+alpha*pi^-5)",
  "(((pi^-4*G^5)*(1-alpha*phi^-2))*(1-alpha*pi^-4))*(1-alpha*e^-6)",
  "((pi^-5*G^4)*(1-alpha*pi^-4))*(1-alpha*e^-6)",
  "(((phi^-5*gamma^-2)*(1-alpha*gamma^2))*(1+alpha*pi^-4))*(1+alpha*pi^-6)",
  "((sqrt(gamma)*G^6)*(1+alpha*phi^-2))*(1+alpha*e^-4)",
  "(sqrt(phi)+phi^-5)*(1+alpha*phi^-4)",
  "(((phi^-3*G^-6)*(1+alpha*phi^-5))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
  "(((sqrt(pi)-sqrt(phi))*(1-alpha*gamma^4))*(1-alpha*e^-5))*(1-alpha*pi^-6)",
  "((gamma^2+G)*(1+alpha*phi^-5))*(1+alpha*pi^-5)",
  "((phi^-3+gamma^5)*(1-alpha*gamma^5))*(1-alpha*pi^-6)",
  "(((e-phi^2)*(1-alpha*gamma^2))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
  "(((pi^6/sqrt(phi))*(1-alpha*gamma))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
  "(((pi^6-e^4)*(1-alpha*G^-5))*(1+alpha*gamma^5))*(1-alpha*pi^-5)",
  "(((e^-4+G^-2)*(1+alpha*phi^-3))*(1-alpha*e^-4))*(1-alpha*e^-6)",
  "((pi^-6+gamma^5)*(1-alpha*phi^-3))*(1-alpha*e^-5)",
  "(pi^-3+G^-3)*(1-alpha*e^-4)",
  "((gamma^4+G^-6)*(1-alpha*gamma^6))*(1-alpha*e^-6)",
  "((((e^-3*sqrt(phi))*(1-alpha*phi^-2))*(1-alpha*e^-4))*(1-alpha*e^-5))*(1+alpha*pi^-6)",
  "((((gamma^4*G^3)*(1+alpha*gamma^3))*(1-alpha*pi^-4))*(1-alpha*pi^-6))*(1-alpha*pi^-6)",
  "((gamma^-3-G^4)*(1+alpha*e^-2))*(1-alpha*pi^-4)",
  "((pi^-2+G^3)*(1+alpha*pi^-3))*(1-alpha*e^-6)",
  "((((e^-2+gamma^5)*(1+alpha*phi^-2))*(1+alpha*e^-4))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
  "((e^-3+gamma^6)*(1+alpha*gamma^3))*(1+alpha*pi^-5)",
  "(((pi^-4*sqrt(G))*(1+alpha*G^4))*(1-alpha*e^-4))*(1+alpha*e^-6)",
  "((pi^-4*G^3)*(1+alpha*e^-4))*(1+alpha*e^-6)",
  "(e+phi^-3)*(1+alpha*e^-3)",
  "((G+sqrt(G))*(1-alpha*phi^-5))*(1+alpha*pi^-6)",
  "(((phi^-5+gamma)*(1+alpha*phi^-2))*(1+alpha*gamma^6))*(1+alpha*pi^-6)",
  "(((e*phi^-6)*(1+alpha*pi^-3))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
  "((e^5*G^3)*(1-alpha*pi^-3))*(1+alpha*pi^-6)",
  "(gamma^-2+gamma^-2)*(1-alpha*gamma^5)",
  "((gamma+gamma^2)*(1-alpha*gamma^3))*(1-alpha*pi^-5)",
  "(((phi^-6*gamma^5)*(1+alpha*e^-4))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
  "(((e^4*G)*(1-alpha*pi^-3))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
  "((((pi^-5*G^2)*(1-alpha*phi^-5))*(1-alpha*e^-5))*(1+alpha*pi^-6))*(1+alpha*pi^-6)",
  "(((e^-6*G)*(1+alpha*sqrt(gamma)))*(1+alpha*pi^-5))*(1+alpha*pi^-6)",
  "((((pi^6*pi^6)*(1+alpha*e^5))*(1+alpha*gamma^-3))*(1+alpha*gamma^4))*(1+alpha*pi^-5)",
  "((phi-G^5)*(1+alpha*e^-5))*(1-alpha*pi^-6)",
  "(((gamma^6+G^4)*(1-alpha*e^-3))*(1+alpha*e^-6))*(1-alpha*pi^-6)",
  "(((phi^-4+G^2)*(1-alpha*phi^-3))*(1+alpha*e^-5))*(1-alpha*pi^-6)",
  "(((sqrt(pi)-G)*(1-alpha*phi^-2))*(1-alpha*gamma^5))*(1+alpha*pi^-6)",
  "((G*G^-2)*(1-alpha*pi^-2))*(1-alpha*pi^-5)",
  "((((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
  "((phi^-2*gamma)*(1-alpha*gamma^2))*(1+alpha*gamma^6)",
  "(((pi^-5-pi^-6)*(1-alpha*G^-6))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
  "(((pi^2-G^-2)*(1+alpha*phi^-5))*(1-alpha*pi^-5))*(1-alpha*pi^-6)",
  "((pi^3-phi^4)*(1-alpha*phi^-2))*(1+alpha*pi^-4)",
  "((e^-5+G^2)*(1+alpha*phi^-3))*(1+alpha*pi^-3)",
  "((pi^-3+G^5)*(1-alpha*phi^-4))*(1-alpha*pi^-6)",
  "((pi^-3+G^5)*(1-alpha*pi^-3))*(1+alpha*e^-5)",
  "((sqrt(e)-G^3)*(1+alpha*pi^-3))*(1-alpha*pi^-6)",
  "(((gamma*gamma)*(1-alpha*e^-2))*(1+alpha*pi^-4))*(1+alpha*pi^-6)",
  "((pi*pi^-5)*(1-alpha*phi))*(1+alpha*e^-6)",
  "((e^2*gamma^-6)*(1+alpha*phi^-4))*(1+alpha*e^-6)"
]

theorem installed_closing_count : installedClosings.length = 79 := rfl

/-- Both routes stay. The stamp does not choose. -/
def heldRoutes : List (String × String × String) := [
  ("(sqrt(pi)-e^-6)*(1-alpha*e^-4)", "((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*e^-6)", "((sqrt(pi)-e^-6)*(1-alpha*e^-4))*(1-alpha*pi^-5)"),
  ("((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4)", "(((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1+alpha*e^-6)", "(((e^-6/G)*(1-alpha*gamma^2))*(1+alpha*e^-4))*(1+alpha*pi^-5)"),
  ("((e^-2*gamma^2)*(1-alpha*phi^-3))*(1-alpha*gamma^6)", "(((e^-2*gamma^2)*(1-alpha*phi^-3))*(1-alpha*gamma^6))*(1-alpha*e^-6)", "(((e^-2*gamma^2)*(1-alpha*phi^-3))*(1-alpha*gamma^6))*(1-alpha*pi^-5)"),
  ("((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3)", "(((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*e^-6)", "(((gamma*sqrt(G))*(1+alpha*phi^-3))*(1+alpha*e^-3))*(1-alpha*pi^-5)")
]

theorem held_route_count : heldRoutes.length = 4 := rfl

/-- Octave discharges this chain. The fourth gap stays wider than the helium match. -/
def yearShareMiss : List String := [
  "(pi^-5*e^-2)*(1+alpha*phi^3)",
  "((pi^-5*e^-2)*(1+alpha*phi^3))*(1+alpha*gamma^3)",
  "(((pi^-5*e^-2)*(1+alpha*phi^3))*(1+alpha*gamma^3))*(1+alpha*e^-4)",
  "((((pi^-5*e^-2)*(1+alpha*phi^3))*(1+alpha*gamma^3))*(1+alpha*e^-4))*(1+alpha*pi^-6)"
]

theorem year_share_step_count : yearShareMiss.length = 4 := rfl

/-
The catalog ceiling, the atomic percent, the burnup ratio, the stack, the higher
burnup, and the stack electric are installed closings that stay unrun.
-/

end FissionStamp
