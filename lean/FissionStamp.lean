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
stamp. The catalog ceiling, the 6.5 percent page, and the 60 GWd/t stack stay
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
  "((((e^-3+G^-2)*(1-alpha*gamma^3))*(1+alpha*e^-4))*(1-alpha*pi^-5))*(1-alpha*pi^-6)"
]

theorem installed_closing_count : installedClosings.length = 68 := rfl

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
