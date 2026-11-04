"""Unit conversions. Factors here are definitional, not measured constants."""

from __future__ import annotations

CENTIMETERS_PER_METER = 100.0
METERS_PER_CENTIMETER = 0.01
CUBIC_METERS_PER_CUBIC_CENTIMETER = 1.0e-6
EV_PER_KEV = 1.0e3
KEV_PER_EV = 1.0e-3
# One barn is 10^-28 m^2 by definition.
BARN_M2 = 1.0e-28


def meters_to_centimeters(length_m: float) -> float:
    return length_m * CENTIMETERS_PER_METER


def cubic_centimeters_to_cubic_meters(volume_cm3: float) -> float:
    return volume_cm3 * CUBIC_METERS_PER_CUBIC_CENTIMETER
