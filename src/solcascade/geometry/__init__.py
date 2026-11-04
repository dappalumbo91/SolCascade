"""Field geometry: axisymmetric coils and the spherical shell."""

from solcascade.geometry.coils import Coil, CoilSet, loop_field
from solcascade.geometry.shell import Shell, ShellCell, ShellLayer

__all__ = [
    "Coil",
    "CoilSet",
    "Shell",
    "ShellCell",
    "ShellLayer",
    "loop_field",
]
