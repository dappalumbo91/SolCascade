"""Command line for the local SolCascade pipeline."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from solcascade import __version__
from solcascade.interfaces import NullProvider, apply_provider, collect_slots
from solcascade.load import load_problem
from solcascade.pipeline import run_pipeline


def provider_scalars(payload: object) -> dict[str, float]:
    """Keep numeric slots from a provider JSON object.

    A null value means the YAML value stands. Keys that start with ``_`` are
    comments. Any other non-number is an error. Unknown slot names are still
    rejected later, when the provider is bound to the problem.
    """
    if not isinstance(payload, dict):
        raise ValueError("provider file must be a JSON object")
    values: dict[str, float] = {}
    for key, raw in payload.items():
        name = str(key)
        if name.startswith("_"):
            continue
        if raw is None:
            continue
        if isinstance(raw, bool) or not isinstance(raw, (int, float)):
            raise ValueError(f"provider slot {name} must be a number or null")
        values[name] = float(raw)
    return values


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="solcascade",
        description="Staged magnetic compression, shell transport, and direct conversion.",
    )
    parser.add_argument("--version", action="version", version=f"solcascade {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Run the pipeline from a YAML config")
    run.add_argument("config", type=Path)
    run.add_argument("--workdir", type=Path, default=None)
    run.add_argument(
        "--provider",
        type=Path,
        default=None,
        help="JSON object of scalar overrides. Keys must be pipeline slots.",
    )

    keys = sub.add_parser("keys", help="List numeric slots a scalar provider may set")
    keys.add_argument("config", type=Path)

    args = parser.parse_args(argv)
    if args.command == "keys":
        config, materials = load_problem(args.config)
        slots, _provenance = apply_provider(collect_slots(config, materials), NullProvider())
        for key in slots:
            print(key)
        return 0
    overrides = None
    if args.provider is not None:
        payload = json.loads(args.provider.read_text(encoding="utf-8"))
        try:
            overrides = provider_scalars(payload)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            return 2
    summary = run_pipeline(args.config, workdir=args.workdir, provider_overrides=overrides)
    conversion = summary["conversion"]
    print(f"workdir: {summary['workdir']}")
    print(f"transport: {summary['transport']['backend']}")
    print(f"fusion energy J: {conversion['fusion_J']:.6g}")
    print(f"net electrical J: {conversion['net_electrical_J']:.6g}")
    print(f"peak triple-product ratio: {summary['plasma']['peak_triple_product_ratio']:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
