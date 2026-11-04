"""Run the checked-in minimal pipeline and print the work directory."""

from pathlib import Path

from solcascade.pipeline import run_pipeline

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    summary = run_pipeline(ROOT / "configs" / "minimal_pipeline.yaml")
    print(summary["workdir"])
    print(summary["conversion"]["net_electrical_J"])


if __name__ == "__main__":
    main()
