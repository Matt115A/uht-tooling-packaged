import argparse
import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import yaml

I7_INDEXES: Dict[str, str] = {
    "701": "TCGCCTTA",
    "702": "CTAGTACG",
    "703": "TTCTGCCT",
    "704": "GCTCAGGA",
    "705": "AGGAGTCC",
    "706": "CATGCCTA",
    "707": "GTAGAGAG",
    "708": "CCTCTCTG",
    "709": "AGCGTAGC",
    "710": "CAGCCTCG",
    "711": "TGCCTCTT",
    "712": "TCCTCTAC",
}

I5_INDEXES: Dict[str, str] = {
    "501": "TAGATCGC",
    "502": "CTCTCTAT",
    "503": "TATCCTCT",
    "504": "AGAGTAGA",
    "505": "GTAAGGAG",
    "506": "ACTGCATA",
    "507": "AAGGAGTA",
    "508": "CTAAGCCT",
    "510": "CGTCTAAT",
    "511": "TCTCTCCG",
    "513": "TCGACTAG",
    "515": "TTCTAGCT",
}

I7_PREFIX = "CAAGCAGAAGACGGCATACGAGAT"
I7_SUFFIX = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
I5_PREFIX = "AATGATACGGCGACCACCGAGATCTACAC"
I5_SUFFIX = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"


def load_binding_sequences(csv_path: Path) -> Tuple[str, str]:
    df = pd.read_csv(csv_path)
    if "binding_region" not in df.columns:
        raise ValueError("CSV must contain a 'binding_region' column")
    if len(df) < 2:
        raise ValueError("CSV must contain at least two rows for i7 and i5 binding regions")
    return df["binding_region"].iloc[0], df["binding_region"].iloc[1]


def load_config(config_path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    if not config_path:
        return {}
    with open(config_path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    return config


def generate_primers(
    template_binding_i7: str,
    template_binding_i5: str,
    i7_indexes: Optional[Dict[str, str]] = None,
    i5_indexes: Optional[Dict[str, str]] = None,
    i7_prefix: str = I7_PREFIX,
    i7_suffix: str = I7_SUFFIX,
    i5_prefix: str = I5_PREFIX,
    i5_suffix: str = I5_SUFFIX,
) -> List[Tuple[str, str]]:
    primers: List[Tuple[str, str]] = []
    i7_map = i7_indexes or I7_INDEXES
    i5_map = i5_indexes or I5_INDEXES

    for idx, seq in i7_map.items():
        name = f"i7_{idx}"
        full_seq = f"{i7_prefix}{seq}{i7_suffix}{template_binding_i7}"
        primers.append((name, full_seq))

    for idx, seq in i5_map.items():
        name = f"i5_{idx}"
        full_seq = f"{i5_prefix}{seq}{i5_suffix}{template_binding_i5}"
        primers.append((name, full_seq))

    return primers


def run_nextera_primer_design(
    binding_csv: Path,
    output_csv: Path,
    log_path: Optional[Path] = None,
    config_path: Optional[Path] = None,
    logger: Optional[logging.Logger] = None,
) -> Path:
    binding_csv = Path(binding_csv)
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    managed_logger = logger is None
    if logger is None:
        logger = logging.getLogger("uht_tooling.nextera")
        logger.setLevel(logging.INFO)
        handler: logging.Handler
        if log_path:
            handler = logging.FileHandler(log_path, mode="w")
        else:
            handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
        logger.handlers = []
        logger.addHandler(handler)
        logger.propagate = False

    try:
        logger.info("Loading binding sequences from %s", binding_csv)
        template_i7, template_i5 = load_binding_sequences(binding_csv)
        logger.info("Loaded binding regions (i7 len=%s, i5 len=%s)", len(template_i7), len(template_i5))

        config = load_config(config_path)
        if config:
            logger.info("Loaded configuration overrides from %s", config_path)

        i7_indexes = config.get("i7_indexes") if config else None
        i5_indexes = config.get("i5_indexes") if config else None
        i7_prefix = config.get("i7_prefix", I7_PREFIX) if config else I7_PREFIX
        i7_suffix = config.get("i7_suffix", I7_SUFFIX) if config else I7_SUFFIX
        i5_prefix = config.get("i5_prefix", I5_PREFIX) if config else I5_PREFIX
        i5_suffix = config.get("i5_suffix", I5_SUFFIX) if config else I5_SUFFIX

        primers = generate_primers(
            template_binding_i7=template_i7,
            template_binding_i5=template_i5,
            i7_indexes=i7_indexes,
            i5_indexes=i5_indexes,
            i7_prefix=i7_prefix,
            i7_suffix=i7_suffix,
            i5_prefix=i5_prefix,
            i5_suffix=i5_suffix,
        )
        logger.info("Generated %s primers", len(primers))

        with output_csv.open("w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["primer_name", "sequence"])
            writer.writerows(primers)

        logger.info("Wrote primers to %s", output_csv)
        return output_csv
    finally:
        if managed_logger and logger:
            for handler in list(logger.handlers):
                handler.close()
                logger.removeHandler(handler)
            logger.propagate = True


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate Nextera XT primers from binding region CSV input.")
    parser.add_argument(
        "--binding-csv",
        required=True,
        type=Path,
        help="CSV file with a 'binding_region' column; first row is i7, second row is i5.",
    )
    parser.add_argument(
        "--output-csv",
        required=True,
        type=Path,
        help="Path to write the generated primer CSV.",
    )
    parser.add_argument(
        "--log-path",
        default=None,
        type=Path,
        help="Optional path to write a log file.",
    )
    parser.add_argument(
        "--config",
        default=None,
        type=Path,
        help="Optional YAML file providing overrides for indexes/prefixes/suffixes.",
    )
    return parser


def main(argv: Optional[List[str]] = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    run_nextera_primer_design(
        binding_csv=args.binding_csv,
        output_csv=args.output_csv,
        log_path=args.log_path,
        config_path=args.config,
    )


if __name__ == "__main__":
    main()
