#!/usr/bin/env python3
"""Prepare or execute the frozen analysis in a new private work directory.

Author: mjwang
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


STAGES = [
    ("20260809_152651/01_reanalysis.R", []),
    ("20260809_152651/03_qc_exploratory.R", []),
    ("20261003_165814/04_correct_mi_pooling.R", []),
    ("20261018_173329/center_comparison.R", []),
    ("20261023_001156/01_steatosis_components.R", ["ROOT"]),
    ("20261016_225206/01_pattern_vs_rest.R", []),
    ("20261027_230711/01_steatosis_probabilities.R", ["--project-root", "ROOT"]),
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-workbook", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--execute", action="store_true", help="Run models; otherwise print the plan only")
    parser.add_argument("--check-environment", action="store_true", help="Report dependencies without running models")
    args = parser.parse_args()
    if args.check_environment:
        expression = ('p<-c("readxl","digest","MASS","ordinal","mice","emmeans","jsonlite");'
                      'missing<-character();for(x in p){ok<-requireNamespace(x,quietly=TRUE);'
                      'cat(x,if(ok)as.character(packageVersion(x))else "MISSING",sep=",",fill=TRUE);'
                      'if(!ok)missing<-c(missing,x)};if(length(missing))quit(status=1)')
        return subprocess.run([args.rscript, "--vanilla", "-e", expression], check=False).returncode
    if args.work_dir is None:
        parser.error("--work-dir is required")
    root = args.work_dir.expanduser().resolve()
    repository = Path(__file__).resolve().parent
    if root == repository or repository in root.parents:
        parser.error("Use a private work directory outside the code repository")
    if root.exists():
        parser.error("Work directory already exists; use a new path to preserve previous runs")
    for relative, trailing in STAGES:
        script = root / "paper1_cmrf_pathology" / "code" / "R" / "versions" / relative
        command = [args.rscript, "--vanilla", str(script), *[str(root) if x == "ROOT" else x for x in trailing]]
        print(" ".join(repr(part) for part in command), flush=True)
    if not args.execute:
        return 0
    if args.input_workbook is None or not args.input_workbook.is_file():
        parser.error("--execute requires an authorised --input-workbook file")
    # This step creates a separate workspace; it does not modify the source workbook.
    root.mkdir(parents=True, exist_ok=False)
    shutil.copytree(repository / "source_scripts" / "paper1_cmrf_pathology", root / "paper1_cmrf_pathology")
    (root / "data").mkdir()
    shutil.copy2(args.input_workbook, root / "data" / "All_patients_SAFE.xlsx")
    results = root / "paper1_cmrf_pathology" / "results" / "versions"
    for version in {relative.split("/", 1)[0] for relative, _ in STAGES}:
        for child in ["source_data", "qa", "inputs"]:
            (results / version / child).mkdir(parents=True, exist_ok=True)
    for relative, trailing in STAGES:
        script = root / "paper1_cmrf_pathology" / "code" / "R" / "versions" / relative
        command = [args.rscript, "--vanilla", str(script), *[str(root) if x == "ROOT" else x for x in trailing]]
        subprocess.run(command, check=True, cwd=root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
