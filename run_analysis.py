#!/usr/bin/env python3
"""Run the CMRF histology analysis in a separate private workspace.

Author: mjwang
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys

# Public step names map to the original internal layout at runtime. This keeps
# the tested scripts unchanged while presenting one ordered analysis pipeline.
SCRIPT_MAP = [('analysis/01_prepare_and_model.R', 'paper1_cmrf_pathology/code/R/versions/20260809_152651/01_reanalysis.R', []), ('analysis/02_diagnostics_and_patterns.R', 'paper1_cmrf_pathology/code/R/versions/20260809_152651/03_qc_exploratory.R', []), ('analysis/03_pool_imputed_models.R', 'paper1_cmrf_pathology/code/R/versions/20261003_165814/04_correct_mi_pooling.R', []), ('analysis/04_compare_centres.R', 'paper1_cmrf_pathology/code/R/versions/20261018_173329/center_comparison.R', []), ('analysis/05_steatosis_components.R', 'paper1_cmrf_pathology/code/R/versions/20261023_001156/01_steatosis_components.R', ['ROOT']), ('analysis/06_common_profiles.R', 'paper1_cmrf_pathology/code/R/versions/20261016_225206/01_pattern_vs_rest.R', []), ('analysis/07_steatosis_probabilities.R', 'paper1_cmrf_pathology/code/R/versions/20261027_230711/01_steatosis_probabilities.R', ['--project-root', 'ROOT']), ('figures/01_primary_figures.R', 'paper1_cmrf_pathology/code/R/versions/20261003_165814/02_figures.R', None), ('figures/02_profile_comparison.R', 'paper1_cmrf_pathology/code/R/versions/20261016_225206/02_pattern_comparison_figure.R', None), ('figures/03_participant_flow.R', 'paper1_cmrf_pathology/code/R/versions/20261021_131343/figure1.R', None), ('figures/04_component_figures.R', 'paper1_cmrf_pathology/code/R/versions/20261022_003505/figures34.R', None), ('figures/05_steatosis_components.R', 'paper1_cmrf_pathology/code/R/versions/20261023_001156/02_figure_steatosis_components.R', None), ('figures/06_complete_steatosis_thresholds.R', 'paper1_cmrf_pathology/code/R/versions/20261024_003114/01_relabel_supplementary_figure.R', None), ('figures/07_fibrosis_count.R', 'paper1_cmrf_pathology/code/R/versions/20261024_094342/figure2_white_labels.R', None), ('figures/08_supplementary_tables.R', 'paper1_cmrf_pathology/code/R/versions/20261025_124331/revise_and_export.R', None), ('figures/09_steatosis_count.R', 'paper1_cmrf_pathology/code/R/versions/20261027_230711/02_figure4_matching_figure2.R', None)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-workbook', type=Path)
    parser.add_argument('--work-dir', type=Path)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--execute', action='store_true', help='Execute the seven statistical steps; default prints the plan')
    parser.add_argument('--check-environment', action='store_true', help='Report dependency versions without running models')
    args = parser.parse_args()
    if args.check_environment:
        expression = ('p<-c("readxl","digest","MASS","ordinal","mice","emmeans","jsonlite");'
                      'missing<-character();for(x in p){ok<-requireNamespace(x,quietly=TRUE);'
                      'cat(x,if(ok)as.character(packageVersion(x))else "MISSING",sep=",",fill=TRUE);'
                      'if(!ok)missing<-c(missing,x)};if(length(missing))quit(status=1)')
        return subprocess.run([args.rscript, '--vanilla', '-e', expression], check=False).returncode
    if args.work_dir is None:
        parser.error('--work-dir is required')
    workspace = args.work_dir.expanduser().resolve()
    repository = Path(__file__).resolve().parent
    if workspace == repository or repository in workspace.parents:
        parser.error('Choose a private work directory outside this repository')
    if workspace.exists():
        parser.error('Work directory already exists; choose a new path')
    analysis_steps = [item for item in SCRIPT_MAP if item[2] is not None]
    for index, (public_path, _, _) in enumerate(analysis_steps, 1):
        print(f'{index:02d}. {public_path}', flush=True)
    if not args.execute:
        return 0
    if args.input_workbook is None or not args.input_workbook.is_file():
        parser.error('--execute requires an authorised --input-workbook')
    workspace.mkdir(parents=True, exist_ok=False)
    for public_path, runtime_path, _ in SCRIPT_MAP:
        target = workspace / runtime_path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repository / public_path, target)
    (workspace / 'data').mkdir()
    shutil.copy2(args.input_workbook, workspace / 'data' / 'All_patients_SAFE.xlsx')
    for _, runtime_path, arguments in analysis_steps:
        version = Path(runtime_path).parts[-2]
        for child in ['source_data', 'qa', 'inputs']:
            (workspace / 'paper1_cmrf_pathology' / 'results' / 'versions' / version / child).mkdir(parents=True, exist_ok=True)
    for public_path, runtime_path, arguments in analysis_steps:
        print(f'Running {public_path}', flush=True)
        command = [args.rscript, '--vanilla', str(workspace / runtime_path),
                   *[str(workspace) if value == 'ROOT' else value for value in arguments]]
        subprocess.run(command, check=True, cwd=workspace)
    return 0


if __name__ == '__main__':
    sys.exit(main())
