#!/bin/bash
#SBATCH --job-name=job_for_barcode31
#SBATCH --time=0-01:00:00
#SBATCH --mem=32GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=clara
#SBATCH --output=/work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/output/batch_30_%j.out

# Environment activation
module purge
module load Anaconda3
eval "$(conda shell.bash hook)"
source ~/.bashrc
conda activate APD

python single_folder_analyzer5.py --sample_dir "/work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/barcode31" --report_path "/work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/report/barcode31_sim_report.xlsx" --blast_db "../DB/ref_viruses_rep_genomes ../DB/ref_prok_rep_genomes ../DB/ref_euk_rep_genomes" --platform nanopore --strictness Fast --file_ext "fastq.gz"
python excel_pathogen_filter.py --report_path "/work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/report/barcode31_sim_report.xlsx"
python excel_pathogen_filter_second_round.py --report_path "/work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/report/barcode31_sim_report.xlsx"
echo 'Finished processing /work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/barcode31'
python -c "import pandas as pd; df = pd.read_excel('/work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/output/analysis_progress.xlsx'); df.loc[df['Folder_Path'] == '/work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/barcode31', 'Progress'] = 'finished'; df.to_excel('/work/ia29ipap-PRJs/PRJS/APD/Data/PREPARETiD_DOMESTICANIMALS_SEPTEMEBER2025_UGANDA/output/analysis_progress.xlsx', index=False)"
