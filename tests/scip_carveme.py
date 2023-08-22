#!/usr/bin/env python

#SBATCH --array=0-5
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1            
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000
#SBATCH --job-name="scip"


from carveme import project_dir
from subprocess import call
import os

data_path = project_dir + '/data/benchmark'
index = int(os.environ['SLURM_ARRAY_TASK_ID'])
org_id = ['bsub', 'ecol', 'mgen', 'paer', 'rsol', 'sone'][index]

organisms = {
    'bsub': 'Bacillus subtilis (168)',
    'ecol': 'Escherichia coli (K-12 MG1655)',
    'mgen': 'Mycoplasma genitalium (G-37)',
    'paer': 'Pseudomonas aeruginosa (PA01)',
    'rsol': 'Ralstonia solenacearum (GMI1000)',
    'sone': 'Shewanella oneidensis (MR-1)'
}

genomes = {
    'bsub': 'Bsubtilis_168.faa',
    'ecol': 'Ecoli_K12_MG1655.faa',
    'mgen': 'M_genitalium_G37.faa',
    'paer': 'Paeruginosa_PAO1.faa',
    'rsol': 'Rsolanacearum_GMI1000.faa',
    'sone': 'Soneidensis_MR1.faa'
}

gram_status = {
    'bsub': 'grampos',
    'ecol': 'gramneg',
    'mgen': 'gramneg',
    'paer': 'gramneg',
    'rsol': 'gramneg',
    'sone': 'gramneg',
}


biolog_media = {
    'bsub': 'M9',
    'ecol': 'M9',
    'paer': 'M9',
    'rsol': 'M9',
    'sone': 'ShewMM'
}

essentiality_media = {
    'bsub': 'LB',
    'ecol': 'M9',
    'mgen': None,
    'paer': 'M9[succ]',
    'sone': 'LB'
}

print(f'Carving model for {organisms[org_id]}')

fasta_file = f"{data_path}/fasta/{genomes[org_id]}"
model_file = f"{data_path}/models/{org_id}.xml"
mediadb = f"{data_path}/media_db.tsv"

media = set()
if org_id in biolog_media and biolog_media[org_id]:
    media.add(biolog_media[org_id])
if org_id in essentiality_media and essentiality_media[org_id]:
    media.add(essentiality_media[org_id])
media = ','.join(media)

gapfill = f'-g "{media}" --mediadb {mediadb}' if media else ''

call(f'carve {fasta_file} -u {gram_status[org_id]} -o {model_file} {gapfill} --fbc2 --solver scip', shell=True)

