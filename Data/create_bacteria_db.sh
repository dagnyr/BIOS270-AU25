#!/bin/bash
#SBATCH --job-name=createdb
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=04:00:00

cd /home/users/dagnyr/bios270/BIOS270-AU25/Data

module load apptainer

rm -f bacteria.db
DATABASE="bacteria.db"
CONTAINER=/home/users/dagnyr/bioinformatics_latest.sif
CLASSDATA=/farmshare/home/classes/bios/270

apptainer exec -B $CLASSDATA:$CLASSDATA $CONTAINER \
    python insert_gff_table.py --database_path "$DATABASE"

apptainer exec -B $CLASSDATA:$CLASSDATA $CONTAINER \
    python insert_protein_cluster_table.py --database_path "$DATABASE"

apptainer exec -B $CLASSDATA:$CLASSDATA $CONTAINER \
    python insert_metadata_table.py --database_path "$DATABASE"

