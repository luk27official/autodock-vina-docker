# Autodock Vina Docker
This is a very simple repository for the Autodock Vina tool performing basic docking a given ligand (SMILES format) to a protein structure (identified by a code, downloads PDB from an online RCSB/AF database).

To compile: `docker-compose build`

To run interactively: `docker-compose run --rm app`

Then, in the interactive session, run the wanted commands (such as `./run_task.py`).