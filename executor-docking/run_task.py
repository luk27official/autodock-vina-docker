#!/usr/bin/env python3

import os
import json
import subprocess

import rdkit
from rdkit.Chem import PandasTools
import pandas as pd
import meeko

def prepare_ligand(identifier: str, inputFile: str):
    smiles = ""
    with open(inputFile) as inp:
        input_json = json.load(inp)
        smiles = input_json['hash']

    lig = rdkit.Chem.MolFromSmiles(smiles)
    protonated_lig = rdkit.Chem.AddHs(lig)
    rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)

    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(protonated_lig)
    lig_pdbqt = meeko_prep.write_pdbqt_string()

    with open(f"{identifier}_ligand.pdbqt", "w") as f:
        f.write(lig_pdbqt)

def prepare_receptor(identifier: str, inputFile: str):
    # TODO: check if this works for AlphaFold DB
    subprocess.run(["wget", f"https://files.rcsb.org/download/{identifier}.pdb", "-O", f"{identifier}.pdb"])
    print(f"Downloaded {identifier}.pdb")

    subprocess.run(["lepro_linux_x86", f"{identifier}.pdb"])
    os.rename('pro.pdb',f'{identifier}_clean_H.pdb') # Output from lepro is pro.pdb
    print(f"Cleaned {identifier}.pdb")

    subprocess.run(["prepare_receptor", "-r", f"./{identifier}_clean_H.pdb", "-o", f"./{identifier}_receptor.pdbqt"])
    print(f"Prepared {identifier}_receptor.pdbqt")

    prepare_ligand(identifier, inputFile)

    #TODO: get sdf from somewhere??
    #subprocess.run(["mk_prepare_ligand.py", "-i", "1iep_ligand.sdf", "-o", f"{identifier}_ligand.pdbqt"])
    #print(f"Prepared {identifier}_ligand.pdbqt")

    # prepare box
    with open(f"{identifier}_receptor_vina_box.txt", "w") as f, open(inputFile) as inp:
        input_json = json.load(inp)
        f.write(f"center_x = {input_json['bounding_box']['center']['x']}\n")
        f.write(f"center_y = {input_json['bounding_box']['center']['y']}\n")
        f.write(f"center_z = {input_json['bounding_box']['center']['z']}\n")

        f.write(f"size_x = {input_json['bounding_box']['size']['x']}\n")
        f.write(f"size_y = {input_json['bounding_box']['size']['y']}\n")
        f.write(f"size_z = {input_json['bounding_box']['size']['z']}\n")

    # docking
    subprocess.run(["vina", "--receptor", f"./{identifier}_receptor.pdbqt", "--ligand", f"./{identifier}_ligand.pdbqt", "--config", f"./{identifier}_receptor_vina_box.txt", "--exhaustiveness=32", "--out", f"./{identifier}_out_vina.pdbqt"])

if __name__ == "__main__":
    prepare_receptor("1iep", "input.json")
    pass