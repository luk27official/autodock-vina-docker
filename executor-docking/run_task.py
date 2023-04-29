#!/usr/bin/env python3

import os
import json
import subprocess

import rdkit
import meeko

def prepare_ligand(identifier: str, inputFile: str):
    smiles = ""
    with open(inputFile) as inp:
        input_json = json.load(inp)
        smiles = input_json['hash']
    lig = rdkit.Chem.MolFromSmiles(smiles)
    protonated_lig = rdkit.Chem.AddHs(lig)
    rdkit.Chem.AllChem.EmbedMolecule(protonated_lig,randomSeed=0xf00d,useRandomCoords=True)

    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(protonated_lig)
    lig_pdbqt = meeko_prep.write_pdbqt_string()

    with open(f"{identifier}_ligand.pdbqt", "w") as f:
        f.write(lig_pdbqt)

def prepare_bounding_box(identifier: str, inputFile: str):
    with open(f"{identifier}_receptor_vina_box.txt", "w") as f, open(inputFile) as inp:
        input_json = json.load(inp)
        f.write(f"center_x = {input_json['bounding_box']['center']['x']}\n")
        f.write(f"center_y = {input_json['bounding_box']['center']['y']}\n")
        f.write(f"center_z = {input_json['bounding_box']['center']['z']}\n")

        f.write(f"size_x = {input_json['bounding_box']['size']['x']}\n")
        f.write(f"size_y = {input_json['bounding_box']['size']['y']}\n")
        f.write(f"size_z = {input_json['bounding_box']['size']['z']}\n")

def dock_molecule(identifier: str, inputFile: str):
    # download pdb
    # TODO: check if this works for AlphaFold DB
    subprocess.run(["wget", f"https://files.rcsb.org/download/{identifier}.pdb", "-O", f"{identifier}.pdb"])
    
    # clean pdb using lePro tool
    subprocess.run(["lepro_linux_x86", f"{identifier}.pdb"])
    os.rename('pro.pdb',f'{identifier}_clean_H.pdb') # Output from lepro is pro.pdb

    # prepare receptor using prepare_receptor tool from ADFR suite
    subprocess.run(["prepare_receptor", "-r", f"./{identifier}_clean_H.pdb", "-o", f"./{identifier}_receptor.pdbqt"])

    prepare_ligand(identifier, inputFile)

    prepare_bounding_box(identifier, inputFile)

    # docking
    subprocess.run(["vina", "--receptor", f"./{identifier}_receptor.pdbqt", "--ligand", f"./{identifier}_ligand.pdbqt", "--config", f"./{identifier}_receptor_vina_box.txt", "--exhaustiveness=32", "--out", f"./{identifier}_out_vina.pdbqt"])

if __name__ == "__main__":
    #dock_molecule("1iep", "input.json")
    #dock_molecule("6WZO", "input_2.json")
    dock_molecule("1c8k", "input_3.json")
    pass