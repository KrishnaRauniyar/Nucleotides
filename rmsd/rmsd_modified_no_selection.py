# #program to calculate the similarity between structures using rmsd

# import csv
# import math
# import Bio.PDB
# from Bio.PDB import PDBParser
# import pandas as pd
# import os
# from scipy.stats import skew
# import numpy as np
# from itertools import combinations
# import os
# from pathlib import Path
# import argparse
# import requests  # Added for PDB downloading
# from urllib.parse import urljoin

# def download_pdb(pdb_id, download_dir="pdb_downloads"):
#     """Download PDB file from RCSB if not already exists"""
#     pdb_id = pdb_id.lower()
#     os.makedirs(download_dir, exist_ok=True)
#     pdb_path = os.path.join(download_dir, f"{pdb_id}.pdb")
    
#     if os.path.exists(pdb_path):
#         return pdb_path
    
#     base_url = "https://files.rcsb.org/download/"
#     url = urljoin(base_url, f"{pdb_id}.pdb")
    
#     try:
#         response = requests.get(url, stream=True, timeout=10)
#         response.raise_for_status()
        
#         with open(pdb_path, 'wb') as f:
#             for chunk in response.iter_content(chunk_size=8192):
#                 f.write(chunk)
#         return pdb_path
#     except Exception as e:
#         print(f"Failed to download {pdb_id}: {str(e)}")
#         return None

# def filter_csv_by_residue(csv_path, residue_type):
#     """Load and filter CSV data by residue type"""
#     df = pd.read_csv(csv_path)
#     df['drug_name'] = df['Protein'].str.split('_').str[-1]
#     df_sorted = df.sort_values('drug_name')
    
#     # Filter by residue type
#     df_filtered = df_sorted[df_sorted['drug_name'] == residue_type].copy()
#     if df_filtered.empty:
#         print(f"No structures found with residue type: {residue_type}")
#         return None
    
#     # Split the Protein column into separate components
#     split_data = df_filtered['Protein'].str.split('_', expand=True)
    
#     # Assign components to separate columns
#     df_filtered['protein'] = split_data[0]    
#     df_filtered['chain'] = split_data[1]       
#     df_filtered['drug_id'] = split_data[2]     
#     df_filtered['drug_name'] = split_data[3]
    
#     return df_filtered

# def calDist(x1, y1, z1, x2, y2, z2):
#     """Calculate distance between two atoms"""
#     return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

# def Cal_Centroid(residue):
#     """Calculates the centroid of structures"""
#     SumX = SumY = SumZ = counter = 0
#     for atom1 in residue:
#         atomCoord = atom1.get_vector()
#         SumX += atomCoord[0]
#         SumY += atomCoord[1]
#         SumZ += atomCoord[2]
#         counter += 1
#     return SumX/counter, SumY/counter, SumZ/counter

# def translation(CenX, CenY, CenZ, residue):
#     """Translates the structure to origin"""
#     residue_coord_dict = {}
#     for counter, atom in enumerate(residue):
#         atomCoord = atom.get_vector()
#         residue_coord_dict[counter] = [
#             CenX - atomCoord[0],
#             CenY - atomCoord[1],
#             CenZ - atomCoord[2]
#         ]
#     return residue_coord_dict

# def kabsch_umeyama(A, B):
#     """Algorithm for alignment by rotation and scaling"""
#     assert A.shape == B.shape, f"Shape mismatch: A.shape = {A.shape}, B.shape = {B.shape}"
#     n, m = A.shape
#     EA = np.mean(A, axis=0)
#     EB = np.mean(B, axis=0)
#     VarA = np.mean(np.linalg.norm(A - EA, axis=1) ** 2)
#     H = ((A - EA).T @ (B - EB)) / n
#     U, D, VT = np.linalg.svd(H)
#     d = np.sign(np.linalg.det(U) * np.linalg.det(VT))
#     S = np.diag([1] * (m - 1) + [d])
#     R = U @ S @ VT
#     c = VarA / np.trace(np.diag(D) @ S)
#     t = EA - c * R @ EB
#     return R, c, t

# def rmsd(A, B):
#     """Calculate RMSD between two structures"""
#     sum_sq_dist = 0
#     for a, b in zip(A, B):
#         sum_sq_dist += calDist(a[0], a[1], a[2], b[0], b[1], b[2]) ** 2
#     return np.sqrt(sum_sq_dist/len(A))

# def main():
#     # Parse command line arguments
#     parser = argparse.ArgumentParser(description='Calculate RMSD between structures')
#     parser.add_argument('--residue', type=str, required=True, help='Residue type to filter (e.g., A, G, C, U)')
#     parser.add_argument('--csv', type=str, help='Path to CSV file with structure data')
#     parser.add_argument('--pdb_dir', type=str, default="pdb_downloads", help='Directory to store/download PDB files')
#     args = parser.parse_args()

#     # Load and filter CSV data
#     df_filtered = filter_csv_by_residue(args.csv, args.residue)
#     if df_filtered is None:
#         return

#     # Extract data for processing
#     PDB_list = df_filtered['protein'].to_list()
#     Chain = df_filtered['chain'].to_list()
#     Drug_name = df_filtered['drug_name'].to_list()
#     Drug_id = df_filtered['drug_id'].to_list()

#     print(f"Processing structures with residue type: {args.residue}")
#     print(df_filtered[['Protein', 'drug_name']])

#     # Process PDB files and collect coordinates
#     drug_coord_dict = {}
#     Dataframe_Index = []
#     missing_pdbs = []

#     for i, (pdb_id, chain_id, drug_name, drug_id) in enumerate(zip(PDB_list, Chain, Drug_name, Drug_id)):
#         Drug = f'{pdb_id}_{chain_id}_{drug_name}_{drug_id}'
#         Dataframe_Index.append(Drug)
        
#         # Download PDB file if needed
#         pdb_path = download_pdb(pdb_id, args.pdb_dir)
#         if pdb_path is None:
#             missing_pdbs.append(pdb_id)
#             continue
            
#         try:
#             p = Bio.PDB.PDBParser()
#             structure = p.get_structure('PrimaryStructureChain', pdb_path)
            
#             for chain in structure[0]:
#                 if chain.id == chain_id:
#                     for residue in chain:
#                         if str(residue.id[1]) == drug_id and residue.get_resname() == args.residue:
#                             drug_coord = []
#                             for atom in residue:
#                                 atomCoord = atom.get_vector()
#                                 drug_coord.append([atomCoord[0], atomCoord[1], atomCoord[2]])
#                             drug_coord_dict[Drug] = drug_coord
#         except Exception as e:
#             print(f"Error processing {pdb_id}: {str(e)}")
#             continue

#     if missing_pdbs:
#         print(f"\nWarning: Could not download these PDB files: {', '.join(missing_pdbs)}")
#         print("These structures will be skipped in the RMSD calculations.")

#     if not drug_coord_dict:
#         print("No valid structures found for RMSD calculation.")
#         return

#     # Calculate RMSD for all pairs
#     output_file = 'similarity_drugs_rmsd.txt'
#     with open(output_file, 'w') as OutputFile:
#         OutputFile.write("Structure1\tStructure2\tRMSD\n")  # Header line
        
#         for (drug1, coords1), (drug2, coords2) in combinations(drug_coord_dict.items(), 2):
#             A = np.array(coords1)
#             B = np.array(coords2)
            
#             if A.shape != B.shape:
#                 print(f"Skipping {drug1} and {drug2} - different number of atoms")
#                 continue
                
#             try:
#                 R, c, t = kabsch_umeyama(A, B)
#                 B_aligned = np.array([t + c * R @ b for b in B])
#                 rmsd_value = rmsd(A, B_aligned)
#                 OutputFile.write(f'{drug1}\t{drug2}\t{rmsd_value:.4f}\n')
#             except Exception as e:
#                 print(f"Error calculating RMSD for {drug1} vs {drug2}: {str(e)}")
#                 continue

#     print(f'\nCompleted successfully. Results saved to {output_file}')

# if __name__ == "__main__":
#     main()


import csv
import math
import Bio.PDB
from Bio.PDB import PDBParser
import pandas as pd
import os
import numpy as np
from itertools import combinations
from pathlib import Path
import argparse
import requests
from urllib.parse import urljoin

# Define common atom subgroups for nucleotides
SUBGROUPS = {
    'backbone': {'P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'C1\''},
    'sugar': {'C5\'', 'C4\'', 'O4\'', 'C3\'', 'C2\'', 'C1\''},
    'phosphate': {'P', 'OP1', 'OP2', 'O5\''},
    'base_A': {'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'},      # Adenine (RNA)
    'base_DA': {'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'},     # Deoxyadenine (DNA)
    'base_G': {'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'}, # Guanine (RNA)
    'base_DG': {'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'},# Deoxyguanine (DNA)
    'base_DT': {'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'},                 # Deoxythymine (DNA, no methyl)
    'base_C': {'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'},                  # Cytosine (RNA)
    'base_U': {'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'},                  # Uracil (RNA)
    'base_DC': {'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'}                  # Deoxycytosine (DNA)
}

def download_pdb(pdb_id, download_dir="pdb_downloads"):
    """Download PDB file from RCSB if not already exists"""
    pdb_id = pdb_id.lower()
    os.makedirs(download_dir, exist_ok=True)
    pdb_path = os.path.join(download_dir, f"{pdb_id}.pdb")
    
    if os.path.exists(pdb_path):
        return pdb_path
    
    base_url = "https://files.rcsb.org/download/"
    url = urljoin(base_url, f"{pdb_id}.pdb")
    
    try:
        response = requests.get(url, stream=True, timeout=10)
        response.raise_for_status()
        with open(pdb_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return pdb_path
    except Exception as e:
        print(f"Failed to download {pdb_id}: {str(e)}")
        return None

def filter_csv_by_residue(csv_path, residue_type=None):
    """Load and filter CSV data by residue type (optional)"""
    df = pd.read_csv(csv_path)
    df['drug_name'] = df['Protein'].str.split('_').str[-1]
    df_sorted = df.sort_values('drug_name')
    
    if residue_type:
        df_filtered = df_sorted[df_sorted['drug_name'] == residue_type].copy()
        if df_filtered.empty:
            print(f"No structures found with residue type: {residue_type}")
            return None
    else:
        df_filtered = df_sorted.copy()  # Use all residues if no filter
    
    split_data = df_filtered['Protein'].str.split('_', expand=True)
    df_filtered['protein'] = split_data[0]
    df_filtered['chain'] = split_data[1]
    df_filtered['drug_id'] = split_data[2]
    df_filtered['drug_name'] = split_data[3]
    
    return df_filtered

def calDist(x1, y1, z1, x2, y2, z2):
    """Calculate distance between two atoms"""
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

def kabsch_umeyama(A, B):
    """Algorithm for alignment by rotation and scaling"""
    assert A.shape == B.shape, f"Shape mismatch: A={A.shape}, B={B.shape}"
    n, m = A.shape
    EA = np.mean(A, axis=0)
    EB = np.mean(B, axis=0)
    VarA = np.mean(np.linalg.norm(A - EA, axis=1) ** 2)
    H = ((A - EA).T @ (B - EB)) / n
    U, D, VT = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U) * np.linalg.det(VT))
    S = np.diag([1] * (m - 1) + [d])
    R = U @ S @ VT
    c = VarA / np.trace(np.diag(D) @ S)
    t = EA - c * R @ EB
    return R, c, t

def rmsd(A, B):
    """Calculate RMSD between two structures"""
    sum_sq_dist = 0
    for a, b in zip(A, B):
        sum_sq_dist += calDist(a[0], a[1], a[2], b[0], b[1], b[2]) ** 2
    return np.sqrt(sum_sq_dist / len(A))

def extract_subgroup_coords(residue, subgroup_name):
    """Extract coordinates for a specified subgroup of atoms"""
    if subgroup_name not in SUBGROUPS:
        raise ValueError(f"Unknown subgroup: {subgroup_name}. Available: {list(SUBGROUPS.keys())}")
    
    subgroup_atoms = SUBGROUPS[subgroup_name]
    coords = []
    atom_names = []
    
    for atom in residue:
        atom_name = atom.get_name()
        if atom_name in subgroup_atoms:
            atomCoord = atom.get_vector()
            coords.append([atomCoord[0], atomCoord[1], atomCoord[2]])
            atom_names.append(atom_name)
    
    # Check if all expected atoms are present
    missing_atoms = subgroup_atoms - set(atom_names)
    if missing_atoms:
        print(f"Warning: Missing atoms {missing_atoms} in {residue} for subgroup {subgroup_name}")
    
    return coords

def process_structures(df, pdb_dir, subgroup_name=None):
    """Process PDB files and extract coordinates, optionally for a subgroup"""
    PDB_list = df['protein'].to_list()
    Chain = df['chain'].to_list()
    Drug_name = df['drug_name'].to_list()
    Drug_id = df['drug_id'].to_list()
    
    drug_coord_dict = {}
    Dataframe_Index = []
    missing_pdbs = []
    
    for pdb_id, chain_id, drug_name, drug_id in zip(PDB_list, Chain, Drug_name, Drug_id):
        Drug = f'{pdb_id}_{chain_id}_{drug_name}_{drug_id}'
        Dataframe_Index.append(Drug)
        
        pdb_path = download_pdb(pdb_id, pdb_dir)
        if pdb_path is None:
            missing_pdbs.append(pdb_id)
            continue
            
        try:
            p = PDBParser()
            structure = p.get_structure('PrimaryStructureChain', pdb_path)
            for chain in structure[0]:
                if chain.id == chain_id:
                    for residue in chain:
                        if str(residue.id[1]) == drug_id and residue.get_resname() == drug_name:
                            if subgroup_name:
                                coords = extract_subgroup_coords(residue, subgroup_name)
                            else:
                                coords = [atom.get_vector()[:] for atom in residue]
                            if coords:  # Only add if coordinates were extracted
                                drug_coord_dict[Drug] = coords
        except Exception as e:
            print(f"Error processing {pdb_id}: {str(e)}")
            continue
    
    return drug_coord_dict, missing_pdbs, Dataframe_Index

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate RMSD between structures')
    parser.add_argument('--csv', type=str, required=True, help='Path to CSV file with structure data')
    parser.add_argument('--residue', type=str, help='Residue type to filter (e.g., A, G, C, U), optional')
    parser.add_argument('--subgroup', type=str, choices=SUBGROUPS.keys(), help='Atom subgroup for alignment (e.g., backbone)')
    parser.add_argument('--pdb_dir', type=str, default="pdb_downloads", help='Directory to store/download PDB files')
    args = parser.parse_args()

    # Load and filter CSV data
    df_filtered = filter_csv_by_residue(args.csv, args.residue)
    if df_filtered is None:
        return

    print(f"Processing structures from {args.csv}")
    if args.residue:
        print(f"Filtered by residue type: {args.residue}")
    if args.subgroup:
        print(f"Using subgroup: {args.subgroup}")
    print(df_filtered[['Protein', 'drug_name']])

    # Process structures
    drug_coord_dict, missing_pdbs, Dataframe_Index = process_structures(df_filtered, args.pdb_dir, args.subgroup)
    
    if missing_pdbs:
        print(f"\nWarning: Could not download these PDB files: {', '.join(missing_pdbs)}")
    
    if not drug_coord_dict:
        print("No valid structures found for RMSD calculation.")
        return

    # Calculate RMSD for all pairs
    output_file = f'similarity_drugs_rmsd{"_" + args.subgroup if args.subgroup else ""}.txt'
    with open(output_file, 'w') as OutputFile:
        OutputFile.write("Structure1\tStructure2\tRMSD\n")
        
        for (drug1, coords1), (drug2, coords2) in combinations(drug_coord_dict.items(), 2):
            A = np.array(coords1)
            B = np.array(coords2)
            
            if A.shape != B.shape:
                print(f"Skipping {drug1} vs {drug2} - different number of atoms ({A.shape} vs {B.shape})")
                continue
                
            try:
                R, c, t = kabsch_umeyama(A, B)
                B_aligned = np.array([t + c * R @ b for b in B])
                rmsd_value = rmsd(A, B_aligned)
                print(f"{drug1} vs {drug2}: RMSD = {rmsd_value:.4f}")
                OutputFile.write(f'{drug1}\t{drug2}\t{rmsd_value:.4f}\n')
            except Exception as e:
                print(f"Error calculating RMSD for {drug1} vs {drug2}: {str(e)}")
                continue

    print(f'\nCompleted successfully. Results saved to {output_file}')

if __name__ == "__main__":
    main()