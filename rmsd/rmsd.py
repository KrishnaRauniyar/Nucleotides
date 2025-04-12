import math
from Bio.PDB import PDBParser
import pandas as pd
import os
import numpy as np
from itertools import combinations
import argparse
import requests
from urllib.parse import urljoin
from multiprocessing import Pool, Manager
from functools import partial
import psutil
import time

# Define common atom subgroups for nucleotides
SUBGROUPS = {
    'backbone': {'P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'C1\''},
    'sugar': {'C5\'', 'C4\'', 'O4\'', 'C3\'', 'C2\'', 'C1\''},
    'phosphate': {'P', 'OP1', 'OP2', 'O5\''},
    'base_A': {'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'},
    'base_DA': {'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'},
    'base_G': {'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'},
    'base_DG': {'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'},
    'base_DT': {'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'},
    'base_C': {'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'},
    'base_U': {'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'},
    'base_DC': {'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'}
}

def download_pdb(pdb_id, download_dir="pdb_downloads"):
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
    df = pd.read_csv(csv_path)
    df['drug_name'] = df['Protein'].str.split('_').str[-1]
    df_sorted = df.sort_values('drug_name')
    if residue_type:
        df_filtered = df_sorted[df_sorted['drug_name'] == residue_type].copy()
        if df_filtered.empty:
            print(f"No structures found with residue type: {residue_type}")
            return None
    else:
        df_filtered = df_sorted.copy()
    split_data = df_filtered['Protein'].str.split('_', expand=True)
    df_filtered['protein'] = split_data[0]
    df_filtered['chain'] = split_data[1]
    df_filtered['drug_id'] = split_data[2]
    df_filtered['drug_name'] = split_data[3]
    return df_filtered

def calDist(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

def kabsch_umeyama(A, B):
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
    sum_sq_dist = 0
    for a, b in zip(A, B):
        sum_sq_dist += calDist(a[0], a[1], a[2], b[0], b[1], b[2]) ** 2
    return np.sqrt(sum_sq_dist / len(A))

def extract_subgroup_coords(residue, subgroup_name):
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
    missing_atoms = subgroup_atoms - set(atom_names)
    if missing_atoms:
        print(f"Warning: Missing atoms {missing_atoms} in {residue} for subgroup {subgroup_name}")
    return coords

def process_structure(args):
    pdb_id, chain_id, drug_name, drug_id, pdb_dir, subgroup_name, shared_dict = args
    Drug = f'{pdb_id}_{chain_id}_{drug_id}_{drug_name}'
    pdb_path = download_pdb(pdb_id, pdb_dir)

    if pdb_path is None:
        return Drug, None
    try:
        p = PDBParser()
        structure = p.get_structure('PrimaryStructureChain', pdb_path)
        for chain in structure[0]:
            if chain.id == chain_id:
                for residue in chain:
                    if str(residue.id[1]) == drug_id and residue.get_resname() == drug_name and residue.get_atoms() != "H":
                        if subgroup_name:
                            coords = extract_subgroup_coords(residue, subgroup_name)
                        else:
                            coords = [atom.get_vector()[:] for atom in residue if not (atom.get_name().startswith('H') or atom.element == 'H')]
                        if coords:
                            shared_dict[Drug] = coords
                            return Drug, coords
    except Exception as e:
        print(f"Error processing {pdb_id}: {str(e)}")
    return Drug, None

def compute_rmsd_pair(pair, drug_coord_dict):
    drug1, drug2 = pair
    A = np.array(drug_coord_dict[drug1])
    B = np.array(drug_coord_dict[drug2])
    if A.shape != B.shape:
        return f"{drug1}\t{drug2}\tSkipped (different atom counts: {A.shape} vs {B.shape})"
    try:
        R, c, t = kabsch_umeyama(A, B)
        B_aligned = np.array([t + c * R @ b for b in B])
        rmsd_value = rmsd(A, B_aligned)
        return f"{drug1}\t{drug2}\t{rmsd_value:.4f}"
    except Exception as e:
        return f"{drug1}\t{drug2}\tError: {str(e)}"

def log_resource_usage(process, label):
    """Log CPU and memory usage for the current process"""
    cpu_percent = process.cpu_percent(interval=1)  # Measure over 1 second
    memory_info = process.memory_info()
    memory_mb = memory_info.rss / 1024 / 1024  # Convert bytes to MB
    print(f"{label} - CPU Usage: {cpu_percent:.2f}%, Memory Usage: {memory_mb:.2f} MB")

def main():
    parser = argparse.ArgumentParser(description='Calculate RMSD between structures')
    parser.add_argument('--csv', type=str, required=True, help='Path to CSV file with structure data')
    parser.add_argument('--residue', type=str, help='Residue type to filter (e.g., A, G, C, U), optional')
    parser.add_argument('--subgroup', type=str, choices=SUBGROUPS.keys(), help='Atom subgroup for alignment (e.g., backbone)')
    parser.add_argument('--pdb_dir', type=str, default="pdb_downloads", help='Directory to store/download PDB files')
    parser.add_argument('--processes', type=int, default=os.cpu_count(), help='Number of parallel processes')
    args = parser.parse_args()

    # Get process object for resource monitoring
    process = psutil.Process(os.getpid())
    
    # Load and filter CSV data
    start_time = time.time()
    df_filtered = filter_csv_by_residue(args.csv, args.residue)
    if df_filtered is None:
        return

    print(f"Processing structures from {args.csv}")
    if args.residue:
        print(f"Filtered by residue type: {args.residue}")
    if args.subgroup:
        print(f"Using subgroup: {args.subgroup}")
    print(df_filtered[['Protein', 'drug_name']])

    # Prepare data for parallel processing
    PDB_list = df_filtered['protein'].to_list()
    Chain = df_filtered['chain'].to_list()
    Drug_name = df_filtered['drug_name'].to_list()
    Drug_id = df_filtered['drug_id'].to_list()
    Dataframe_Index = [f'{p}_{c}_{n}_{i}' for p, c, n, i in zip(PDB_list, Chain, Drug_name, Drug_id)]

    # Parallelize structure processing
    manager = Manager()
    drug_coord_dict = manager.dict()
    missing_pdbs = set()

    process_args = [(pdb_id, chain_id, drug_name, drug_id, args.pdb_dir, args.subgroup, drug_coord_dict)
                    for pdb_id, chain_id, drug_name, drug_id in zip(PDB_list, Chain, Drug_name, Drug_id)]
    
    print("\nStarting structure processing...")
    log_resource_usage(process, "Before structure processing")
    with Pool(processes=args.processes) as pool:
        results = pool.map(process_structure, process_args)
    log_resource_usage(process, "After structure processing")
    
    for drug, coords in results:
        if coords is None:
            missing_pdbs.add(drug.split('_')[0])

    if missing_pdbs:
        print(f"\nWarning: Could not process these PDB files: {', '.join(sorted(missing_pdbs))}")
    
    if not drug_coord_dict:
        print("No valid structures found for RMSD calculation.")
        return

    # Parallelize RMSD calculations
    output_file = f'similarity_drugs_rmsd{"_" + args.subgroup if args.subgroup else args.residue}.txt'
    pairs = list(combinations(drug_coord_dict.keys(), 2))
    
    print(f"\nStarting RMSD calculations for {len(pairs)} pairs...")
    log_resource_usage(process, "Before RMSD calculation")
    with Pool(processes=args.processes) as pool:
        rmsd_results = pool.map(partial(compute_rmsd_pair, drug_coord_dict=dict(drug_coord_dict)), pairs)
    log_resource_usage(process, "After RMSD calculation")
    
    # Write results to file
    with open(output_file, 'w') as OutputFile:
        OutputFile.write("Structure1\tStructure2\tRMSD\n")
        for result in rmsd_results:
            OutputFile.write(result + "\n")
            print(result.replace("\t", " "))

###### This section is for generating input file for clustermap
    # Generate clustermap_input.txt
    print("\nGenerating clustermap_input.txt...")
    clustermap_file = 'clustermap_input.txt'
    
    # Get sorted list of structures
    structures = sorted(drug_coord_dict.keys())
    n = len(structures)
    
    # Initialize RMSD matrix with zeros
    rmsd_matrix = [[0.0] * n for _ in range(n)]
    struct_to_idx = {struct: idx for idx, struct in enumerate(structures)}
    
    # Fill matrix with RMSD values
    for result in rmsd_results:
        parts = result.split('\t')
        if len(parts) == 3 and parts[2].replace('.', '').replace('-', '').isdigit():
            s1, s2, rmsd = parts
            i, j = struct_to_idx[s1], struct_to_idx[s2]
            rmsd_value = float(rmsd)
            rmsd_matrix[i][j] = rmsd_value
            rmsd_matrix[j][i] = rmsd_value  # Assume symmetry
    
    # Write clustermap_input.txt
    with open(clustermap_file, 'w') as f:
        for i, struct in enumerate(structures):
            rmsd_values = ','.join(f"{val:.4f}" for val in rmsd_matrix[i])
            f.write(f"{struct};{rmsd_values}\n")
    
    print(f"Clustermap input saved to {clustermap_file}")

    # Final resource usage and timing
    end_time = time.time()
    execution_time = end_time - start_time
    log_resource_usage(process, "Final")
    print(f"\nExecution time: {execution_time:.2f} seconds")
    print(f"Completed successfully. Results saved to {output_file}")

if __name__ == "__main__":
    main()