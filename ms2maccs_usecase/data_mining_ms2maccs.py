import os
import numpy as np
import pickle
from concurrent.futures import ProcessPoolExecutor, as_completed
from matchms.importing import load_from_mgf
from matchms.filtering import reduce_to_number_of_peaks
from tqdm import tqdm
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
import numpy as np
from rdkit import Chem

import ms2x
from adapted_BertzCT import BertzCT
ms = ms2x.MS2X()

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
from multiprocessing import Lock

import re
from collections import Counter

import sys
sys.path.append("../")

unsanitized_mols = []  # collect (frag, submol) tuples globally
unsanitized_lock = Lock()


ELEMENTS = ["C", "H", "O", "N", "Cl", "F", "Br", "I", "S", "P"]

def parse_formula(formula):
    parsed = Counter()
    for elem, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula):
        parsed[elem] += int(count) if count else 1
    return parsed

def formula_to_vector(formula):
    vec = np.zeros(len(ELEMENTS), dtype=np.float32)
    if not formula:
        return vec
    counts = parse_formula(str(formula))
    for i, elem in enumerate(ELEMENTS):
        vec[i] = counts.get(elem, 0)
    mx = vec.max()
    if mx > 1e-8:
        vec /= mx
    return vec


def process_spec_simple(spec, charge):
    local_mapping = {}
    local_unsanitized = []  # collect locally, merge after
    spec = reduce_to_number_of_peaks(spec, n_max=60)
    smiles = spec.get("smiles")
    if not smiles:
        return local_mapping, local_unsanitized
    for frag in spec.peaks.mz:
        try:
            submol, subformula = ms.approximate_substructure_by_fragments([frag], smiles, charge=charge)
        except ValueError:
            continue
        if not submol:
            continue
        submol = submol[0]
        subformula = subformula[0]

        local_unsanitized.append((round(frag, 2), submol))
        try:
            Chem.GetSymmSSSR(submol)
            Chem.SanitizeMol(submol)
            maccs = np.array(MACCSkeys.GenMACCSKeys(submol)) 
        except:
            Chem.GetSymmSSSR(submol)
            maccs = np.array(MACCSkeys.GenMACCSKeys(submol))
        frag_rounded = round(frag, 2)
        formula_vector = formula_to_vector(subformula)
        if frag_rounded not in local_mapping:
            local_mapping[frag_rounded] = np.concatenate([maccs,formula_vector]) # with formula vector
        else:
            local_mapping[frag_rounded] += np.concatenate([maccs, formula_vector])
    return local_mapping, local_unsanitized

if __name__ == "__main__":
    from threading import Lock

    # variables
    path_name = "test_specs_H_p_mode.mgf" # input.mgf file
    save_name = "fp_bit_map_test_H_p_mode.pkl" # output.pkl file
    charge = 1 # 1 for positive, -1 for negative

    print("Loading spectra...")
    
    spectra = list(tqdm(load_from_mgf(path_name)))
    print(f"Loaded {len(spectra)} spectra")

    fp_bit_mapping = {}
    unsanitized_mols = []

    with ProcessPoolExecutor(max_workers=os.cpu_count() - 1) as executor:
        futures = {executor.submit(process_spec_simple, spec, charge): spec for spec in spectra}
        for future in tqdm(as_completed(futures), total=len(spectra), desc="Building fp_bit_map"):
            try:
                local_mapping, local_unsanitized = future.result()
                unsanitized_mols.extend(local_unsanitized)  # merge unsanitized
                for frag, maccs in local_mapping.items():
                    if frag not in fp_bit_mapping:
                        fp_bit_mapping[frag] = maccs
                    else:
                        fp_bit_mapping[frag] += maccs
            except Exception as e:
                print(f"Worker failed: {e}")
                continue

    print(f"Built mapping with {len(fp_bit_mapping)} unique fragments")
    print(f"Collected {len(unsanitized_mols)} unsanitized mols")

    with open(save_name, "wb") as f:
        pickle.dump(fp_bit_mapping, f)
    print("Saved")