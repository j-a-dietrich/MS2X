from collections import Counter
import re

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
#from rdkit.Chem.GraphDescriptors import BertzCT
from ms2x.adapted_BertzCT import BertzCT
#from adapted_BertzCT import BertzCT
from rdkit.Chem.rdmolfiles import CanonicalRankAtoms

from find_mfs import FormulaFinder
import time
import threading


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

def parse_formula(formula: str) -> Counter:
    """Parse molecular formula -> element counter (without H)."""
    pattern = r"([A-Z][a-z]?)(\d*)"

    parsed = Counter()

    for element, count in re.findall(pattern, formula):
        n = int(count) if count else 1
        if element != "H":
            parsed[element] += n

    return parsed


def parse_formula_H(formula: str) -> Counter:
    """Parse molecular formula -> element counter (without H)."""
    pattern = r"([A-Z][a-z]?)(\d*)"

    parsed = Counter()

    for element, count in re.findall(pattern, formula):
        n = int(count) if count else 1
        parsed[element] += n

    return parsed


def mol_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    formula = rdMolDescriptors.CalcMolFormula(mol)

    return formula, mol


# ---------------------------------------------------------
# Main class
# ---------------------------------------------------------

class MS2X:

    def __init__(self):

        self.ff = FormulaFinder(
            "CHNOPSFClBrISiBAsSeSb"
        )
        
        self.predefined_substructures = {
            "[C5H5]+": "ccccc"
        }

        self.max_time = 0.2
    # -----------------------------------------------------
    # Public API
    # -----------------------------------------------------

    def extract_info(self, spec):
        """
        Extract:
        - molecular formula
        - mol object
        - fragment formulas
        """

        formula, mol = self._get_molecule(spec)

        peaks = spec.peaks.mz

        frag_formulas = [
            self._get_fragment_formula(mz, formula)
            for mz in peaks
        ]

        frag_formulas = [f for f in frag_formulas if f is not None]

        return {
            "formula": formula,
            "mol": mol,
            "frag_formulas": frag_formulas,
        }

    # -----------------------------------------------------
    # Molecule
    # -----------------------------------------------------

    def _get_molecule(self, spec):

        smiles = spec.get("smiles")

        formula, mol = mol_from_smiles(smiles)

        return formula, mol

    # -----------------------------------------------------
    # Fragment formula
    # -----------------------------------------------------

    def _get_fragment_formula(
        self,
        frag_mz,
        mol_formula,
    ):

        df = (
            self.ff.find_formulae(
                frag_mz,
                charge=self.charge,
                error_ppm=50,
                filter_rdbe=(0, 20),
                max_counts=mol_formula,
            )
            .to_dataframe()
        )

        if df.empty:
            return None

        # try first few candidates

        for i in range(min(5, len(df))):

            f = df.formula.iloc[i]

            parsed = parse_formula_H(f)

            if not self._formula_is_valid(parsed):
                continue

            return f

        return None

    # -----------------------------------------------------
    # Filters
    # -----------------------------------------------------

    def _formula_is_valid(self, parsed):

        h = parsed.get("H", 0)
        c = parsed.get("C", 1)

        if h / c <= 0.25:
            return False

        return True

    # -----------------------------------------------------
    # Subgraph search
    # -----------------------------------------------------
    
    def complete_dfs(self, mol, target_formula):

        target = parse_formula(target_formula)

        memory = set()
        results = set()

        self.start_time = time.time()
        stop_event = threading.Event() 

        mol_counts = Counter(
            atom.GetSymbol() for atom in mol.GetAtoms()
        )

        start_element = min(
            target.keys(),
            key=lambda el: target[el] * mol_counts.get(el, 999),
        )

        nC = target.get("C", 0)

        def run():
            # =========================================
            # TEMPLATE MODE
            # =========================================
            
            if nC >= 10:

                ring_info = mol.GetRingInfo()

                for ring in ring_info.AtomRings():
                    if stop_event.is_set():
                        return

                    atom_set = set(ring)

                    count = Counter(
                        mol.GetAtomWithIdx(i).GetSymbol()
                        for i in atom_set
                        if mol.GetAtomWithIdx(i).GetSymbol() != "H"
                    )

                    # skip too big
                    if any(count[k] > target.get(k, 0) for k in count):
                        continue

                    self._dfs(
                        mol,
                        atom_set,
                        count,
                        target,
                        memory,
                        results,
                        stop_event
                    )

                return

            # =========================================
            # NORMAL MODE
            # =========================================
            for atom in mol.GetAtoms():
                if stop_event.is_set():
                    return
                if atom.GetSymbol() != start_element:
                    continue

                start = {atom.GetIdx()}
                count = Counter([start_element])

                self._dfs(
                    mol,
                    start,
                    count,
                    target,
                    memory,
                    results,
                    stop_event
                )

        thread = threading.Thread(target=run)
        thread.start()
        thread.join(timeout=self.max_time)
        stop_event.set() 
        thread.join()   

        #print(f"DFS took {time.time() - self.start_time:.2f}s, found {len(results)} results")
        return results
    

    def _dfs(
        self,
        mol,
        atom_set,
        element_count,
        target,
        memory,
        results,
        stop_event
    ):

        if stop_event.is_set():
            return

        key = (
            frozenset(atom_set),
            frozenset(element_count.items()),
        )

        if key in memory:
            return

        memory.add(key)

        # end condition

        if element_count == target:
            results.add(frozenset(atom_set))
            return
        
        #if not self._heuristics(
        #    mol,
        #    atom_set,
        #    results
        #    
        #):
        #    return

        # expand

        neighbors = set()

        for idx in atom_set:

            atom = mol.GetAtomWithIdx(idx)

            for n in atom.GetNeighbors():

                n_idx = n.GetIdx()

                if n_idx not in atom_set:
                    neighbors.add(n_idx)

        for n_idx in neighbors:

            atom = mol.GetAtomWithIdx(n_idx)

            el = atom.GetSymbol()

            if el == "H":
                continue

            if el not in target:
                continue

            new_count = element_count.copy()
            new_count[el] += 1

            if new_count[el] > target[el]:
                continue

            self._dfs(
                mol,
                atom_set | {n_idx},
                new_count,
                target,
                memory,
                results,
                stop_event
            )

    # -------------------------
    # Early stopping criteria for DFS
    # -------------------------

    def _heuristics(
            self,
            mol,
            atom_set,
            results,
        ):

            # not enough info yet
            if not results or len(atom_set) < 6:
                return True

            len_atom_set = len(atom_set)

            best_score = None
            best_overlap = 0

            for sg in results:

                overlap = len(atom_set.intersection(sg))
                overlap_ratio = overlap / len_atom_set

                if overlap_ratio < 0.7:
                    continue

                # get score
                if sg not in self._score_cache:
                    frag = self.build_fragment(mol, sg)
                    losses = self.build_losses(mol, sg)
                    self._score_cache[sg] = self.compute_complexity(
                        frag,
                        losses,
                    )

                score = self._score_cache[sg]

                if best_score is None or score > best_score:
                    best_score = score
                    best_overlap = overlap_ratio

            # prune if very similar to a good fragment
            if best_score is not None and best_overlap > 0.85:
                return False

            return True

    
    # -----------------------------------------------------
    # Build substructures
    # -----------------------------------------------------

    def build_fragment(self, mol, atom_set):

        rw = Chem.RWMol()

        idx_map = {}

        for old_idx in atom_set:

            atom = mol.GetAtomWithIdx(old_idx)

            new_idx = rw.AddAtom(
                Chem.Atom(atom.GetAtomicNum())
            )

            idx_map[old_idx] = new_idx

        for old_idx in atom_set:

            atom = mol.GetAtomWithIdx(old_idx)

            for n in atom.GetNeighbors():

                j = n.GetIdx()

                if j not in atom_set:
                    continue

                if idx_map[old_idx] < idx_map[j]:

                    bond = mol.GetBondBetweenAtoms(
                        old_idx,
                        j,
                    )

                    rw.AddBond(
                        idx_map[old_idx],
                        idx_map[j],
                        bond.GetBondType(),
                    )

        return rw.GetMol()


    def build_losses(self, mol, atom_set):

        rw = Chem.RWMol(mol)

        for idx in sorted(atom_set, reverse=True):
            rw.RemoveAtom(idx)

        frags = Chem.GetMolFrags(
            rw.GetMol(),
            asMols=True,
            sanitizeFrags=False,
        )

        return frags

    # -----------------------------------------------------
    # Complexity
    # -----------------------------------------------------

    def compute_complexity(
        self,
        fragment,
        losses,
    ):

        score = BertzCT(fragment)

        for m in losses:
            score += BertzCT(m)

        return score
    

    # -----------------------------------------------------
    # Checking predefined pattern
    # -----------------------------------------------------

    def check_predefined_substructures(self, mol, frag_formula):
        for formula, smiles in self.predefined_substructures.items():
            if formula == frag_formula:
                substructure = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol.HasSubstructMatch(substructure):
                    return substructure
        return False

    # -----------------------------------------------------
    # Pipeline
    # -----------------------------------------------------

    def approximate_substructure(self, spec, charge=1):

        results = []
        formulas = []
        complexities = []
        self.charge = charge

        self._score_cache = set()

        info = self.extract_info(spec)

        mol = info["mol"]
        mol_formula = info["formula"]
        frag_formulas = info["frag_formulas"]

        for f in frag_formulas:

            substructure = self.check_predefined_substructures(mol, f)
            if substructure:
                results.append(substructure)
                formulas.append(f)
                continue

            subgraphs = self.complete_dfs(
                mol,
                f,
            )

            best = None
            best_score = 0

            for sg in subgraphs:

                frag = self.build_fragment(
                    mol,
                    sg,
                )

                losses = self.build_losses(
                    mol,
                    sg,
                )

                score = self.compute_complexity(
                    frag,
                    losses,
                )

                score = score / BertzCT(mol)

                if score > best_score:
                    best_score = score
                    best = frag

            if best is not None:
                results.append(best)
                formulas.append(f)
                complexities.append(best_score)
                
        return results, formulas, complexities
    

    def approximate_substructure_by_fragments(self, fragments, smiles, charge=1):

        results = []
        formulas = []
        complexities = []
        self.charge = charge

        self._score_cache = set()


        mol = Chem.MolFromSmiles(smiles)
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        frag_formulas = [self._get_fragment_formula(frag_mz, mol_formula) for frag_mz in fragments]

        for f in frag_formulas:
            if not f: # this needs to be pushed
                results.append(None)
                formulas.append(None)
                complexities.append(0)
                continue

            substructure = self.check_predefined_substructures(mol, f)
            if substructure:
                results.append(substructure)
                formulas.append(f)
                continue

            subgraphs = self.complete_dfs(
                mol,
                f,
            )

            best = None
            best_score = 0

            for sg in subgraphs:

                frag = self.build_fragment(
                    mol,
                    sg,
                )

                losses = self.build_losses(
                    mol,
                    sg,
                )

                score = self.compute_complexity(
                    frag,
                    losses,
                )

                score = score / BertzCT(mol)

                if score > best_score:
                    best_score = score
                    best = frag

            if best is not None:
                results.append(best)
                formulas.append(f)
                complexities.append(best_score)

        return results, formulas, complexities
