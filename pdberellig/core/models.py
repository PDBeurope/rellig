from dataclasses import dataclass

import rdkit
from pdbeccdutils.computations import parity_method
from pdbeccdutils.core.models import ParityResult


@dataclass
class Similarity:
    target_id: str
    query_id: str
    result: ParityResult = None


@dataclass
class CompareObj:
    id: str
    mol: rdkit.Chem.rdchem.Mol

    def similarity_to(self, other, threshold=0.01):
        result = parity_method.compare_molecules(self.mol, other.mol, threshold)
        return Similarity(self.id, other.id, result)


@dataclass
class CofactorSim:
    query_id: str
    template_sim: Similarity
    representative_sim: Similarity
