from rdkit import Chem
from rdkit.Chem import AllChem

class QCMolecule:
    def __init__(self, smiles: str, charge: int = 0, spin: int = 0):
        self.smiles = smiles
        self.charge = charge
        self.spin = spin
        self.xyz_string = self._generate_3d_coordinates()

    def _generate_3d_coordinates(self) -> str:
        mol = Chem.MolFromSmiles(self.smiles)
        mol = Chem.AddHs(mol)
        
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)  # Quick classical pre-optimization
        
        conf = mol.GetConformer()
        xyz_lines = []
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
            
        return "; ".join(xyz_lines)