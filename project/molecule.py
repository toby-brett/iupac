import random
from functools import partial

import requests
from rdkit import Chem

def generate_acid_anhydride(max_length: int = 10):

    mol = Chem.RWMol()

    length1 = random.randint(2, max_length)
    length2 = random.randint(2, max_length)

    atoms1 = [mol.AddAtom(Chem.Atom("C")) for _ in range(length1)]
    O1 = mol.AddAtom(Chem.Atom("O"))
    O2 = mol.AddAtom(Chem.Atom("O"))
    O3 = mol.AddAtom(Chem.Atom("O"))
    atoms2 = [mol.AddAtom(Chem.Atom("C")) for _ in range(length2)]

    i = 0
    for i in range(length1 - 1):
        mol.AddBond(atoms1[i], atoms1[i + 1], Chem.BondType.SINGLE)

    mol.AddBond(atoms1[i + 1], O1, Chem.BondType.DOUBLE)
    mol.AddBond(atoms1[i + 1], O2, Chem.BondType.SINGLE)
    mol.AddBond(O2, atoms2[0], Chem.BondType.SINGLE)
    mol.AddBond(atoms2[0], O3, Chem.BondType.DOUBLE)

    for i in range(length2 - 1):
        mol.AddBond(atoms2[i], atoms2[i + 1], Chem.BondType.SINGLE)

    return mol

def generate_ester(max_length):

    mol = Chem.RWMol()

    length1 = random.randint(2, max_length)
    length2 = random.randint(2, max_length)

    atoms1 = [mol.AddAtom(Chem.Atom("C")) for _ in range(length1)]
    O1 = mol.AddAtom(Chem.Atom("O"))
    O2 = mol.AddAtom(Chem.Atom("O"))
    atoms2 = [mol.AddAtom(Chem.Atom("C")) for _ in range(length2)]

    i = 0
    for i in range(length1 - 1):
        mol.AddBond(atoms1[i], atoms1[i + 1], Chem.BondType.SINGLE)

    mol.AddBond(atoms1[i + 1], O1, Chem.BondType.DOUBLE)
    mol.AddBond(atoms1[i + 1], O2, Chem.BondType.SINGLE)
    mol.AddBond(O2, atoms2[0], Chem.BondType.SINGLE)

    for i in range(length2 - 1):
        mol.AddBond(atoms2[i], atoms2[i + 1], Chem.BondType.SINGLE)

    return mol

def generate_amide(max_length):

    mol = Chem.RWMol()

    length1 = random.randint(2, max_length)
    length2 = random.randint(1, max_length)

    atoms1 = [mol.AddAtom(Chem.Atom("C")) for _ in range(length1)]
    O1 = mol.AddAtom(Chem.Atom("O"))
    O2 = mol.AddAtom(Chem.Atom("N"))
    atoms2 = [mol.AddAtom(Chem.Atom("C")) for _ in range(length2)]

    i = -1
    for i in range(length1 - 1):
        mol.AddBond(atoms1[i], atoms1[i + 1], Chem.BondType.SINGLE)

    mol.AddBond(atoms1[i + 1], O1, Chem.BondType.DOUBLE)
    mol.AddBond(atoms1[i + 1], O2, Chem.BondType.SINGLE)
    if len(atoms2) > 0:
        mol.AddBond(O2, atoms2[0], Chem.BondType.SINGLE)

    for i in range(length2 - 1):
        mol.AddBond(atoms2[i], atoms2[i + 1], Chem.BondType.SINGLE)

    return mol

def generate_alkyl(length: int):
    mol = Chem.RWMol()
    atoms = [mol.AddAtom(Chem.Atom("C")) for _ in range(length)]
    for i in range(length - 1):
        mol.AddBond(atoms[i], atoms[i + 1], Chem.BondType.SINGLE)

    return mol

def generate_benzene():
    # Build from SMILES — guaranteed valid aromatic ring
    mol = Chem.RWMol(Chem.MolFromSmiles("c1ccccc1"))
    return mol

def get_valid_alkene_bonds(mol: Chem.RWMol):
    valid = []

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        if a1.GetSymbol() == "C" and a2.GetSymbol() == "C":

            bonds_a1 = [bond.GetBondTypeAsDouble() for bond in a1.GetBonds()]
            bonds_a2 = [bond.GetBondTypeAsDouble() for bond in a2.GetBonds()]

            if 2.0 in bonds_a1 or 2.0 in bonds_a2:
                continue

            if sum(bonds_a1) <= 3 and sum(bonds_a2) <= 3:
                valid.append(bond.GetIdx())

    return valid

def get_available_carbons(mol: Chem.RWMol, extra_bonds=1):
    valid = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C":
            bonds = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
            remaining_valence = 4 - bonds
            if remaining_valence >= extra_bonds:
                valid.append(atom.GetIdx())
    return valid

def add_alkyl_branch(mol: Chem.RWMol, max_length: int = 5):

    carbon_indices = get_available_carbons(mol, 1)
    if not carbon_indices:
        return mol

    target = random.choice(carbon_indices)

    for i in range(max_length):
        new_c = mol.AddAtom(Chem.Atom("C"))
        mol.AddBond(target, new_c, Chem.BondType.SINGLE)
        target = new_c

    return mol

def add_alcohol(mol: Chem.RWMol):

    carbon_indices = get_available_carbons(mol, 1)
    if not carbon_indices:
        return mol

    target = random.choice(carbon_indices)

    o_idx = mol.AddAtom(Chem.Atom("O"))
    mol.AddBond(target, o_idx, Chem.BondType.SINGLE)

    return mol

def add_carboxylic(mol: Chem.RWMol):

    carbon_indices = get_available_carbons(mol, 3)
    if not carbon_indices:
        return mol

    target = random.choice(carbon_indices)

    o1_idx = mol.AddAtom(Chem.Atom("O"))
    o2_index = mol.AddAtom(Chem.Atom("O"))

    mol.AddBond(target, o1_idx, Chem.BondType.DOUBLE)
    mol.AddBond(target, o2_index, Chem.BondType.SINGLE)

    return mol

def add_alkene(mol: Chem.RWMol):

    bonds = get_valid_alkene_bonds(mol)
    if not bonds:
        return mol

    bond_idx = random.choice(bonds)
    bond = mol.GetBondWithIdx(bond_idx)

    bond.SetBondType(Chem.BondType.DOUBLE)

    return mol

def add_aldehyde(mol: Chem.RWMol):
    carbon_indices = get_available_carbons(mol, 1)
    if not carbon_indices:
        return mol

    target = random.choice(carbon_indices)

    c_idx = mol.AddAtom(Chem.Atom("C"))  # aldehyde carbon
    o_idx = mol.AddAtom(Chem.Atom("O"))

    mol.AddBond(target, c_idx, Chem.BondType.SINGLE)
    mol.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)

    return mol

def add_ketone(mol: Chem.RWMol):
    carbon_indices = get_available_carbons(mol, 2)
    if not carbon_indices:
        return mol
    target = random.choice(carbon_indices)
    o_idx = mol.AddAtom(Chem.Atom("O"))
    mol.AddBond(target, o_idx, Chem.BondType.DOUBLE)
    return mol

def add_benzene(mol: Chem.RWMol):
    carbon_indices = get_available_carbons(mol, 1)
    if not carbon_indices:
        return mol

    target = random.choice(carbon_indices)
    c_idx = mol.AddAtom(Chem.Atom("C"))
    mol.AddBond(target, c_idx, Chem.BondType.SINGLE)

    atoms = [mol.AddAtom(Chem.Atom("C")) for _ in range(6)]
    for i in range(6):
        mol.AddBond(atoms[i], atoms[(i+1) % 6], Chem.BondType.DOUBLE if i%2==0 else Chem.BondType.SINGLE)

    mol.AddBond(c_idx, atoms[0], Chem.BondType.SINGLE)

    return mol

def add_amine(mol: Chem.RWMol):
    carbon_indices = get_available_carbons(mol, 1)
    if not carbon_indices:
        return mol
    target = random.choice(carbon_indices)
    n_idx = mol.AddAtom(Chem.Atom("N"))
    mol.AddBond(target, n_idx, Chem.BondType.SINGLE)
    return mol

def add_nitrile(mol: Chem.RWMol):
    carbon_indices = get_available_carbons(mol, 3)
    if not carbon_indices:
        return mol
    target = random.choice(carbon_indices)
    n_idx = mol.AddAtom(Chem.Atom("N"))
    mol.AddBond(target, n_idx, Chem.BondType.TRIPLE)
    return mol

def add_acyl_chloride(mol: Chem.RWMol):
    carbon_indices = get_available_carbons(mol, 1)
    if not carbon_indices:
        return mol

    target = random.choice(carbon_indices)

    c_idx = mol.AddAtom(Chem.Atom("C"))  # acyl carbon
    o_idx = mol.AddAtom(Chem.Atom("O"))
    cl_idx = mol.AddAtom(Chem.Atom("Cl"))

    mol.AddBond(target, c_idx, Chem.BondType.SINGLE)
    mol.AddBond(c_idx, o_idx, Chem.BondType.DOUBLE)
    mol.AddBond(c_idx, cl_idx, Chem.BondType.SINGLE)

    return mol

def add_halo(mol: Chem.RWMol):

    code = random.choice(["F", "Cl", "Br", "I"])

    carbon_indices = get_available_carbons(mol, 1)
    if not carbon_indices:
        return mol
    target = random.choice(carbon_indices)
    n_idx = mol.AddAtom(Chem.Atom(code))
    mol.AddBond(target, n_idx, Chem.BondType.SINGLE)
    return mol

def add_nothing(mol: Chem.RWMol):
    return mol

def smiles_to_name(smiles):
    result = smiles_to_name_pubchem(smiles)
    if result is None:
        result = smiles_to_name_cactus(smiles)

    return result
def smiles_to_name_pubchem(smiles):
    encoded = requests.utils.quote(smiles)

    # Try PreferredIUPACName first — more systematic than IUPACName
    for prop in ["PreferredIUPACName", "IUPACName"]:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
            f"{encoded}/property/{prop}/JSON"
        )
        r = requests.get(url, timeout=8)
        if r.status_code != 200:
            continue
        props = r.json().get("PropertyTable", {}).get("Properties", [{}])
        name = props[0].get(prop)
        if name is None:
            continue
        return name.lower()

    return None
def smiles_to_name_cactus(smiles):
    encoded = requests.utils.quote(smiles)
    r = requests.get(
        f"https://cactus.nci.nih.gov/chemical/structure/{encoded}/iupac_name", timeout=8
    )
    if r.status_code == 200 and "Page not found" not in r.text:
        name = r.text.strip()
        return name.lower()

    return None
import subprocess

def opsin_name(smiles):
    try:
        out = subprocess.run(
            ["java", "-jar", "opsin-cli.jar", "-osmi", smiles],
            capture_output=True,
            text=True,
            timeout=5
        )
        return out.stdout.strip().lower()
    except:
        return None

BASE_GROUP_RULES_EASY = {
    "alkyl": [add_alcohol, add_alkene, add_aldehyde, add_ketone, add_amine, add_nitrile, add_halo, add_carboxylic],
    "benzene": [add_alkyl_branch, add_halo, add_amine, add_alcohol],
    "ester": [add_nothing],
    "amide": [add_nothing],
    "acid_anhydride": [add_nothing],
}

BASE_GROUP_RULES_MEDIUM = {
    "alkyl": [add_alcohol, add_alkene, add_amine, add_halo],  # fewer than easy
    "benzene": [add_alkyl_branch, add_halo, add_amine],       # drop alcohol for medium
    "ester": [add_alkyl_branch, add_halo],                    # can add simple branches or halogens
    "amide": [add_alkyl_branch, add_halo, add_alcohol],       # slightly more complex than easy
    "acid_anhydride": [add_alkyl_branch, add_halo],          # only simple additions
}

BASE_GROUP_RULES_HARD = {
    "alkyl": [
        add_alcohol,
        add_alkene,
        add_aldehyde,
        add_ketone,
        add_amine,
        add_nitrile,
        add_halo,
        add_carboxylic,
        add_acyl_chloride,
    ],  # almost all groups allowed
    "benzene": [
        add_alkyl_branch,
        add_halo,
        add_amine,
        add_alcohol,
        add_alkene,
        add_ketone,
    ],  # can decorate benzene heavily
    "ester": [
        add_alkyl_branch,
        add_halo,
        add_alcohol,
        add_alkene,
    ],  # functionalize esters
    "amide": [
        add_alkyl_branch,
        add_halo,
        add_alcohol,
        add_ketone,
        add_alkene,
    ],  # allow multiple decorations
    "acid_anhydride": [
        add_alkyl_branch,
        add_halo,
        add_alkene,
        add_alcohol,
    ],  # slightly more complex than medium
}


def hard():

    lengths = list(range(1, 10))
    length = random.choice(lengths)

    base_functions = {
        "benzene": generate_benzene,
        "ester": partial(generate_ester, max(3, length)),
        "amide": partial(generate_amide, max(3, length)),
        "acid_anhydride": partial(generate_acid_anhydride, max(3, length)),
        "alkyl": partial(generate_alkyl, max(3, length)),
    }
    weights = [
        5,
        5,
        5,
        5,
        5
    ]

    base = random.choices(list(base_functions.keys()), weights=weights)[0]
    mol = base_functions[base]()

    for _ in range(3):
        fn = random.choice(BASE_GROUP_RULES_HARD[base])
        mol = fn(mol)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)

    return mol

def medium():

    lengths = list(range(1, 9))
    length = random.choice(lengths)

    base_functions = {
        "benzene": generate_benzene,
        "ester": partial(generate_ester, max(3, length)),
        "amide": partial(generate_amide, max(3, length)),
        "acid_anhydride": partial(generate_acid_anhydride, max(3, length)),
        "alkyl": partial(generate_alkyl, max(3, length)),
    }
    weights = [
        5,
        4,
        4,
        3,
        8
    ]

    base = random.choices(list(base_functions.keys()), weights=weights)[0]
    mol = base_functions[base]()

    for _ in range(2):
        fn = random.choice(BASE_GROUP_RULES_MEDIUM[base])
        mol = fn(mol)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)

    return mol

def easy():

    lengths = list(range(1, 8))
    length = random.choice(lengths)

    base_functions = {
        "benzene": generate_benzene,
        "ester": partial(generate_ester, max(3, length)),
        "amide": partial(generate_amide, max(3, length)),
        "acid_anhydride": partial(generate_acid_anhydride, max(3, length)),
        "alkyl": partial(generate_alkyl, max(3, length)),
    }
    weights = [
        5,
        2,
        3,
        1,
        10
    ]

    base = random.choices(list(base_functions.keys()), weights=weights)[0]
    mol = base_functions[base]()

    for _ in range(1):
        fn = random.choice(BASE_GROUP_RULES_EASY[base])
        mol = fn(mol)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)

    return mol

