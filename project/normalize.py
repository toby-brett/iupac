import random
import re
from inspect import CORO_RUNNING

COMMON_TO_IUPAC = {
    # aldehydes
    "formaldehyde": "methanal",
    "acetaldehyde": "ethanal",

    # ketones
    "acetone": "propanone",

    # acids
    "formicacid": "methanoicacid",
    "aceticacid": "ethanoicacid",

    # alcohols
    "glycerol": "propane-1,2,3-triol",

    # aromatics
    "toluene": "methylbenzene",
    "xylene": "dimethylbenzene",
    "styrene": "ethenylbenzene",
    "aniline": "aminobenzene",

    # s
    "formonitrile": "methanenitrile",
    "acetonitrile": "ethanenitrile",

    # amides
    "acetamide": "ethanamide",

    # others
    "chloroform": "trichloromethane",
    "urea": "carbamide",
    "acetylchloride": "ethanoyl chloride",
    "cyanicacid": "3-hydroxy propane nitrile",
    "1-chloroethylacetate": "1-chloroethyl ethanoate",
    "oxaldehyde": "ethanedial",
    "ethenyl2-aminoacetate": "ethenyl 2-aminoethanoate"
}

FUNCTIONAL_GROUP_FIXES = [
    # turn "methane nitrile" → "methanenitrile"
    (r'(\w+)anenitrile', r'\1anenitrile'),  # already correct safeguard
    (r'(\w+)ane?nitrile', r'\1anenitrile'),

    # "ethane amide" → "ethanamide"
    (r'(\w+)ane?amide', r'\1anamide'),

    # "ethane oic acid" → "ethanoicacid"
    (r'(\w+)ane?oicacid', r'\1anoicacid'),

    # aldehyde spacing issues
    (r'(\w+)anal', r'\1anal'),

    # ketone spacing issues
    (r'(\w+)anone', r'\1anone'),
]

DEF_STARTERS = [
    "meth",
    "eth",
    "prop",
    "but",
    "pent",
    "hex",
    "hept",
    "oct",
    "non",
    "dec"
]

DEF_PREFIX_SUFFIX = {
    "oxo": "one",
    "hydroxy": "ol",
    "amino": "amine",
    "carboxy": "oic acid",
    "cyano": "nitrile",
    "phenyl": "benzene"
}

vowels = ["a", "e", "i", "o", "u"]


def collect(pattern, string, label, kind):
    matches = []
    for m in re.finditer(pattern, string):
        matches.append({
            "type": kind,
            "group": label,
            "match": m.group(0),
            "groups": m.groups(),
            "start": m.start(),
            "end": m.end()
        })
    return matches

def tokenize(string):

    starters = ["meth", "eth", "prop", "but", "pent", "hex", "hept", "oct", "non", "dec", "benzene", "phen"]
    starters = sorted(starters, key=len, reverse=True)
    starters_pattern = "|".join(starters)

    halogens = ["fluoro", "chloro", "bromo", "iodo"]
    halogens = sorted(halogens, key=len, reverse=True)
    halogens_pattern = "|".join(halogens)

    branch_suffix_pattern = rf'({starters_pattern})(ane|ene)$'
    branch_suffix_groups = collect(branch_suffix_pattern, string, "branch", "suffix")

    alcohol_suffix_pattern = (rf"({starters_pattern})"                          # prop
                rf"(an|ane)?"                                           # an
                rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*)(en|ene))?"     # optional -3,4-en
                rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*))?"             # -4,6-
                rf"(di|tri|tetra|penta)?(ol)$")
    alcohol_suffix_groups = collect(alcohol_suffix_pattern, string, "alcohol", "suffix")

    amine_suffix_pattern = ( rf"({starters_pattern})"                   # prop
        rf"(an)?"                                          # an
        rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*)(en|ene))?"   # optional -3,4-en
        rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*))?"         # -4,6-
        rf"(di|tri|tetra|penta)?(amine)$" )
    amine_suffix_groups = collect(amine_suffix_pattern, string, "amine", "suffix")

    aldehyde_suffix_pattern = (rf"({starters_pattern})"  # prop
        rf"(an|ane)?"  # an
        rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*)(en|ene))?"  # optional -3,4-en|ene
        rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*))?"  # -4,6-
        rf"(di|tri|tetra|penta)?(al)$")
    aldehyde_suffix_groups = collect(aldehyde_suffix_pattern, string, "aldehyde", "suffix")

    carboxylic_suffix_pattern = (rf"({starters_pattern})"  # prop
        rf"(an|ane)?"  # an
        rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*)(en|ene))?"  # optional -3,4-en|ene
        rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*))?"     # -4,6-
        rf"(di|tri|tetra|penta)?(oic acid)$")
    carboxylic_suffix_groups = collect(carboxylic_suffix_pattern, string, "carboxylic", "suffix")

    nitrile_suffix_pattern = (rf"({starters_pattern})"  # prop
        rf"(an|ane)?"  # an
        rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*)(en|ene))?"  # optional -3,4-en|ene
        rf"(?:(\s*-?\s*)(\d+(?:,\d+)*)(\s*-?\s*))?"  # -4,6-
        rf"(di|tri|tetra|penta)?(nitrile)$")
    nitrile_suffix_groups = collect(nitrile_suffix_pattern, string, "nitrile", "suffix")

    benzene_suffix_pattern = (r'(?:-(\d+(?:,\d+)*))?(\s*-?\s*)(di|tri|tetra|penta)?(benzene)')   # -5,3-dibenzene
    benzene_suffix_groups = collect(benzene_suffix_pattern, string, "benzene", "suffix")



    branch_prefix_pattern = rf'(?:([\d]+(?:,\d+)*)(\s*-?\s*)?)?(di|tri|tetra|penta)?({starters_pattern})yl'
    branch_prefix_groups = collect(branch_prefix_pattern, string, "methyl", "prefix")

    alcohol_prefix_pattern = ( r'(\d+(?:,\d+)*)(\s*-?\s*)(di|tri|tetra|penta)?(hydroxy)')
    alcohol_prefix_groups = collect(alcohol_prefix_pattern, string, "alcohol", "prefix")

    aldehyde_prefix_pattern = (r'(\d+(?:,\d+)*)(\s*-?\s*)(di|tri|tetra|penta)?(formyl)')
    aldehyde_prefix_groups = collect(aldehyde_prefix_pattern, string, "aldehyde", "prefix")

    benzene_prefix_pattern = (r'(\d+(?:,\d+)*)?(\s*-?\s*)?(di|tri|tetra|penta)?(phenyl)')
    benzene_prefix_groups = collect(benzene_prefix_pattern, string, "benzene", "prefix")

    amine_prefix_pattern = (r'(\d+(?:,\d+)*)?(\s*-?\s*)?(di|tri|tetra|penta)?(amino)')
    amine_prefix_groups = collect(amine_prefix_pattern, string, "amine", "prefix")

    carboxylic_prefix_pattern = (r'(\d+(?:,\d+)*)?(\s*-?\s*)?(di|tri|tetra|penta)?(carboxy)')
    carboxylic_prefix_groups = collect(carboxylic_prefix_pattern, string, "carboxylic acid", "prefix")

    nitrile_prefix_pattern = (r'(\d+(?:,\d+)*)?(-)?(di|tri|tetra|penta)?cyano')
    nitrile_prefix_groups = collect(nitrile_prefix_pattern, string, "nitrile", "prefix")

    nitro_prefix_pattern = (r'(\d+(?:,\d+)*)?(-)?(di|tri|tetra|penta)?nitro')
    nitro_prefix_groups = collect(nitro_prefix_pattern, string, "nitro", "prefix")

    halo_prefix_pattern = (fr'(\d+(?:,\d+)*)?(-)?(di|tri|tetra|penta)?({halogens_pattern})')
    halo_prefix_groups = collect(halo_prefix_pattern, string, "halo", "prefix")

    return {
        "alkyl": {"suffix": branch_suffix_groups, "prefix": branch_prefix_groups},
        "alcohol": {"suffix": alcohol_suffix_groups, "prefix": alcohol_prefix_groups},
        "amine": {"suffix": amine_suffix_groups, "prefix": amine_prefix_groups},
        "aldehyde": {"suffix": aldehyde_suffix_groups, "prefix": aldehyde_prefix_groups},
        "carboxylic": {"suffix": carboxylic_suffix_groups, "prefix": carboxylic_prefix_groups},
        "nitrile": {"suffix": nitrile_suffix_groups, "prefix": nitrile_prefix_groups},
        "benzene": {"suffix": benzene_suffix_groups, "prefix": benzene_prefix_groups},
        "halo": {"suffix": [], "prefix": halo_prefix_groups},
        "nitro": {"suffix": [], "prefix": nitro_prefix_groups}
    }


def get_orders(tokenised):

    i = 0
    for group, parts in sorted([item for item in tokenised.items() if item[1]["prefix"]], key=lambda x: x[1]["prefix"][0]["start"]):
        i += 1
        parts["prefix"][0]["start"] = i

    return tokenised


def generate_cross_hints(answer_tokenised, correct_tokenised):
    """
    Generates hints by comparing groups to themselves and each other

    :param answer_tokenised:
    :param correct_tokenised:
    :return:
    """
    hints = []

    for group, correct_parts in correct_tokenised.items():
        answer_parts = answer_tokenised[group]
        # iterating through each group

        answered_suffix = answer_parts["suffix"]
        correct_suffix = correct_parts["suffix"]
        answered_prefix = answer_parts["prefix"]
        correct_prefix = correct_parts["prefix"]

        # checks for suffix when prefix should have been used
        if answered_suffix and correct_prefix and not correct_suffix:
            hints.append((f"Are you sure {group} is the main group?", 2))

        # checks for prefix when suffix should have been used
        if answered_prefix and correct_suffix and not correct_prefix:
            hints.append((f"Are you sure {group} isn't the main group?", 2))

        # checks the suffix hints
        if answered_suffix and correct_suffix and answered_suffix != correct_suffix:

            answered_suffix_tokens = answered_suffix[0]["groups"]
            correct_suffix_tokens = correct_suffix[0]["groups"]

            i = 1
            for token_answered, token_correct in zip(answered_suffix_tokens, correct_suffix_tokens):
                if token_answered == token_correct:
                    i += 1
                    continue
                if i == 1:
                    hints.append((f"Double count your carbons", 1))
                elif i == 2:
                    hints.append((f"You need \"{token_correct}\" not \"{token_answered}\" on your {group}", 3))
                elif i in (3, 5, 7, 9):
                    hints.append((f"Don't forget dashes/spaces", 4))
                elif i == 6:
                    hints.append((f"Is \"{token_answered}\" correct?", 5))
                elif i in (4, 8):
                    hints.append((f"Are your carbon numbers correct?", 3))
                elif i == 10:
                    hints.append((f"How many {group} groups are there?", 4))
                elif i == 11:
                    hints.append((f"Is that the correct suffix for {"an" if group[0].lower() in vowels else "a"} {group} group?", 5))
                i += 1

        # checks the prefix hints
        if answered_prefix and correct_prefix and answered_prefix != correct_prefix:
            answered_prefix_tokens = answered_prefix[0]["groups"]
            correct_prefix_tokens = correct_prefix[0]["groups"]

            answered_prefix_position = answered_prefix[0]["start"]
            correct_prefix_position = correct_prefix[0]["start"]

            if answered_prefix_position != correct_prefix_position:
                hints.append(("Are your prefixes in alphabetical order?", 5))

            i = 1
            for token_answered, token_correct in zip(answered_prefix_tokens, correct_prefix_tokens):
                if token_answered == token_correct:
                    i += 1
                    continue

                if i == 1:
                    hints.append((f"Which carbon is you {group} group on?", 2))
                elif i == 2:
                    hints.append(("Don't forget dashes/spaces", 4))
                elif i == 3:
                    hints.append((f"How many {group} groups are there?", 4))
                elif i == 4:
                    hints.append((f"Is {token_answered} correct?", 5))

                i += 1

        # if groups exist that arn't present at all
        if (answered_prefix or answered_suffix) and (not correct_prefix and not correct_suffix):
            hints.append((f"Is {group} present here?", 1))

    if not hints:
        hints.append(("Double check everything", 0))

    return hints


def pick_hint(answer, correct):

    tokenized_answer = get_orders(tokenize(answer))
    tokenized_correct = get_orders(tokenize(correct))
    hints = generate_cross_hints(tokenized_answer, tokenized_correct)
    options = {}
    lowest_score = 100
    for hint in hints:
        if int(hint[1]) < lowest_score:
            lowest_score = hint[1]
        if int(hint[1]) not in options.keys():
            options[hint[1]] = []
        options[hint[1]].append(hint[0])

    return random.choice(options[lowest_score])


def normalize(name: str) -> str:

    print("before:", name)

    s = name.lower().strip()

    # --- basic cleanup ---
    s = s.replace('–', '-')

    # --- replace common names FIRST ---
    for k, v in COMMON_TO_IUPAC.items():
        if k in s:
            s = s.replace(k, v)

    # --- fix functional group formatting ---
    for pattern, repl in FUNCTIONAL_GROUP_FIXES:
        s = re.sub(pattern, repl, s)

    # --- enforce nitrile rule explicitly ---
    # methane nitrile → methanenitrile
    s = re.sub(r'(\w+)ane?nitrile', lambda m: m.group(1) + 'anenitrile', s)

    # --- remove duplicate hyphens ---
    s = re.sub(r'-+', '-', s)

    print("after:", s)

    return s

