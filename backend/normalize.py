import re

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

def which_starters(string):
    starters = []
    for starter in DEF_STARTERS:
        if starter in string:
            starters.append(starter)

    return starters

def same_starters(answer, correct):
    answer_starters = which_starters(answer)
    correct_starters = which_starters(correct)

    wrong_in_answer = []
    for starter in answer_starters:
        if starter not in correct_starters:
            wrong_in_answer.append(starter)

    if len(wrong_in_answer) == 0:
        return None

    return wrong_in_answer

def list_to_string(items):
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    elif len(items) == 2:
        return f"{items[0]} and {items[1]}"
    else:
        string = ""
        for i in range(len(items) - 1):
            string += f"{items[i]}-, "
        string += f"and {items[-1]}-"

    return string


def check_number(answer, correct):
    digits_answer = ''.join(c for c in answer if c.isdigit()).split()
    digits_correct = ''.join(c for c in correct if c.isdigit()).split()

    digits_correct.sort()
    digits_answer.sort()

    incorrect = []
    for a, c in zip(digits_answer, digits_correct):
        if a != c:
            incorrect.append(a)

    return incorrect if len(incorrect) > 0 else None

def check_prefix_suffix(answer, correct):
    for prefix in DEF_PREFIX_SUFFIX:
        suffix = DEF_PREFIX_SUFFIX[prefix]

        if suffix in answer and suffix not in correct and prefix in correct and prefix not in answer:
            return "prefix", suffix
        elif prefix in answer and prefix not in correct and suffix in correct and suffix not in answer:
            return "suffix", prefix

    return None

def normalize(name: str) -> str:
    s = name.lower().strip()

    # --- basic cleanup ---
    s = s.replace('–', '-')
    s = re.sub(r'[()\s]', '', s)

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

    return s

def get_message(answer, correct):

    message = None

    if re.search(r'-(\d+)-', correct) and not re.search(r'-(\d+)-', answer):
        message = "Which carbon is the functional group on?"
    elif same_starters(answer, correct):
        message = f"Are you sure that the starters \"{list_to_string(same_starters(answer, correct))}\" is correct?"
    elif check_prefix_suffix(answer, correct) is not None:
        message = f"Should you use a {check_prefix_suffix(answer, correct)[0]} for {check_prefix_suffix(answer, correct)[1]}?"
    elif check_number(answer, correct) is not None:
        message = f"Are the carbon numbers {list_to_string(check_number(answer, correct))} correct?"

    return message