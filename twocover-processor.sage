#! /usr/bin/env sage

"""
Here's the JSON data scheme that this script produces:
{
    "label": <str>, # the curve's label
    "conductor": <int>, # the curve's conductor,
    "disc": <int>, # the absolute value of the discriminant of the curve
    "coeffs": [<int>], # a list of the coefficients of the sextic polynomial f such that y^2 = f(x) defines the curve
    "pts": [[x, y, z]], # a list of the (projective) points on the curve
    "count": <int>, # the number of points
    "root": <int>, # a chosen root of the above sextic
    "g": [[<int>]], # a list whose entries are lists of the coefficients of the factors of f(x)/(x - root) over Q
    "x-coords": [[numerator, denominator]], # list of x-coordinates of points on the curve
    "verified": <bool>, # is the above list verified complete (conditional on GRH)?
    "obstruction_found": <bool>, # have we found an obstruction to the algorithm completing?
    "search_bound": <int>, # how far have we searched before for rational points on the twists?
    "hasse_principle": <bool>, # do the curve's twists all appear to satisfy the Hasse principle?
    "class_group_fields": [[[<int>], <int>, <bool>]], # list of (field, class number, unconditional?) triples needed to prove rank bounds
    "twists": [
        {
            "coeffs": [str], # list of the coefficients of the twist parameter in the form "num/denom"
            "found_pts": <int>, # number of points that have been found
            "base_pt": [<int>], # chosen base point of the twist
            "loc_solv": <bool>, # is the twist locally solvable?
            "pts": [[<int>]], # list of coordinates of points
            "verified": <bool>, # have we verified that the list of points is complete (conditional on GRH)?
            "g1": [
                {
                    "g": [<int>], # list of coefficients of the corresponding factor of f
                    "aInv": [[<str>]], # a-invariant data (the strings are in the form "num/denom")
                    "Ecov": [[[[<str>]], [[<int>]]]], # equation data for map Ecov: E -> P^1
                    "rank": <int>, # computed rank of MW group
                    "MW_proven": <bool>, # has the MW group been provably computed (conditional on GRH)?
                    "MW_orders": [<int>], # orders of generators of MW group
                    "gens": [[[str]]] # generators of MW group (str in the form "num/denom")
                    "gens_reduced": bool # whether the MW generators have already been reduced
                    "x-coords": [[numerator, denominator]] # list of x-coordinates found by elliptic Chabauty
                }
            ]
        }
    ],
    "runtime": <float>, # estimated runtime in seconds taken so far to compute this data
    "exception": False # did an unhandled exception occur?
}
"""

import argparse
import json
import logging
import resource
import time

# Set up command line options
parser = argparse.ArgumentParser()
label_group = parser.add_mutually_exclusive_group()
label_group.add_argument("--index", type=int, help="the index of the LMFDB label of the curve to process")
label_group.add_argument("--label", help="the LMFDB label of the curve to process")
parser.add_argument("--database", help="the database file to read from", default="data/g2c_curves-r2-w1.json")
parser.add_argument("--label_list", help="the list of labels the index is based on", default="data/labels.json")
parser.add_argument("--output_directory", help="directory for output files", default="./results")
parser.add_argument("--stages", help="'all' or comma-separated list from: 'search', 'locsolv', 'ainv', 'map', 'class', 'mw', 'reduce', 'chabauty'", default="all")
parser.add_argument("--missing", help="use list of labels with no preexisting file", action="store_true")
parser.add_argument("--force", help="Run regardless of whether points or obstructions already found", action="store_true")
parser.add_argument("--unconditional", help="Generate Magma file that, if run, makes the rank bounds unconditional", action="store_true")
args = parser.parse_args()

# Identify label of curve
if args.label is not None:
    LABEL = args.label
elif args.index is not None:
    with open(args.label_list, "r") as f:
        LABEL_LIST = json.load(f)
    if args.missing:
        # This mode is for filling in missing curves that were skipped for whatever reason.
        files_in_output_dir = os.listdir(args.output_directory)
        LABEL_LIST = [label for label in LABEL_LIST if "curve-{}.json".format(label) not in files_in_output_dir]
    if args.index < 0 or args.index >= len(LABEL_LIST):
        print("Index {} out of bounds; exiting.".format(args.index))
        exit()
    LABEL = LABEL_LIST[args.index]
else:
    raise ValueError("Must provide a label or label index.")

DATA_FILE = args.database
if args.force:
    HALT_ON_OBSTRUCTION = False
else:
    HALT_ON_OBSTRUCTION = True # stop immediately if obstruction found

class Obstruction(Exception):
    """Exception thrown when an obstruction to success of the algorithm is found"""
    pass

class PointCountError(Exception):
    """Exception thrown when the computed point count doesn't agree with the known point count"""
    pass

# The computation is broken up into stages, which can be selected using the command-line option --stages
ALL_STAGES = {"search", "locsolv", "ainv", "map", "class", "mw", "reduce", "chabauty"}
if args.stages.lower() == "all":
    STAGES = ALL_STAGES
else:
    STAGES = frozenset(args.stages.lower().split(","))
    if not STAGES.issubset(ALL_STAGES):
        raise ValueError("Invalid stage label.")

# Load the Magma code that implements the constructions of the paper
magma.load("twocovers.m")

# Limit memory usage to 8 GB
MEMORY_LIMIT = 8 * 1024 * 1024 * 1024
resource.setrlimit(resource.RLIMIT_AS, (MEMORY_LIMIT, MEMORY_LIMIT))
resource.setrlimit(resource.RLIMIT_RSS, (MEMORY_LIMIT, MEMORY_LIMIT))
magma.eval("SetMemoryLimit({})".format(2 * MEMORY_LIMIT))

# Bounds on how far we should search for points on the main curve and on the twists
CURVE_SEARCH_BOUND = 10000
TWIST_SEARCH_BOUND = 100000

# Filenames where the output will be written
OUTPUT_FILE_PREFIX = "{}/curve-{}".format(args.output_directory, LABEL)
OUTPUT_FILE = OUTPUT_FILE_PREFIX + ".json"
LOG_FILE = OUTPUT_FILE_PREFIX + ".log"
CLASS_GROUP_FILE = OUTPUT_FILE_PREFIX + "-class-groups.m"

# Set up log file
logging.basicConfig(
    filename=LOG_FILE,
    encoding="utf-8",
    format="%(asctime)s %(levelname)s: %(message)s",
    level=logging.DEBUG
)

def integral_proj_pt(P):
    """Takes a list of rational numbers, clears denominators, returns list of ints"""
    denom = LCM([QQ(c).denominator() for c in P])
    return [int(ZZ(c * denom)) for c in P]

def build_curve_data(label, poly_coeffs):
    """Initializes curve data from the given label and coefficients"""
    R.<x> = QQ[]
    C_mag = magma.HyperellipticCurve(R(poly_coeffs[0]), R(poly_coeffs[1]))
    C_even_mag = magma.function_call("even_weierstrass_model", C_mag)
    sextic = magma.HyperellipticPolynomials(C_even_mag)
    coeffs = [int(ZZ(c)) for c in magma.Eltseq(sextic)]
    f = R(coeffs)
    roots = [r[0] for r in f.roots()]
    # integer root that's smallest in absolute value
    int_roots = [r for r in roots if r.is_integer()]
    if len(int_roots) > 0:
        root = int(ZZ(min(int_roots, key=lambda x: abs(x))))
    else:
        # handle case where there's no integral Weierstrass point
        denom = min(r.denominator() for r in roots)
        f = denom^6 * f(x / denom)
        f /= gcd(list(f))
        coeffs = [int(ZZ(c)) for c in list(f)]
        root = int(ZZ(min((r[0] for r in f.roots() if r[0].is_integer()), key=lambda x: abs(x))))
    # factors of f(x)/(x - root), ordered by degree
    g_factors = sorted([p[0] for p in (f / (x - QQ(root))).factor()], key=lambda p: p.degree())
    return {
        "label": label,
        "conductor": int(label.split(".")[0]),
        "disc": int(label.split(".")[2]),
        "coeffs": coeffs,
        "pts": None,
        "count": None,
        "root": root,
        "g": [integral_proj_pt(list(p)) for p in g_factors],
        "x-coords": None,
        "verified": False,
        "hasse_principle": None,
        "obstruction_found": False,
        "search_bound": int(0),
        "class_group_fields": None,
        "twists": None,
        "runtime": 0,
        "exception": False
    }

def twist_coeffs(curve):
    """Coefficients of the twist parameters associated to the curve"""
    R.<x> = QQ[]
    f = R(curve["coeffs"])
    gens = magma.function_call("twist_param_generators", f)
    twists = list(magma.function_call("products_of_subsets", gens))
    return [[str(QQ(c)) for c in magma.Eltseq(d)] for d in twists]

def twist_data(curve):
    """Builds initial data structure for the twists"""
    if curve["twists"] is not None:
        # already computed
        return curve["twists"]
    coeff_data = twist_coeffs(curve)
    return [{
        "coeffs": d,
        "found_pts": None,
        "base_pt": None,
        "loc_solv": None,
        "pts": None,
        "verified": False,
        "g1": [{
                "g": g,
                "aInv": None,
                "Ecov": None,
                "rank": None,
                "MW_proven": None,
                "gens": None,
                "gens_reduced": None,
                "x-coords": None
            } for g in curve["g"]]
    } for d in coeff_data]

def twist_from_coeffs(coeff_data):
    """Turns the coefficient data for a delta-value back into a polynomial"""
    R.<x> = QQ[]
    return R([QQ(c) for c in coeff_data])

def twist_point_search(curve, twist_index, bound=TWIST_SEARCH_BOUND):
    """Search for points on the twist. Skips search if point already found or if not locally solvable."""
    R.<x> = QQ[]
    f = R(curve["coeffs"])
    root = R(curve["root"])
    assert len(curve["twists"]) >= 1
    twist = curve["twists"][twist_index]
    if twist["found_pts"] is not None:
        found_pts = twist["found_pts"]
        base_pt = twist["base_pt"]
    elif twist["loc_solv"] is False:
        found_pts = 0
        base_pt = None
    else:
        delta = twist_from_coeffs(twist["coeffs"])
        Z = magma.function_call("genus5_canonical", [f, root, delta])
        pts = magma.PointSearch(Z, bound)
        found_pts = len(pts)
        if len(pts) == 0:
            base_pt = None
        else:
            base_pt_mag = sorted(list(pts), key=lambda P: int(ZZ(magma.HeightOnAmbient(P))))[0]
            base_pt = integral_proj_pt(base_pt_mag.Eltseq())
    return found_pts, base_pt

def twists_locally_solvable(curve, twist_index):
    """Wrapper for Magma function testing local solvability"""
    R.<x> = QQ[]
    f = R(curve["coeffs"])
    root = R(curve["root"])
    bad_primes = [p[0] for p in ZZ(curve["conductor"]).factor()]
    twist = curve["twists"][twist_index]
    if twist["loc_solv"] is not None:
        return twist["loc_solv"]
    elif twist["base_pt"] is not None:
        return True
    else:
        delta = twist_from_coeffs(twist["coeffs"])
        Z = magma.function_call("genus5_canonical", [f, root, delta])
        # we test 2 last because this causes problems for some twists
        ls = magma.function_call("is_genus5_locally_solvable", [Z, bad_primes], {"Skip": [2]})
        # hack because conversion of Magma booleans is broken
        ls = (ls == magma(True))
        if ls:
            ls2 = magma.function_call("is_genus5_locally_solvable_at_2", Z)
            return (ls2 == magma(True))
        else:
            return False

def hasse_principle(curve):
    """Do all the twists satisfy the Hasse principle? Returns None if unknown."""
    for twist in curve["twists"]:
        if twist["found_pts"] is None or twist["loc_solv"] is None:
            return None
        elif twist["found_pts"] == 0 and twist["loc_solv"]:
            return False
    return True

def get_aInv_data(curve, twist_index, g_index):
    """Extracts elliptic curve a-invariants and returns them in JSON-compatible format"""
    twist = curve["twists"][twist_index]
    R.<x> = QQ[]
    f = R(curve["coeffs"])
    root = QQ(curve["root"])
    delta = twist_from_coeffs(twist["coeffs"])
    D = twist["g1"][g_index]
    g = R(D["g"])
    E = magma.function_call("twist_ell_curve", [f, root, g, delta])
    aInvs = E.aInvariants()
    return [[str(QQ(c)) for c in a.Eltseq()] for a in aInvs]

def aInv_data_to_ell_curve(g, aInv_data):
    """Builds elliptic curve in Magma with given a-invariant data"""
    K = magma.NumberField(g)
    aInvs = [K([QQ(c) for c in a]) for a in aInv_data]
    E = magma.EllipticCurve(aInvs)
    return E

def get_class_group_fields(D):
    """Computes (field, class number) pairs needed for unconditional MW computation."""
    R.<x> = QQ[]
    g = R(D["g"])
    E = aInv_data_to_ell_curve(g, D["aInv"])
    cgf = [(c[1], len(c[2])) for c in magma.function_call("two_selmer_class_groups", [E])]
    return cgf

def extract_MW_gens(mw, A):
    """Extract MW generators from map mw: A -> E(K)"""
    gens_mag = [magma("{}({})".format(mw.name(), a.name())) for a in A.gens()]
    return gens_mag

def parse_MW_gens(gens_mag):
    """Parse generators of MW group into JSON-compatible format"""
    gens = []
    for P in gens_mag:
        coords = [[str(QQ(c)) for c in a.Eltseq()] for a in P.Eltseq()]
        gens.append(coords)
    return gens

def build_MW_gens(D):
    """Given generator data in JSON-compatible format, build MW generators as Magma objects"""
    R.<x> = QQ[]
    g = R(D["g"])
    E = aInv_data_to_ell_curve(g, D["aInv"])
    K = E.BaseRing()
    gens_mag = [E([K(coord) for coord in gen]) for gen in D["gens"]]
    return gens_mag

def reduce_MW_gens(gens_mag):
    """Given a list of MW generators, apply LLL reduction to get smaller generators"""
    indep_gens = []
    torsion_gens = []
    for P in gens_mag:
        if P.Order() == 0:
            indep_gens.append(P)
        else:
            torsion_gens.append(P)
    r = len(indep_gens)
    if r <= 1:
        return gens_mag
    M = magma.HeightPairingMatrix(indep_gens)
    I = magma.ScalarMatrix(r, M.BaseRing()(1))
    L = magma.Lattice(I, M)
    L_, T = L.BasisReduction(nvals=2)
    new_indep_gens = []
    for i in range(r):
        P = indep_gens[0].Curve().Identity()
        for j in range(r):
            P += T[i+1][j+1] * indep_gens[j]
        new_indep_gens.append(P)
    return torsion_gens + new_indep_gens

def twist_MW_group(D):
    """Wrapper for Mordell-Weil group computation"""
    R.<x> = QQ[]
    g = R(D["g"])
    assert D["aInv"] is not None
    E = aInv_data_to_ell_curve(g, D["aInv"])
    magma.eval("SetClassGroupBounds(\"GRH\");")
    A, mw, bool_rank, bool_index = magma.MordellWeilGroup(E, nvals=4)
    rank = int(magma.TorsionFreeRank(A))
    MW_proven = (bool_rank == bool_index == magma(True))
    orders = [int(gen.Order()) for gen in A.gens()]
    gens_mag = extract_MW_gens(mw, A)
    gens = parse_MW_gens(gens_mag)
    return rank, MW_proven, orders, gens

def rank_degree_ordering(D):
    """Sorting key for ordering first by rank, then by degree"""
    if D["rank"] is not None:
        rank = D["rank"]
    else:
        rank = 999 # arbitrary large number to make this sorted last
    degree = len(D["g"]) - 1
    return (rank, degree)

def twist_chabauty_map(curve, twist_index, g_index):
    """Wrapper for Magma function computing the map Ecov: E -> P^1"""
    twist = curve["twists"][twist_index]
    R.<x> = QQ[]
    f = R(curve["coeffs"])
    root = QQ(curve["root"])
    delta = twist_from_coeffs(twist["coeffs"])
    D = twist["g1"][g_index]
    g = R(D["g"])
    Ecov = magma.function_call("twist_chabauty_map", [f, root, g, delta, twist["base_pt"], D["aInv"]])
    return magma.function_call("pack_map", [Ecov]).sage()

def twist_chabauty(curve, twist_index, g_index):
    """Wrapper for Magma function that runs elliptic Chabauty on the twist"""
    twist = curve["twists"][twist_index]
    R.<x> = QQ[]
    f = R(curve["coeffs"])
    root = QQ(curve["root"])
    delta = twist_from_coeffs(twist["coeffs"])
    D = twist["g1"][g_index]
    g = R(D["g"])
    assert D["MW_proven"] and D["rank"] < g.degree()
    x_coords, pts = magma.function_call("twist_chabauty",
            [f, root, g, delta, D["aInv"], D["Ecov"], D["MW_orders"], D["gens"]], nvals=2)
    return tuple([integral_proj_pt(P.Eltseq()) for P in pt_set] for pt_set in (x_coords, pts))

def x_coords_of_twist_pts(curve):
    """Finds all points on P^1 associated with points on twists whose point-sets have been verified"""
    R.<x> = QQ[]
    f = R(curve["coeffs"])
    root = curve["root"]
    coords = []
    verified = True
    for twist in curve["twists"]:
        if not twist["verified"]:
            verified = False
            continue
        delta = twist_from_coeffs(twist["coeffs"])
        dup = magma.function_call("twisted_duplication_map", [f, root, delta])
        dup_x, dup_z = dup.DefiningEquations().sage()
        for P in twist["pts"]:
            x, z = dup_x(P), dup_z(P)
            if z < 0:
                x, z = -x, -z # normalize signs
            coords.append([int(ZZ(x / gcd(x, z))), int(ZZ(z / gcd(x, z)))])
    return list(set([tuple(c) for c in coords])), verified

def lift_to_hyperelliptic_curve(f_coeffs, x, z):
    """Return all (x : y : z) such that y^2 = f(x, z)"""
    R.<t> = QQ[]
    f = R(f_coeffs)
    c = f.leading_coefficient()
    x = ZZ(x)
    z = ZZ(z)
    weight = ceil(f.degree() / 2)
    if z == 0:
        if f.degree() % 2 == 1:
            lifts = [[1, 0, 0]]
        elif c.is_square():
            lifts = [[1, sqrt(c), 0], [1, -sqrt(c), 0]]
        else:
            lifts = []
    else:
        y2 = f(x / z)
        if y2 == 0:
            lifts = [[x, 0, z]]
        elif y2.is_square():
            lifts = [[x, sqrt(y2) * z^weight, z], [x, -sqrt(y2) * z^weight, z]]
        else:
            lifts = []
    return [[int(ZZ(n)) for n in lift] for lift in lifts]

def known_point_count(curve, bound=CURVE_SEARCH_BOUND):
    """Number of points on the curve up to the given bound"""
    R.<x> = QQ[]
    f = R(curve["coeffs"])
    return len(magma.HyperellipticCurve(f).Points(Bound=bound))

def build_class_group_code(class_group_fields):
    """Generate Magma file to unconditionally verify class groups"""
    # Generate Magma file needed to make class group computations unconditional
    code_lines = [
            'SetVerbose("ClassGroup", 3);',
            'R<x> := PolynomialRing(Rationals());'
        ]
    for c in class_group_fields:
        if c[2]:
            # Skip class groups already computed unconditionally
            continue
        code_lines += [
                'K := NumberField(R!{}); K;'.format([QQ(n) for n in c[0]]),
                'C := ClassGroup(K : Proof := "Full"); C;',
                'assert #C eq {};'.format(int(c[1])),
                'printf "Verified class number of %o is %o.\\n\\n", K, #C;'
            ]
    code_lines.append('exit;')
    return "\n".join(code_lines)

class SageJSONEncoder(json.JSONEncoder):
    """Extends the default JSON encoder to handle Sage's Integer and Rational types"""
    def default(self, o):
        if type(o) is Integer:
            return int(o)
        elif type(o) is Rational:
            return str(o)
        else:
            return json.JSONEncoder.default(self, o)

def record_data(data, filename, timestamp, final=False):
    """Write data as JSON to the given file, updating the runtime with time since the provided timestamp.
    Returns a new timestamp so this can be called repeatedly to provide good estimates of total runtime."""
    t = time.time()
    data["runtime"] += t - timestamp
    s = json.dumps(data, cls=SageJSONEncoder)
    with open(filename, "w") as f:
        f.write(s)
    if HALT_ON_OBSTRUCTION and data["obstruction_found"] and not final:
        raise Obstruction
    return t

# Start recording runtime
t = time.time()

# Try to load existing curve data, and build basic curve data if not found
try:
    with open(OUTPUT_FILE, "r") as f:
        curve = json.load(f)
except FileNotFoundError:
    with open(DATA_FILE, "r") as database:
        curve_database = json.load(database)
    # Build the basic data structure recording information about the curve
    curve = build_curve_data(LABEL, curve_database[LABEL])
    logging.info("Initialized curve data.")
    t = record_data(curve, OUTPUT_FILE, t)
else:
    assert curve["label"] == LABEL
    logging.info("Loaded curve data from file.")
    if not args.force:
        if curve["verified"]:
            logging.info("Rational points already verified; exiting.")
            exit()
        elif HALT_ON_OBSTRUCTION and curve["obstruction_found"]:
            logging.info("Obstruction already found. Exiting.")
            exit()

try:
    if curve["twists"] is None:
        # Compute coefficients of twist parameters
        logging.info("Setting up initial twist data...")
        curve["twists"] = twist_data(curve)
        t = record_data(curve, OUTPUT_FILE, t)
        logging.info("Twist setup complete.")

    if "search" in STAGES and ("search_bound" not in curve or curve["search_bound"] < TWIST_SEARCH_BOUND):
        # Search for points on each twist, and choose a base point
        logging.info("Searching for rational points on each twist...")
        for i in range(len(curve["twists"])):
            twist = curve["twists"][i]
            if twist["base_pt"] is not None or twist["verified"]:
                continue
            found_pts, base_pt = twist_point_search(curve, twist_index=i, bound=TWIST_SEARCH_BOUND)
            twist["found_pts"] = found_pts
            twist["base_pt"] = base_pt
        curve["search_bound"] = int(TWIST_SEARCH_BOUND)
        t = record_data(curve, OUTPUT_FILE, t)
        logging.info("Finished searching for rational points on each twist.")

    if "locsolv" in STAGES and any(twist["loc_solv"] is None for twist in curve["twists"]):
        # Test whether the twists are locally solvable
        logging.info("Testing local solvability of twists...")
        for i in range(len(curve["twists"])):
            twist = curve["twists"][i]
            if twist["loc_solv"] is not None:
                continue
            loc_solv = twists_locally_solvable(curve, twist_index=i)
            twist["loc_solv"] = loc_solv
            if not loc_solv:
                twist["pts"] = []
                twist["verified"] = True
            t = record_data(curve, OUTPUT_FILE, t)
        logging.info("Finished local solvability testing.")

        # Record whether all the twists appear to satisfy the Hasse principle
        curve["hasse_principle"] = hasse_principle(curve)
        if not curve["hasse_principle"]:
            curve["obstruction_found"] = True
        t = record_data(curve, OUTPUT_FILE, t)

    if "ainv" in STAGES and any(twist["base_pt"] is not None and any(D["aInv"] is None for D in twist["g1"]) for twist in curve["twists"]):
        # Compute a-invariants of the elliptic curve associated to each twist with found points
        logging.info("Computing elliptic curve a-invariants...")
        for i in range(len(curve["twists"])):
            twist = curve["twists"][i]
            if twist["base_pt"] is None:
                continue
            for j in range(len(curve["g"])):
                D = twist["g1"][j]
                if D["aInv"] is None:
                    D["aInv"] = get_aInv_data(curve, twist_index=i, g_index=j)
                    t = record_data(curve, OUTPUT_FILE, t)
        logging.info("Finished computation of a-invariants.")
        t = record_data(curve, OUTPUT_FILE, t)

    if "map" in STAGES:
        # Compute defining equations of elliptic Chabauty map E -> P^1
        for i in range(len(curve["twists"])):
            twist = curve["twists"][i]
            if twist["base_pt"] is None:
                continue
            for j in range(len(twist["g1"])):
                D = twist["g1"][j]
                if D["aInv"] is not None and ("Ecov" not in D or D["Ecov"] is None):
                    logging.info("Computing elliptic Chabauty map... (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                    D["Ecov"] = twist_chabauty_map(curve, twist_index=i, g_index=j)
                    logging.info("Elliptic Chabauty map computed. (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                    t = record_data(curve, OUTPUT_FILE, t)

    if "class" in STAGES and ("class_group_fields" not in curve or curve["class_group_fields"] is None):
        # Compute (field, class number) pairs needed to make MW computation unconditional
        logging.info("Conditionally computing class numbers needed for MW computation...")
        class_group_fields = []
        for twist in curve["twists"]:
            if twist["base_pt"] is None:
                continue
            for D in twist["g1"]:
                for c in get_class_group_fields(D):
                    # again have to use a hack for conversion of Magma booleans--I think this is fixed in Sage 9.4
                    if not any((magma.IsIsomorphic(c[0], cgf[0]) == magma(True)) for cgf in class_group_fields):
                        class_group_fields.append(c)
        curve["class_group_fields"] = [[c[0].DefiningPolynomial().Eltseq().sage(), c[1], False] for c in class_group_fields]
        t = record_data(curve, OUTPUT_FILE, t)
        logging.info("Class numbers computed and recorded.")

    if "mw" in STAGES:
        # Compute the Mordell-Weil group of each twist where we found a base point
        for i in range(len(curve["twists"])):
            twist = curve["twists"][i]
            if twist["base_pt"] is None:
                logging.info("Skipped MW computation because no base point (delta = {})".format(twist["coeffs"]))
                continue
            for j in range(len(curve["g"])):
                D = twist["g1"][j]
                if D["gens"] is not None:
                    logging.info("MW generators already computed (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                    continue
                logging.info("Computing MW group (delta = {}, g = {})...".format(twist["coeffs"], D["g"]))
                rank, MW_proven, orders, gens = twist_MW_group(D)
                logging.info("Finished MW computation (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                D["rank"] = rank
                D["MW_proven"] = MW_proven
                D["MW_orders"] = orders
                D["gens"] = gens
                D["gens_reduced"] = False
                t = record_data(curve, OUTPUT_FILE, t)
            if all(D["rank"] is not None and D["rank"] >= len(D["g"]) - 1 for D in twist["g1"]):
                logging.info("Unable to use elliptic Chabauty on twist (delta = {})".format(twist["coeffs"]))
                curve["obstruction_found"] = True
                t = record_data(curve, OUTPUT_FILE, t)
        t = record_data(curve, OUTPUT_FILE, t)

    if "reduce" in STAGES and any(
            twist["base_pt"] is not None
            and not twist["verified"]
            and not all(D["gens_reduced"] for D in twist["g1"])
            for twist in curve["twists"]):
        # Use LLL reduction to find smaller MW generators
        logging.info("Reducing MW generators...")
        for twist in curve["twists"]:
            if twist["base_pt"] is None or twist["verified"]:
                continue
            for D in twist["g1"]:
                if D["gens"] is None:
                    continue
                elif "gens_reduced" not in D:
                    D["gens_reduced"] = False
                elif D["gens_reduced"]:
                    continue
                gens_mag = build_MW_gens(D)
                reduced_gens = reduce_MW_gens(gens_mag)
                D["gens"] = parse_MW_gens(reduced_gens)
                D["gens_reduced"] = True
        logging.info("Finished reducing MW generators.")
        t = record_data(curve, OUTPUT_FILE, t)

    if "chabauty" in STAGES:
        # Check for obstructions due to MW ranks being too high
        for twist in curve["twists"]:
            if twist["base_pt"] is not None and all(D["rank"] is not None and D["rank"] >= len(D["g"]) - 1 for D in twist["g1"]):
                logging.info("Unable to use elliptic Chabauty on twist (delta = {})".format(twist["coeffs"]))
                curve["obstruction_found"] = True
                t = record_data(curve, OUTPUT_FILE, t)
        # Run elliptic Chabauty on the twists, where possible
        for i in range(len(curve["twists"])):
            twist = curve["twists"][i]
            if twist["base_pt"] is None:
                logging.info("Skipped Chabauty because no base point (delta = {})".format(twist["coeffs"]))
                continue
            g_index_order = sorted(range(len(curve["g"])), key=lambda j: rank_degree_ordering(twist["g1"][j]))
            for j in g_index_order:
                D = twist["g1"][j]
                if twist["verified"]:
                    logging.info("Skipped Chabauty because twist points already verified. (delta = {})".format(twist["coeffs"]))
                    break
                elif D["gens"] is None:
                    logging.info("Skipped Chabauty because MW generators haven't been computed. (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                    continue
                elif not D["MW_proven"]:
                    logging.info("Skipped Chabauty because MW group hasn't been provably computed. (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                    continue
                elif D["rank"] >= len(D["g"]) - 1:
                    logging.info("Skipped Chabauty because rank >= degree(g). (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                    continue
                if "Ecov" not in D or D["Ecov"] is None:
                    logging.info("Computing elliptic Chabauty map... (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                    D["Ecov"] = twist_chabauty_map(curve, twist_index=i, g_index=j)
                    logging.info("Elliptic Chabauty map computed. (delta = {}, g = {})".format(twist["coeffs"], D["g"]))
                    t = record_data(curve, OUTPUT_FILE, t)
                logging.info("Running elliptic Chabauty (delta = {}, g = {})...".format(twist["coeffs"], D["g"]))
                x_coords, pts = twist_chabauty(curve, twist_index=i, g_index=j)
                D["x-coords"] = x_coords
                twist["pts"] = pts
                if len(pts) < twist["found_pts"]:
                    raise PointCountError("Twist point count less than known point count! (delta = {})".format(twist["coeffs"]))
                twist["verified"] = True
                logging.info("Chabauty completed successfully (delta = {}, g = {}).".format(twist["coeffs"], D["g"]))
                t = record_data(curve, OUTPUT_FILE, t)

        # Extract x-coordinates for points on the curve
        if not curve["verified"]:
            coords, verified = x_coords_of_twist_pts(curve)
            curve["x-coords"] = coords
            curve["verified"] = verified
            t = record_data(curve, OUTPUT_FILE, t)

        if curve["pts"] is None and curve["verified"]:
            pts = []
            for P in curve["x-coords"]:
                pts += lift_to_hyperelliptic_curve(curve["coeffs"], P[0], P[1])
            curve["pts"] = pts
            curve["count"] = len(pts)
            curve["exception"] = False
            if len(pts) != known_point_count(curve, bound=CURVE_SEARCH_BOUND):
                raise PointCountError("Computed point count on genus 2 curve does not match known point count!")
            logging.info("Rational points successfully computed.")
            t = record_data(curve, OUTPUT_FILE, t)

    if "class" in STAGES and args.unconditional and not all(c[2] for c in curve["class_group_fields"]):
        code = build_class_group_code(curve["class_group_fields"])
        with open(CLASS_GROUP_FILE, "w") as f:
            f.write(code)
        logging.info("To make class group computations unconditional, run the following Magma file: {}".format(CLASS_GROUP_FILE.split("/")[-1]))
except Obstruction:
    curve["exception"] = False
    logging.info("Obstruction found to provably computing rational points. Exiting.")
except:
    # If an uncaught exception happens at any point, record that it happened first
    curve["exception"] = True
    logging.exception("An exception occurred. Recording data and exiting.")
else:
    logging.info("Tasks complete. Recording data and exiting.")
finally:
    record_data(curve, OUTPUT_FILE, t, final=True)

