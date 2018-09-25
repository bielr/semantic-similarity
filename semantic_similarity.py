from itertools import groupby
from math import inf, log
import csv
import networkx as nx

# import go_tools
# import TCSS.main
# from TCSS.ontology import GOGraph

# async def load_tcss_objs(
#         ontology,
#         evidence_codes = go_tools.get_curated_evidence_codes(),
#         ontology_file = go_tools.get_obo_path()):

#     objs = {}
#     ontology = ontology.split(",")
#     g = GOGraph()

#     g._obo_parser(ontology_file)

#     with stringdb_connect() as stringdb_conn:
#         with stringdb_conn.cursor() as cursor:
#             for term in g.go_annotations.keys():
#                 prots = go_tools.string_get_explicit_annotations(string_cursor, term, evidence_codes)
#                 g.go_annotations[term]['gene'] = set(prots)

#     run = {'C':g._cellular_component, 'P':g._biological_process, 'F':g._molecular_function}
#     ont = {'C':"Cellular Component", 'P':"Biological Process", 'F':"Molecular Function"}

#     for i in ontology:
#         ns, cutoff = i.split(":")
#         cutoff = float(cutoff)

#         objs[ns] = run[ns]()
#         objs[ns]._species()
#         objs[ns]._clustering(cutoff)

#     return objs


def init_ic(freqs_file):
    with open(freqs_file) as f:
        reader = csv.reader(f, delimiter='\t')
        freqs = {go: int(freq) for go, freq in reader if int(freq) > 0}

    total_anns = sum(freqs.values())
    freqs = {go: -log(freq/total_anns) for go, freq in freqs.items()}
    return lambda go: freqs[go] #.get(go, inf)

def get_root(rel_g, ic, term):
    desc = nx.descendants(rel_g, term).union({term})
    return min(desc, key=ic, default=None)

def get_mica(rel_g, ic, term1, term2):
    desc1 = nx.descendants(rel_g, term1).union({term1})
    desc2 = nx.descendants(rel_g, term2).union({term2})
    ca = desc1.intersection(desc2)

    return max(ca, key=ic, default=None)

def get_mil(rel_g, ic, term):
    anc = nx.descendants(rel_g, term).union({term})
    return max(anc, key=ic)

def ic_dist(ic, u, v):
    ic_u = ic(u)
    ic_v = ic(v)
    return abs(ic(v) - ic(u))

def get_mica_dissim(rel_g, ic, term1, term2):
    mica = get_mica(rel_g, ic, term1, term2)
    return ic_dist(ic, term1, mica) + ic_dist(ic, mica, term2)

def get_hrss_sim(rel_g, ic, term1, term2):
    mica = get_mica(rel_g, ic, term1, term2)
    root = get_root(rel_g, ic, mica)

    alpha = ic_dist(ic, root, mica)
    gamma = ic_dist(ic, term1, mica) + ic_dist(ic, mica, term2)
    beta = ( ic(get_mil(rel_g, ic, term1)) + ic(get_mil(rel_g, ic, term2)) ) / 2

    return 1/(1+gamma) * alpha/(alpha+beta)

def agg_bma_min(f, gos1, gos2):
    num = sum(min(f(go1, go2) for go1 in gos1) for go2 in gos2) + sum(min(f(go1, go2) for go2 in gos2) for go1 in gos1)
    denom = len(gos1) + len(gos2)
    return num/denom

def agg_bma_max(f, gos1, gos2):
    mat = [[f(go1, go2) for go2 in gos2] for go1 in gos1]

    num = sum(max(mat[i][j] for j in range(len(gos2))) for i in range(len(gos1))) \
        + sum(max(mat[i][j] for i in range(len(gos1))) for j in range(len(gos2)))
    # num = sum(max(f(go1, go2) for go1 in gos1) for go2 in gos2) + sum(max(f(go1, go2) for go2 in gos2) for go1 in gos1)
    denom = len(gos1) + len(gos2)
    return num/denom

def agg_min(f, gos1, gos2):
    return min(f(go1, go2) for go1 in gos1 for go2 in gos2)

def agg_max(f, gos1, gos2):
    return max(f(go1, go2) for go1 in gos1 for go2 in gos2)

def namespace_wise_comparisons(onto, agg_f, gos1, gos2):
    def namespace(go):
        return onto[go].other['namespace']

    grouped_gos1 = [(ns1, list(g_gos1)) for ns1, g_gos1 in groupby(sorted(gos1, key=namespace), key=namespace)]
    grouped_gos2 = [(ns2, list(g_gos2)) for ns2, g_gos2 in groupby(sorted(gos2, key=namespace), key=namespace)]

    for ns1, g_gos1 in grouped_gos1:
        for ns2, g_gos2 in grouped_gos2:
            if ns1 == ns2:
                #assert g_gos1 and g_gos2, str((grouped_gos1, grouped_gos2))
                yield agg_f(g_gos1, g_gos2)
