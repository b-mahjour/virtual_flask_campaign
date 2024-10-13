import flask
import webkit
from shared.hpc.commands import get_nodes, get_edges, get_network, connect_to_rds
from shared.hpc.reset_tables import connect_to_local_psql
import networkx as nx
from shared.util import get_state_img, find_convex_hull
from psycopg2.extras import execute_batch
import boto3
import botocore
from botocore.exceptions import ClientError
import json
from rdkit import Chem
from rdkit.Chem import Descriptors

# from shared.aimnet2sph import AIMNet2Calculator
# from shared.energy_util import calc2
# import torch

# def load_calculator(model_file):
#     """
#     Load a pysisyphus calculator from a model file

#     Returns:
#         calc: pysisyphus.Calculator
#     """
#     model = torch.jit.load(model_file)
#     calc = AIMNet2Calculator(model)
#     return calc


# models_dir = webkit.app.config["MODELS_FOLDER"]

# model_file = f"{models_dir}/aimnet2_wb97m-d3_0.jpt"
# calculator = load_calculator(model_file)


def link_hits_from_experimental_to_vf():

    s3_client = boto3.client(
        "s3",
        aws_access_key_id=webkit.app.config["AWS_ACCESS_KEY_ID"],
        aws_secret_access_key=webkit.app.config["AWS_SECRET_ACCESS_KEY"],
    )

    bucket = "coley-hte-data"

    response = s3_client.list_objects_v2(
        Bucket=bucket, Prefix="bmahjour/rxrange/", Delimiter="/"
    )
    folders = [prefix.get("Prefix") for prefix in response.get("CommonPrefixes", [])]

    skippers = [
        "bmahjour/rxrange/init_random_96/",
        "bmahjour/rxrange/initial_24_test/",
        "bmahjour/rxrange/initial_24_test_tecan/",
    ]

    mass_to_wells = {}

    for f in folders:
        if f in skippers:
            continue
        object_key = f + "metadata.json"
        try:

            response = s3_client.get_object(Bucket=bucket, Key=object_key)

            metadata = json.loads(response["Body"].read().decode("utf-8"))
            res = metadata["reactionMetaData"]["results"]
            for k in res:
                for kk in res[k]:
                    if kk["description"] == "mz":
                        if kk["value"] not in mass_to_wells:
                            mass_to_wells[kk["value"]] = []
                        mass_to_wells[kk["value"]].append((f, k))
        except botocore.exceptions.ClientError as e:
            print(e)
            return
        except Exception as e:
            print(e)
            return
    return mass_to_wells


from collections import defaultdict


def assign_xy(nodes):
    nodes_by_propagations = defaultdict(list)
    for node in nodes:
        nodes_by_propagations[node["propagations"]].append(node)
    poss = {}
    for y_value, nodes_in_group in nodes_by_propagations.items():
        num_nodes = len(nodes_in_group)
        if num_nodes == 1:
            x_values = [0.5]
        else:
            x_values = [i / (num_nodes - 1) for i in range(num_nodes)]
        for node, x_value in zip(nodes_in_group, x_values):
            poss[node["sm"]] = (x_value, y_value)

    return poss


def create_condensed_network_on_db_if_not_yet_made(reset=False):

    conn = connect_to_rds(host="18.222.178.6", user="rxrange", password="pass")
    cur = conn.cursor()

    if reset:
        cur.execute("DROP TABLE IF EXISTS rxrange_condensed_nodes CASCADE;")
        cur.execute("DROP TABLE IF EXISTS rxrange_condensed_edges CASCADE;")
        conn.commit()
    cur.execute(
        f"""
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.tables
            WHERE table_name = 'rxrange_condensed_nodes'
        );
        """
    )
    if cur.fetchone()[0]:
        return

    net = get_network(cur, 1, "rxrange")
    nodes = get_nodes(cur, 1, name="rxrange")
    edges = get_edges(cur, 1, name="rxrange")
    net.load_data(nodes, edges)

    hit_nodes = []
    hit_edges = []
    root = net.get_origin_node()
    nx_network = nx.DiGraph()
    for k in nodes:
        if k.other_data["structural_failure"] == False:
            path = net.get_shortest_path(root, k.node_id)
            for idx, jj in enumerate(path):
                sm = net.get_node(jj).unmapped_smiles
                product = False
                if sm not in nx_network:
                    if idx == len(path) - 1:
                        product = True
                    pro = 0
                    if idx > 0:
                        edge = net.get_edge(path[idx - 1], jj)
                        hit_edges.append(
                            (
                                net.get_node(path[idx - 1]).unmapped_smiles,
                                sm,
                                edge["edge_id"],
                            )
                        )
                        pro = edge["other_data"]["propagations"]

                    if k.other_data["target_molecule"] is not None:
                        target_sm = k.other_data["target_molecule"]
                    else:
                        target_sm = None
                    img = get_state_img(sm)
                    hull = find_convex_hull(img)
                    node_data = {
                        "hull": hull,
                        "img": img,
                        "sm": sm,
                        "target_sm": target_sm,
                        "x": 0,
                        "y": 0,
                        "product": product,
                        "propagations": pro,
                        "parent_node_id": jj,
                    }

                    nx_network.add_node(sm, **node_data)
                    hit_nodes.append(node_data)

    print(len(hit_nodes), len(hit_edges))
    for kk in hit_edges:
        prev_sm, sm, idx = kk
        nx_network.add_edge(prev_sm, sm, idx=idx)

    pos = assign_xy(hit_nodes)

    normalized_positions = normalize_positions(pos, padding=0.025)

    for k in hit_nodes:
        k["x"], k["y"] = normalized_positions[k["sm"]]

    table_name_nodes = "rxrange_condensed_nodes"
    table_name_edges = "rxrange_condensed_edges"
    create_table_command = f"""
        CREATE TABLE IF NOT EXISTS {table_name_nodes} (
            id SERIAL PRIMARY KEY,
            hull INT[],
            img TEXT,
            sm TEXT,
            target_sm TEXT,
            x FLOAT,
            y FLOAT,
            product BOOLEAN,
            propagations INT,
            parent_node_id INT
        );    
        """
    cur.execute(create_table_command)
    create_table_command = f"""
    CREATE TABLE IF NOT EXISTS {table_name_edges} (
        edge_id SERIAL PRIMARY KEY,
        source_node_smiles TEXT,
        destination_node_smiles TEXT,
        parent_link_id INT
    );
    """
    cur.execute(create_table_command)

    data_to_insert = []
    for node in hit_nodes:
        data_to_insert.append(
            (
                node["hull"],
                node["img"],
                node["sm"],
                node["target_sm"],
                node["x"],
                node["y"],
                node["product"],
                node["propagations"],
                node["parent_node_id"],
            )
        )
    execute_batch(
        cur,
        f"INSERT INTO rxrange_condensed_nodes (hull, img, sm, target_sm, x, y, product, propagations, parent_node_id) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)",
        data_to_insert,
    )

    l2 = []
    for link in hit_edges:
        l2.append((link[0], link[1], link[2]))
    execute_batch(
        cur,
        f"INSERT INTO rxrange_condensed_edges (source_node_smiles, destination_node_smiles, parent_link_id) VALUES (%s, %s, %s)",
        l2,
    )

    conn.commit()


class CondensedNode:
    def __init__(self, id, hull, img, sm, target_sm, x, y, product, propagations, parent_node_id):
        self.id = id
        self.hull = hull
        self.img = img
        self.sm = sm
        self.target_sm = target_sm
        self.x = x
        self.y = y
        self.product = product
        self.propagations = propagations
        self.parent_node_id = parent_node_id

    def to_dict(self):
        return {
            "id": self.id,
            "hull": self.hull,
            "img": self.img,
            "sm": self.sm,
            "target_sm": self.target_sm,
            "x": self.x,
            "y": self.y,
            "product": self.product,
            "propagations": self.propagations,
            "parent_node_id": self.parent_node_id,
        }


class CondensedEdge:
    def __init__(self, id, source_node_smiles, destination_node_smiles, parent_link_id):
        self.id = id
        self.source_node_smiles = source_node_smiles
        self.destination_node_smiles = destination_node_smiles
        self.parent_link_id = parent_link_id

    def to_dict(self):
        return {
            "source_node_smiles": self.source_node_smiles,
            "destination_node_smiles": self.destination_node_smiles,
            "parent_link_id": self.parent_link_id,
        }


def remove_atom_mapping(sm_in):
    mol = Chem.MolFromSmiles(sm_in)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)

@webkit.app.route("/api/get_network", methods=["GET"])
def serve_network():
    context = {}

    conn = connect_to_rds(host="18.222.178.6", user="rxrange", password="pass")

    cur = conn.cursor()

    create_condensed_network_on_db_if_not_yet_made()

    mass_to_wells = link_hits_from_experimental_to_vf()
    # for k in mass_to_wells:
        # print(k, mass_to_wells[k])

    # print(1)
    table_name_nodes = "rxrange_condensed_nodes"
    table_name_edges = "rxrange_condensed_edges"

    cur.execute(f"SELECT * FROM {table_name_nodes}")
    nodes = cur.fetchall()
    hits2 = [CondensedNode(*node).to_dict() for node in nodes]

    sm_to_exact_mass = {}
    masses_found = []
    for h in hits2:
        h["reaction_hits"] = []
        sms = h["sm"].split(".")
        for sm in sms:
            mol = Chem.MolFromSmiles(sm)
            sm_to_exact_mass[sm] = Descriptors.ExactMolWt(mol)
            if sm_to_exact_mass[sm] < 30:
                continue
            if sm == "NCCc1ccc(O)cc1" or sm == "CC(C)(C)OC(=O)N[C@@H](CCC(=O)O)C(=O)O":
                continue
            # print(sm_to_exact_mass[sm], sm)
            target = False
            if sm == remove_atom_mapping(h["target_sm"]):
                target = True
            if "target_found" not in h:
                h["target_found"] = False
            for mass in mass_to_wells:
                if abs(float(mass) - (sm_to_exact_mass[sm]+1)) < 0.1:
                    wells = [k[1] for k in mass_to_wells[mass]]
                    plates = [k[0][0:-1].split("/")[-1] for k in mass_to_wells[mass]]
                    well_plate = [f"{plates[i]}-{wells[i]}" for i in range(len(wells))]
                    h["reaction_hits"].append(
                        {
                            "compound": sm,
                            "exact_mass": sm_to_exact_mass[sm],
                            "mass_found": mass,
                            "wells": well_plate,
                            "adduct": "M+H",
                        }
                    )
                    if mass not in masses_found:
                        masses_found.append(mass)
                    if target:
                        # print(remove_atom_mapping(h["target_sm"]), sm)
                        h["target_found"] = True
                    # print(well_plate)
                if abs(float(mass) - (sm_to_exact_mass[sm] + 23)) < 0.1:
                    wells = [k[1] for k in mass_to_wells[mass]]
                    plates = [k[0][0:-1].split("/")[-1] for k in mass_to_wells[mass]]
                    well_plate = [f"{plates[i]}-{wells[i]}" for i in range(len(wells))]
                    h["reaction_hits"].append(
                        {
                            "compound": sm,
                            "exact_mass": sm_to_exact_mass[sm],
                            "mass_found": mass,
                            "wells": well_plate,
                            "adduct": "M+Na",
                        }
                    )
                    if mass not in masses_found:
                        masses_found.append(mass)
                    if target:
                        h["target_found"] = True
                    # print(well_plate)

    for mass in mass_to_wells:
        if mass not in masses_found:
            if mass > 200:
                print(mass, mass_to_wells[mass])

    cur.execute(f"SELECT * FROM {table_name_edges}")
    edges = cur.fetchall()
    # print(len(edges))
    seen_links = []
    links = []
    link_ids = []
    for edge in edges:
        # print(edge)
        if (edge[1], edge[2]) not in seen_links:
            link_ids.append(str(edge[3]))
            # print(data)
            seen_links.append((edge[1], edge[2]))
            links.append(CondensedEdge(*edge).to_dict())
            # links[-1]["other_data"] = data[3]

    placeholders = ", ".join(["%s"] * len(link_ids))
    query = f"SELECT * FROM rxrange_edges WHERE edge_id IN ({placeholders})"
    cur.execute(query, link_ids)
    out = cur.fetchall()
    for i, edge in enumerate(edges):
        links[i]["other_data"] = out[i][3]

    # print(3)

    context["network"] = {"nodes": hits2, "links": links}
    return flask.jsonify(**context)

    # net = get_network(cur, 1, "rxrange")
    # print("debug1")
    # nodes = get_nodes(cur, 1, name="rxrange")
    # print("debug2")
    # edges = get_edges(cur, 1, name="rxrange")
    # print("debug3")
    # net.load_data(nodes, edges)
    # hits = []
    # neighbor_hits = []
    # paths = []
    # edges = []
    # all_nodes_collect_ids = []
    # root = net.get_origin_node()
    # for k in nodes:
    #     if k.other_data["structural_failure"] == False:
    #         hits.append(k)
    #         all_nodes_collect_ids.append(k.node_id)
    #         path = net.get_shortest_path(root, k.node_id)
    #         for i, p in enumerate(path):
    #             if i > 0:
    #                 prev_sm = net.get_node(path[i - 1]).unmapped_smiles
    #                 edges.append((prev_sm, net.get_node(p).unmapped_smiles))
    #         # if len(hits) > 5:
    #         #     break

    #             # neighbor = net.find_neighbors(p)
    #             # for kk in neighbor:
    #             #     if kk == k:
    #             #         continue
    #             #     if kk in path:
    #             #         continue
    #             #     edges.append(
    #             #         (
    #             #             net.get_node(kk).unmapped_smiles,
    #             #             net.get_node(p).unmapped_smiles,
    #             #         )
    #             #     )
    #             #     neighbor_hits.append(net.get_node(kk))
    #             #     all_nodes_collect_ids.append(kk)

    #         paths.append(path)

    # print(len(hits), len(paths))
    # hits2 = []
    # nx_network = nx.DiGraph()
    # seen_node = []
    # for i, k in enumerate(hits):
    #     product = False
    #     path = paths[i]
    #     for idx, jj in enumerate(path):
    #         sm = net.get_node(jj).unmapped_smiles
    #         if idx == len(path) - 1:
    #             product = True
    #         if sm not in nx_network:

    #             neighbor_edge_ids = net.nx_graph.in_edges(jj)
    #             props = []
    #             for id in neighbor_edge_ids:
    #                 ed_data = net.get_edge(id[0], id[1])["other_data"]
    #                 props.append(ed_data["propagations"])

    #             pro = min(props)
    #             img = get_state_img(sm)
    #             hull = find_convex_hull(img)
    #             node_data = {
    #                 "id": sm,
    #                 "hull": hull,
    #                 "img": img,
    #                 "sm": sm,
    #                 "x": 0,
    #                 "y": 0,
    #                 "product": product,
    #                 "propagations": pro,
    #             }

    #             nx_network.add_node(sm, **node_data)

    # links = []
    # for kk in edges:
    #     prev_sm, sm = kk
    #     nx_network.add_edge(prev_sm, sm)
    #     links.append((prev_sm, sm))

    # for hit in hits2:
    #     if hit["sm"] == net.get_node(root).unmapped_smiles:
    #         hit["propagations"] = 0
    # print(len(hits2))
    # pos = assign_xy(hits2)

    # normalized_positions = normalize_positions(pos)

    # for k in hits2:
    #     k["x"], k["y"] = normalized_positions[k["sm"]]


def normalize_positions(pos, padding=0.1):
    # Determine the bounds of the current positions
    min_x = min(k[0] for k in pos.values())
    max_x = max(k[0] for k in pos.values())
    min_y = min(k[1] for k in pos.values())
    max_y = max(k[1] for k in pos.values())

    # Calculate the range considering padding
    width = max_x - min_x
    height = max_y - min_y
    # print(width, height)
    scale_x = (1 - padding) / width
    scale_y = (1 - padding) / height

    # Adjust positions with padding offset
    padded_min_x = min_x - width * (padding / 2)
    padded_min_y = min_y - height * (padding / 2)

    # Normalize and scale coordinates
    for k in pos:
        pos[k] = (
            (pos[k][0] - padded_min_x) * scale_x,
            (pos[k][1] - padded_min_y) * scale_y,
        )

    return pos
