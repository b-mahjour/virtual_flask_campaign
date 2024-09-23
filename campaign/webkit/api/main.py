import flask
import webkit
from shared.hpc.commands import get_nodes, get_edges,get_network
from shared.hpc.reset_tables import connect_to_local_psql
import networkx as nx
from shared.util import get_state_img, find_convex_hull

@webkit.app.route("/api/get_network", methods=["GET"])
def serve_network():
    context = {}


    conn, cur = connect_to_local_psql()

    net = get_network(cur,2, "rxrange")
    nodes = get_nodes(cur,2, name="rxrange")
    edges = get_edges(cur,2, name="rxrange")
    net.load_data(nodes, edges)

    hits = []
    for k in nodes:
        if k.other_data["forward_predictor_augmented_transformer"]:
            hits.append(k)

    root = net.get_origin_node()
    print(root)
    links = []
    hits2 = []
    nx_network = nx.DiGraph()
    for k in hits:
        path = net.get_shortest_path(root, k.node_id)
        product = False
        for idx,jj in enumerate(path):
            sm = net.get_node(jj).unmapped_smiles
            if idx == len(path)-1:
                product = True
            if sm not in nx_network:
                img = get_state_img(sm)
                hull = find_convex_hull(img)
                node_data = {
                    "id": sm,
                    "hull": hull,
                    "img": img,
                    "sm": sm,
                    "x": 0,
                    "y": 0,
                    "product": product
                }

                nx_network.add_node(sm, **node_data)
                if idx > 0:
                    prev_sm = net.get_node(path[idx-1]).unmapped_smiles
                    nx_network.add_edge(prev_sm, sm)
                    links.append((prev_sm, sm))
                hits2.append(node_data)
    pos = nx.nx_agraph.graphviz_layout(nx_network, prog="dot")

    # max_x = max([k[0] for k in pos.values()]) * .9
    # max_y = max([k[1] for k in pos.values()]) * .9
    # for k in hits2:
    #     k["x"] = pos[k["sm"]][0]/max_x - .05
    #     k["y"] = pos[k["sm"]][1]/max_y - .05
    #     print(k["x"], k["y"])

    normalized_positions = normalize_positions(pos)

    for k in hits2:
        k['x'], k['y'] = normalized_positions[k['sm']]


    context["network"] = {
        "nodes": hits2,
        "links": links
    }
    return flask.jsonify(**context)



def normalize_positions(pos, padding=0.1):
    # Determine the bounds of the current positions
    min_x = min(k[0] for k in pos.values())
    max_x = max(k[0] for k in pos.values())
    min_y = min(k[1] for k in pos.values())
    max_y = max(k[1] for k in pos.values())

    # Calculate the range considering padding
    width = max_x - min_x
    height = max_y - min_y
    scale_x = (1 - padding) / width
    scale_y = (1 - padding) / height

    # Adjust positions with padding offset
    padded_min_x = min_x - width * (padding / 2)
    padded_min_y = min_y - height * (padding / 2)

    # Normalize and scale coordinates
    for k in pos:
        pos[k] = ((pos[k][0] - padded_min_x) * scale_x, 
                  (pos[k][1] - padded_min_y) * scale_y)

    return pos