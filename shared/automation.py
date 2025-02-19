from collections import deque
import numpy as np


def batch_liquid_handler(actions):

    tips = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}
    tip_content = {
        0: None,
        1: None,
        2: None,
        3: None,
        4: None,
        5: None,
        6: None,
        7: None,
    }
    current_tip = 0
    pending_dispenses = []
    pending_aspirates = []
    new_order = []
    qu = deque(actions)
    while qu:
        k = qu.popleft()
        if current_tip == 8:
            current_tip = 0
            total = 0
            t = pending_aspirates[0]["arg4"]
            new_tips = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}
            new_tips2 = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: []}
            # for i in pending_aspirates:
            # print(i["arg3"], i["arg4"])
            for i in pending_aspirates:
                # print(i["arg3"], i["arg4"])
                if i["arg4"] == t:
                    total = total + i["arg3"]
                    # print(total)
                else:
                    total = round(total, 2)
                    # print("a", total)
                    new_order.append(
                        {
                            "action": "aspirate",
                            "arg1": i["arg1"],
                            "arg2": i["arg2"],
                            "arg3": total,
                            "arg4": t,
                            "arg5": i["arg5"],
                        }
                    )
                    new_tips2[t].append(total)
                    new_tips[t] = total
                    t = i["arg4"]
                    total = i["arg3"]

            new_order.append(
                {
                    "action": "aspirate",
                    "arg1": i["arg1"],
                    "arg2": i["arg2"],
                    "arg3": round(total, 2),
                    "arg4": t,
                    "arg5": i["arg5"],
                }
            )

            t = 0
            cur = pending_dispenses[0]["arg5"]
            # print(new_tips2)
            for i in pending_dispenses:
                disp = round(i["arg3"], 2)
                if i["arg5"] != cur or np.abs(new_tips[t]) <= 0.1:
                    t = t + 1
                    cur = i["arg5"]
                new_tips[t] = new_tips[t] - disp
                i["arg3"] = disp
                i["arg4"] = t
                new_order.append(i)
                # print(disp, new_tips)
            # print(new_tips)
            new_order.append(
                {
                    "action": "wash",
                    "arg1": "",
                    "arg2": "",
                    "arg3": "",
                    "arg5": "",
                    "arg4": "",
                }
            )

            pending_aspirates = []
            pending_dispenses = []
            tips = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}
            tip_content = {
                0: None,
                1: None,
                2: None,
                3: None,
                4: None,
                5: None,
                6: None,
                7: None,
            }

        if k["action"] == "dispense":
            pending_dispenses.append(k)
            continue
        if k["action"] == "aspirate":
            if tip_content[current_tip] != None:
                if tip_content[current_tip] != k["arg5"]:
                    current_tip = current_tip + 1
                    qu.appendleft(k)
                    continue
                else:
                    new_vol = tips[current_tip] + k["arg3"]
                    if new_vol >= 400:
                        current_tip = current_tip + 1
                        qu.appendleft(k)
                        continue
                    else:
                        tips[current_tip] = new_vol
                        k["arg4"] = current_tip
                        pending_aspirates.append(k)
                        continue
            else:
                tip_content[current_tip] = k["arg5"]
                tips[current_tip] = k["arg3"]
                k["arg4"] = current_tip
                pending_aspirates.append(k)
                continue
    return new_order


def generate_worklist(solv, stock_info):
    actions = []
    for k in solv:
        actions.append(
            {
                "action": "aspirate",
                "arg1": "solventA",
                "arg2": "A1",
                "arg3": round(solv[k], 2),
                "arg5": "solventA_A1",
            }
        )
        actions.append(
            {
                "action": "dispense",
                "arg1": "rA",
                "arg2": k,
                "arg3": round(solv[k], 2),
                "arg5": "solventA_A1",
            }
        )
    for k in stock_info:
        if stock_info[k]["library"] == "mcr_redux25_s2":
            src = "iB"
        else:
            src = "iA"
        src_well = stock_info[k]["well"]
        if len(stock_info[k]["doses"]) == 0:
            continue
        for j in stock_info[k]["doses"]:
            # if j[0] == "D5":
            # print(k, j)
            actions.append(
                {
                    "action": "aspirate",
                    "arg1": src,
                    "arg2": src_well,
                    "arg3": round(j[1], 2),
                    "arg5": f"{src}_{src_well}",
                }
            )
            actions.append(
                {
                    "action": "dispense",
                    "arg1": "rA",
                    "arg2": j[0],
                    "arg3": round(j[1], 2),
                    "arg5": f"{src}_{src_well}",
                }
            )
    return actions


from shared.worklist_batcher import batch_tips


def generate_experiment(solv, stock_info):
    actions = generate_worklist(solv, stock_info)
    # for k in actions:
    # if k["arg2"] == "D5" and k["arg1"]=="rA":
    # print(k)
    batched_actions = batch_tips(actions)
    # batched_actions = batch_liquid_handler(actions)

    initial_seq = [
        {"action": "load_plate", "arg1": "iA", "arg2": "B1", "arg3": "mcr_redux25_s1"},
        {"action": "load_plate", "arg1": "iB", "arg2": "C1", "arg3": "mcr_redux25_s2"},
        {
            "action": "load_plate",
            "arg1": "rA",
            "arg2": "C3",
            "arg3": "mcr_redux25_1_24_scale_up_reaction_plate2",
        },
        {
            "action": "load_plate",
            "arg1": "solventA",
            "arg2": "C2",
            "arg3": "mcr_redux25_dmso_solvent",
        },
        *batched_actions,
    ]

    return initial_seq


def generate_worklist2(stock_info):
    actions = []
    src = "iA"
    for k in stock_info:
        src_well = stock_info[k]["well"]
        for j in stock_info[k]["doses"]:
            actions.append(
                {
                    "action": "aspirate",
                    "arg1": src,
                    "arg2": src_well,
                    "arg3": round(j[1], 2),
                    "arg5": f"{src}_{src_well}",
                }
            )
            actions.append(
                {
                    "action": "dispense",
                    "arg1": "rA",
                    "arg2": j[0],
                    "arg3": round(j[1], 2),
                    "arg5": f"{src}_{src_well}",
                }
            )
    return actions


def generate_experiment2(stock_info):
    actions = generate_worklist2(stock_info)
    batched_actions = batch_tips(actions)

    initial_seq = [
        {
            "action": "load_plate",
            "arg1": "iA",
            "arg2": "B1",
            "arg3": "scale_up_mcr_redux_inventory_dmso_nmp",
        },
        {
            "action": "load_plate",
            "arg1": "rA",
            "arg2": "C3",
            "arg3": "mcr_redux25_2_SU_rxn",
        },
        *batched_actions,
    ]

    return initial_seq


def generate_quench(solv, stock_info):
    # actions = generate_worklist(solv, stock_info)

    initial_seq = [
        {"action": "load_plate", "arg1": "iA", "arg2": "B1", "arg3": "mcr_redux25_s1"},
        {"action": "load_plate", "arg1": "iB", "arg2": "C1", "arg3": "mcr_redux25_s2"},
        {
            "action": "load_plate",
            "arg1": "rA",
            "arg2": "E1",
            "arg3": "mcr_redux25_96_reaction_plate",
        },
        {
            "action": "load_plate",
            "arg1": "solventA",
            "arg2": "C2",
            "arg3": "mcr_redux25_dmso_solvent",  # Water
        },
        {
            "action": "load_plate",
            "arg1": "solventB",
            "arg2": "C1",
            "arg3": "rxrange_mecn_solvent",  # MeCN
        },
    ]

    return initial_seq


import psycopg2
from psycopg2.extras import Json


def connect_to_rds(
    host="10.63.5.43",
    port=5432,
    dbname="postgres",
    user="postgres",
    password="pass",
):
    """Create a connection to the PostgreSQL database"""
    try:
        conn = psycopg2.connect(
            host=host, port=port, dbname=dbname, user=user, password=password
        )
        print("Connection established")
        return conn
    except Exception as e:
        print(f"Unable to connect to the database: {e}")
        print("No connection to RDS")
        return None


def submit_to_scheduler_db(campaign_name, experiment_name, userName, stage):
    conn = connect_to_rds()
    cur = conn.cursor()

    cur.execute(
        """SELECT * FROM campaign_experiments WHERE campaign_name = %s;""",
        (campaign_name,),
    )
    nid = cur.fetchall()
    been_done = False
    for jj in nid:
        if jj[3] == experiment_name:
            been_done = True
            break
    nid = len(nid) + 1
    if experiment_name == "" or been_done == False:
        if experiment_name == "":
            exp = f"{userName}-{campaign_name}-{nid}"
        else:
            exp = experiment_name
        cur.execute(
            "INSERT INTO campaign_experiments (campaign_name, experiment_name, variable_name, variable_value) VALUES (%s, %s, %s, %s);",
            (campaign_name, exp, "protocol", Json([stage])),
        )
        cur.execute(
            "INSERT INTO campaign_experiments (campaign_name, experiment_name, variable_name, variable_value) VALUES (%s, %s, %s, %s);",
            (campaign_name, exp, "output", Json(None)),
        )

        cur.execute(
            "INSERT INTO campaign_experiments (campaign_name, experiment_name, variable_name, variable_value) VALUES (%s, %s, %s, %s);",
            (campaign_name, exp, "state", Json("initialized")),
        )
    else:
        exp = experiment_name
        cur.execute(
            "UPDATE campaign_experiments SET variable_value = %s WHERE campaign_name = %s AND experiment_name = %s AND variable_name = %s;",
            (Json([stage]), campaign_name, exp, "protocol"),
        )
        cur.execute(
            "UPDATE campaign_experiments SET variable_value = %s WHERE campaign_name = %s AND experiment_name = %s AND variable_name = %s;",
            (Json("initialized"), campaign_name, exp, "state"),
        )

    conn.commit()
    cur.close()
    conn.close()
    return exp


import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


def draw_hypergraph_sub_no_filter(
    network, node_size=0.3, fig_size=(2, 2), name=None, colorm="Blues", shift=0
):
    hypx = network.nx
    sub = list(network.nodes.keys())

    pos = nx.nx_agraph.graphviz_layout(hypx, prog="dot")
    propagation_values = [network.nodes[node].propagations for node in sub]
    norm = mcolors.Normalize(
        vmin=min(propagation_values) - 1, vmax=max(propagation_values) + shift
    )
    cmap = plt.cm._colormaps.get_cmap(colorm)
    node_colors = [cmap(norm(value)) for value in propagation_values]

    node_colors_final = []
    s = []
    # idxes = []
    for i, k in enumerate(sub):
        # if k in highlight:
        #     node_colors_final.append("gold")
        #     s.append(node_size * 2)
        #     idxes.append(i)
        #     continue
        node_colors_final.append(node_colors[i])
        s.append(node_size)

    fig, ax = plt.subplots(dpi=300)
    fig.set_size_inches(fig_size[0], fig_size[1])

    nodes_x, nodes_y = zip(*[pos[n] for n in hypx.nodes()])
    nodes_y = [-y for y in nodes_y]
    ax.scatter(
        nodes_x,
        nodes_y,
        s=s,
        color=node_colors_final,
        alpha=0.6,
        zorder=2,
    )

    # for i in idxes:
    #     ax.scatter(
    #         nodes_x[i],
    #         nodes_y[i],
    #         s=s[i],
    #         color="gold",
    #         alpha=1,
    #         zorder=4,
    #     )

    for edge in hypx.edges():
        start_x, start_y = pos[edge[0]]
        start_y = -start_y
        end_x, end_y = pos[edge[1]]
        end_y = -end_y

        edge_vec = np.array([end_x - start_x, end_y - start_y])
        edge_length = np.linalg.norm(edge_vec)
        edge_direction = edge_vec / edge_length

        start_x, start_y = (
            start_x + edge_direction[0] * node_size,
            start_y + edge_direction[1] * node_size,
        )
        end_x, end_y = (
            end_x - edge_direction[0] * node_size,
            end_y - edge_direction[1] * node_size,
        )

        # if edge[0] in highlight and edge[1] in highlight:
        #     ax.plot(
        #         [start_x, end_x],
        #         [start_y, end_y],
        #         color="gold",
        #         alpha=1,
        #         zorder=3,
        #         linewidth=0.5,
        #     )
        # else:

        ax.plot(
            [start_x, end_x],
            [start_y, end_y],
            color="black",
            alpha=0.1,
            zorder=1,
            linewidth=0.3,
        )

    plt.axis("off")

    if name == None:
        plt.show()
    else:
        plt.savefig(
            f"{name}.png",
            dpi=900,
            transparent=True,
            bbox_inches="tight",
            pad_inches=0.01,
        )
