from shared.hpc.batch import SingleNetworkDistributedPipeline
from rdkit import RDLogger
import sys
from shared.hpc.reset_tables import cancel_all_transactions, reset_psql
RDLogger.DisableLog("rdApp.*")
from shared.hpc.commands import get_nodes, get_edges,get_network
import psycopg2

def connect_to_rds(
    host="18.222.178.6",
    port=6432,
    dbname="testdb",
    user="rxrange",
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
        return None


conn = connect_to_rds()
cur = conn.cursor()

net = get_network(cur,1, "rxrange")
nodes = get_nodes(cur,1, name="rxrange")
edges = get_edges(cur,1, name="rxrange")
net.load_data(nodes, edges)


origin = net.get_origin_node()
print(origin)
all_nodes_to_be_evald = []
for node in nodes:
    if node.other_data["structural_failure"] == False:
        path = net.get_all_shortest_paths(origin, node.node_id)
        valid_path = True
        for p in path:
            for pp in p:
                path_node = net.get_node(pp)
                if "single_point_energy" in path_node.other_data:
                    if path_node.other_data["single_point_energy"] == -999999:
                        valid_path = False
                        break
            if valid_path:
                for pp in p:
                    if pp not in all_nodes_to_be_evald:
                        all_nodes_to_be_evald.append(pp)
        node.other_data["path_from_origin"] = path
print(len(all_nodes_to_be_evald))


x = SingleNetworkDistributedPipeline("ten", all_nodes_to_be_evald[400:], 500)
import time 
time0 = time.time()
x.pipe("g16")

print(time.time()-time0, "seconds")
