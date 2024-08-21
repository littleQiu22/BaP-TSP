from gurobipy import Model, GRB, quicksum, Var, Constr
import networkx as nx
from typing import Dict, List, Tuple, Set
import time
import math
import os

######### functions to read tsp instance from file


def deg_to_radians(x):
    deg = int(x)
    mm = 100 * (x - deg)
    return (deg + mm / 60) / 180 * math.pi


def get_graph(data_path):
    N = None
    weight_graph = nx.Graph()
    with open(data_path, "r") as f:
        r = 0
        c = 0
        earth_radius = 6378.388
        line = f.readline()
        N = None
        weight_type = None
        weight_format = None
        while line != "":
            if line.startswith("DIMENSION"):
                tokens = line.split()
                N = int(tokens[1])
            elif line.startswith("EDGE_WEIGHT_TYPE"):
                tokens = line.split()
                weight_type = tokens[1]
            elif line.startswith("EDGE_WEIGHT_FORMAT"):
                tokens = line.split()
                weight_format = tokens[1]
            elif line.startswith("EOF"):
                break
            elif line.startswith("EDGE_WEIGHT_SECTION"):
                if weight_format == "UPPER_ROW":
                    c = 1
                    while True:
                        finished = False
                        line = f.readline()
                        weight_strs = line.split()
                        weight_ints = [int(weight_str) for weight_str in weight_strs]
                        for weight in weight_ints:
                            if c >= N:
                                r += 1
                                c = r + 1
                            weight_graph.add_edge(r, c, weight=weight)
                            c += 1
                            if c == N and r == N - 2:
                                finished = True
                                break
                        if finished:
                            break
                if weight_format == "LOWER_DIAG_ROW":
                    while True:
                        finished = False
                        line = f.readline()
                        weight_strs = line.split()
                        weight_ints = [int(weight_str) for weight_str in weight_strs]
                        for weight in weight_ints:
                            if c == N - 1 and r == N - 1:
                                finished = True
                                break
                            if c > r:
                                r += 1
                                c = 0
                            if c < r:
                                weight_graph.add_edge(r, c, weight=weight)
                            c += 1
                        if finished:
                            break
                elif weight_format == "FULL_MATRIX":
                    while True:
                        finished = False
                        line = f.readline()
                        weight_strs = line.split()
                        weight_ints = [int(weight_str) for weight_str in weight_strs]
                        for weight in weight_ints:
                            if c == N - 1 and r == N - 1:
                                finished = True
                                break
                            if c >= N:
                                r += 1
                                c = 0
                            if c < r:
                                weight_graph.add_edge(r, c, weight=weight)
                            c += 1
                        if finished:
                            break
            elif line.startswith("NODE_COORD_SECTION"):
                positions = []
                for _ in range(N):
                    line = f.readline()
                    tokens = line.split()
                    positions.append((float(tokens[1]), float(tokens[2])))
                for i in range(N):
                    for j in range(i):
                        weight = None
                        if weight_type == "ATT":
                            dx = abs(positions[i][0] - positions[j][0])
                            dy = abs(positions[i][1] - positions[j][1])
                            weight = math.ceil(math.sqrt((dx * dx + dy * dy) / 10))
                        elif weight_type == "GEO":
                            lat_i = deg_to_radians(positions[i][0])
                            long_i = deg_to_radians(positions[i][1])
                            lat_j = deg_to_radians(positions[j][0])
                            long_j = deg_to_radians(positions[j][1])
                            q1 = math.cos(long_i - long_j)
                            q2 = math.cos(lat_i - lat_j)
                            q3 = math.cos(lat_i + lat_j)
                            weight = math.floor(
                                earth_radius
                                * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3))
                                + 1.0
                            )
                        elif weight_type == "EUC_2D":
                            dx = abs(positions[i][0] - positions[j][0])
                            dy = abs(positions[i][1] - positions[j][1])
                            weight = round(math.sqrt(dx * dx + dy * dy))
                        weight_graph.add_edge(i, j, weight=weight)
            line = f.readline()
    if N % 2 != 0:
        new_loc_idx = N
        overlap_loc_idx = 0
        for i in range(N):
            if i == overlap_loc_idx:
                weight = 0
            else:
                weight = weight_graph[i][overlap_loc_idx]["weight"]
            weight_graph.add_edge(i, new_loc_idx, weight=weight)
        N += 1
    return N, weight_graph


class Partition:
    """
    Record partition of two subset of vertexs which violate DFJ subtour constraints
    """

    NEXT_PARTITION_ID = 1

    def __init__(self, sub_vertexs1, sub_vertexs2, connect_edges) -> None:
        self.sub_vertexs1 = sub_vertexs1
        self.sub_vertexs2 = sub_vertexs2
        self.connect_edges: Set[Tuple] = connect_edges
        self.partiton_id = Partition.NEXT_PARTITION_ID
        Partition.NEXT_PARTITION_ID += 1

    def __repr__(self) -> str:
        return f"SubVertexs1: {self.sub_vertexs1}; SubVertexs2: {self.sub_vertexs2}; Connected edges: {self.connect_edges}"


class Matching:
    """
    Solution of TSP with even vertexs consists of two disjoint perfect matching
    """

    NEXT_MATCHING_ID = 1

    def __init__(self, sum_weight, edge_set: Set[Tuple]) -> None:
        self.sum_weight = sum_weight
        self.edge_set: Set[Tuple] = edge_set
        self.matching_id = self.NEXT_MATCHING_ID
        Matching.NEXT_MATCHING_ID += 1

    def contains(self, edge):
        return (edge[0], edge[1]) in self.edge_set or (
            edge[1],
            edge[0],
        ) in self.edge_set

    def __repr__(self) -> str:
        return f"Sum Weight: {self.sum_weight}; Edges: {self.edge_set}"


class ColumnGeneration:
    def __init__(
        self,
        N,
        weight_graph: nx.Graph,
    ) -> None:
        self.N = N
        self.weight_graph: nx.Graph = weight_graph  # weight between vertexs
        self.flow_graph: nx.Graph = (
            self.weight_graph.copy()
        )  # to check whether DFJ constraints are violated
        self.pricing_graph: nx.Graph = (
            self.weight_graph.copy()
        )  # to solve pricing problem
        self.iter = 0
        self.rmp_obj = 1e20  # current objective value of restricted master problem
        self.mp_lb = -1e20  # lower bound of master problem
        self.eps = 1e-5
        self.gap_tol = 1e-5
        self.log_template = "{iter:<10d}{rmp_obj:<15.6e}{mp_lb:<15.6e}{min_rc:<15.6e}{gap:<10.2f}{time:<10.2f}"
        self.rmp = Model()
        self.matching_var_map = (
            {}
        )  # key: Matching object --> value: Gurobi continuous variable
        self.edge_constr_map = (
            {}
        )  # key: edge like (u, v) tuple format --> value: EdgeFlow constr
        self.subtour_constr_map = (
            {}
        )  # key: Partition object --> value: DFJ subtour constr
        self.rmp.setParam(GRB.Param.LogToConsole, 0)
        # add initial variables to rmp
        init_machings = self.get_init_matchings(self.N, self.weight_graph)
        for matching in init_machings:
            self.matching_var_map[matching] = self.rmp.addVar(
                lb=0,
                ub=1,
                obj=matching.sum_weight,
                vtype=GRB.CONTINUOUS,
                name=f"Matching-{matching.matching_id}",
            )
        for edge in self.weight_graph.edges():
            constr = self.rmp.addLConstr(
                quicksum(
                    [
                        var if matching.contains(edge) else 0
                        for matching, var in self.matching_var_map.items()
                    ]
                )
                <= 1,
                name=f"EdgeFlow-{edge}",
            )
            self.edge_constr_map[edge] = constr
        self.two_deg_constr = self.rmp.addLConstr(
            quicksum(self.matching_var_map.values()) >= 2, name="TwoDegree"
        )

    def cal_gap(self, obj, bound):
        return 100 * abs(obj - bound) / max((abs(obj), abs(bound), self.eps))

    def reset_graph(self, graph: nx.Graph, attr="weight", attr_val=0):
        """
        Reset `graph` attribute `attr` to value `attr_val`
        """
        for _u, _v, data in graph.edges(data=True):
            data[attr] = attr_val

    def copy_weight(self, dest_graph: nx.Graph, src_graph: nx.Graph, opposite=False):
        """
        Copy `src_graph`'s (optional `opposite`) weight to `dest_graph`.
        """
        sign = -1 if opposite else 1
        for u, v, weight in src_graph.edges(data="weight"):
            dest_graph[u][v]["weight"] = sign * weight

    def get_init_matchings(self, N, graph: nx.Graph):
        seq = [i for i in range(N)] + [0]
        edge_set1 = set()
        sum_weight_1 = 0
        edge_set2 = set()
        sum_weight_2 = 0
        for idx in range(0, N, 2):
            edge_set1.add((seq[idx], seq[idx + 1]))
            sum_weight_1 += graph[seq[idx]][seq[idx + 1]]["weight"]
        for idx in range(1, N, 2):
            edge_set2.add((seq[idx], seq[idx + 1]))
            sum_weight_2 += graph[seq[idx]][seq[idx + 1]]["weight"]
        return [
            Matching(sum_weight_1, edge_set1),
            Matching(sum_weight_2, edge_set2),
        ]

    def solve(self):
        st = time.time()
        print(
            "{iter:<10}{rmp_obj:<15}{mp_lb:<15}{min_rc:<15}{gap:<10}{time:<10}".format(
                iter="Iter.",
                rmp_obj="RMP Obj.",
                mp_lb="MP LB",
                min_rc="Min RC",
                gap="Gap(%)",
                time="Time(s)",
            )
        )
        while True:
            self.iter += 1
            has_new_column = False
            self.rmp.optimize()
            if self.rmp.Status != GRB.Status.OPTIMAL:
                print("Error: RMP is not solved to optimal")
                return
            self.rmp_obj = self.rmp.ObjVal
            gap = self.cal_gap(self.rmp.ObjVal, self.mp_lb)
            reduced_cost = 0
            if gap > self.gap_tol:
                # since networkx only supports max weight perfect matching, we should taking the opposite of the weights
                self.copy_weight(self.pricing_graph, self.weight_graph, opposite=True)
                two_deg_pi = self.two_deg_constr.Pi
                for (u, v), edge_constr in self.edge_constr_map.items():
                    edge_constr_pi = edge_constr.Pi
                    self.pricing_graph[u][v]["weight"] += edge_constr_pi
                for partition, subtour_constr in self.subtour_constr_map.items():
                    subtour_constr_pi = subtour_constr.Pi
                    for u, v in partition.connect_edges:
                        if self.pricing_graph.has_edge(u, v):
                            self.pricing_graph[u][v]["weight"] += subtour_constr_pi
                max_weight_matching = nx.max_weight_matching(
                    self.pricing_graph, maxcardinality=True
                )
                total_pricing_weight = 0
                total_weight = 0
                for u, v in max_weight_matching:
                    total_pricing_weight += self.pricing_graph[u][v]["weight"]
                    total_weight += self.weight_graph[u][v]["weight"]
                total_pricing_weight *= -1
                reduced_cost = total_pricing_weight - two_deg_pi
                cur_mp_lb = self.rmp.ObjVal + reduced_cost * 2
                if cur_mp_lb > self.mp_lb:
                    self.mp_lb = cur_mp_lb
                new_machting = None
                if reduced_cost < -self.eps:
                    new_machting = Matching(total_weight, max_weight_matching)
                if new_machting is not None:
                    has_new_column = True
                    x = self.rmp.addVar(
                        lb=0,
                        ub=1,
                        obj=new_machting.sum_weight,
                        vtype=GRB.CONTINUOUS,
                        name=f"Matching-{new_machting.matching_id}",
                    )
                    self.matching_var_map[new_machting] = x
                    self.rmp.chgCoeff(self.two_deg_constr, x, 1)
                    for edge, edge_constr in self.edge_constr_map.items():
                        if new_machting.contains(edge):
                            self.rmp.chgCoeff(edge_constr, x, 1)
                    for partition, subtour_constr in self.subtour_constr_map.items():
                        connect_edges = partition.connect_edges
                        edge_num_in_new_matching = 0
                        for edge in connect_edges:
                            if new_machting.contains(edge):
                                edge_num_in_new_matching += 1
                        if edge_num_in_new_matching > 0:
                            self.rmp.chgCoeff(
                                subtour_constr, x, edge_num_in_new_matching
                            )
            print(
                self.log_template.format(
                    iter=self.iter,
                    rmp_obj=self.rmp.ObjVal,
                    mp_lb=self.mp_lb,
                    min_rc=reduced_cost,
                    gap=gap,
                    time=time.time() - st,
                )
            )
            if has_new_column:
                continue
            self.reset_graph(self.flow_graph)
            for matching, var in self.matching_var_map.items():
                val = var.X
                if val > self.eps:
                    for u, v in matching.edge_set:
                        self.flow_graph[u][v]["weight"] += val
            cut_val, partition = nx.stoer_wagner(self.flow_graph)
            if cut_val < 2 - self.eps:
                print(
                    "Find vertex subset that violates DFJ constraint",
                )
                connect_edges = set(
                    [
                        (vertex1, vertex2)
                        for vertex1 in partition[0]
                        for vertex2 in partition[1]
                        if self.flow_graph.has_edge(vertex1, vertex2)
                    ]
                )
                p = Partition(partition[0], partition[1], connect_edges)
                subtour_constr = self.rmp.addLConstr(
                    quicksum(
                        [
                            var if matching.contains(edge) else 0
                            for edge in connect_edges
                            for matching, var in self.matching_var_map.items()
                        ]
                    )
                    >= 2,
                    name=f"Subtour-{p.partiton_id}",
                )
                self.subtour_constr_map[p] = subtour_constr
                continue
            print(f"Solving finished. Solution type is {self.get_sol_type()}")
            break

    def get_sol_type(self):
        if self.rmp.SolCount <= 0:
            return None
        self.reset_graph(self.flow_graph, "weight", 0)
        for matching, var in self.matching_var_map.items():
            val = var.X
            if val > self.eps:
                for u, v in matching.edge_set:
                    self.flow_graph[u][v]["weight"] += val
        sol_type = "Integer"
        for u, v, data in self.flow_graph.edges(data=True):
            weight = data["weight"]
            if weight < self.eps or weight > 1 - self.eps:
                continue
            else:
                sol_type = "Fractional"
                break
        return sol_type


if __name__ == "__main__":
    _m = Model()  # just show Gurobi license in advance
    src_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(os.path.dirname(src_dir), "data")
    data_path = os.path.join(data_dir, "gr24.tsp")
    # data_path = os.path.join(data_dir, "fri26.tsp")
    # data_path = os.path.join(data_dir, "burma14.tsp")
    # data_path = os.path.join(data_dir, "ulysses16.tsp")
    # data_path = os.path.join(data_dir, "ulysses22.tsp")
    # data_path = os.path.join(data_dir, "gr17.tsp")
    # data_path = os.path.join(data_dir, "gr21.tsp")
    # data_path = os.path.join(data_dir, "bays29.tsp")
    # data_path = os.path.join(data_dir, "bayg29.tsp")
    # data_path = os.path.join(data_dir, "hk48.tsp")
    # data_path = os.path.join(data_dir, "att48.tsp")
    # data_path = os.path.join(data_dir, "berlin52.tsp")
    # data_path = os.path.join(data_dir, "swiss42.tsp")
    # data_path = os.path.join(data_dir, "dantzig42.tsp")
    # data_path = os.path.join(data_dir, "gr48.tsp")
    N, graph = get_graph(data_path)
    cg = ColumnGeneration(N, graph)
    cg.solve()
