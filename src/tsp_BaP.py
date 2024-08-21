from gurobipy import Model, GRB, quicksum, Var, Constr
import networkx as nx
from typing import Dict, List, Tuple, Set
import time
import math
import os
import heapq

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

    def __init__(self, sum_weight, edge_set: Set[Tuple], is_artificial=False) -> None:
        self.sum_weight = sum_weight
        self.edge_set: Set[Tuple] = edge_set
        self.matching_id = self.NEXT_MATCHING_ID
        Matching.NEXT_MATCHING_ID += 1
        # artificial matching will not be deleted after any branching
        self.is_artificial = is_artificial

    def contains(self, edge):
        return (edge[0], edge[1]) in self.edge_set or (
            edge[1],
            edge[0],
        ) in self.edge_set

    def __repr__(self) -> str:
        return f"Sum Weight: {self.sum_weight}; Edges: {self.edge_set}"


class Node:
    NEXT_NODE_ID = 1

    SOLVE_OPTIMALITY = 0
    SOLVE_TIMELIMIT_ABORT = 1
    SOLVE_INFEASIBLE = 2
    SOLVE_NOT_START = 3
    SOLVE_INVALID = 4
    SOLVE_SOLVING = 5

    NODE_FRACTIONAL = 0
    NODE_INTEGER = 1
    NODE_PRUNED = 2
    NODE_INFEASIBLE = 3
    NODE_NOT_DEFINE = 4
    NODE_INVALID = 5

    SOLVE_STATUS_REPR_MAP = {
        SOLVE_OPTIMALITY: "Optimal",
        SOLVE_TIMELIMIT_ABORT: "TimeLimit",
        SOLVE_INFEASIBLE: "Infeasible",
        SOLVE_NOT_START: "NotStart",
        SOLVE_INVALID: "Invalid",
        SOLVE_SOLVING: "Solving",
    }

    NODE_TYPE_REPR_MAP = {
        NODE_FRACTIONAL: "Fractional",
        NODE_INTEGER: "Integer",
        NODE_PRUNED: "Pruned",
        NODE_INFEASIBLE: "Infeasible",
        NODE_NOT_DEFINE: "NotDefine",
        NODE_INVALID: "Invalid",
    }

    NODE_COMPARE_BOUND = 0  # bound first search
    NODE_COMPARE_DEPTH = 1  # integer first search

    class NodeBound:
        def __init__(self, lb, is_valid=True, is_feasible=True) -> None:
            self.lb = lb
            self.is_valid = is_valid
            self.is_feasible = is_feasible

        def __lt__(self, other: "Node.NodeBound"):
            if self.is_feasible and not other.is_feasible:
                return True
            if not self.is_feasible and other.is_feasible:
                return False
            return self.lb < other.lb

    def __init__(
        self,
        parent_node: "Node",
        N,
        weight_graph: nx.Graph,
        rmp: Model = None,
        matching_var_map={},
        two_deg_constr=None,
        edge_constr_map={},
        subtour_constr_map={},
        node_compare=0,
    ) -> None:
        self.node_id = Node.NEXT_NODE_ID
        Node.NEXT_NODE_ID += 1
        self.N = N
        self.weight_graph: nx.Graph = weight_graph  # weight between vertexs
        self.flow_graph: nx.Graph = (
            self.weight_graph.copy()
        )  # to check whether DFJ constraints are violated
        self.pricing_graph: nx.Graph = (
            self.weight_graph.copy()
        )  # to solve pricing problem
        self.rmp: Model = rmp
        # key: Matching object --> value: Gurobi continuous variable
        self.matching_var_map: Dict[Matching, Var] = matching_var_map
        self.two_deg_constr: Constr = two_deg_constr
        # key: edge like (u, v) tuple format --> value: EdgeFlow constr
        self.edge_constr_map: Dict[Tuple, Constr] = edge_constr_map
        # key: Partition object --> value: DFJ subtour constraint
        self.subtour_constr_map: Dict[Partition, Constr] = subtour_constr_map
        self.iter = 0
        self.rmp_obj = 1e20
        self.mp_lb = parent_node.mp_lb if parent_node is not None else -1e20
        self.eps = 1e-5
        self.gap_tol = 1e-5
        self.solve_status = Node.SOLVE_NOT_START
        self.node_type = Node.NODE_NOT_DEFINE
        self.parent_node: "Node" = parent_node
        self.node_bound: "Node.NodeBound" = None
        self.depth = parent_node.depth + 1 if parent_node is not None else 1
        self.node_compare = node_compare

    def gen_node_bound_to_BaP(self, BaP: "BranchandPrice"):
        if self.node_bound is not None:
            self.node_bound.is_valid = False
        if (
            self.node_type == Node.NODE_INFEASIBLE
            or self.node_type == Node.NODE_INVALID
        ):
            self.node_bound = Node.NodeBound(self.mp_lb, is_feasible=False)
        else:
            self.node_bound = Node.NodeBound(self.mp_lb)
        BaP.add_node_bound(self.node_bound)

    def lazy_delete_node_bound(self):
        if self.node_bound is not None:
            self.node_bound.is_valid = False

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
            Matching(10 * sum_weight_1, edge_set1, True),
            Matching(10 * sum_weight_2, edge_set2, True),
        ]

    def build_root_model(self):
        self.rmp = Model(f"Node-{self.node_id}")
        self.rmp.setParam(GRB.Param.LogToConsole, 0)
        init_machings = self.get_init_matchings(self.N, self.weight_graph)
        for matching in init_machings:
            self.matching_var_map[matching] = self.rmp.addVar(
                lb=0,
                ub=1,
                obj=matching.sum_weight,
                vtype=GRB.CONTINUOUS,
                name=f"Matching-{matching.matching_id}",
            )
        self.two_deg_constr = self.rmp.addLConstr(
            quicksum(self.matching_var_map.values()) >= 2, name="TwoDegree"
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

    def solve(self, incumbent_obj):
        self.solve_status = Node.SOLVE_SOLVING
        if self.mp_lb >= incumbent_obj:
            self.node_type = Node.NODE_PRUNED
            return
        while True:
            self.iter += 1
            has_new_column = False
            self.rmp.optimize()
            if self.rmp.Status != GRB.Status.OPTIMAL:
                if self.rmp.Status == GRB.Status.TIME_LIMIT:
                    self.solve_status = Node.SOLVE_TIMELIMIT_ABORT
                    sol_type = self.get_sol_type()
                    if sol_type is not None:
                        self.node_type = sol_type
                    else:
                        self.node_type = Node.NODE_INVALID
                elif self.rmp.Status == GRB.Status.INFEASIBLE:
                    self.solve_status = Node.SOLVE_INFEASIBLE
                    self.node_type = Node.NODE_INFEASIBLE
                else:
                    self.solve_status = Node.SOLVE_INVALID
                    self.node_type = Node.NODE_INVALID
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
            self.solve_status = Node.SOLVE_OPTIMALITY
            if self.mp_lb >= incumbent_obj:
                self.node_type = Node.NODE_PRUNED
            else:
                self.node_type = self.get_sol_type()
            break

    def get_sol_type(self):
        if self.rmp.SolCount <= 0:
            return None
        self.reset_graph(self.flow_graph, "weight", 0)
        self.reset_graph(self.flow_graph, "matching_num", 0)
        for matching, var in self.matching_var_map.items():
            val = var.X
            if val > self.eps:
                for u, v in matching.edge_set:
                    self.flow_graph[u][v]["weight"] += val
                    self.flow_graph[u][v]["matching_num"] += 1
        sol_type = Node.NODE_INTEGER
        for u, v, data in self.flow_graph.edges(data=True):
            weight = data["weight"]
            if weight < self.eps or weight > 1 - self.eps:
                continue
            else:
                sol_type = Node.NODE_FRACTIONAL
                break
        return sol_type

    def get_sol(self):
        if self.rmp.SolCount <= 0:
            return None
        prec_map = [[] for _ in range(self.N)]
        for u, v, data in self.flow_graph.edges(data=True):
            weight = data["weight"]
            if weight > 1 - self.eps:
                prec_map[u].append(v)
                prec_map[v].append(u)
        route = []
        prev_loc = None
        cur_loc = 0
        route.append(cur_loc + 1)
        while len(route) <= self.N:
            next_locs = prec_map[cur_loc]
            for next_loc in next_locs:
                if next_loc == prev_loc:
                    continue
                route.append(next_loc + 1)
                prev_loc = cur_loc
                cur_loc = next_loc
                break
        return route

    def branch(self) -> List["Node"]:
        if self.node_type != Node.NODE_FRACTIONAL:
            return []
        branch_edge = None
        branch_edge_flow_gap = None
        branch_edge_matching_num = None
        for u, v, data in self.flow_graph.edges(data=True):
            edge = (u, v)
            weight = data["weight"]
            matching_num = data["matching_num"]
            if weight < self.eps or weight > 1 - self.eps:
                continue
            flow_gap = abs(weight - 0.5)
            if branch_edge_flow_gap is None:
                branch_edge = edge
                branch_edge_flow_gap = flow_gap
                branch_edge_matching_num = matching_num
            else:
                if flow_gap < branch_edge_flow_gap:
                    branch_edge = edge
                    branch_edge_flow_gap = flow_gap
                    branch_edge_matching_num = matching_num
                elif flow_gap == branch_edge_flow_gap:
                    if matching_num < branch_edge_matching_num:
                        branch_edge = edge
                        branch_edge_flow_gap = flow_gap
                        branch_edge_matching_num = matching_num
        nodes = [
            self.gen_remove_edge_node(branch_edge),
            self.gen_fix_edge_node(branch_edge),
        ]
        return nodes

    def gen_remove_edge_node(self, branch_edge) -> "Node":
        # remove the edge from graph
        # delete matching variables containing the edge
        # delete edge flow constraint regarding the edge

        new_weight_graph = self.weight_graph.copy()
        new_model = self.rmp.copy()
        new_weight_graph.remove_edge(branch_edge[0], branch_edge[1])

        new_matching_var_map = {}
        new_two_deg_constr = new_model.getConstrByName("TwoDegree")
        new_edge_constr_map = {}
        new_subtour_constr_map = {}
        for matching in self.matching_var_map:
            var = new_model.getVarByName(f"Matching-{matching.matching_id}")
            new_matching_var_map[matching] = var
            if matching.contains(branch_edge) and not matching.is_artificial:
                var.UB = 0
        for edge in new_weight_graph.edges():
            edge_constr = new_model.getConstrByName(f"EdgeFlow-{edge}")
            if edge_constr is None:
                edge_constr = new_model.getConstrByName(
                    f"EdgeFlow-{(branch_edge[1], branch_edge[0])}"
                )
            if (edge[0] == branch_edge[0] and edge[1] == branch_edge[1]) or (
                edge[0] == branch_edge[1] and edge[1] == branch_edge[0]
            ):
                new_model.remove(edge_constr)
            else:
                new_edge_constr_map[edge] = edge_constr
        for partition in self.subtour_constr_map:
            subtour_constr = new_model.getConstrByName(
                f"Subtour-{partition.partiton_id}"
            )
            new_subtour_constr_map[partition] = subtour_constr
        new_model.update()
        new_node = Node(
            self,
            self.N,
            new_weight_graph,
            new_model,
            new_matching_var_map,
            new_two_deg_constr,
            new_edge_constr_map,
            new_subtour_constr_map,
            self.node_compare,
        )
        return new_node

    def gen_fix_edge_node(self, branch_edge) -> "Node":
        # just modify the edge flow constraint from <= 1 to >= 1
        new_weight_graph = self.weight_graph.copy()
        new_model = self.rmp.copy()
        edge_constr = new_model.getConstrByName(f"EdgeFlow-{branch_edge}")
        if edge_constr is None:
            edge_constr = new_model.getConstrByName(
                f"EdgeFlow-{(branch_edge[1], branch_edge[0])}"
            )
        expr = new_model.getRow(edge_constr)
        new_model.remove(edge_constr)
        new_model.addLConstr(
            expr,
            GRB.GREATER_EQUAL,
            1,
            name=f"EdgeFlow-{branch_edge}",
        )
        new_model.update()

        new_matching_var_map = {}
        new_two_deg_constr = new_model.getConstrByName("TwoDegree")
        new_edge_constr_map = {}
        new_subtour_constr_map = {}
        for matching in self.matching_var_map:
            var = new_model.getVarByName(f"Matching-{matching.matching_id}")
            new_matching_var_map[matching] = var
        for edge in new_weight_graph.edges():
            edge_constr = new_model.getConstrByName(f"EdgeFlow-{edge}")
            if edge_constr is None:
                edge_constr = new_model.getConstrByName(
                    f"EdgeFlow-{(edge[1], edge[0])}"
                )
            new_edge_constr_map[edge] = edge_constr
        for partition in self.subtour_constr_map:
            subtour_constr = new_model.getConstrByName(
                f"Subtour-{partition.partiton_id}"
            )
            new_subtour_constr_map[partition] = subtour_constr
        new_model.update()
        new_node = Node(
            self,
            self.N,
            new_weight_graph,
            new_model,
            new_matching_var_map,
            new_two_deg_constr,
            new_edge_constr_map,
            new_subtour_constr_map,
            self.node_compare,
        )
        return new_node

    def __lt__(self, other: "Node"):
        # determine the order in which nodes are searched
        if self.node_compare == Node.NODE_COMPARE_BOUND:
            if self.parent_node is None and other.parent_node is not None:
                return True
            if self.parent_node is not None and other.parent_node is None:
                return False
            if self.parent_node is not None and other.parent_node is not None:
                return self.parent_node.mp_lb < other.parent_node.mp_lb
            return False
        elif self.node_compare == Node.NODE_COMPARE_DEPTH:
            if self.node_id > other.node_id:
                return True
            return False


class BranchandPrice:
    def __init__(
        self,
        N,
        graph: nx.Graph,
        node_compare=Node.NODE_COMPARE_BOUND,
    ) -> None:
        self.N = N
        self.graph = graph
        self.incumbent_obj = 1e20
        self.incumbent_solution = None
        self.bound = -1e20  # bound of BaP = lowest bound among all leaf nodes
        self.infeasible = False
        self.node_bounds: List[Node.NodeBound] = []
        self.to_search_nodes: List[Node] = []
        self.searched_nodes: List[Node] = []
        self.eps = 1e-5
        self.gap_tol = 1e-5
        self.iter = 0
        self.node_compare = node_compare
        self.log_template = "{iter:<10d}{node_id:<10d}{pnode_id:<10}{incum_obj:<15.6e}{bn:<15.6e}{gap:<10.2f}{node_obj:<15.6e}{node_lb:<15.6e}{node_type:<15}{time:<10.2f}{remain_node:<15d}{searched_node:<15d}"

    def add_node_bound(self, node_bound):
        heapq.heappush(self.node_bounds, node_bound)

    def cal_gap(self, obj, bound):
        return 100 * abs(obj - bound) / max((abs(obj), abs(bound), self.eps))

    def search(self):
        log_header = "{iter:<10}{node_id:<10}{pnode_id:<10}{incum_obj:<15}{bn:<15}{gap:<10}{node_obj:<15}{node_lb:<15}{node_type:<15}{time:<10}{remain_node:<15}{searched_node:<15}".format(
            iter="Iter.",
            node_id="Node Id",
            pnode_id="PNode Id",
            incum_obj="Incum Obj.",
            bn="Bound",
            gap="Gap(%)",
            node_obj="Node Obj.",
            node_lb="Node LB",
            node_type="Node Type",
            time="Time(s)",
            remain_node="#Remain. Node",
            searched_node="#Searched Node",
        )
        print(log_header)
        st = time.time()
        root_node = Node(None, N, graph, node_compare=self.node_compare)
        root_node.build_root_model()
        self.to_search_nodes.append(root_node)
        root_node.gen_node_bound_to_BaP(self)
        while self.to_search_nodes:
            self.iter += 1
            node: Node = heapq.heappop(self.to_search_nodes)
            node.solve(self.incumbent_obj)
            node.gen_node_bound_to_BaP(self)
            while self.node_bounds:
                top_node_bound = self.node_bounds[0]
                if not top_node_bound.is_valid:
                    heapq.heappop(self.node_bounds)
                    continue
                self.bound = max(top_node_bound.lb, self.bound)
                if not top_node_bound.is_feasible:
                    self.infeasible = True
                break
            if self.infeasible:
                print("Problem is infeasible")
                return
            if node.node_type == Node.NODE_INTEGER:
                if self.incumbent_obj > node.rmp_obj:
                    self.incumbent_obj = node.rmp_obj
                    self.incumbent_solution = node.get_sol()
            elif node.node_type == Node.NODE_FRACTIONAL:
                new_nodes = node.branch()
                for new_node in new_nodes:
                    heapq.heappush(self.to_search_nodes, new_node)
                    new_node.gen_node_bound_to_BaP(self)
            if node.node_type != Node.NODE_INTEGER:
                node.lazy_delete_node_bound()
            self.searched_nodes.append(node)

            cur_node_id = node.node_id
            parent_node_id = (
                node.parent_node.node_id if node.parent_node is not None else 0
            )
            gap = self.cal_gap(self.incumbent_obj, self.bound)
            node_obj = node.rmp_obj
            node_lb = node.mp_lb
            node_type_repr = Node.NODE_TYPE_REPR_MAP[node.node_type]
            t = time.time() - st
            remain_node_num = len(self.to_search_nodes)
            search_node_num = len(self.searched_nodes)
            iter_info = self.log_template.format(
                iter=self.iter,
                node_id=cur_node_id,
                pnode_id=parent_node_id,
                incum_obj=self.incumbent_obj,
                bn=self.bound,
                gap=gap,
                node_obj=node_obj,
                node_lb=node_lb,
                node_type=node_type_repr,
                time=t,
                remain_node=remain_node_num,
                searched_node=search_node_num,
            )
            print(iter_info)
            if gap < self.gap_tol:
                break
        if self.incumbent_solution is not None:
            print(f"Search finished. Incumbent solution: {self.incumbent_solution}")
        else:
            print("Search finished. Not found feasible solution")


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
    node_compare = Node.NODE_COMPARE_BOUND
    # node_compare = Node.NODE_COMPARE_DEPTH
    BaP = BranchandPrice(N, graph, node_compare=node_compare)
    BaP.search()
