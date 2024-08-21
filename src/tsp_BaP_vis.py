from gurobipy import Model, GRB, quicksum, Var, Constr
import networkx as nx
from typing import Dict, List, Tuple, Set, Union
import time
import heapq
import math
import tkinter as tk
from tkinter import ttk, Canvas, messagebox
import threading
import platform
import os
import traceback

OS_NAME = platform.system()

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
        evt_queue: List["BranchandPrice.Event"],
        rmp: Model = None,
        matching_var_map={},
        two_deg_constr=None,
        edge_constr_map={},
        subtour_constr_map={},
        preferred_col_offset=0,
        node_compare=0,
    ) -> None:
        self.evt_queue = evt_queue
        self.node_id = Node.NEXT_NODE_ID
        Node.NEXT_NODE_ID += 1
        self.N = N
        self.weight_graph: nx.Graph = weight_graph
        self.flow_graph: nx.Graph = self.weight_graph.copy()
        self.pricing_graph: nx.Graph = self.weight_graph.copy()
        self.rmp: Model = rmp
        # key: Matching object --> value: Continuous Var
        self.matching_var_map: Dict[Matching, Var] = matching_var_map
        self.two_deg_constr: Constr = two_deg_constr
        # key: edge --> value: EdgeFlow constraint
        self.edge_constr_map: Dict[Tuple, Constr] = edge_constr_map
        # key: Partition object --> value: Subtour constraint
        self.subtour_constr_map: Dict[Partition, Constr] = subtour_constr_map
        self.iter = 0
        self.iter_log_size_map = [0]
        self.init_column_num = len(self.matching_var_map)
        self.added_column_num = 0
        self.init_subtour_cut_num = len(self.subtour_constr_map)
        self.added_subtour_cut_num = 0
        self.rmp_obj = 1e20
        self.mp_lb = parent_node.mp_lb if parent_node is not None else -1e20
        self.t_master = 0
        self.t_pricing = 0
        self.t_detect_subtour = 0
        self.t_total = 0
        self.eps = 1e-5
        self.gap_tol = 1e-5
        self.log_history = []
        self.log_template = "{iter:<10d}{rmp_obj:<15.6e}{mp_lb:<15.6e}{min_rc:<15.6e}{gap:<10.2f}{time:<10.2f}\n"
        self.solve_status = Node.SOLVE_NOT_START
        self.node_type = Node.NODE_NOT_DEFINE
        self.parent_node: "Node" = parent_node
        self.node_bound: "Node.NodeBound" = None
        self.depth = parent_node.depth + 1 if parent_node is not None else 1
        self.preferred_col_offset = preferred_col_offset
        self.node_compare = node_compare
        BranchandPrice.Event.merge_evt_queue(
            self.evt_queue,
            BranchandPrice.Event(
                BranchandPrice.EVT_NEW_NODE,
                BranchandPrice.EventNodeData(self),
            ),
        )

    def gen_node_bound_to_BaP(self, bbs: "BranchandPrice"):
        if self.node_bound is not None:
            self.node_bound.is_valid = False
        if self.node_type == Node.NODE_INFEASIBLE:
            self.node_bound = Node.NodeBound(self.mp_lb, is_feasible=False)
        else:
            self.node_bound = Node.NodeBound(self.mp_lb)
        bbs.add_node_bound(self.node_bound)

    def lazy_delete_node_bound(self):
        if self.node_bound is not None:
            self.node_bound.is_valid = False

    def cal_gap(self, obj, bound):
        return 100 * abs(obj - bound) / max((abs(obj), abs(bound), self.eps))

    def reset_graph(self, graph: nx.Graph, attr="weight", attr_val=0):
        for _u, _v, data in graph.edges(data=True):
            data[attr] = attr_val

    def copy_weight(self, dest_graph: nx.Graph, src_graph: nx.Graph, opposite=False):
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
        self.init_column_num = len(init_machings)
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
        self.init_subtour_cut_num = 0

    def solve(self, incumbent_obj):
        init_st = time.time()
        self.solve_status = Node.SOLVE_SOLVING
        BranchandPrice.Event.merge_evt_queue(
            self.evt_queue,
            BranchandPrice.Event(
                BranchandPrice.EVT_NODE_UPDATE,
                BranchandPrice.EventNodeData(self),
            ),
        )
        if self.mp_lb >= incumbent_obj:
            self.node_type = Node.NODE_PRUNED
            self.t_total = time.time() - init_st
            BranchandPrice.Event.merge_evt_queue(
                self.evt_queue,
                BranchandPrice.Event(
                    BranchandPrice.EVT_NODE_UPDATE,
                    BranchandPrice.EventNodeData(self),
                ),
            )
            return
        while True:
            self.iter += 1
            has_new_column = False
            st = time.time()
            self.rmp.optimize()
            self.t_master += time.time() - st
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
                self.t_total = time.time() - init_st
                BranchandPrice.Event.merge_evt_queue(
                    self.evt_queue,
                    BranchandPrice.Event(
                        BranchandPrice.EVT_NODE_UPDATE,
                        BranchandPrice.EventNodeData(self),
                    ),
                )
                return
            self.rmp_obj = self.rmp.ObjVal
            gap_rel = self.cal_gap(self.rmp.ObjVal, self.mp_lb)
            cur_mp_lb = self.mp_lb
            reduced_cost = 0
            if gap_rel > self.gap_tol:
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
                st = time.time()
                max_weight_matching = nx.max_weight_matching(
                    self.pricing_graph, maxcardinality=True
                )
                self.t_pricing += time.time() - st
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
                    self.added_column_num += 1
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
            et = time.time()
            self.log_history.append(
                self.log_template.format(
                    iter=self.iter,
                    rmp_obj=self.rmp.ObjVal,
                    mp_lb=self.mp_lb,
                    min_rc=reduced_cost,
                    gap=gap_rel,
                    time=et - init_st,
                )
            )
            self.iter_log_size_map.append(len(self.log_history))
            self.t_total = et - init_st
            BranchandPrice.Event.merge_evt_queue(
                self.evt_queue,
                BranchandPrice.Event(
                    BranchandPrice.EVT_NODE_UPDATE,
                    BranchandPrice.EventNodeData(self),
                ),
            )
            if has_new_column:
                continue
            self.reset_graph(self.flow_graph)
            for matching, var in self.matching_var_map.items():
                val = var.X
                if val > self.eps:
                    for u, v in matching.edge_set:
                        self.flow_graph[u][v]["weight"] += val
            st = time.time()
            cut_val, partition = nx.stoer_wagner(self.flow_graph)
            et = time.time()
            self.t_detect_subtour += et - st
            if cut_val < 2 - self.eps:
                self.log_history.append(
                    "Find vertex subset that violates DFJ constraint\n",
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
                self.added_subtour_cut_num += 1
                BranchandPrice.Event.merge_evt_queue(
                    self.evt_queue,
                    BranchandPrice.Event(
                        BranchandPrice.EVT_NODE_UPDATE,
                        BranchandPrice.EventNodeData(self),
                    ),
                )
                continue
            self.solve_status = Node.SOLVE_OPTIMALITY
            if self.mp_lb >= incumbent_obj:
                self.node_type = Node.NODE_PRUNED
            else:
                self.node_type = self.get_sol_type()
            self.t_total = time.time() - init_st
            BranchandPrice.Event.merge_evt_queue(
                self.evt_queue,
                BranchandPrice.Event(
                    BranchandPrice.EVT_NODE_UPDATE,
                    BranchandPrice.EventNodeData(self),
                ),
            )
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
            self.evt_queue,
            new_model,
            new_matching_var_map,
            new_two_deg_constr,
            new_edge_constr_map,
            new_subtour_constr_map,
            -1,
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
            self.evt_queue,
            new_model,
            new_matching_var_map,
            new_two_deg_constr,
            new_edge_constr_map,
            new_subtour_constr_map,
            1,
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
    EVT_NEW_NODE = 0
    EVT_NODE_UPDATE = 1
    EVT_BAP_NEW_LOG = 2

    SEARCH_STRA_BFS = 0
    SEARCH_STRA_DFS = 1

    class EventNodeData:
        def __init__(self, node: Node) -> None:
            self.node_id = node.node_id
            self.parent_node_id = (
                node.parent_node.node_id if node.parent_node is not None else None
            )
            self.node_type = node.node_type
            self.solve_status = node.solve_status
            self.depth = node.depth
            self.node_obj = node.rmp_obj
            self.node_lb = node.mp_lb
            self.gap = node.cal_gap(node.rmp_obj, node.mp_lb)
            self.init_cols = node.init_column_num
            self.added_cols = node.added_column_num
            self.init_cuts = node.init_subtour_cut_num
            self.added_cuts = node.added_subtour_cut_num
            self.t_master = node.t_master
            self.t_pricing = node.t_pricing
            self.t_detect_cut = node.t_detect_subtour
            self.t_total = node.t_total
            self.iter = node.iter
            self.log_history = node.log_history
            self.iter_log_size_map = node.iter_log_size_map
            self.preferred_col_offset = node.preferred_col_offset

    class Event:
        NEXT_EVENT_ID = 1

        def __init__(self, evt_code, evt_data: Union[str, Node]) -> None:
            self.evt_code = evt_code
            self.evt_data = evt_data
            self.evt_id = BranchandPrice.Event.NEXT_EVENT_ID
            BranchandPrice.Event.NEXT_EVENT_ID += 1

        def merge_evt_queue(
            evt_queue: List["BranchandPrice.Event"],
            evt: "BranchandPrice.Event",
        ):
            if not evt_queue:
                evt_queue.append(evt)
            else:
                prev_evt = evt_queue[-1]
                if (
                    prev_evt.evt_code < BranchandPrice.EVT_BAP_NEW_LOG
                    and evt.evt_code < BranchandPrice.EVT_BAP_NEW_LOG
                    and prev_evt.evt_data.node_id == evt.evt_data.node_id
                ):
                    evt_queue[-1] = evt
                else:
                    evt_queue.append(evt)

    def __init__(
        self,
        N,
        graph: nx.Graph,
        evt_queue: List["BranchandPrice.Event"],
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
        self.node_map: Dict[str, Node] = {}
        self.eps = 1e-5
        self.gap_tol = 1e-5
        self.iter = 0
        self.evt_queue = evt_queue
        self.node_compare = node_compare
        self.log_template = "{iter:<10d}{node_id:<10d}{pnode_id:<10}{incum_obj:<15.6e}{bn:<15.6e}{gap:<10.2f}{node_obj:<15.6e}{node_lb:<15.6e}{node_type:<15}{time:<10.2f}{remain_node:<15d}{searched_node:<15d}\n"

    def add_node_bound(self, node_bound):
        heapq.heappush(self.node_bounds, node_bound)

    def cal_gap(self, obj, bound):
        return 100 * abs(obj - bound) / max((abs(obj), abs(bound), self.eps))

    def search(self):
        st = time.time()
        root_node = Node(None, N, graph, evt_queue, node_compare=self.node_compare)
        root_node.build_root_model()
        self.to_search_nodes.append(root_node)
        self.node_map[root_node.node_id] = root_node
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
                messagebox.showinfo("Search finished", "Problem is infeasible")
                return
            if node.node_type == Node.NODE_INTEGER:
                if self.incumbent_obj > node.rmp_obj:
                    self.incumbent_obj = node.rmp_obj
                    self.incumbent_solution = node.get_sol()
            elif node.node_type == Node.NODE_FRACTIONAL:
                new_nodes = node.branch()
                for new_node in new_nodes:
                    heapq.heappush(self.to_search_nodes, new_node)
                    self.node_map[new_node.node_id] = new_node
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
            BranchandPrice.Event.merge_evt_queue(
                self.evt_queue,
                BranchandPrice.Event(BranchandPrice.EVT_BAP_NEW_LOG, iter_info),
            )
            print(iter_info, end="")
            if gap < self.gap_tol:
                break
        if self.incumbent_solution:
            messagebox.showinfo(
                "Search finished", f"Incumbent solution: {self.incumbent_solution}"
            )
            print(f"Search finished. Incumbent solution: {self.incumbent_solution}")
        else:
            messagebox.showinfo("Search finished", "Not found feasible solution")
            print("Search finished. Not found feasible solution")


class BranchandPriceTreeCanvas(Canvas):
    def __init__(self, root, BaP_vis: "BranchandPriceVisual", *args, **kwargs):
        super().__init__(root, *args, **kwargs)
        self.BaP_vis: "BranchandPriceVisual" = BaP_vis
        self.node_square_length = 230
        self.original_font_size = 12
        self.current_font_size = self.original_font_size
        self.unit_dy = 300
        self.unit_dx = 240
        self.default_node_row = 1
        self.default_node_col = 3
        self.default_node_x1 = self.node_square_length * (self.default_node_col - 1)
        self.default_ndoe_y1 = self.node_square_length * self.default_node_row
        self.cur_scale = 1
        self.scale_factor = 1.0

        self.COLOR_FRACTIONAL = ["grey", "white"]
        self.COLOR_INTEGER = ["green", "white"]
        self.COLOR_PRUNED = ["red", "white"]
        self.COLOR_INFEAS_INVALID = ["black", "white"]
        self.COLOR_SOLVING = ["blue", "white"]
        self.COLOR_NOT_START = ["white", "black"]

        self.node_row_col_map = {}
        self.row_free_cols_map = {}
        self.label_template = (
            "Node Id: {node_id:<10d}\n"
            "Node Depth: {node_depth:<10d}\n"
            "Node Obj.: {node_obj:<10.6e}\n"
            "Node LB: {node_lb:<10.6e}\n"
            "Gap(%): {gap:>5.2f}\n"
            "Solve Stat: {solve_stat:<10}\n"
            "Node Type: {node_type:<10}\n"
            "Time(s): {time:<10.2f}"
        )

        self.bind("<MouseWheel>", self.zoom)
        self.pan_start_x = 0
        self.pan_start_y = 0
        self.bind("<ButtonPress-3>", self.start_pan)
        self.bind("<ButtonRelease-3>", self.end_pan)
        self.tag_bind("square", "<Button-1>", self.on_node_click)
        self.tag_bind("label", "<Button-1>", self.on_node_click)

        self.is_lock_solving = True
        self.lock_button = tk.Button(
            self,
            text="To Untrack Solving Node",
            command=self.switch_lock,
            font=("Courier", 12),
        )
        self.create_window(10, 10, anchor="nw", window=self.lock_button)

        self.locate_frame = tk.Frame(self)
        self.locate_button = tk.Button(
            self.locate_frame,
            text="Locate Node Id",
            command=self.locate_node,
            font=("Courier", 12),
        )
        self.locate_entry = tk.Entry(self.locate_frame, font=("Courier", 12))
        self.locate_button.pack(side=tk.LEFT)
        self.locate_entry.pack(side=tk.LEFT)
        self.create_window(10, 50, anchor="nw", window=self.locate_frame)

        self.cur_solving_node_id = None
        self.cur_focus_node_id = None

    def locate_node(self):
        try:
            node_id = int(self.locate_entry.get().strip())
        except Exception:
            messagebox.showerror("Invalid Input", "Please enter a valid integer")
            return
        if node_id not in self.node_row_col_map:
            messagebox.showerror("Invalid Input", "No such node id")
            return
        self.is_lock_solving = False
        if self.is_lock_solving:
            self.lock_button.config(text="To Untrack Solving Node")
        else:
            self.lock_button.config(text="To Track Solving Node")
        if node_id != self.cur_focus_node_id:
            self.cur_focus_node_id = node_id
            self.focus_node(node_id)
            self.BaP_vis.update_node_detail(node_id=node_id, is_node_chg=True)

    def switch_lock(self):
        self.is_lock_solving = not self.is_lock_solving
        if self.is_lock_solving:
            self.lock_button.config(text="To Untrack Solving Node")
        else:
            self.lock_button.config(text="To Track Solving Node")
        if self.is_lock_solving and self.cur_focus_node_id != self.cur_solving_node_id:
            self.cur_focus_node_id = self.cur_solving_node_id
            self.focus_node(self.cur_solving_node_id)
            self.BaP_vis.update_node_detail(
                node_id=self.cur_solving_node_id, is_node_chg=True
            )

    def find_col_with_inplace_chg(
        self, preferred_col, free_cols: List[List], break_tie_to_left
    ):
        left_idx = 0
        right_idx = len(free_cols) - 1
        while left_idx <= right_idx:
            mid_idx = (left_idx + right_idx) // 2
            mid_col_intv = free_cols[mid_idx]
            if mid_col_intv[1] < preferred_col:
                left_idx = mid_idx + 1
            else:
                if left_idx == right_idx:
                    break
                right_idx = mid_idx
        if left_idx >= len(free_cols):
            raise ValueError("Cannot find free canvas col")
        col_intv_idx = left_idx
        col_intv = free_cols[col_intv_idx]
        if left_idx > 0:
            alt_col_intv = free_cols[left_idx - 1]
            dis_to_col_intv = (
                col_intv[0] - preferred_col if preferred_col < col_intv[0] else 0
            )
            dis_to_alt_col_intv = preferred_col - alt_col_intv[1]
            if dis_to_alt_col_intv < dis_to_col_intv or (
                dis_to_alt_col_intv == dis_to_col_intv and break_tie_to_left
            ):
                col_intv = alt_col_intv
                col_intv_idx = left_idx - 1
        found_col = preferred_col
        if preferred_col < col_intv[0]:
            found_col = col_intv[0]
        if preferred_col > col_intv[1]:
            found_col = col_intv[1]
        free_cols.pop(col_intv_idx)
        if found_col < col_intv[1]:
            free_cols.insert(col_intv_idx, [found_col + 1, col_intv[1]])
        if found_col > col_intv[0]:
            free_cols.insert(col_intv_idx, [col_intv[0], found_col - 1])
        return found_col

    def get_square_label_color(self, node_type, solve_status):
        if node_type == Node.NODE_FRACTIONAL:
            return self.COLOR_FRACTIONAL
        elif node_type == Node.NODE_INTEGER:
            return self.COLOR_INTEGER
        elif node_type == Node.NODE_PRUNED:
            return self.COLOR_PRUNED
        elif node_type == Node.NODE_INFEASIBLE or node_type == Node.NODE_INVALID:
            return self.COLOR_INFEAS_INVALID
        else:
            if solve_status == Node.SOLVE_NOT_START:
                return self.COLOR_NOT_START
            elif solve_status == Node.SOLVE_SOLVING:
                return self.COLOR_SOLVING
        return self.COLOR_NOT_START

    def update_node(self, evt_data: BranchandPrice.EventNodeData):
        node_id = evt_data.node_id
        square_tag = f"square-{node_id}"
        label_tag = f"label-{node_id}"
        if node_id not in self.node_row_col_map:
            self.create_node(evt_data)
        else:
            square_color, label_color = self.get_square_label_color(
                evt_data.node_type, evt_data.solve_status
            )
            label_text = self.label_template.format(
                node_id=node_id,
                node_depth=evt_data.depth,
                node_obj=evt_data.node_obj,
                node_lb=evt_data.node_lb,
                gap=evt_data.gap,
                solve_stat=Node.SOLVE_STATUS_REPR_MAP[evt_data.solve_status],
                node_type=Node.NODE_TYPE_REPR_MAP[evt_data.node_type],
                time=evt_data.t_total,
            )
            self.itemconfig(square_tag, fill=square_color)
            self.itemconfig(label_tag, fill=label_color, text=label_text)
        self.cur_solving_node_id = node_id
        is_cur_focus_node_changed = False
        if self.is_lock_solving and self.cur_focus_node_id != node_id:
            is_cur_focus_node_changed = True
            self.cur_focus_node_id = node_id
            self.focus_node(node_id)
        if is_cur_focus_node_changed:
            self.BaP_vis.update_node_detail(evt_data=evt_data, is_node_chg=True)
        elif self.cur_focus_node_id == node_id:
            self.BaP_vis.update_node_detail(evt_data=evt_data, is_node_chg=False)

    def create_node(self, evt_data: BranchandPrice.EventNodeData):
        if self.cur_solving_node_id is None:
            self.cur_solving_node_id = evt_data.node_id
            self.cur_focus_node_id = evt_data.node_id
        node_id = evt_data.node_id
        parent_node_id = evt_data.parent_node_id
        if parent_node_id is None:
            node_row = self.default_node_row
            node_col = self.default_node_col
            child_x1 = self.default_node_x1
            child_y1 = self.default_ndoe_y1
        else:
            parent_row, parent_col = self.node_row_col_map[parent_node_id]
            node_row = parent_row + 1
            preferred_col = parent_col + evt_data.preferred_col_offset
            if node_row not in self.row_free_cols_map:
                self.row_free_cols_map[node_row] = [[-1e20, 1e20]]
            free_cols = self.row_free_cols_map[node_row]
            node_col = self.find_col_with_inplace_chg(
                preferred_col, free_cols, evt_data.preferred_col_offset > 0
            )
            parent_square_tag = f"square-{parent_node_id}"
            x1, y1, x2, y2 = self.coords(parent_square_tag)
            dx = self.unit_dx * (node_col - parent_col)
            dy = self.unit_dy
            dx *= self.cur_scale
            dy *= self.cur_scale
            child_x1 = x1 + dx
            child_y1 = y1 + dy
        child_x2, child_y2 = (
            child_x1 + self.node_square_length * self.cur_scale,
            child_y1 + self.node_square_length * self.cur_scale,
        )
        square_color, label_color = self.get_square_label_color(
            evt_data.node_type, evt_data.solve_status
        )
        self.create_rectangle(
            child_x1,
            child_y1,
            child_x2,
            child_y2,
            outline="black",
            fill=square_color,
            tags=["node", "square", f"square-{node_id}"],
        )

        label_text = self.label_template.format(
            node_id=node_id,
            node_depth=evt_data.depth,
            node_obj=evt_data.node_obj,
            node_lb=evt_data.node_lb,
            gap=evt_data.gap,
            solve_stat=Node.SOLVE_STATUS_REPR_MAP[evt_data.solve_status],
            node_type=Node.NODE_TYPE_REPR_MAP[evt_data.node_type],
            time=evt_data.t_total,
        )
        self.create_text(
            (child_x1 + child_x2) / 2,
            (child_y1 + child_y2) / 2,
            text=label_text,
            font=("Courier", self.current_font_size),
            fill=label_color,
            tags=["node", "label", f"label-{node_id}"],
        )
        if parent_node_id is not None:
            parent_bottom_center_x = (x1 + x2) // 2
            parent_bottom_center_y = y2
            child_upper_center_x = parent_bottom_center_x + dx
            child_upper_center_y = child_y1
            self.create_line(
                parent_bottom_center_x,
                parent_bottom_center_y,
                child_upper_center_x,
                child_upper_center_y,
                fill="black",
                tags=["node", "line", f"line-{parent_node_id}-{node_id}"],
            )
        self.node_row_col_map[node_id] = (node_row, node_col)

    def zoom(self, event):
        scale = 1.0
        if event.delta > 0:
            scale = 1.1
        elif event.delta < 0:
            scale = 0.9
        self.cur_scale *= scale

        self.scale("node", event.x, event.y, scale, scale)

        self.current_font_size = max(1, int(self.original_font_size * self.cur_scale))
        self.itemconfig("label", font=("Courier New", self.current_font_size))

    def start_pan(self, event):
        self.is_lock_solving = False
        self.lock_button.config(text="To Track Solving Node")
        self.pan_start_x = event.x
        self.pan_start_y = event.y
        self.config(cursor="fleur")

    def end_pan(self, event):
        dx = event.x - self.pan_start_x
        dy = event.y - self.pan_start_y
        self.move("node", dx, dy)
        self.config(cursor="")

    def on_node_click(self, event):
        item_id = event.widget.find_closest(event.x, event.y)

        tags = self.gettags(item_id)

        node_id = None
        for tag in tags:
            if tag.startswith("label-"):
                node_id = int(tag[6:])
                break
            if tag.startswith("square-"):
                node_id = int(tag[7:])
                break
        if node_id is not None:
            self.is_lock_solving = False
            self.lock_button.config(text="To Track Solving Node")
            if node_id != self.cur_focus_node_id:
                self.cur_focus_node_id = node_id
                self.BaP_vis.update_node_detail(node_id=node_id, is_node_chg=True)
            else:
                self.BaP_vis.update_node_detail(node_id=node_id, is_node_chg=False)

        dot_size = 5
        dot_id = self.create_oval(
            event.x - dot_size,
            event.y - dot_size,
            event.x + dot_size,
            event.y + dot_size,
            fill="red",
            outline="red",
            tags="click_dot",
        )

        self.after(300, self.delete, dot_id)

    def focus_node(self, node_id):
        square_tag = f"square-{node_id}"
        x1, y1, x2, y2 = self.coords(square_tag)
        dx = self.default_node_x1 - x1
        dy = self.default_ndoe_y1 - y1
        self.move("node", dx, dy)


class CustomTk(tk.Tk):
    def report_callback_exception(self, exc_type, exc_value, exc_traceback):
        print("Exception in Custom Tkinter callback")
        traceback.print_exception(exc_type, exc_value, exc_traceback)
        os._exit(1)


class BranchandPriceVisual:
    def __init__(
        self, BaP: BranchandPrice, evt_queue: List[BranchandPrice.Event]
    ) -> None:
        self.BaP = BaP
        self.root = CustomTk()
        self.root.title("Branch-and-Price Procedure")
        if OS_NAME in ["Windows", "Darwin"]:
            self.root.state("zoomed")
        elif OS_NAME == "Linux":
            self.root.wm_attributes("-zoomed", True)
        else:
            raise ValueError(f"Unknown OS {OS_NAME}")
        self.panedwindow = ttk.PanedWindow(self.root, orient="horizontal")
        self.panedwindow.pack(fill=tk.BOTH, expand=True)
        self.BaP_frame = ttk.Frame(self.panedwindow, relief="ridge")
        self.panedwindow.add(self.BaP_frame, weight=3)
        self.BaP_panedwindow = ttk.PanedWindow(self.BaP_frame, orient="vertical")
        self.BaP_panedwindow.pack(fill=tk.BOTH, expand=True)
        self.BaP_graph_frame = ttk.Frame(self.BaP_panedwindow, relief="groove")
        self.BaP_panedwindow.add(self.BaP_graph_frame, weight=2)
        self.BaP_log_frame = ttk.Frame(self.BaP_panedwindow, relief="groove")
        self.BaP_log_frame.pack(expand=True, fill=tk.BOTH)
        self.BaP_panedwindow.add(self.BaP_log_frame, weight=1)
        self.BaP_log_header_text = tk.Text(
            self.BaP_log_frame,
            wrap=tk.NONE,
            width=0,
            height=0,
            state="disabled",
            font=("Courier", 12),
            bg="#F5DEB3",
        )
        BaP_log_header = "{iter:<10}{node_id:<10}{pnode_id:<10}{incum_obj:<15}{bn:<15}{gap:<10}{node_obj:<15}{node_lb:<15}{node_type:<15}{time:<10}{remain_node:<15}{searched_node:<15}".format(
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
        print(BaP_log_header)
        self.add_text(self.BaP_log_header_text, BaP_log_header)
        self.BaP_log_text = tk.Text(
            self.BaP_log_frame,
            wrap=tk.NONE,
            width=0,
            height=0,
            state="disabled",
            font=("Courier", 12),
            bg="#F5DEB3",
        )
        self.BaP_log_vsb = tk.Scrollbar(
            self.BaP_log_frame, orient="vertical", command=self.BaP_log_text.yview
        )
        self.BaP_log_hsb = tk.Scrollbar(
            self.BaP_log_frame,
            orient="horizontal",
            command=self.on_BaP_log_horizon_scroll,
        )
        self.BaP_log_header_text.configure(xscrollcommand=self.BaP_log_hsb.set)
        self.BaP_log_text.configure(
            yscrollcommand=self.BaP_log_vsb.set, xscrollcommand=self.BaP_log_hsb.set
        )
        self.BaP_log_header_text.grid(row=0, column=0, sticky="nsew")
        self.BaP_log_text.grid(row=1, column=0, sticky="nsew")
        self.BaP_log_vsb.grid(row=1, column=1, sticky="ns")
        self.BaP_log_hsb.grid(row=2, column=0, sticky="ew")
        self.BaP_log_frame.grid_rowconfigure(0, weight=1)
        self.BaP_log_frame.grid_rowconfigure(1, weight=20)
        self.BaP_log_frame.grid_columnconfigure(0, weight=1)

        self.node_frame = ttk.Frame(self.panedwindow, relief="ridge")
        self.panedwindow.add(self.node_frame, weight=2)
        self.node_frame_panedwindow = ttk.PanedWindow(
            self.node_frame, orient="vertical"
        )
        self.node_frame_panedwindow.pack(fill=tk.BOTH, expand=True)
        self.node_info_frame = ttk.Frame(self.node_frame_panedwindow, relief="groove")
        self.node_frame_panedwindow.add(self.node_info_frame, weight=2)
        self.node_info_text = tk.Text(
            self.node_info_frame,
            wrap=tk.NONE,
            width=0,
            height=0,
            state="disabled",
            font=("Courier", 12),
            bg="#F5DEB3",
        )
        self.node_info_vsb = tk.Scrollbar(
            self.node_info_frame, orient="vertical", command=self.node_info_text.yview
        )
        self.node_info_hsb = tk.Scrollbar(
            self.node_info_frame, orient="horizontal", command=self.node_info_text.xview
        )
        self.node_info_text.configure(
            xscrollcommand=self.node_info_hsb.set, yscrollcommand=self.node_info_vsb.set
        )
        self.node_info_text.grid(row=0, column=0, sticky="nsew")
        self.node_info_vsb.grid(row=0, column=1, sticky="ns")
        self.node_info_hsb.grid(row=1, column=0, sticky="ew")
        self.node_info_frame.grid_rowconfigure(0, weight=1)
        self.node_info_frame.grid_columnconfigure(0, weight=1)
        self.node_detail_text_template = (
            "Node Id: {node_id:<10d}\n"
            "Node Depth: {node_depth:<10d}\n"
            "Node Obj.: {node_obj:<10.6e}\n"
            "Node LB: {node_lb:<10.6e}\n"
            "Gap(%): {gap:>5.2f}\n"
            "Solve Stat: {solve_stat:<10}\n"
            "Node Type: {node_type:<10}\n"
            "Iter: {iter:<10}\n"
            "#Init Cols: {init_cols:<10d}\n"
            "#Added Cols: {added_cols:<10d}\n"
            "#Init Cuts: {init_cuts:<10d}\n"
            "#Added Cuts: {added_cuts:<10d}\n"
            "Master Time(s): {t_master:<10.2f}\n"
            "Pricing Time(s): {t_pricing:<10.2f}\n"
            "Detect Cut Time(s): {t_detect_cut:<10.2f}\n"
            "Total Time(s): {t_total:<10.2f}"
        )
        self.node_log_frame = ttk.Frame(self.node_frame_panedwindow, relief="groove")
        self.node_frame_panedwindow.add(self.node_log_frame, weight=2)
        self.node_log_header_text = tk.Text(
            self.node_log_frame,
            wrap=tk.NONE,
            width=0,
            height=0,
            state="disabled",
            font=("Courier", 12),
            bg="#F5DEB3",
        )
        node_log_header = (
            "{iter:<10}{rmp_obj:<15}{mp_lb:<15}{min_rc:<15}{gap:<10}{time:<10}".format(
                iter="Iter.",
                rmp_obj="RMP Obj.",
                mp_lb="MP LB",
                min_rc="Min RC",
                gap="Gap(%)",
                time="Time(s)",
            )
        )
        self.node_log_text = tk.Text(
            self.node_log_frame,
            wrap=tk.NONE,
            width=0,
            height=0,
            state="disabled",
            font=("Courier", 12),
            bg="#F5DEB3",
        )
        self.add_text(self.node_log_header_text, node_log_header)
        self.node_log_vsb = tk.Scrollbar(
            self.node_log_frame, orient="vertical", command=self.node_log_text.yview
        )
        self.node_log_hsb = tk.Scrollbar(
            self.node_log_frame,
            orient="horizontal",
            command=self.on_node_log_horizon_scroll,
        )
        self.node_log_header_text.configure(xscrollcommand=self.node_log_hsb.set)
        self.node_log_text.configure(
            yscrollcommand=self.node_log_vsb.set, xscrollcommand=self.node_log_hsb.set
        )
        self.node_log_header_text.grid(row=0, column=0, sticky="nsew")
        self.node_log_text.grid(row=1, column=0, sticky="nsew")
        self.node_log_vsb.grid(row=1, column=1, sticky="ns")
        self.node_log_hsb.grid(row=2, column=0, sticky="ew")
        self.node_log_frame.grid_rowconfigure(0, weight=1)
        self.node_log_frame.grid_rowconfigure(1, weight=20)
        self.node_log_frame.grid_columnconfigure(0, weight=1)

        self.canvas = BranchandPriceTreeCanvas(self.BaP_graph_frame, self, bg="#F5DEB3")
        self.canvas.pack(expand=True, fill=tk.BOTH)

        self.evt_queue = evt_queue
        self.cur_evt_id = 0
        self.cur_evt_idx = 0

    def on_node_log_horizon_scroll(self, *args):
        self.node_log_header_text.xview(*args)
        self.node_log_text.xview(*args)

    def on_BaP_log_horizon_scroll(self, *args):
        self.BaP_log_header_text.xview(*args)
        self.BaP_log_text.xview(*args)

    def add_text(self, tk_text: tk.Text, text_str: str):
        cur_text_yview = tk_text.yview()[1]
        scroll_to_bottom = False
        if cur_text_yview > 0.99:
            scroll_to_bottom = True
        tk_text.config(state="normal")
        tk_text.insert(tk.END, text_str)
        tk_text.config(state="disabled")
        if scroll_to_bottom:
            tk_text.yview(tk.MOVETO, 1)

    def update_node_detail(
        self,
        evt_data: BranchandPrice.EventNodeData = None,
        node_id=None,
        is_node_chg=False,
    ):
        if node_id is not None:
            evt_data = BranchandPrice.EventNodeData(self.BaP.node_map[node_id])
        new_text = self.node_detail_text_template.format(
            node_id=evt_data.node_id,
            node_depth=evt_data.depth,
            node_obj=evt_data.node_obj,
            node_lb=evt_data.node_lb,
            gap=evt_data.gap,
            solve_stat=Node.SOLVE_STATUS_REPR_MAP[evt_data.solve_status],
            node_type=Node.NODE_TYPE_REPR_MAP[evt_data.node_type],
            iter=evt_data.iter,
            init_cols=evt_data.init_cols,
            added_cols=evt_data.added_cols,
            init_cuts=evt_data.init_cuts,
            added_cuts=evt_data.added_cuts,
            t_master=evt_data.t_master,
            t_pricing=evt_data.t_pricing,
            t_detect_cut=evt_data.t_detect_cut,
            t_total=evt_data.t_total,
        )
        self.node_info_text.config(state="normal")
        self.node_info_text.delete("1.0", "end")
        self.node_info_text.insert(tk.END, new_text)
        self.node_info_text.config(state="disabled")
        self.node_log_text.config(state="normal")
        cur_log_text_yview = self.node_log_text.yview()[1]
        scroll_to_bottom = False
        if cur_log_text_yview > 0.99:
            scroll_to_bottom = True
        if is_node_chg:
            self.node_log_text.delete("1.0", "end")
            log = "".join(evt_data.log_history)
            if log:
                self.node_log_text.insert(tk.END, log)
        else:
            line_count = int(self.node_log_text.index("end-1c").split(".")[0]) - 1
            log_size = (
                evt_data.iter_log_size_map[evt_data.iter]
                if evt_data.iter < len(evt_data.iter_log_size_map)
                else evt_data.iter_log_size_map[-1]
            )
            log = "".join(evt_data.log_history[line_count:log_size])
            if log:
                self.node_log_text.insert(tk.END, log)
        if scroll_to_bottom:
            self.node_log_text.yview(tk.MOVETO, 1)
        self.node_log_text.config(state="disabled")

    def update(self):
        if self.evt_queue:
            update_evt = None
            cur_evt = self.evt_queue[self.cur_evt_idx]
            if cur_evt.evt_id > self.cur_evt_id:
                update_evt = cur_evt
                self.cur_evt_id = cur_evt.evt_id
            else:
                if self.cur_evt_idx < len(self.evt_queue) - 1:
                    self.cur_evt_idx += 1
                    self.cur_evt_id = self.evt_queue[self.cur_evt_idx].evt_id
                    update_evt = self.evt_queue[self.cur_evt_idx]
            if update_evt is not None:
                evt_code = update_evt.evt_code
                evt_data = update_evt.evt_data
                if evt_code == BranchandPrice.EVT_NEW_NODE:
                    self.canvas.create_node(evt_data)
                elif evt_code == BranchandPrice.EVT_NODE_UPDATE:
                    self.canvas.update_node(evt_data)
                elif evt_code == BranchandPrice.EVT_BAP_NEW_LOG:
                    self.add_text(self.BaP_log_text, evt_data)
        self.root.after(100, self.update)

    def run(self):
        t = threading.Thread(target=self.BaP.search)
        t.start()
        self.root.after(1000, self.update)
        self.root.mainloop()


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
    evt_queue = []
    node_compare = Node.NODE_COMPARE_BOUND
    # node_compare = Node.NODE_COMPARE_DEPTH
    BaP = BranchandPrice(N, graph, evt_queue, node_compare=node_compare)
    vis = BranchandPriceVisual(BaP, evt_queue)
    try:
        vis.run()
        os._exit(0)
    except Exception as e:
        print(e)
        os._exit(1)
