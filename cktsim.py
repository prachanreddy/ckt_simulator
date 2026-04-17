import networkx as nx
import re
import numpy as np
import pandas as pd
import sympy as sp


class Component:
    def __init__(self, name, n1, n2, value_str):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = self._parse_value(value_str)

    def _parse_value(self, value_str):
        unit_prefix = {"k": 1e3, "M": 1e6, "m": 1e-3, "u": 1e-6, "n": 1e-9}
        match = re.match(r"([-+]?\d*\.?\d+)([a-zA-Z]*)", value_str)
        if not match:
            raise ValueError(f"Invalid value for '{self.name}': '{value_str}'")

        mag, suffix = match.groups()
        scale = unit_prefix.get(suffix, 1.0)
        return sp.sympify(mag) * sp.sympify(scale)

    def get_impedance(self, s):
        raise NotImplementedError

    def get_voltage_source(self, s):
        return sp.Integer(0)


class Resistor(Component):
    def get_impedance(self, s):
        return self.value


class Capacitor(Component):
    def get_impedance(self, s):
        return 1 / (s * self.value)


class Inductor(Component):
    def get_impedance(self, s):
        return s * self.value


class VoltageSource(Component):
    def get_impedance(self, s):
        return sp.Integer(0)

    def get_voltage_source(self, s):
        return self.value / s


class BB(Component):
    def get_impedance(self, s):
        return s**2


class ACVoltageSource(Component):
    def get_impedance(self, s):
        return sp.Integer(0)

    def get_voltage_source(self, s):
        omega = 2 * sp.pi * 60  # 60 Hz
        return self.value * s / (s**2 + omega**2)


COMPONENT_CLASSES = {
    'R': Resistor,
    'C': Capacitor,
    'L': Inductor,
    'V': VoltageSource,
    'B': BB
}


def component_factory(name, n1, n2, value_str):
    name_upper = name.upper()

    if name_upper.startswith("AC"):
        return ACVoltageSource(name, n1, n2, value_str)

    key = name_upper[0]
    if key in COMPONENT_CLASSES:
        return COMPONENT_CLASSES[key](name, n1, n2, value_str)

    raise ValueError(f"Unknown component type for '{name}'")


class Circuit:
    def __init__(self, netlist_string):
        self.netlist = netlist_string
        self.graph = nx.MultiDiGraph()
        self.components = {}

        self._parse_netlist()
        self.twigs, self.links = self._get_spanning_tree_and_branches()

        self.s = sp.Symbol('s', complex=True)
        self.currents = None
        self.voltages = None

    def _parse_netlist(self):
        for line in self.netlist.strip().split('\n'):
            parts = line.split()
            if not parts or parts[0].startswith('#'):
                continue
            self._parse_component(parts)

    def _parse_component(self, parts):
        name, n1, n2 = parts[0], parts[1], parts[2]

        if len(parts) < 4:
            raise ValueError(f"Component '{name}' is missing a value.")

        self.graph.add_node(n1)
        self.graph.add_node(n2)
        self.graph.add_edge(n1, n2, key=name, name=name)

        self.components[name] = component_factory(name, n1, n2, parts[3])

    def _get_spanning_tree_and_branches(self):
        if len(self.graph.nodes) == 0 or not nx.is_connected(self.graph.to_undirected()):
            raise ValueError("Graph not connected.")

        simple_topo = nx.Graph()
        simple_topo.add_nodes_from(self.graph.nodes())

        for u, v in self.graph.edges():
            simple_topo.add_edge(u, v)

        root = sorted(list(self.graph.nodes()))[0]
        tree_edges = {tuple(sorted(e)) for e in nx.bfs_edges(simple_topo, source=root)}

        all_branches = list(self.graph.edges(keys=True))
        twigs, links = [], []

        for branch in all_branches:
            u, v, k = branch
            edge = tuple(sorted((u, v)))

            if edge in tree_edges:
                twigs.append(branch)
                tree_edges.remove(edge)
            else:
                links.append(branch)

        return twigs, links

    def analyze(self):
        cutset = self._generate_cutset_matrix()
        tieset = self._generate_tieset_matrix()

        branch_names = [b[2] for b in self.links + self.twigs]

        Z = self._build_impedance_matrix(branch_names)
        Vs = self._build_source_vector(branch_names)

        B = sp.Matrix(tieset.values.tolist()) if not tieset.empty else sp.zeros(0, len(branch_names))
        A = sp.Matrix(cutset.values.tolist())

        I_mat = sp.eye(Z.shape[0])

        top = sp.Matrix.hstack(A, sp.zeros(A.shape[0], Z.shape[0]))
        mid = sp.Matrix.hstack(Z, -I_mat)
        bot = sp.Matrix.hstack(sp.zeros(B.shape[0], Z.shape[0]), B)

        M = sp.Matrix.vstack(top, mid, bot)
        rhs = sp.Matrix.vstack(
            sp.zeros(A.shape[0], 1),
            -Vs,
            sp.zeros(B.shape[0], 1)
        )

        try:
            sol = M.gauss_jordan_solve(rhs)[0]
        except:
            raise np.linalg.LinAlgError("System matrix is singular.")

        n = len(branch_names)
        self.currents = {}
        self.voltages = {}

        for i, b in enumerate(branch_names):
            self.currents[f"I_{b}"] = sp.simplify(sol[i])
            self.voltages[f"V_{b}"] = sp.simplify(sol[n + i])

        print("--- Circuit Structure ---")
        print(f"Twigs: {[t[2] for t in self.twigs]}")
        print(f"Links: {[l[2] for l in self.links]}")

        print("\n--- Laplace Domain Results ---")
        for k, v in self.currents.items():
            print(f"{k}(s) = {v}")
        for k, v in self.voltages.items():
            print(f"{k}(s) = {v}")

    def _build_impedance_matrix(self, branch_names):
        Z = sp.zeros(len(branch_names))
        for i, b in enumerate(branch_names):
            Z[i, i] = self.components[b].get_impedance(self.s)
        return Z

    def _build_source_vector(self, branch_names):
        Vs = sp.zeros(len(branch_names), 1)
        for i, b in enumerate(branch_names):
            Vs[i] = self.components[b].get_voltage_source(self.s)
        return Vs

    def _generate_cutset_matrix(self):
        all_branches = self.links + self.twigs
        branch_names = [b[2] for b in all_branches]
        twig_names = [t[2] for t in self.twigs]

        cutset = pd.DataFrame(0, index=twig_names, columns=branch_names)

        tree = nx.MultiDiGraph()
        tree.add_nodes_from(self.graph.nodes())
        tree.add_edges_from(self.twigs)

        for twig in self.twigs:
            u, v, key = twig
            tree.remove_edge(u, v, key=key)

            part = set(nx.node_connected_component(tree.to_undirected(), u))

            for b in all_branches:
                bu, bv, bk = b
                if bu in part and bv not in part:
                    cutset.loc[key, bk] = 1
                elif bv in part and bu not in part:
                    cutset.loc[key, bk] = -1

            tree.add_edge(u, v, key=key)

        return cutset

    def _generate_tieset_matrix(self):
        all_branches = self.links + self.twigs
        link_names = [l[2] for l in self.links]
        branch_names = [b[2] for b in all_branches]

        tieset = pd.DataFrame(0, index=link_names, columns=branch_names)

        if not self.links:
            return tieset

        for link in self.links:
            tieset.loc[link[2], link[2]] = 1

        tree = nx.Graph()
        for u, v, k in self.twigs:
            tree.add_edge(u, v, key=k)

        for link in self.links:
            u, v, key = link
            path = nx.shortest_path(tree, source=v, target=u)

            for i in range(len(path) - 1):
                pu, pv = path[i], path[i + 1]

                for twig in self.twigs:
                    tu, tv, tk = twig

                    if (tu, tv) == (pu, pv):
                        tieset.loc[key, tk] = 1
                        break
                    elif (tu, tv) == (pv, pu):
                        tieset.loc[key, tk] = -1
                        break

        return tieset


if __name__ == "__main__":
    netlist = """
        AC1 1 0 10
        R1 1 2 10
        L1 2 3 25m
        C1 3 0 100u
    """

    circuit = Circuit(netlist)
    circuit.analyze()
