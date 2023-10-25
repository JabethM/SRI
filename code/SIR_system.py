import copy

import numpy as np
import random
import networkx as nx
from Variant import Variant


class SIR:

    def __init__(self, nodes, initialisation=(0, 0.1), variants=2, probabilities=None):

        if probabilities is None:
            probabilities = ((0.5, 0.5), (0.5, 0.5))
        assert (variants == len(probabilities))

        self.num_of_variants = variants
        # Increments everytime a new variant is added to the system
        self.variant_count = 0
        self.variants = []
        self.relation_matrix = None

        self.single_variant = True

        self.generate_variants()

        self.dt = 0.1
        self.time = 0

        self.P = np.asarray(probabilities)

        self.num_nodes = nodes
        self.G = None

        """
        There's two approaches I could take, one which iterates through the entire graph and computes the next state
        based on neighbours but if the graph doesnt percolate and has large number of nodes, this could become wildly 
        inefficient.
        """
        self.susceptible_set = [set() for v in range(self.num_of_variants)]
        self.infected_set = [set() for v in range(self.num_of_variants)]
        self.recovered_set = [set() for v in range(self.num_of_variants)]

        self.initialise_graph(initialisation)
        self.set_disease()

    def initialise_graph(self, init_tuple):
        # Random
        if init_tuple[0] == 0:
            self.G = nx.erdos_renyi_graph(self.num_nodes, init_tuple[1])
        elif init_tuple[0] == 1:
            self.G = nx.watts_strogatz_graph(self.num_nodes, init_tuple[1], init_tuple[2])
        # Scale Free
        elif init_tuple[0] == 2:
            self.G = nx.barabasi_albert_graph(self.num_nodes, init_tuple[1])
        # High Modularity
        elif init_tuple[0] == 3:
            assert (sum(init_tuple[1]) == self.num_nodes)
            assert (len(init_tuple[2]) == len(init_tuple[1]))
            self.G = nx.stochastic_block_model(init_tuple[1], init_tuple[2])

        state_dict = {i: np.zeros(self.num_of_variants) for i in range(self.num_nodes)}
        nx.set_node_attributes(self.G, state_dict, name="state")
        return self.G

    def generate_variants(self):
        # arbitrary high number to pick from
        large_picking_set = 1000
        if self.num_of_variants > 1:
            self.single_variant = False

            # Data that is used to sort binary tree
            variant_data = np.random.randint(1, self.num_of_variants * large_picking_set, self.num_of_variants)
            range_of_variant_data = max(variant_data) - min(variant_data)

            # e^(-1 *(Delta(data) / range(data)) * (number of variants between))

            variant_relation = [[np.exp(-1 * (abs(variant_data[j] - variant_data[i]) / range_of_variant_data)
                                        * abs(j - i)) if j != i else 1
                                 for j in range(len(variant_data))]
                                for i in range(len(variant_data))]
            Variant.relation_matrix = np.array(variant_relation)
            self.relation_matrix = Variant.relation_matrix

            name = ord('A')
            root = Variant(variant_data[0], name=chr(name))
            for vr in range(len(variant_data)):
                if vr == 0:
                    disease = root
                else:
                    disease = root.insert(variant_data[vr], name=chr(name + vr))

                self.variants.append(disease)
            return self.variants
        else:
            self.single_variant = True
            return None

    def potential_patient_zero_set(self):
        if self.variant_count == 0:
            potentially_infectious = range(self.num_nodes)
        else:
            state_data = self.G.nodes(data=True)[1].get("state")
            potentially_infectious = np.where(state_data[self.variant_count - 1] >= 1)[0]

        return potentially_infectious

    def set_disease(self):
        number_of_initial_outbreaks = 1
        if self.variant_count >= self.num_of_variants:
            return False

        potential_set = self.potential_patient_zero_set()
        patient_zero = np.random.choice(potential_set, number_of_initial_outbreaks)

        for case in patient_zero:
            self.G.nodes[case]["state"][self.variant_count] = 1
            self.node_change([self.variant_count], case, True, True)

        self.variant_count += 1
        return True

    def node_change(self, variant: [int], node, is_infected=False, update_system=True,
                    infected_set=None, recovered_set=None):

        assert (update_system or ((not update_system) and (infected_set is not None) and (recovered_set is not None)))
        assert(isinstance(variant, np.ndarray) or isinstance(variant, list))

        s_set = copy.deepcopy(self.susceptible_set)
        if update_system:
            i_set = copy.deepcopy(self.infected_set)
            r_set = copy.deepcopy(self.recovered_set)
        else:
            i_set = infected_set
            r_set = recovered_set

        sets = [i_set, r_set]
        s_size = len(sets)
        for v in variant:
            try:
                s_set[v].remove(node)
            except KeyError:
                pass

            try:
                sets[int(is_infected)][v].remove(node)
            except KeyError:
                pass

            sets[(int(is_infected) + 1) % len(sets)][v].add(node)

        if update_system:
            self.susceptible_set = s_set
            self.infected_set = i_set
            self.recovered_set = r_set
        return s_set, i_set, r_set

    def iterate(self):
        G_copy = self.G.copy()
        transition_probabilities = np.random.random((self.num_nodes, self.num_of_variants + 1))
        current_states = np.array([attribute.get("state") for node, attribute in G_copy.nodes(data=True)])
        neighbours = [[n for n in G_copy.neighbors(i)] for i in range(self.num_nodes)]

        potentially_infected_set = copy.deepcopy(self.infected_set)
        potentially_recovered_set = copy.deepcopy(self.recovered_set)

        for variant in range(self.num_of_variants):
            for i_node in self.infected_set[variant]:
                nbs = neighbours[i_node]
                nbs.append(i_node)

                nbs_state = current_states[nbs]
                nbs_var_state = nbs_state[:, variant]

                a = transition_probabilities[nbs]  # [:, variant]
                is_infected = np.where(a < self.P[variant][0], True, False)
                is_recovered = np.where(a < self.P[variant][1], True, False)
                for n in range(len(nbs)):

                    if nbs_var_state[n] == 0 and is_infected[n][variant]:
                        placeholder, potentially_infected_set, potentially_recovered_set \
                            = self.node_change([variant], nbs[n], True, False, potentially_infected_set,
                                               potentially_recovered_set)

                    if nbs_var_state[n] == 1 and is_recovered[n][variant]:
                        immune_evolution_chance = np.random.random(self.num_of_variants)
                        immune_evolution_chance[variant] = 0

                        variant_digits = np.arange(self.num_of_variants)
                        has_immunity_evolved = np.where(immune_evolution_chance < self.relation_matrix[variant],
                                                        variant_digits, -1)
                        HIE = has_immunity_evolved[has_immunity_evolved != -1]
                        placeholder, potentially_infected_set, potentially_recovered_set \
                            = self.node_change(HIE, nbs[n], False, False, potentially_infected_set,
                                               potentially_recovered_set)

                    else:
                        pass

        # TODO: check if a node is infected with nore than one disease
        self.infected_set = potentially_infected_set
        self.recovered_set = potentially_recovered_set
        return


def main():
    execute = SIR(50, initialisation=(0, 0.1), variants=7, probabilities=tuple([(random.random(), random.random())
                                                                                for _ in range(7)]))

    for i in range(10):
        execute.iterate()


main()
