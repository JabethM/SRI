import copy

import numpy as np
import random
import networkx as nx
from Variant import Variant
from itertools import chain


class SIR:
    NEW_DISEASE_CHANCE = 0.01
    probability_precision = 7

    def __init__(self, nodes, initialisation=(0, 0.1), variants=2, probabilities=None, end_time=np.inf, seed=None):

        if probabilities is None:
            probabilities = ((0.5, 0.5), (0.5, 0.5))
        assert (variants == len(probabilities))

        self.num_of_variants = variants
        # Increments everytime a new variant is added to the system
        self.variant_count = 0
        self.variants = []
        self.relation_matrix = None
        self.root = None

        self.single_variant = True

        self.dt = 0.1
        self.time = 0
        self.end = False
        self.end_time = end_time
        np.random.seed(seed)

        self.P = np.asarray(probabilities)

        self.num_nodes = nodes
        self.G = None

        """
        There's two approaches I could take, one which iterates through the entire graph and computes the next state
        based on neighbours but if the graph doesnt percolate and has large number of nodes, this could become wildly 
        inefficient.
        """
        self.infected_set = [set() for v in range(self.num_of_variants)]
        self.recovered_set = [set() for v in range(self.num_of_variants)]
        self.possible_outcomes = []

        self.initialise_graph(initialisation)
        self.set_initial_outbreak()

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

    def generate_new_disease(self, parent=None):
        if self.variant_count == 0:
            variant_data = random.randint(1, 1000)
            self.root = Variant(variant_data)
            Variant.current_data_set.append(variant_data)

        else:
            large_picking_set = 30 * self.num_of_variants

            repeated = True  # Prevents the same number being chosen
            variant_data = 0
            while repeated is True:
                variant_data = random.randint(1, large_picking_set)
                if variant_data not in Variant.current_data_set:
                    Variant.current_data_set.append(variant_data)
                    repeated = False

        if parent is None:
            parent = self.root

        name = ord('A') + self.variant_count
        disease = parent.insert(variant_data, name=chr(name))
        Variant.add_to_relation_matrix()

        self.variants.append(disease)
        self.variant_count += 1
        if self.variant_count > self.num_of_variants:
            self.end = True

        return self.variant_count

    def set_initial_outbreak(self):
        number_of_initial_outbreaks = 1
        if self.variant_count >= self.num_of_variants:
            return False

        potential_set = range(self.num_nodes)
        patient_zero = np.random.choice(potential_set, number_of_initial_outbreaks)

        for case in patient_zero:
            self.infected_set[self.variant_count].update({case})
            self.G.nodes(data=True)[case]["state"][self.variant_count] = 1

            self.G = self.double_infection_check(self.G, case)

        self.generate_new_disease()

        return True

    def choose_outcome(self):
        temp_infected_set = copy.deepcopy(self.infected_set)[:self.variant_count]

        nbs, gillespe_line, possible_outcomes = self.choice_probability_setup(temp_infected_set)

        if self.is_list_empty(gillespe_line):
            self.end = True
            return None

        gillespe_hold = gillespe_line  # temporary array which will be used to examine parts of the line
        choices = []
        for i in range(3):
            if i == 0:
                shortened_gillespe = [np.sum(list(chain(*variant))) for variant in gillespe_hold]
            else:
                shortened_gillespe = [np.sum(choice) for choice in gillespe_hold]

            shortened_gillespe = np.array(shortened_gillespe) / np.sum(shortened_gillespe)
            chosen_option = np.random.choice(np.arange(np.shape(gillespe_hold)[0]), size=1, p=shortened_gillespe)[0]
            choices.append(chosen_option)
            gillespe_hold = gillespe_hold[chosen_option]

        if choices[1] == 0:
            cumulative_lengths = np.cumsum(list(map(len, nbs[choices[0]])))
            infector_idx = np.searchsorted(cumulative_lengths, choices[2], side='right')
            victim_idx = choices[2] - (cumulative_lengths[infector_idx - 1] if infector_idx > 0 else 0)
            victim = nbs[choices[0]][infector_idx][victim_idx]

        else:
            victim = list(temp_infected_set[choices[0]])[choices[2]]

        # True when infecting, False when recovering
        infect_others = not bool(choices[1])

        return victim, choices[0], infect_others

    def choice_probability_setup(self, temp_infected_set):
        initial_buffer = 10

        nbs = [[list(np.setdiff1d(nbs, list(variant))) for nbs in self.node_neighbours(variant)]
               for variant in temp_infected_set]

        # Randomly Choose one of 2 x N outcomes where N represents the different variants
        possible_outcomes = []
        gillespe_line = []

        for variant in range(self.variant_count):
            inf_nbs_flat = list(np.array(nbs[variant]).flatten())
            inf_nbs_flat_test = [nb for inf_node in nbs[variant] for nb in inf_node]

            gillespe_line.append([list(np.full(len(inf_nbs_flat_test), initial_buffer * self.P[variant, 0])),
                                  list(np.full(len(temp_infected_set[variant]), initial_buffer * self.P[variant, 1]))])

            possible_outcomes.append([inf_nbs_flat, list(temp_infected_set[variant])])

        return nbs, gillespe_line, possible_outcomes

    def set_node_infected(self, node, variant, graph):
        # CHANCE OF NEW INFECTION
        new_disease = np.random.choice([True, False], size=1,
                                       p=[self.NEW_DISEASE_CHANCE, 1 - self.NEW_DISEASE_CHANCE])
        if new_disease:
            variant = self.generate_new_disease(self.variants[variant]) - 1

        if self.end or node in self.recovered_set[variant]:
            return graph

        self.infected_set[variant].update({node})
        graph.nodes(data=True)[node]["state"][variant] = 1
        return graph

    def set_node_recovery(self, node, variant, graph):

        relation_array = Variant.get_relation()[variant]
        assert relation_array[variant] == 1
        probabilities = np.column_stack((relation_array, 1 - relation_array))
        immunity_choices = np.array([np.random.choice([True, False], size=1, p=probs)[0]
                                     for probs in probabilities])
        immunity_developed = np.where(immunity_choices)[0]

        for i in immunity_developed:
            self.infected_set[i].discard(node)
            self.recovered_set[i].update({node})
            graph.nodes(data=True)[node]["state"][i] = 2
        return graph

    # ###### HELPERS ##### #
    def double_infection_check(self, graph, node):
        node_state = graph.nodes(data=True)[node]["state"]
        infections = np.where(node_state == 1)
        if np.shape(infections)[-1] > 1:
            predominant_infection = np.random.choice(infections)
            node_state[infections] = 0
            node_state[predominant_infection] = 1
        return graph

    def is_list_empty(self, in_list):
        if isinstance(in_list, list):  # Is a list
            return all(map(self.is_list_empty, in_list))
        return False

    def node_neighbours(self, list_of_nodes):
        nbs = [[n for n in self.G.neighbors(node)] for node in list_of_nodes]
        return nbs

    # ###### HELPERS ##### #

    def iterate(self):
        G_copy = self.G.copy()
        outcomes = self.choose_outcome()
        if outcomes is not None:
            node, disease, infection_status = outcomes
            if infection_status:
                G_copy = self.set_node_infected(node, disease, G_copy)
            else:
                G_copy = self.set_node_recovery(node, disease, G_copy)

            self.G = G_copy

        return

    def run(self):
        while not self.end:
            self.iterate()
            self.time = round(self.time + self.dt, 2)
            print(self.time)
            if self.time >= self.end_time:
                self.end = True
        print(self.P)
        print(self.variant_count)

    def step_run(self, dt=0.1):
        print(self.P)
        print(self.variant_count)

        if not self.end:
            self.iterate()
            self.time = round(self.time + dt, 2)
            print(self.time)

            if self.time >= self.end_time:
                self.end = True
        return self.end

def main():
    seed = 345
    np.random.seed(seed)
    probs = tuple([(np.random.randint(0, 1000), np.random.randint(0, 600)) for _ in range(7)])
    execute = SIR(50, initialisation=(0, 0.1), variants=7, probabilities=probs, seed=seed)
    execute.run()
