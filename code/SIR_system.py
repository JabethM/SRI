import copy

import numpy as np
import random
import networkx as nx
from Variant import Variant


class SIR:
    NEW_DISEASE_CHANCE = 0.01
    probability_precision = 7

    def __init__(self, nodes, initialisation=(0, 0.1), variants=2, probabilities=None):

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

        self.P = np.asarray(probabilities)
        # Normalising probabilities w.r.t themselves
        self.P = np.apply_along_axis(lambda x: x / np.sum(x), axis=1, arr=self.P)

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

    def generate_new_disease(self):
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

        name = ord('A') + self.variant_count
        disease = self.root.insert(variant_data, name=chr(name))
        disease.add_to_relation_matrix()

        self.variants.append(disease)
        self.variant_count += 1
        if self.variant_count > self.num_of_variants:
            self.end = True

        return len(self.variants)

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

        nbs, outcome_distribution, outcome_distribution_flat = self.choice_probability_setup(temp_infected_set)

        cumulative_lengths = np.cumsum(list(map(len, outcome_distribution)))

        chosen_node = np.random.choice(np.arange(len(outcome_distribution_flat)), size=1, p=outcome_distribution_flat)
        chosen_action = np.searchsorted(cumulative_lengths, chosen_node, side='right')[0]

        node_idx = chosen_node[0] - (cumulative_lengths[chosen_action - 1] if chosen_action > 0 else 0)
        disease_chosen = (chosen_action // 2)

        infect_others = chosen_action % 2 == 0
        if infect_others:
            # INFECT NEIGHBOUR
            potential_victims = nbs[disease_chosen][node_idx]
            victim = np.random.choice(potential_victims)
            changing_node = victim
        else:
            # RECOVER
            infected_node = list(temp_infected_set[disease_chosen])[node_idx]
            changing_node = infected_node

        return changing_node, disease_chosen, infect_others

    def choice_probability_setup(self, temp_infected_set):
        nbs = [[list(np.setdiff1d(nbs, list(variant))) for nbs in self.node_neighbours(variant)]
               for variant in temp_infected_set]
        infected_count = np.array([len(nb) for nb in nbs])
        infected_count = np.round(infected_count / (np.sum(infected_count) if np.sum(infected_count) > 0 else 1),
                                  self.probability_precision)

        outcome_probabilities = np.round(self.P[:self.variant_count], self.probability_precision)
        outcome_probabilities *= infected_count[:, np.newaxis]

        # Randomly Choose one of 2 x N outcomes where N represents the different variants
        possible_outcomes = []
        outcome_distribution = []

        for variant in range(self.variant_count):
            inf_nbs_flat = list(np.array(nbs[variant]).flatten())
            possible_outcomes.extend([inf_nbs_flat, list(temp_infected_set[variant])])

            nb_lengths = np.array([len(neighbours) for neighbours in nbs[variant]])
            nb_weights = nb_lengths / (np.sum(nb_lengths) if np.sum(nb_lengths) > 0 else 1)

            infected_nodes = list(temp_infected_set[variant])
            infected_weights = np.ones_like(infected_nodes) / len(infected_nodes)
            outcome_distribution.extend([list(outcome_probabilities[variant][0] * nb_weights),
                                         list(outcome_probabilities[variant][1] * infected_weights)])

        outcome_distribution_flat = np.concatenate(outcome_distribution)

        # Account for deleted probabilities due to redundant repetitions (infecting an infected node)
        outcome_distribution_flat = outcome_distribution_flat / np.sum(outcome_distribution_flat)

        return nbs, outcome_distribution, outcome_distribution_flat

    def node_neighbours(self, list_of_nodes):
        nbs = [[n for n in self.G.neighbors(node)] for node in list_of_nodes]
        return nbs

    def set_node_infected(self, node, variant, graph):
        # CHANCE OF NEW INFECTION
        new_disease = np.random.choice([True, False], size=1,
                                       p=[self.NEW_DISEASE_CHANCE, 1 - self.NEW_DISEASE_CHANCE])
        if new_disease:
            variant = self.generate_new_disease() - 1

        a = self.infected_set[variant]

        self.infected_set[variant].update({node})
        graph.nodes(data=True)[node]["state"][variant] = 1
        return graph

    def set_node_recovery(self, node, variant, graph):

        relation_array = Variant.get_relation()[variant]
        assert relation_array[variant] == 1
        probabilities = np.column_stack((relation_array, 1-relation_array))
        immunity_choices = np.array([np.random.choice([True, False], size=1, p=probs)[0]
                                     for probs in probabilities])
        immunity_developed = np.where(immunity_choices)[0]

        for i in immunity_developed:
            self.infected_set[i].discard({node})
            self.recovered_set[i].update({node})
            graph.nodes(data=True)[node]["state"][i] = 2
        return graph

    def iterate(self):
        G_copy = self.G.copy()
        node, disease, infection_status = self.choose_outcome()

        if infection_status:
            G_copy = self.set_node_infected(node, disease, G_copy)
        else:
            G_copy = self.set_node_recovery(node, disease, G_copy)

        self.G = G_copy
        return

    def double_infection_check(self, graph, node):
        node_state = graph.nodes(data=True)[node]["state"]
        infections = np.where(node_state == 1)
        if np.shape(infections)[-1] > 1:
            predominant_infection = np.random.choice(infections)
            node_state[infections] = 0
            node_state[predominant_infection] = 1
        return graph

    def run(self, end_time=np.inf):
        while not self.end:
            self.iterate()
            self.time = round(self.time + self.dt, 2)
            print(self.time)
            if self.time >= end_time:
                self.end = True



def main():
    probs = tuple([(random.randint(0, 10000), random.randint(0, 1000)) for _ in range(7)])
    execute = SIR(50, initialisation=(0, 0.1), variants=7, probabilities=probs)
    execute.run()


main()
