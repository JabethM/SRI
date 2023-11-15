import copy

import numpy as np
import random
import networkx as nx
from Variant import Variant


class SIR:
    NEW_DISEASE_CHANCE = 0.01

    def __init__(self, nodes, initialisation=(0, 0.1), variants=2, probabilities=None, precision=1000):

        if probabilities is None:
            probabilities = ((0.5, 0.5), (0.5, 0.5))
        assert (variants == len(probabilities))

        self.num_of_variants = variants
        # Increments everytime a new variant is added to the system
        self.variant_count = 0
        self.variants = []
        self.relation_matrix = None
        self.root = Variant(random.randint(1, 1000))

        self.single_variant = True

        self.generate_variants()

        self.dt = 0.1
        self.time = 0
        self.end = 1000
        self.probability_precision = precision

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
        large_picking_set = 30 * self.num_of_variants
        if self.single_variant > 1:
            repeated = True  # Prevents the same number being chosen
            variant_data = 0
            while repeated is True:
                variant_data = random.randint(1, large_picking_set)
                if variant_data not in Variant.current_data_set:
                    Variant.current_data_set.append(variant_data)
                    repeated = False

            name = ord('A') + self.variant_count
            disease = self.root.insert(variant_data, name=name)
            self.root.add_to_relation_matrix()

            self.variants.append(disease)
            self.variant_count += 1

            return self.variants
        else:
            self.single_variant = True
            return None

    def set_disease(self):
        number_of_initial_outbreaks = 1
        if self.variant_count >= self.num_of_variants:
            return False

        potential_set = range(self.num_nodes)
        patient_zero = np.random.choice(potential_set, number_of_initial_outbreaks)

        for case in patient_zero:
            self.node_change(self.variant_count, np.array([case]), np.array([True]), True, True)

            self.G = self.double_infection_check(self.G, case)
        self.variant_count += 1
        return True

    def choose_outcome(self):
        nbs = [[] for _ in range(self.variant_count)]
        outcome_probabilities = self.P[:self.variant_count]
        for i in range(self.variant_count):
            inf_nbs = self.node_neighbours(self.infected_set[i])
            nbs[i] = inf_nbs

        # Randomly Choose one of 2 x N outcomes where N represents the different variants
        outcome_distribution = outcome_probabilities.flatten() / self.variant_count
        chosen_action = np.random.choice(np.arange(len(outcome_distribution)), size=1, p=outcome_distribution)
        action_row, action_column = divmod(chosen_action[0], 2)  # 2 is the number of different outcomes for infectious

        disease_variant = action_row
        choosing_set = nbs[disease_variant]

        # length of choosing set might change but max distinct choices will remain the same
        distinct_num_node_choices = len(choosing_set)

        infect_others = None
        if action_column == 0:
            # Infect Random Neighbour
            new_infectious_prob = self.NEW_DISEASE_CHANCE / distinct_num_node_choices
            old_infectious_prob = (1 - self.NEW_DISEASE_CHANCE) / distinct_num_node_choices

            new_infection_dist = np.full(distinct_num_node_choices, new_infectious_prob)
            old_infection_dist = np.full(distinct_num_node_choices, old_infectious_prob)
            infection_dist = np.concatenate((old_infection_dist, new_infection_dist))

            test = np.sum(infection_dist)
            choosing_set = choosing_set + choosing_set
            infect_others = True

        else:
            # Recover node
            infection_dist = np.ones_like(choosing_set) / distinct_num_node_choices
            infection_dist = np.array(infection_dist).astype(float)
            infect_others = False

        assert (abs(np.sum(infection_dist) - 1) < 1e-9)

        infected_node_idx = np.random.choice(range(len(choosing_set)), size=1, p=infection_dist)[0]
        infected_node_idx = infected_node_idx % distinct_num_node_choices

        infected_node = list(self.infected_set[disease_variant])[infected_node_idx]
        victim_node = np.random.choice(nbs[disease_variant][infected_node_idx])

        if infect_others:
            return (infected_node, victim_node), disease_variant, infect_others
        else:
            return infected_node, disease_variant, infect_others

    def node_neighbours(self, list_of_nodes):
        nbs = [[n for n in self.G.neighbors(node)] for node in list_of_nodes]
        return nbs

    def update_state_sets(self):

        return

    def node_change(self, variant, node, transition_success: np.ndarray, infected=False, update_system=True,
                    infected_set=None, recovered_set=None, graph=None):

        assert (update_system or ((not update_system) and (infected_set is not None) and (recovered_set is not None)))

        trans_nodes = np.where(transition_success == True)
        trans_nodes = node[trans_nodes]
        trans_nodes_set = set(trans_nodes)

        if update_system:
            i_set = copy.deepcopy(self.infected_set)
            r_set = copy.deepcopy(self.recovered_set)
        else:
            i_set = infected_set
            r_set = recovered_set

        if graph is None:
            graph = self.G

        if len(trans_nodes) == 0:
            return i_set, r_set, graph

        sets = [i_set, r_set]
        s_size = len(sets)

        try:
            
            sets[int(infected)][variant] -= trans_nodes_set
        except KeyError:
            pass

        # Zero for infected, One for Recovered
        infected_or_recovered = (int(infected) + 1) % len(sets)
        sets[infected_or_recovered][variant].update(trans_nodes_set)

        for n in trans_nodes:
            graph.nodes(data=True)[n]["state"][variant] = infected_or_recovered + 1

        if update_system:
            self.infected_set = i_set
            self.recovered_set = r_set
        return i_set, r_set, graph

    def iterate(self):
        G_copy = self.G.copy()
        nodes, disease, infection_status = self.choose_outcome()

        if infection_status:
            a = 2
            infector = nodes[0]
            victim = nodes[1]
            self.infected_set
            # Infect victim in nodes
            # update self.infected with victim at diseased row
            # remove victim from other rows
        else:
            infector = nodes
            # Recover node in nodes
            # update self.recovered with infector at diseased row
            # remove infector from other rows

        """G_copy = self.G.copy()

        transition_probabilities = np.random.random((self.num_nodes, self.num_of_variants + 1))
        current_states = np.array([attribute.get("state") for node, attribute in G_copy.nodes(data=True)])
        neighbours = [[n for n in G_copy.neighbors(i)] for i in range(self.num_nodes)]

        potentially_infected_set = copy.deepcopy(self.infected_set)
        potentially_recovered_set = copy.deepcopy(self.recovered_set)

        for variant in range(self.num_of_variants):

            is_infected = np.where(transition_probabilities < self.P[variant][0], True, False)
            is_recovered = np.where(transition_probabilities < self.P[variant][1], True, False)
            for i_node in self.infected_set[variant]:
                nbs = list(neighbours[i_node])
                nbs.append(i_node)
                nbs = np.array(nbs)

                nbs_state = np.array(current_states[nbs])
                nbs_var_state = nbs_state[:, variant]

                should_infect = nbs_var_state == 0
                should_infect &= is_infected[nbs][:, variant]

                should_recover = nbs_var_state == 1
                should_recover &= is_recovered[nbs][:, variant]

                immune_evolution_chance = np.random.random(self.num_of_variants)
                immune_evolution_chance[variant] = 0

                has_immunity_evolved = np.where(immune_evolution_chance < self.relation_matrix[variant])[0]

                potentially_infected_set, potentially_recovered_set, G_copy \
                    = self.node_change(variant, nbs, should_infect, True, False, potentially_infected_set,
                                       potentially_recovered_set, G_copy)

                potentially_infected_set, potentially_recovered_set, G_copy \
                    = self.node_change(has_immunity_evolved, nbs, should_recover, False, False,
                                       potentially_infected_set, potentially_recovered_set, G_copy)"""
        self.G = G_copy
        self.infected_set = potentially_infected_set
        self.recovered_set = potentially_recovered_set
        return

    def double_infection_check(self, graph, node):
        node_state = graph.nodes(data=True)[node]["state"]
        infections = np.where(node_state == 1)
        if np.shape(infections)[-1] > 1:
            predominant_infection = np.random.choice(infections)
            node_state[infections] = 0
            node_state[predominant_infection] = 1
        return graph

    def run(self, time_between_diseases=1):
        while self.time < self.end:
            self.iterate()

            self.time = round(self.time + self.dt, 2)
            if self.time % time_between_diseases == 0:
                next_stage = self.set_disease()
                if next_stage is False:
                    self.time = self.end


def main():
    execute = SIR(50, initialisation=(0, 0.1), variants=7, probabilities=tuple([(random.random() / 2, random.random())
                                                                                for _ in range(7)]))
    execute.run()


main()
