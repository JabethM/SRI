import copy
import numpy as np
import networkx as nx
from .Variant import Variant
from itertools import chain
import json


class SIR:

    def __init__(self, nodes=150, initialisation=(0, 0.1), variants=2, probabilities=None, end_time=np.inf, seed=None,
                 epsilon=300, config_file=None, delta=50, picking_set=500, rho=((-np.log(0.4)/(100^2)))):

        if config_file is not None:
            # Simulation Initialised with Configuration JSON File
            with open(config_file, 'r') as file:
                config_data = json.load(file)
            # Internal Params
            self.new_disease_chance = config_data["internal_params"]["new_disease_rate"]
            self.deterministic_spawning = \
                config_data["internal_params"]["deterministic_spawning"]["deterministic_mode"]

            self.epsilon = config_data["internal_params"]["new_infection_std_dev"]
            self.delta = config_data["internal_params"]["variant_info"]["data_std_dev"]
            self.picking_set = config_data["internal_params"]["variant_info"]["data_range"]
            self.rho = config_data["internal_params"]["variant_info"]["relationship_dropoff_coefficient"]
            # Input Params
            self.num_of_variants = config_data["input_params"]["maximum_num_of_variants"]
            self.end_time = config_data["input_params"]["end_time"]

            self.num_nodes = config_data["input_params"]["number_of_nodes"]
            infection_rate = config_data["input_params"]["transmission_time"]
            recovery_rate = config_data["input_params"]["recovery_time"]
            probabilities = (infection_rate, recovery_rate)

            network_type = config_data["input_params"]["network_info"]["type_of_network"]
            network_properties = config_data["input_params"]["network_info"]["network_properties"]
            network_properties.insert(0, network_type)
            initialisation = tuple(network_properties)
            seed = config_data["input_params"]["seed"]
        else:
            # Simulation Initialised with Constructor arguments
            if probabilities is None:
                probabilities = (0.5, 0.5)
            self.new_disease_chance = 0.0025
            self.deterministic_spawning = False

            self.epsilon = epsilon
            self.delta = delta
            self.picking_set = picking_set
            self.rho = rho

            self.num_nodes = nodes
            self.num_of_variants = variants
            self.end_time = end_time

        np.random.seed(seed)
        self.dt = 0.1
        self.time = 0
        self.old_time = 0
        self.end = False
        self.steps = 0

        # Increments everytime a new variant is added to the system
        self.variant_count = 0
        self.variants = []
        self.relation_matrix = None
        self.root = None

        self.trans_time = np.zeros((self.num_of_variants, 2))
        self.trans_time[0] = np.asarray(probabilities)

        self.rates = np.zeros((self.num_of_variants, 2))
        self.rates[0] = 1 / self.trans_time[0]

        if self.deterministic_spawning:
            self.deterministic_count = 0
            self.new_disease_spawn_time = \
                config_data["internal_params"]["deterministic_spawning"]["time_separation"]

            self.descendant_transmissions = \
                config_data["internal_params"]["deterministic_spawning"]["child_transmission_T"]

            self.descendant_recovery = \
                config_data["internal_params"]["deterministic_spawning"]["child_recovery_T"]

            self.determined_relationship_matrix = \
                np.array(config_data["internal_params"]["deterministic_spawning"]["relation_matrix"])

            assert (len(self.determined_relationship_matrix) == self.num_of_variants)
        self.G = None
        Variant.data_range = self.picking_set
        Variant.rho = self.rho
        
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

    def generate_new_disease(self, parent=None, exception=None):
        if self.variant_count >= self.num_of_variants:
            self.variant_count += 1
            return self.variant_count

        if self.variant_count == 0:
            variant_data = float(np.random.randint(1, self.picking_set))
            self.root = Variant(variant_data)
            Variant.current_data_set.append(variant_data)
            if parent is None:
                parent = self.root

        else:
            if parent is None:
                parent = self.root

            if not self.deterministic_spawning:
                self.trans_time[self.variant_count] = abs(np.random.normal(self.trans_time[self.variant_count - 1],
                                                                       scale=self.epsilon))
            else:
                self.trans_time[self.variant_count] = [self.descendant_transmissions[self.variant_count - 1],
                                                       self.descendant_recovery[self.variant_count - 1]]

            variant_data = 0
            self.rates[self.variant_count] = 1 / self.trans_time[self.variant_count]

            repeated = True  # Prevents the same number being chosen
            while repeated is True:
                # Ensures disease is considered 'close' to its parent
                variant_data = np.random.normal(parent.data, scale=self.delta)
                if variant_data not in Variant.current_data_set:
                    Variant.current_data_set.append(variant_data)
                    repeated = False

        name = ord('A') + self.variant_count
        disease = parent.insert(variant_data, name=chr(name))
        disease.infectious_rate = self.rates[self.variant_count, 0]
        disease.recovery_rate = self.rates[self.variant_count, 1]

        if self.deterministic_spawning:
            Variant.relation_matrix = self.determined_relationship_matrix[:self.variant_count+1, :self.variant_count+1]
        else:
            Variant.add_to_relation_matrix()

        self.variants.append(disease)
        self.variant_count += 1

        if self.variant_count > 1 and not self.deterministic_spawning:
            if exception is None:
                exception = False
            self.backdate_immunity(exception)

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
        temp_infected_set = self.infected_set[:self.variant_count]
        temp_rec_set = self.recovered_set[:self.variant_count]

        nbs, gillespe_line, possible_outcomes, rate_total = self.choice_probability_setup(temp_infected_set,
                                                                                          temp_rec_set)

        # No more actions? end
        numpy_gillespe = np.array(gillespe_line)
        if np.all(numpy_gillespe[..., 0] == 0):
            self.end = True
            return None
        self.dt = (1 / rate_total) * np.log(1 / (np.random.rand()))
        gillespe_hold = gillespe_line  # temporary array which will be used to examine parts of the line
        choices = []
        for i in range(3):
            if i == 0:
                shortened_gillespe = np.prod(gillespe_hold, axis=-1).sum(axis=-1)

            elif i == 1:
                shortened_gillespe = np.prod(gillespe_hold, axis=-1)

            elif i == 2:
                shortened_gillespe = np.full(gillespe_hold[0], gillespe_hold[1])
                gillespe_hold = shortened_gillespe

            shortened_gillespe = np.array(shortened_gillespe) / np.sum(shortened_gillespe)
            chosen_option = np.random.choice(np.arange(len(gillespe_hold)), size=1, p=shortened_gillespe)[0]
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

    def choice_probability_setup(self, temp_infected_set, temp_rec_set):
        iterator = range(len(temp_rec_set))

        nbs = [[node for node in
                self.node_neighbours(temp_infected_set[variant], temp_rec_set[variant])]
               for variant in iterator]

        # Randomly Choose one of 2 x N outcomes where N represents the different variants
        possible_outcomes = []
        gillespe_line = []
        rate_total = 0
        for variant in range(self.variant_count):
            inf_nbs_flat = list(chain(*nbs[variant]))

            rate_total += len(inf_nbs_flat) * self.rates[variant, 0] + \
                          len(temp_infected_set[variant]) + self.rates[variant, 1]

            gillespe_line.append([[len(inf_nbs_flat), self.rates[variant, 0]],
                                  [len(temp_infected_set[variant]), self.rates[variant, 1]]])
            possible_outcomes.append([inf_nbs_flat, list(temp_infected_set[variant])])

        return nbs, gillespe_line, possible_outcomes, rate_total

    def set_node_infected(self, node, variant, graph):
        # CHANCE OF NEW INFECTION
        if not self.deterministic_spawning:
            new_disease = np.random.choice([True, False], size=1,
                                           p=[self.new_disease_chance, 1 - self.new_disease_chance])
        else:
            new_disease = not((self.time + self.dt - self.deterministic_count) < self.new_disease_spawn_time)
            self.deterministic_count += int(new_disease) * self.new_disease_spawn_time

        if new_disease and self.variant_count < self.num_of_variants:
            variant = self.generate_new_disease(self.variants[variant], exception=node) - 1

        if self.end or node in self.recovered_set[variant]:
            return graph

        self.infected_set[variant].update({node})
        graph.nodes(data=True)[node]["state"][variant] = 1
        return graph

    def set_node_recovery(self, node, variant, graph):
        if self.deterministic_spawning:
            relation_array = self.determined_relationship_matrix[variant]
        else:
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

    def backdate_immunity(self, exception):
        relation_array = Variant.get_relation()[self.variant_count - 1]
        for variant in range(self.variant_count - 1):
            probability = [relation_array[variant], 1 - relation_array[variant]]
            nodes = np.array(list(self.recovered_set[variant]))

            immunity_developed = np.random.choice([True, False], size=len(nodes), p=probability)
            immune_nodes = nodes[np.where(immunity_developed)[0]]
            immune_nodes = immune_nodes[immune_nodes != exception]
            for node in immune_nodes:
                self.G.nodes(data=True)[node]["state"][self.variant_count - 1] = 2

            self.recovered_set[self.variant_count - 1].update(set(immune_nodes))

        return self.G

    # ###### START OF HELPERS ##### #
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

    def node_neighbours(self, set_of_nodes, removal_set):
        removal_set = set_of_nodes & removal_set
        nbs = [[n for n in self.G.neighbors(node) if n not in removal_set] for node in set_of_nodes]
        return nbs

    # ###### END OF HELPERS ##### #

    def iterate(self):
        G_copy = self.G
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
        print(self.rates)
        print(self.variant_count)

    def step_run(self):
        if not self.end:
            self.iterate()

            self.old_time = self.time
            self.time = round(self.time + self.dt, 7)
            self.steps += 1

            if (self.steps % 1000) == 0:
                print(".", end='', flush=True)
            if (self.steps % 10000) == 0:
                print()

            if self.time >= self.end_time:
                self.end = True
        #print(self.variant_count)
        return self.time


def main():
    seed = 39
    np.random.seed(seed)
    probs = tuple([(np.random.randint(0, 1000), np.random.randint(0, 1000)) for _ in range(1)])
    execute = SIR(50, initialisation=(0, 0.1), variants=7, probabilities=probs, seed=seed)
    execute.run()
