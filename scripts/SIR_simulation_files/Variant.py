import numpy as np
import random


class Variant:
    relation_matrix = [[1]]
    current_data_set = []

    def __init__(self, data, name=None):
        self.left = None
        self.right = None
        self.parent = None
        self.data = data
        self.name = name
        self.children = []

    def insert(self, data, name=None):
        if data != self.data:
            self.children.append(Variant(data, name))
        else:
            self.data = data
            if name is not None:
                self.name = name


    @classmethod
    def add_to_relation_matrix(cls):

        range_of_variant_data = max(Variant.current_data_set) - min(Variant.current_data_set)
        num_of_variants = len(Variant.current_data_set)

        if num_of_variants <= 1:
            Variant.relation_matrix = np.array([[1]])
        else:
            old_matrix = Variant.relation_matrix
            new_length = num_of_variants  # name change for clarity

            # e^(-1 *(Delta(data2) / range(data2)) * (number of variants between))
            variant_relation = np.array([
                [np.exp(-1 * (abs(Variant.current_data_set[j] - Variant.current_data_set[i]) / range_of_variant_data)
                        * abs(j - i)) if j != i else 1
                 for j in range(new_length)]
                for i in range(new_length)])

            variant_relation[:num_of_variants - 1, :num_of_variants - 1] = old_matrix
            Variant.relation_matrix = variant_relation
            return Variant.relation_matrix

    @classmethod
    def get_relation(cls):
        a = Variant.relation_matrix
        return a

    def print_tree(self, root, level=0, prefix="Root: "):
        if root is not None:
            print(" " * (level * 4) + prefix + str(root.name))
            for child in root.children:
                self.print_tree(child, level + 1, "Child: ")
