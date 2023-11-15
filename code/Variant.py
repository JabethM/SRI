import numpy as np
import random


class Variant:
    relation_matrix = None
    current_data_set = []

    def __init__(self, data, name=None):
        self.left = None
        self.right = None
        self.parent = None
        self.data = data
        self.name = name
        self.relation = {}

    def insert(self, data, name=None):
        if self.data:
            if data < self.data:
                if self.left is None:
                    self.left = Variant(data, name)
                    self.left.parent = self
                    return self.left
                else:
                    return self.left.insert(data, name)
            elif data > self.data:
                if self.right is None:
                    self.right = Variant(data, name)
                    self.right.parent = self
                    return self.right
                else:
                    return self.right.insert(data, name)
            else:
                self.data = data
                return self

    def add_to_relation_matrix(self):

        range_of_variant_data = max(Variant.current_data_set) - min(Variant.current_data_set)
        num_of_variants = len(self.current_data_set)

        if num_of_variants <= 1:
            self.relation_matrix = np.array([[1]])
        else:
            old_matrix = self.relation_matrix
            new_length = num_of_variants  # name change for clarity

            # e^(-1 *(Delta(data) / range(data)) * (number of variants between))
            variant_relation = [
                [np.exp(-1 * (abs(self.current_data_set[j] - self.current_data_set[i]) / range_of_variant_data)
                        * abs(j - i)) if j != i else 1
                 for j in range(new_length)]
                for i in range(new_length)]

            variant_relation[:num_of_variants, :num_of_variants] = old_matrix
            self.relation_matrix = variant_relation
            return self.relation_matrix
