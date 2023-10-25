import numpy as np
import random


class Variant:
    relation_matrix = None

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
