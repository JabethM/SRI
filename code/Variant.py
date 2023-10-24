import numpy as np
import random


class Variant:

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
                else:
                    self.left.insert(data)
            elif data > self.data:
                if self.right is None:
                    self.right = Variant(data, name)
                    self.right.parent = self
                else:
                    self.right.insert(data)
            else:
                self.data = data
