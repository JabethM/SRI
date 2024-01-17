
class node_colour:
    def __init__(self, four_digit):
        assert(len(four_digit) == 4)
        assert(isinstance(four_digit, tuple))
        assert(isinstance(four_digit[0], float))
        self.four_digit = four_digit

    def set_four_digit(self, four_digit):
        self.four_digit = four_digit
        return

    def get_four_digit(self):
        four_digit = self.four_digit
        return  four_digit
