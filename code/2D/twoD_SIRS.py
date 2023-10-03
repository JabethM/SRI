import numpy as np


class sim:

    def __init__(self, dt, dim, sim_type, p1=0.5, p2=0.1, p3=0.015): # Pulsating Probs (0.5, 0.1, 0.015)
        self.dt = dt
        self.time = 0

        if type(dim) != tuple:
            raise TypeError("dim parameter not given as tuple of integers")
        self.dim = dim
        self.dim_size = len(self.dim)

        self.type = sim_type  # 0 = Batch update; 1 = random sequential update

        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

        self.board = np.random.randint(3, size=self.dim)
        self.board = np.zeros(self.dim)
        self.board[25, 25] = 1

    def iterate(self):

        if self.type == 0:
            s, i, r = self.neighbours()  # recovered, potential, susceptible

            probs = np.random.rand(*self.dim)

            check = s + i + r
            """
            s_probs = s_probs * np.where(self.board == 0, 1, 0)
            i_probs = i_probs * np.where(self.board == 1, 1, 0)
            r_probs = r_probs * np.where(self.board == 2, 1, 0)
            """

            updated_board = np.copy(self.board)
            updated_board[(self.board == 0) & (i >= 1) & (probs < self.p1)] = 1
            updated_board[(self.board == 1) & (probs < self.p2)] = 2
            updated_board[(self.board == 2) & (probs < self.p3)] = 0
            self.board = updated_board
        else:
            selection = np.random.randint(0, self.dim, size=(self.dim, 2))
            for coordinate in selection:
                s, i, r = self.neighbours()
                s_single = s[coordinate[0], coordinate[1]]
                i_single = i[coordinate[0], coordinate[1]]
                r_single = r[coordinate[0], coordinate[1]]
                b0 = coordinate[0]
                b1 = coordinate[1]

                b = self.board[coordinate[0], coordinate[1]]

                probs = np.random.random(3)

                if b == 0 and (i_single > 0) and (probs[0] > self.p1):
                    self.board[coordinate[0], coordinate[1]] = 1

                elif b == 1 and (r_single > 0) and (probs[1] > self.p2):
                    self.board[coordinate[0], coordinate[1]] = 2

                elif b == 2 and (s_single > 0) and (probs[2] > self.p3):
                    self.board[coordinate[0], coordinate[1]] = 0

        self.time += self.dt
        return self.board

    def neighbours(self):
        S_board = np.where(self.board == 0, 1, 0)
        R_board = np.where(self.board == 1, 1, 0)
        I_board = np.where(self.board == 2, 1, 0)

        nbs = lambda x: (np.roll(x, 1, axis=0) + np.roll(x, -1, axis=0) +
                         np.roll(x, 1, axis=1) + np.roll(x, -1, axis=1))

        S_board = nbs(S_board)
        R_board = nbs(R_board)
        I_board = nbs(I_board)

        return S_board, R_board, I_board

    def minority(self):
        R_board = np.sum(np.where(self.board == 0, 1, 0))
        S_board = np.sum(np.where(self.board == 1, 1, 0))
        P_board = np.sum(np.where(self.board == 2, 1, 0))

        minority = np.min([R_board, S_board, P_board])
        return minority


def main():
    dt = 0.5
    L = 10
    sim_type = 0
    b = sim(dt, L, sim_type)
