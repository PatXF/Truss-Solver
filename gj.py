import numpy as np


# Checks if valid matrices are given or not
def check_valid(A: np.array, B: np.array):
    if (A.shape[0] == 0) | (A.shape[0] != B.shape[0]) | (A.shape[0] != A.shape[1]):
        return False
    return True


# finds the nearest pivot (both functions together)
def find_nearest(i, n, arr):
    for j in range(i + 1, n):
        if arr[j] >= i:
            return j
    return -1


def update_least(A: np.array, s: int, n: int):
    arr = np.arange(n)
    for i in range(s + 1, n):
        for j in range(i, n):
            if A[i][j] != 0:
                arr[i] = j
                break
    return arr


# gaussian elimination function to solve a system of linear equations
def solve(A: np.array, n: int, least_non_zero: np.array):
    i = 0
    while i < n:
        if A[i][i] == 0:
            j = find_nearest(i, n, least_non_zero)
            if j == -1:
                return np.array([np.nan for i in range(n)])
            temp = A[i]
            A[i] = A[j]
            A[j] = temp
        else:
            A[i] = A[i] / A[i][i]
            for j in range(n):
                if i == j:
                    continue
                A[j] = A[j] - A[i]*A[j][i]
            least_non_zero = update_least(A, i, n)
        i += 1

    return A


def gj(A: np.array, B: np.array):
    if check_valid(A, B):
        n = A.shape[0]
        least_non_zero = update_least(A, -1, n)
        adj = np.zeros((n, n + 1))
        for i in range(n):
            for j in range(n):
                adj[i][j] = A[i][j]
        for i in range(n):
            adj[i][n] = B[i]
        return solve(adj, n, least_non_zero)[0:n, n]
