from typing import List


def U(x: float, t: float) -> float:
    return x ** 4 + t ** 4


def F(x: float, t: float, v: float) -> float:
    return 4 * t ** 3 - 12 * v * x ** 2


def solve(L: float, T: float, n: int, m: int, v: float) -> List[List[float]]:
    dx: float = L / (n - 1)
    dt: float = T / (m - 1)
    u: List[List[float]] = [[0 for _ in range(n)] for _ in range(m)]
    u[0] = [U(dx * i, 0) for i in range(n)]
    for j in range(1, m):
        alpha: List[float] = [0 for _ in range(n - 1)]
        beta: List[float] = [0 for _ in range(n - 1)]
        beta[0] = U(0, j * dt)
        for i in range(1, n - 1):
            alpha[i] = v * dt / (2 * v * dt + 2 * dx ** 2 - alpha[i - 1] * v * dt)
            beta[i] = (beta[i - 1] * v * dt + v * dt * (u[j - 1][i + 1] + u[j - 1][i - 1]) + (
                    2 * dx ** 2 - 2 * v * dt) * u[j - 1][i]
                       + 2 * dx ** 2 * dt * F(i * dx, j * dt, v)) / (2 * v * dt + 2 * dx ** 2 - alpha[i - 1] * v * dt)
        u[j][-1], u[j][0] = U(dx * (n - 1), j * dt), U(0, j * dt)
        for k in range(n - 2, -1, -1):
            u[j][k] = alpha[k] * u[j][k + 1] + beta[k]
    return u


def deviation(u: List[List[float]], u_real: List[List[float]]) -> float:
    return max(abs(u[i][j] - u_real[i][j]) for j in range(len(u[0])) for i in range(len(u)))


def main():
    L: float = 1  # rod length
    n: int = 400  # number of x steps
    T: float = 1  # period of time
    m: int = 200  # number of T steps
    v: float = 0.019  # heat conduction coefficient

    u: List[List[float]] = solve(L, T, n, m, v)  # solution for dx and dt


if __name__ == "__main__":
    main()
