# problems.py
import math
import numpy as np
import matlab.engine

ROUND_DECIMALS = 12


def round_array(arr, decimals=ROUND_DECIMALS):
    return np.round(np.asarray(arr, dtype=float), decimals)


def to_2d_rounded(x, decimals=ROUND_DECIMALS):
    xv = np.atleast_2d(x).astype(float)
    return np.round(xv, decimals)


# Start MATLAB engine once here for the simulator problem
eng = matlab.engine.start_matlab()
eng.addpath(r'C:\Documents\Symbolicregression\Codes', nargout=0)  #Simulator file path


# ====================================================
# Problem evaluation functions
# ====================================================

def evaluate_bnh(x):
    xv = to_2d_rounded(x)
    F_list, G_list = [], []

    for xi in xv:
        f1 = 4*xi[0]**2 + 4*xi[1]**2
        f2 = (xi[0]-5)**2 + (xi[1]-5)**2
        g1 = (xi[0]-5)**2 + xi[1]**2 - 25
        g2 = 7.7 - ((xi[0]-8)**2 + (xi[1]+3)**2)
        F_list.append([f1, f2])
        G_list.append([g1, g2])

    return round_array(F_list), round_array(G_list)


def evaluate_tnk(x):
    xv = to_2d_rounded(x)
    F_list, G_list = [], []

    for xi in xv:
        f1 = xi[0]
        f2 = xi[1]
        theta = np.arctan2(xi[0], xi[1])
        g1 = 1 - xi[0]**2 - xi[1]**2 + 0.1 * np.cos(16*theta)
        g2 = (xi[0]-0.5)**2 + (xi[1]-0.5)**2 - 0.5
        F_list.append([f1, f2])
        G_list.append([g1, g2])

    return round_array(F_list), round_array(G_list)

def evaluate_simulator(x):
    xv = to_2d_rounded(x)
    F_list = []

    for xi in xv:
        d1, d2, d3, mat_code = xi

        if 0 <= mat_code < 1/3:
            mat = 1
        elif 1/3 <= mat_code < 2/3:
            mat = 2
        else:
            mat = 3

        md, mass = eng.simulation(float(d1), float(d2), float(d3), float(mat), nargout=2)
        F_list.append([float(md), float(mass)])

    F = round_array(F_list)
    n = len(F_list)
    G = np.zeros((n, 0))
    return F, G

def evaluate_zdt3(x):
    xv = to_2d_rounded(x)
    F_list = []

    for xi in xv:
        n = len(xi)
        f1 = xi[0]
        g = 1.0 + (9.0 / (n - 1)) * sum(xi[1:])
        h = 1.0 - math.sqrt(f1 / g) - (f1 / g) * math.sin(10.0 * math.pi * f1)
        f2 = g * h
        F_list.append([f1, f2])

    F = round_array(F_list)
    n_samples = len(F_list)
    G = np.zeros((n_samples, 0))
    return F, G


def evaluate_dtlz2(x):
    xv = to_2d_rounded(x)
    F_list = []

    for xi in xv:
        g = sum((xj - 0.5)**2 for xj in xi[1:])
        f1 = (1.0 + g) * math.cos(xi[0] * math.pi / 2.0)
        f2 = (1.0 + g) * math.sin(xi[0] * math.pi / 2.0)
        F_list.append([f1, f2])

    F = round_array(F_list)
    n_samples = len(F_list)
    G = np.zeros((n_samples, 0))
    return F, G


def evaluate_welded_beam(x):
    xv = to_2d_rounded(x)
    F_list, G_list = [], []

    for xi in xv:
        x1, x2, x3, x4 = xi
        f1 = 1.10471 * x1**2 * x2 + 0.04811 * x3 * x4 * (14.0 + x2)
        f2 = 2.1952 / (x4 * x3**3)

        P = 6000.0
        L = 14.0
        t_max = 13600.0
        s_max = 30000.0

        R = np.sqrt(0.25 * (x2**2 + (x1 + x3)**2))
        M = P * (L + x2 / 2.0)
        J = 2.0 * np.sqrt(0.5) * x1 * x2 * (x2**2 / 12.0 + 0.25 * (x1 + x3)**2)
        t1 = P / (np.sqrt(2.0) * x1 * x2)
        t2 = M * R / J
        t = np.sqrt(t1**2 + t2**2 + t1 * t2 * x2 / R)
        s = 6.0 * P * L / (x4 * x3**2)
        P_c = 64746.022 * (1.0 - 0.0282346 * x3) * x3 * x4**3

        g1 = t - t_max
        g2 = s - s_max
        g3 = x1 - x4
        g4 = P - P_c

        F_list.append([f1, f2])
        G_list.append([g1, g2, g3, g4])

    return round_array(F_list), round_array(G_list)


# ====================================================
# Problem dictionaries
# ====================================================

problem_dtlz2 = {
    "concepts": {
        "dtlz2": {
            "eval_func": evaluate_dtlz2,
            "D": 5,
            "xl": [0, 0, 0, 0, 0],
            "xu": [1, 1, 1, 1, 1],
        }
    }
}

problem_zdt3 = {
    "concepts": {
        "zdt3": {
            "eval_func": evaluate_zdt3,
            "D": 5,
            "xl": [0, 0, 0, 0, 0],
            "xu": [1, 1, 1, 1, 1],
        }
    }
}

problem_bnh = {
    "concepts": {
        "bnh": {"eval_func": evaluate_bnh, "D": 2, "xl": [0, 0], "xu": [5, 3]}
    }
}

problem_tnk = {
    "concepts": {
        "tnk": {
            "eval_func": evaluate_tnk,
            "D": 2,
            "xl": [0, 0],
            "xu": [np.pi, np.pi],
        }
    }
}

problem_simulator = {
    "concepts": {
        "simulator": {
            "eval_func": evaluate_simulator,
            "D": 4,
            "xl": [0.5, 0.5, 0.5, 0],
            "xu": [3.5, 7.0, 2.5, 1],
        }
    }
}

problem_welded_beam = {
    "concepts": {
        "welded_beam": {
            "eval_func": evaluate_welded_beam,
            "D": 4,
            "xl": [0.125, 0.1, 0.1, 0.125],
            "xu": [5.0, 10.0, 10.0, 5.0],
        }
    }
}
