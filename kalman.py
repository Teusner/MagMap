import numpy as np


def kalman_predict(xup, Gup, u, Γα, A):
    Γ1 = A @ Gup @ A.T + Γα
    x1 = A @ xup + u
    return x1, Γ1

def kalman_correc(x0, Γ0, y, Γβ, C):
    S = C @ Γ0 @ C.T + Γβ
    K = Γ0 @ C.T @ np.linalg.inv(S)
    ytilde = y - C @ x0
    Gup = (np.eye(len(x0)) - K @ C) @ Γ0
    xup = x0 + K @ ytilde
    return xup, Gup
    
def kalman(x0, Γ0, u, y, Γα, Γβ, A, C):
    xup, Gup = kalman_correc(x0, Γ0, y, Γβ, C)
    x1, Γ1 = kalman_predict(xup, Gup, u, Γα, A)
    return x1, Γ1