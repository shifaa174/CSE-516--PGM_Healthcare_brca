import numpy as np
from scipy.optimize import minimize

__all__ = ["CoxPH", "c_index"]


def c_index(times, events, risk):


    times, events, risk = map(np.asarray, (times, events, risk))
    num = den = 0.0
    for i in np.where(events == 1)[0]:
        comparable = times > times[i]
        den += comparable.sum()
        num += (risk[comparable] < risk[i]).sum() \
             + 0.5 * (risk[comparable] == risk[i]).sum()
    return num / den if den > 0 else np.nan


class CoxPH:


    def __init__(self, penalizer=0.1, max_iter=500):
        self.penalizer = penalizer
        self.max_iter = max_iter

    def fit(self, Zf, times, events):
        Zf = np.asarray(Zf, float)
        self.mean_, self.std_ = Zf.mean(0), Zf.std(0) + 1e-8
        Z = (Zf - self.mean_) / self.std_
        order = np.argsort(-np.asarray(times, float))
        Z, ev = Z[order], np.asarray(events, float)[order]

        def nll(b):
            eta = Z @ b
            m = eta.max()
            logcum = np.log(np.cumsum(np.exp(eta - m)) + 1e-300) + m
            ll = np.sum(ev * (eta - logcum))
            return -ll / len(ev) + self.penalizer * np.sum(b ** 2)

        res = minimize(nll, np.zeros(Z.shape[1]), method="L-BFGS-B",
                       options={"maxiter": self.max_iter})
        self.beta_ = res.x
        self.converged_ = bool(res.success)
        return self

    def risk(self, Zf):

        return ((np.asarray(Zf, float) - self.mean_) / self.std_) @ self.beta_
