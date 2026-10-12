"""
survival.py -- Downstream survival head and evaluation metric.
===============================================================================
Dependency-free on purpose: `lifelines` does not build in every environment,
and writing the partial likelihood out explicitly makes the objective
auditable for the report.

    CoxPH    L2-penalised Cox proportional hazards (Breslow tie handling),
             fitted by L-BFGS on the partial log-likelihood.
    c_index  Harrell's concordance index, the primary downstream KPI.

Why C-index and not accuracy: ~90% of patients in this cohort are
right-censored (95 events in 974 patients), so accuracy is undefined.
C-index measures whether the model RANKS riskier patients above less risky
ones, and handles censoring correctly.
"""
import numpy as np
from scipy.optimize import minimize

__all__ = ["CoxPH", "c_index"]


def c_index(times, events, risk):
    """Harrell's concordance index. `risk` is higher = worse prognosis.

    Counts comparable pairs: patient i had an event at t_i, patient j was
    still under observation after t_i. Concordant if risk_i > risk_j.
    """
    times, events, risk = map(np.asarray, (times, events, risk))
    num = den = 0.0
    for i in np.where(events == 1)[0]:
        comparable = times > times[i]
        den += comparable.sum()
        num += (risk[comparable] < risk[i]).sum() \
             + 0.5 * (risk[comparable] == risk[i]).sum()
    return num / den if den > 0 else np.nan


class CoxPH:
    """Cox PH fitted by maximising the L2-penalised Breslow partial likelihood.

        log L(beta) = sum_{i: event} [ eta_i - log sum_{j in R(t_i)} exp(eta_j) ]
        eta = Z beta,  R(t) = risk set = patients still observed at time t
    """

    def __init__(self, penalizer=0.1, max_iter=500):
        self.penalizer = penalizer
        self.max_iter = max_iter

    def fit(self, Zf, times, events):
        Zf = np.asarray(Zf, float)
        self.mean_, self.std_ = Zf.mean(0), Zf.std(0) + 1e-8
        Z = (Zf - self.mean_) / self.std_
        order = np.argsort(-np.asarray(times, float))     # descending time
        Z, ev = Z[order], np.asarray(events, float)[order]

        def nll(b):
            eta = Z @ b
            m = eta.max()
            # risk set in descending-time order is a running prefix sum
            logcum = np.log(np.cumsum(np.exp(eta - m)) + 1e-300) + m
            ll = np.sum(ev * (eta - logcum))
            return -ll / len(ev) + self.penalizer * np.sum(b ** 2)

        res = minimize(nll, np.zeros(Z.shape[1]), method="L-BFGS-B",
                       options={"maxiter": self.max_iter})
        self.beta_ = res.x
        self.converged_ = bool(res.success)
        return self

    def risk(self, Zf):
        """Linear predictor eta = z'beta. Higher = higher hazard."""
        return ((np.asarray(Zf, float) - self.mean_) / self.std_) @ self.beta_
