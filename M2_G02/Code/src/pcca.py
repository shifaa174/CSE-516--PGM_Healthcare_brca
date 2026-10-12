"""
pcca.py -- The graphical model.
===============================================================================
Two-view latent-variable Gaussian PGM (probabilistic CCA; Bach & Jordan, 2005).

    GRAPH            Z
                   /   \
                  v     v
                 X       Y          X = genomics, Y = histology, Z = latent

    GENERATIVE MODEL
        z   ~ N(0, I_d)
        x|z ~ N(W_x z + mu_x, Psi_x)
        y|z ~ N(W_y z + mu_y, Psi_y)

    FACTORISATION      p(z, x, y) = p(z) p(x|z) p(y|z)
    CONDITIONAL IND.   X  _||_  Y  |  Z

Psi_x and Psi_y are FULL within-modality covariance blocks. Overall Psi is
block-diagonal, and that block-diagonal constraint is exactly how the graph's
conditional independence is imposed: within-modality correlation lives in the
noise blocks, so only CROSS-modality correlation is forced through Z.

PARAMETER LEARNING -- two routes on the SAME model:
    fit_cca_mle()  closed-form MLE via classical CCA   [base paper / baseline]
    fit_em()       Expectation-Maximization            [our extension]

STRUCTURE LEARNING -- not performed. The graph is a modelling commitment, not
something searched for; see the M2 report. The one structural hyperparameter
is the latent dimension d, selected on the validation split.

INFERENCE -- exact, closed form (Gaussian conditioning). No approximation.
    posterior(X=..., Y=...)   p(z | x, y)
    posterior(X=...)          p(z | x)     [histology missing]
    posterior(Y=...)          p(z | y)     [genomics missing]
    map_estimate(...)         argmax_z p(z | .)  == posterior mean
    score(X, Y)               log p(x, y)  -- the evidence / learning KPI

Group G02, CSE516 Probabilistic Graphical Models.
"""
import numpy as np
from scipy.linalg import eigh, cholesky, solve_triangular, LinAlgError

__all__ = ["ProbabilisticCCA"]


def _ridge(S, eps):
    """Ridge-regularise a covariance so it stays invertible when n < p."""
    return S + eps * np.trace(S) / S.shape[0] * np.eye(S.shape[0])


def _pd(S, floor=1e-6):
    """Project a symmetric matrix onto the positive-definite cone."""
    S = 0.5 * (S + S.T)
    w, V = eigh(S)
    w = np.maximum(w, floor * max(1.0, w.max()))
    return V @ np.diag(w) @ V.T


def _mvn_loglik(Oc, Sigma):
    """sum_i log N(o_i | 0, Sigma). Sigma is formed explicitly (D is small)."""
    D = Sigma.shape[0]
    try:
        L = cholesky(Sigma, lower=True)
    except LinAlgError:
        L = cholesky(_pd(Sigma), lower=True)
    logdet = 2.0 * np.sum(np.log(np.diag(L)))
    A = solve_triangular(L, Oc.T, lower=True)
    quad = np.einsum("dn,dn->n", A, A)
    return float(np.sum(-0.5 * (D * np.log(2 * np.pi) + logdet + quad)))


class ProbabilisticCCA:
    """The PGM. See module docstring for the model definition.

    Parameters
    ----------
    d : int       latent dimension (structural hyperparameter)
    reg : float   ridge on sample covariances in the CCA route
    random_state : int   seed for EM initialisation (reproducibility)
    """

    def __init__(self, d=10, reg=1e-2, random_state=0, noise="block"):
        self.d = d
        self.reg = reg
        self.random_state = random_state
        self.noise = noise          # "block" (full Psi_x, Psi_y) or "diag"
        self.rng = np.random.default_rng(random_state)

    def _noise_block(self, S):
        """Apply the chosen noise structure to a within-modality covariance.

        "block" keeps the full within-modality covariance -- the faithful
        probabilistic-CCA model, valid when n > p.
        "diag" keeps only the variances -- needed when p >> n, where a full
        within-modality covariance is not identifiable from the sample.
        Either way X _||_ Y | Z still holds; only the richness of the noise
        model changes.
        """
        if self.noise == "diag":
            return np.diag(np.maximum(np.diag(S), 1e-6))
        return _pd(S)

    # ------------------------------------------------------------------ #
    # helpers
    # ------------------------------------------------------------------ #
    def _assemble(self):
        """Stack the per-modality parameters into joint (W, Psi)."""
        W = np.vstack([self.W_x_, self.W_y_])
        p = self.W_x_.shape[0]
        Psi = np.zeros((W.shape[0], W.shape[0]))
        Psi[:p, :p] = self.Psi_x_
        Psi[p:, p:] = self.Psi_y_          # <- block-diagonal == X _||_ Y | Z
        return W, Psi

    def _center(self, X, Y):
        return np.hstack([X - self.mu_x_, Y - self.mu_y_])

    # ------------------------------------------------------------------ #
    # LEARNING ROUTE A -- closed-form MLE via classical CCA  (BASELINE)
    # ------------------------------------------------------------------ #
    def fit_cca_mle(self, X, Y):
        """Maximum-likelihood parameters via classical CCA (Bach & Jordan).

        Solve the canonical correlation problem, then map canonical directions
        and correlations to (W_x, W_y, Psi_x, Psi_y):

            W_x = Sxx A M,  W_y = Syy B M,  M M^T = diag(rho_1..rho_d)
            Psi_x = Sxx - W_x W_x^T,  Psi_y = Syy - W_y W_y^T
        """
        n = X.shape[0]
        self.mu_x_, self.mu_y_ = X.mean(0), Y.mean(0)
        Xc, Yc = X - self.mu_x_, Y - self.mu_y_

        Sxx = _ridge(Xc.T @ Xc / (n - 1), self.reg)
        Syy = _ridge(Yc.T @ Yc / (n - 1), self.reg)
        Sxy = Xc.T @ Yc / (n - 1)

        def invsqrt(S):
            w, V = eigh(S)
            w = np.maximum(w, 1e-10)
            return V @ np.diag(w ** -0.5) @ V.T

        Sxx_i, Syy_i = invsqrt(Sxx), invsqrt(Syy)
        U, s, Vt = np.linalg.svd(Sxx_i @ Sxy @ Syy_i, full_matrices=False)
        d = min(self.d, len(s))
        rho = np.clip(s[:d], 1e-6, 1 - 1e-6)            # canonical correlations
        A, B = Sxx_i @ U[:, :d], Syy_i @ Vt[:d].T       # canonical directions
        M = np.diag(np.sqrt(rho))

        self.W_x_, self.W_y_ = Sxx @ A @ M, Syy @ B @ M
        self.Psi_x_ = self._noise_block(Sxx - self.W_x_ @ self.W_x_.T)
        self.Psi_y_ = self._noise_block(Syy - self.W_y_ @ self.W_y_.T)
        self.rho_, self.d_eff_, self.method_ = rho, d, "cca-mle"
        self.n_iter_ = 1
        self.loglik_trace_ = [self.score(X, Y)]
        return self

    # ------------------------------------------------------------------ #
    # LEARNING ROUTE B -- Expectation-Maximization  (OUR EXTENSION)
    # ------------------------------------------------------------------ #
    def fit_em(self, X, Y, max_iter=300, tol=1e-7, init="random", verbose=False):
        """EM for the same model.

        E-step   q(z_i) = N( G W' Psi^-1 o_i , G ),  G = (I + W' Psi^-1 W)^-1
        M-step   W   <- C (E[zz'])^-1      with C = (1/n) sum_i o_i E[z_i]'
                 Psi <- blockdiag( S - W C' )   <- the CI constraint lives here
        """
        n, p = X.shape
        q = Y.shape[1]
        D, d = p + q, self.d
        self.mu_x_, self.mu_y_ = X.mean(0), Y.mean(0)
        O = self._center(X, Y)
        S = O.T @ O / n

        if init == "cca":
            self.fit_cca_mle(X, Y)
            W, Psi = self._assemble()
        else:
            W = self.rng.normal(scale=0.1, size=(D, d))
            Psi = np.zeros((D, D))
            Psi[:p, :p] = self._noise_block(S[:p, :p])
            Psi[p:, p:] = self._noise_block(S[p:, p:])

        trace, prev = [], -np.inf
        for it in range(max_iter):
            # ---- E-step ----
            Pinv = np.linalg.inv(_pd(Psi))
            WtP = W.T @ Pinv
            G = np.linalg.inv(np.eye(d) + WtP @ W)      # posterior covariance
            Ez = O @ WtP.T @ G.T                        # (n, d) posterior means
            Ezz = G + Ez.T @ Ez / n                     # (1/n) sum E[z z']

            # ---- M-step ----
            C = O.T @ Ez / n
            W = C @ np.linalg.inv(Ezz)
            Full = S - W @ C.T
            Psi = np.zeros((D, D))
            Psi[:p, :p] = self._noise_block(Full[:p, :p])   # block-diagonal
            Psi[p:, p:] = self._noise_block(Full[p:, p:])   # projection

            ll = _mvn_loglik(O, W @ W.T + Psi)
            trace.append(ll)
            if verbose and it % 25 == 0:
                print(f"    EM iter {it:4d}  loglik = {ll:.2f}")
            if abs(ll - prev) < tol * max(1.0, abs(prev)):
                break
            prev = ll

        self.W_x_, self.W_y_ = W[:p], W[p:]
        self.Psi_x_, self.Psi_y_ = Psi[:p, :p], Psi[p:, p:]
        self.d_eff_, self.method_, self.n_iter_ = d, "em", it + 1
        self.loglik_trace_ = trace
        return self

    # ------------------------------------------------------------------ #
    # INFERENCE -- exact Gaussian conditioning
    # ------------------------------------------------------------------ #
    def posterior(self, X=None, Y=None):
        """(mean, cov) of p(z | observed modalities).

        A MISSING MODALITY IS MARGINALISED OUT EXACTLY: its term is simply
        absent from the posterior precision and mean. Nothing is imputed.

        NOTE for the report: in a linear-Gaussian model the posterior
        COVARIANCE does not depend on the observed values -- only on WHICH
        modalities are present. Posterior spread therefore measures modality
        availability, not per-patient difficulty.
        """
        parts = []
        if X is not None:
            parts.append((self.W_x_, self.Psi_x_, X - self.mu_x_))
        if Y is not None:
            parts.append((self.W_y_, self.Psi_y_, Y - self.mu_y_))
        if not parts:
            raise ValueError("at least one modality must be observed")
        S = np.eye(self.d_eff_)
        acc = 0.0
        for W, Psi, Oc in parts:
            WtP = W.T @ np.linalg.inv(_pd(Psi))
            S = S + WtP @ W
            acc = acc + Oc @ WtP.T
        G = np.linalg.inv(S)
        return acc @ G.T, G

    def map_estimate(self, X=None, Y=None):
        """argmax_z p(z | .) -- equals the posterior mean for a Gaussian."""
        return self.posterior(X, Y)[0]

    def posterior_entropy(self, X=None, Y=None):
        """Differential entropy (nats) of the latent posterior."""
        _, G = self.posterior(X, Y)
        _, logdet = np.linalg.slogdet(G)
        return 0.5 * (self.d_eff_ * np.log(2 * np.pi * np.e) + logdet)

    def score(self, X, Y):
        """Joint log-likelihood log p(x, y) of (held-out) data."""
        W, Psi = self._assemble()
        return _mvn_loglik(self._center(X, Y), W @ W.T + Psi)

    def score_per_sample(self, X, Y):
        return self.score(X, Y) / X.shape[0]
