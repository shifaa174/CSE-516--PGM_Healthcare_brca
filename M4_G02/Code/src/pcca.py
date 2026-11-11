import numpy as np
from scipy.linalg import eigh, cholesky, solve_triangular, LinAlgError

__all__ = ["ProbabilisticCCA"]


def _ridge(S, eps):

    return S + eps * np.trace(S) / S.shape[0] * np.eye(S.shape[0])


def _pd(S, floor=1e-6):

    S = 0.5 * (S + S.T)
    w, V = eigh(S)
    w = np.maximum(w, floor * max(1.0, w.max()))
    return V @ np.diag(w) @ V.T


def _mvn_loglik(Oc, Sigma):

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


    def __init__(self, d=10, reg=1e-2, random_state=0, noise="block"):
        self.d = d
        self.reg = reg
        self.random_state = random_state
        self.noise = noise
        self.rng = np.random.default_rng(random_state)

    def _noise_block(self, S):


        if self.noise == "diag":
            return np.diag(np.maximum(np.diag(S), 1e-6))
        return _pd(S)

    def _assemble(self):

        W = np.vstack([self.W_x_, self.W_y_])
        p = self.W_x_.shape[0]
        Psi = np.zeros((W.shape[0], W.shape[0]))
        Psi[:p, :p] = self.Psi_x_
        Psi[p:, p:] = self.Psi_y_
        return W, Psi

    def _center(self, X, Y):
        return np.hstack([X - self.mu_x_, Y - self.mu_y_])

    def fit_cca_mle(self, X, Y):


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
        rho = np.clip(s[:d], 1e-6, 1 - 1e-6)
        A, B = Sxx_i @ U[:, :d], Syy_i @ Vt[:d].T
        M = np.diag(np.sqrt(rho))

        self.W_x_, self.W_y_ = Sxx @ A @ M, Syy @ B @ M
        self.Psi_x_ = self._noise_block(Sxx - self.W_x_ @ self.W_x_.T)
        self.Psi_y_ = self._noise_block(Syy - self.W_y_ @ self.W_y_.T)
        self.rho_, self.d_eff_, self.method_ = rho, d, "cca-mle"
        self.n_iter_ = 1
        self.loglik_trace_ = [self.score(X, Y)]
        return self

    def fit_em(self, X, Y, max_iter=300, tol=1e-7, init="random", verbose=False):


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
            Pinv = np.linalg.inv(_pd(Psi))
            WtP = W.T @ Pinv
            G = np.linalg.inv(np.eye(d) + WtP @ W)
            Ez = O @ WtP.T @ G.T
            Ezz = G + Ez.T @ Ez / n

            C = O.T @ Ez / n
            W = C @ np.linalg.inv(Ezz)
            Full = S - W @ C.T
            Psi = np.zeros((D, D))
            Psi[:p, :p] = self._noise_block(Full[:p, :p])
            Psi[p:, p:] = self._noise_block(Full[p:, p:])

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

    def posterior(self, X=None, Y=None):


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

        return self.posterior(X, Y)[0]

    def posterior_entropy(self, X=None, Y=None):

        _, G = self.posterior(X, Y)
        _, logdet = np.linalg.slogdet(G)
        return 0.5 * (self.d_eff_ * np.log(2 * np.pi * np.e) + logdet)

    def score(self, X, Y):

        W, Psi = self._assemble()
        return _mvn_loglik(self._center(X, Y), W @ W.T + Psi)

    def score_per_sample(self, X, Y):
        return self.score(X, Y) / X.shape[0]
