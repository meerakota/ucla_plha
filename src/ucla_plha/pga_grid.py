"""Magnitude x ln(PGA) grid of the ground motion rates of a site.

The liquefaction hazard of a soil profile is a sum over the events (ruptures, with the weights of
the source models, ground motion models, and their branches) of rate x P(FS < fsl | event). For a
model whose conditional probability depends on the event only through its magnitude m and its
peak acceleration (NGL SMT 2024 with the PGA shared by all layers of the profile, Boulanger and
Idriss 2016, and Moss et al. 2006 conditioned on PGA),

    sum_i rate_i P(FS < fsl | event i)
        = sum_i rate_i integral P(FS < fsl | m_i, ln PGA) f(ln PGA | mu_i, sigma_i) d ln PGA
        ~ sum_j sum_k R[j, k] P(FS < fsl | m_j, x_k)

where R is a grid of rates over magnitude nodes m_j and ln(PGA) nodes x_k (RateGrid). R does not
depend on the liquefaction model or the soil profile, so it is built once per site and shared by
the liquefaction models, and the conditional probabilities of a model are needed only at the grid
nodes (ConditionalTable).

Quadrature in ln(PGA): the nodes x_k are the centres of bins of width d_ln_pga between pga_min
and pga_max (the first and last bins extend to -inf and +inf). The integral over a bin is
P(x_k) times the probability of the bin, plus first and second moment corrections with central
differences of P (P' and P''), so that the rule is exact for P quadratic over three bins. The
ground motion distribution is log-normal, truncated at mu + truncation_level sigma and
renormalized as in the PSHA (plha.get_exceedance, nshmp-lib TRUNCATION_UPPER_ONLY), or
untruncated.

Magnitude: linear interpolation of P between magnitude nodes of spacing d_m (the rate of an event
is shared between its two nodes).

Ground motion distributions: the rates of the events are accumulated in a histogram over nodes of
mu (spacing d_mu) and sigma (spacing d_sigma), and the histogram is converted to R with the
quadrature weights of the nodes. An event is shared between its two mu nodes (linear weights),
which preserves the mean and adds a(1 - a) d_mu^2 to the variance; this is subtracted from
sigma^2, and the event is then shared between its two sigma nodes with weights linear in sigma^2,
so that the mixture has the mean and variance of the event (the higher moments differ by terms of
order d_mu^4 and d_sigma^2). Means below ln(pga_min) - 5 or above ln(pga_max) are moved to these
limits (their ground motions are in the first or last bins).
"""

import numpy as np
from scipy.special import ndtr

# Default grid parameters (config "pga_grid" of a liquefaction model). With these, the hazard
# curves of NGL SMT 2024 at the benchmark sites agree with an exact per-event integration within
# about 2e-4 (relative) at the rates of interest; see the CHANGELOG.
DEFAULTS = {
    "pga_min": 1.0e-4,
    "pga_max": 10.0,
    "d_ln_pga": 0.1,
    "d_m": 0.05,
    "d_mu": 0.05,
    "d_sigma": 0.005,
}

# events per chunk in the histogram and the per-event weights (bounds the memory)
_CHUNK = 20000


def quadrature_cdfs(
    conditional_cdfs, m, mu_ln_pga, sigma_ln_pga, truncation_level=None, n_quadrature=128,
    chunk=2000,
):
    """P(FS < fsl | event) by Gauss-Legendre quadrature of conditional probabilities over the
    (truncated) normal distribution of ln(PGA) of each event, from -9 standard deviations to the
    truncation level (9 if untruncated), renormalized to the probability between the limits.

    Inputs:
        conditional_cdfs (callable): conditional_cdfs(m, ln_pga), arrays of length K, returns
            P(FS < fsl | m, PGA), shape = K x L
        m, mu_ln_pga, sigma_ln_pga (Numpy arrays): magnitude, mean and standard deviation of
            ln(PGA) of each event, length = N (m and sigma_ln_pga may be scalars)
        truncation_level (float): upper truncation in standard deviations (None: untruncated)
        n_quadrature (int): number of quadrature points

    Returns:
        fsl_cdfs (Numpy ndarray): shape = N x L
    """
    mu_ln_pga = np.asarray(mu_ln_pga, dtype=float)
    sigma_ln_pga = np.broadcast_to(np.asarray(sigma_ln_pga, dtype=float), mu_ln_pga.shape)
    m = np.broadcast_to(np.asarray(m, dtype=float), mu_ln_pga.shape)
    z_lo = -9.0
    z_hi = 9.0 if truncation_level is None else float(truncation_level)
    x, w = np.polynomial.legendre.leggauss(int(n_quadrature))
    z = 0.5 * (z_hi - z_lo) * x + 0.5 * (z_hi + z_lo)
    w = w * 0.5 * (z_hi - z_lo) * np.exp(-0.5 * z * z) / np.sqrt(2.0 * np.pi)
    # renormalized to the probability between the limits (the truncated distribution)
    w = w / (ndtr(z_hi) - ndtr(z_lo))
    out = None
    for start in range(0, len(mu_ln_pga), chunk):
        sl = slice(start, start + chunk)
        n = len(mu_ln_pga[sl])
        ln_pga = mu_ln_pga[sl, np.newaxis] + sigma_ln_pga[sl, np.newaxis] * z[np.newaxis, :]
        p = conditional_cdfs(np.repeat(m[sl], len(z)), ln_pga.ravel())
        p = p.reshape(n, len(z), p.shape[-1])
        if out is None:
            out = np.zeros((len(mu_ln_pga), p.shape[-1]))
        out[sl] = np.einsum("q,nql->nl", w, p)
    if out is None:
        out = np.zeros((0, conditional_cdfs(np.array([6.0]), np.array([0.0])).shape[-1]))
    return out


def grid_parameters(pga_grid=None):
    """Grid parameters: DEFAULTS updated with the given dict (config "pga_grid")."""
    params = dict(DEFAULTS)
    if pga_grid:
        unknown = set(pga_grid) - set(DEFAULTS)
        if unknown:
            raise ValueError(f"unknown pga_grid parameters: {sorted(unknown)}")
        params.update({k: float(v) for k, v in pga_grid.items()})
    if not 0.0 < params["pga_min"] < params["pga_max"]:
        raise ValueError("pga_grid: pga_min must be positive and smaller than pga_max")
    for key in ("d_ln_pga", "d_m", "d_mu", "d_sigma"):
        if params[key] <= 0.0:
            raise ValueError(f"pga_grid: {key} must be positive")
    return params


def _phi(z):
    with np.errstate(invalid="ignore"):
        out = np.exp(-0.5 * z * z) * (1.0 / np.sqrt(2.0 * np.pi))
    return np.where(np.isfinite(z), out, 0.0)


def bin_weights(mu, sigma, edges, truncation_level=None):
    """Quadrature weights of the ln(PGA) bins for normal distributions of ln(PGA).

    Inputs:
        mu (Numpy array): means of ln(PGA), length = N
        sigma (Numpy array): standard deviations of ln(PGA), length = N
        edges (Numpy array): bin edges of ln(PGA) (uniform spacing), length = n + 1; the first
            and last bins extend to -inf and +inf
        truncation_level (float): upper truncation in standard deviations (None: untruncated)

    Returns:
        W (Numpy ndarray): weights, shape = N x n, such that the expected value of a smooth
            function P of ln(PGA) is approximately W @ P(centres of the bins)
    """
    mu = np.asarray(mu, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    h = edges[1] - edges[0]
    centres = 0.5 * (edges[:-1] + edges[1:])
    z = (edges[np.newaxis, :] - mu[:, np.newaxis]) / sigma[:, np.newaxis]
    z[:, 0] = -np.inf
    z[:, -1] = np.inf
    norm = 1.0
    if truncation_level is not None:
        t = float(truncation_level)
        z = np.minimum(z, t)
        norm = ndtr(t)
    F = ndtr(z)
    ph = _phi(z)
    zph = np.where(np.isfinite(z), z, 0.0) * ph
    # moments of each bin about the mean, then about the bin centre
    m0 = np.diff(F, axis=1) / norm
    a1 = sigma[:, np.newaxis] * (ph[:, :-1] - ph[:, 1:]) / norm
    a2 = sigma[:, np.newaxis] ** 2 * (m0 - np.diff(zph, axis=1) / norm)
    d = mu[:, np.newaxis] - centres[np.newaxis, :]
    m1 = a1 + d * m0
    m2 = a2 + 2.0 * d * a1 + d * d * m0
    # no corrections in the first and last bins (semi-infinite)
    c1 = m1[:, 1:-1] / (2.0 * h)
    c2 = 0.5 * m2[:, 1:-1] / h**2
    W = m0
    W[:, 2:] += c1 + c2
    W[:, :-2] += c2 - c1
    W[:, 1:-1] -= 2.0 * c2
    return W


class RateGrid:
    """Rates of the ground motions of a site on a grid of magnitude and ln(PGA) (see the module
    docstring).

    Inputs:
        truncation_level (float): upper truncation of the ground motion distribution in standard
            deviations, as in the PSHA (None: untruncated)
        **params: grid parameters (DEFAULTS)
    """

    def __init__(self, truncation_level=None, **params):
        params = grid_parameters(params)
        self.params = params
        self.truncation_level = truncation_level
        lo, hi = np.log(params["pga_min"]), np.log(params["pga_max"])
        h = params["d_ln_pga"]
        n = max(int(np.ceil((hi - lo) / h - 1e-9)), 3)
        self.edges = lo + h * np.arange(n + 1)
        self.ln_pga = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.d_m = params["d_m"]
        self.d_mu = params["d_mu"]
        self.d_sigma = params["d_sigma"]
        self.mu_lo = lo - 5.0
        self.n_mu = int(np.ceil((hi - self.mu_lo) / self.d_mu)) + 1
        self.mu_hi = self.mu_lo + self.d_mu * (self.n_mu - 1)
        # histogram over (magnitude node, sigma node, mu node), magnitude nodes j0.. and sigma
        # nodes s0.. (bounds grown as needed)
        self._hist = None
        self._j0 = 0
        self._s0 = 0
        self._rates = None

    def magnitude_nodes(self, m):
        """Lower magnitude node (integer, node magnitude = j d_m) and the linear weight of the
        upper node, for each magnitude."""
        x = np.asarray(m, dtype=float) / self.d_m
        j = np.floor(x + 1e-9).astype(np.int64)
        return j, np.clip(x - j, 0.0, 1.0)

    def _grow(self, j_lo, j_hi, s_lo, s_hi):
        """Extends the histogram to magnitude nodes j_lo..j_hi and sigma nodes s_lo..s_hi."""
        if self._hist is None:
            self._j0, self._s0 = j_lo, s_lo
            self._hist = np.zeros((j_hi - j_lo + 1, s_hi - s_lo + 1, self.n_mu))
            return
        nj, ns, _ = self._hist.shape
        j0, s0 = min(j_lo, self._j0), min(s_lo, self._s0)
        j1, s1 = max(j_hi, self._j0 + nj - 1), max(s_hi, self._s0 + ns - 1)
        if (j0, s0, j1, s1) == (self._j0, self._s0, self._j0 + nj - 1, self._s0 + ns - 1):
            return
        hist = np.zeros((j1 - j0 + 1, s1 - s0 + 1, self.n_mu))
        hist[self._j0 - j0 : self._j0 - j0 + nj, self._s0 - s0 : self._s0 - s0 + ns] = self._hist
        self._hist, self._j0, self._s0 = hist, j0, s0

    def add(self, m, mu_ln_pga, sigma_ln_pga, weight):
        """Adds events: magnitudes, means and standard deviations of ln(PGA), and weights (rate
        times the logic tree weights), each of length N (sigma_ln_pga may be a scalar)."""
        m = np.asarray(m, dtype=float)
        mu = np.asarray(mu_ln_pga, dtype=float)
        sigma = np.broadcast_to(np.asarray(sigma_ln_pga, dtype=float), mu.shape)
        weight = np.broadcast_to(np.asarray(weight, dtype=float), mu.shape)
        keep = weight != 0.0
        if not np.all(keep):
            m, mu, sigma, weight = m[keep], mu[keep], sigma[keep], weight[keep]
        if len(m) == 0:
            return
        self._rates = None
        for start in range(0, len(m), _CHUNK):
            sl = slice(start, start + _CHUNK)
            self._add(m[sl], mu[sl], sigma[sl], weight[sl])

    def _add(self, m, mu, sigma, weight):
        u = (np.clip(mu, self.mu_lo, self.mu_hi) - self.mu_lo) / self.d_mu
        ui = np.minimum(np.floor(u).astype(np.int64), self.n_mu - 2)
        uf = u - ui
        j, jf = self.magnitude_nodes(m)
        # variance of the event less that added by the split between the mu nodes
        s2 = np.maximum(sigma**2 - uf * (1.0 - uf) * self.d_mu**2, self.d_sigma**2)
        si = np.maximum(np.floor(np.sqrt(s2) / self.d_sigma).astype(np.int64), 1)
        sf = np.clip((s2 / self.d_sigma**2 - si**2) / (2 * si + 1), 0.0, 1.0)
        self._grow(int(j.min()), int(j.max()) + 1, int(si.min()), int(si.max()) + 1)
        nj, ns, nu = self._hist.shape
        base = ((j - self._j0) * ns + (si - self._s0)) * nu + ui
        keys = []
        weights = []
        for dj, wj in ((0, 1.0 - jf), (ns * nu, jf)):
            for ds, ws in ((0, 1.0 - sf), (nu, sf)):
                for du, wu in ((0, 1.0 - uf), (1, uf)):
                    keys.append(base + (dj + ds + du))
                    weights.append(weight * wj * ws * wu)
        keys = np.concatenate(keys)
        k0 = int(keys.min())
        counts = np.bincount(keys - k0, weights=np.concatenate(weights))
        flat = self._hist.reshape(-1)
        flat[k0 : k0 + len(counts)] += counts

    def rates(self):
        """Magnitude nodes (j, magnitude j d_m) and the rate grid R, shape = len(j) x n, with
        rows only for nodes with events."""
        if self._rates is None:
            if self._hist is None:
                self._rates = (np.zeros(0, dtype=np.int64), np.zeros((0, len(self.ln_pga))))
                return self._rates
            nj, ns, nu = self._hist.shape
            mu_nodes = self.mu_lo + self.d_mu * np.arange(nu)
            R = np.zeros((nj, len(self.ln_pga)))
            for s in range(ns):
                block = self._hist[:, s, :]
                rows = np.flatnonzero(block.any(axis=1))
                if len(rows) == 0:
                    continue
                cols = np.flatnonzero(block[rows].any(axis=0))
                sigma = (self._s0 + s) * self.d_sigma
                K = bin_weights(
                    mu_nodes[cols], np.full(len(cols), sigma), self.edges, self.truncation_level
                )
                R[rows] += block[np.ix_(rows, cols)] @ K
            used = np.flatnonzero(self._hist.any(axis=(1, 2)))
            self._rates = (self._j0 + used, R[used])
        return self._rates

    def event_weights(self, mu_ln_pga, sigma_ln_pga):
        """Quadrature weights (bin_weights) of the ln(PGA) bins for each event, shape = N x n."""
        mu = np.asarray(mu_ln_pga, dtype=float)
        sigma = np.broadcast_to(np.asarray(sigma_ln_pga, dtype=float), mu.shape)
        return bin_weights(mu, sigma, self.edges, self.truncation_level)

    def nbytes(self):
        """Memory of the histogram and rate grid in bytes."""
        n = 0 if self._hist is None else self._hist.nbytes
        if self._rates is not None:
            n += self._rates[1].nbytes
        return n


class ConditionalTable:
    """Conditional probabilities P(FS < fsl | m_j, x_k) of a soil profile at the nodes of a
    RateGrid, computed when first needed.

    Inputs:
        grid (RateGrid): the grid
        conditional_cdfs (callable): conditional_cdfs(m, ln_pga), with arrays m and ln_pga of
            the same length K, returns the probabilities, shape = K x L
    """

    def __init__(self, grid, conditional_cdfs):
        self.grid = grid
        self.conditional_cdfs = conditional_cdfs
        self._rows = {}

    def rows(self, nodes):
        """Probabilities at the magnitude nodes (integers), shape = len(nodes) x n x L."""
        missing = [j for j in np.unique(nodes) if j not in self._rows]
        if missing:
            x = self.grid.ln_pga
            mm = np.repeat(np.asarray(missing, dtype=float) * self.grid.d_m, len(x))
            xx = np.tile(x, len(missing))
            p = self.conditional_cdfs(mm, xx).reshape(len(missing), len(x), -1)
            for j, pj in zip(missing, p):
                self._rows[int(j)] = pj
        return np.stack([self._rows[int(j)] for j in nodes])

    def curve(self):
        """Sum of R x P over the grid (the hazard of the events added to the grid), length L."""
        nodes, R = self.grid.rates()
        if len(nodes) == 0:
            return 0.0
        return np.maximum(np.einsum("jk,jkl->l", R, self.rows(nodes)), 0.0)

    def event_cdfs(self, m, mu_ln_pga, sigma_ln_pga, weights=None):
        """P(FS < fsl | event) of events (magnitudes, means and standard deviations of ln(PGA),
        each of length N; sigma_ln_pga may be a scalar), from the probabilities at the grid
        nodes: the quadrature weights of the ln(PGA) bins of each event (RateGrid.event_weights,
        or the given weights, N x n, e.g. shared by several tables of the grid) and linear
        interpolation between the magnitude nodes. Returns shape = N x L."""
        m = np.asarray(m, dtype=float)
        mu = np.asarray(mu_ln_pga, dtype=float)
        sigma = np.broadcast_to(np.asarray(sigma_ln_pga, dtype=float), mu.shape)
        j, f = self.grid.magnitude_nodes(m)
        if len(j) == 0:
            n_values = self.conditional_cdfs(np.array([6.0]), np.array([0.0])).shape[1]
            return np.zeros((0, n_values))
        out = np.empty((len(j), self.rows(j[:1]).shape[2]))
        for start in range(0, len(j), _CHUNK):
            sl = slice(start, start + _CHUNK)
            w_chunk = (
                self.grid.event_weights(mu[sl], sigma[sl]) if weights is None else weights[sl]
            )
            jc, fc = j[sl], f[sl]
            part = out[sl]
            for jj in np.unique(jc):
                sel = np.flatnonzero(jc == jj)
                lower, upper = self.rows([jj, jj + 1])
                w = w_chunk[sel]
                part[sel] = (1.0 - fc[sel])[:, np.newaxis] * (w @ lower) + fc[sel][
                    :, np.newaxis
                ] * (w @ upper)
        return np.clip(out, 0.0, 1.0)


def shared_event_cdfs(tables, m, mu_ln_pga, sigma_ln_pga):
    """P(FS < fsl | event) of the events for several ConditionalTables (ConditionalTable.event_cdfs),
    with the quadrature weights of the ln(PGA) bins of each event computed once for each grid
    and shared by the tables of that grid. Processed in chunks of events (bounds the memory).

    Inputs:
        tables (list of ConditionalTable): the tables
        m, mu_ln_pga, sigma_ln_pga (Numpy arrays): magnitude, mean and standard deviation of
            ln(PGA) of each event, length = N (sigma_ln_pga may be a scalar)

    Returns:
        list of Numpy ndarrays, one for each table, shape = N x L
    """
    m = np.asarray(m, dtype=float)
    mu = np.asarray(mu_ln_pga, dtype=float)
    sigma = np.broadcast_to(np.asarray(sigma_ln_pga, dtype=float), mu.shape)
    if len(m) == 0:
        return [table.event_cdfs(m, mu, sigma) for table in tables]
    out = [None] * len(tables)
    for start in range(0, len(m), _CHUNK):
        sl = slice(start, start + _CHUNK)
        weights = {}
        for i, table in enumerate(tables):
            key = id(table.grid)
            if key not in weights:
                weights[key] = table.grid.event_weights(mu[sl], sigma[sl])
            p = table.event_cdfs(m[sl], mu[sl], sigma[sl], weights=weights[key])
            if out[i] is None:
                out[i] = np.empty((len(m), p.shape[1]))
            out[i][sl] = p
    return out
