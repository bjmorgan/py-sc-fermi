"""Scale-aware numerical solve of element-pool chemical potentials.

Given site-exclusion groups of variable charge states, this module solves
for one chemical potential per constrained element such that each
element's total content matches its target. The working representation is
`GroupData`: one ``(n_free, log_w, s)`` tuple per group, holding the free
site count, the per-state log Boltzmann weights (shape ``(n,)``), and the
per-state stoichiometry vectors over the constrained elements (shape
``(n, K)``).

The target equations are the stationarity condition of the convex grand
potential ``F(mu) = sum_g n_free_g * log(Z_g(mu)) - mu . target``, whose
gradient is the element content minus the target and whose Hessian is,
group by group, ``n_free`` times the per-site covariance of the
stoichiometry vectors -- positive semi-definite, so `F` has a single
minimum whenever the targets are jointly achievable.

Solver strategy
---------------
A single off-the-shelf root-find does not suffice, because defect
concentrations span ~1e-30..1 per cell and every general-purpose
convergence test is absolute or norm-based, so it declares success while
dilute elements are still far from their targets: at ``mu = 0`` the
residual of a 1e-18 target is already ~1e-18, below any absolute
threshold, so the unconstrained populations are accepted as a solution.
The solve therefore runs in up to three stages, each engaged only when
the previous stage's *independently measured* relative deviation misses
tolerance:

1. A Newton root-find on the per-element-scaled residual
   ``content_X(mu) / target_X - 1`` (scipy ``hybr``), seeded out of the
   underflow plateau at ``mu = 0``. Scaling makes every residual O(1)
   however dilute the targets; this settles the common case.
2. A Newton root-find on the log residual ``ln(content_X / target_X)``,
   whose full Jacobian resolves shared-species coupling (one species
   carrying most of several pools' content) that the first stage's seed
   cannot. Second because the log flattens at saturation.
3. A derivative-free bracketing sweep: cyclic exact line solves of `F`
   that need only the sign of the residual, so the flat plateaux that
   stall Newton -- saturated states, extreme dilution, and zero/negative
   net-content targets, which have no scaled residual or log seed at all
   -- cannot defeat it. Globally convergent for feasible targets, so it
   is also the entry point whenever any target is non-positive.

Correctness does not rest on any stage's own verdict: the achieved
content is recomputed and checked against every target before `mu` is
returned, fails closed on NaN, and raises `ElementPoolError` otherwise.
The guarantee is therefore asymmetric -- a solve that meets tolerance is
correct, and a feasible solve the stages happen to miss raises loudly
rather than returning a wrong number.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import bisect, root
from scipy.special import logsumexp

GroupData = list[tuple[float, np.ndarray, np.ndarray]]

# Maximum acceptable relative deviation of each element's variable-state
# content from its remaining target (the pool target less fixed
# contributions); a solve outside this is raised as an error rather than
# returned as silently unconstrained concentrations.
_element_pool_tolerance = 1e-6


class ElementPoolError(ValueError):
    """Raised when element-pool constraints cannot be satisfied: targets
    inconsistent with fixed concentrations, infeasible against site
    capacities, or unmet by the chemical-potential solve."""


def _sanitised_mu(result_x: np.ndarray) -> np.ndarray:
    """A root-find's returned `mu`, with non-finite components mapped back
    into the representable range (NaN to 0, infinities to +/-700, where
    700 < log(largest double) keeps exp(mu) finite)."""
    return np.nan_to_num(
        np.asarray(result_x, dtype=float), nan=0.0, posinf=700.0, neginf=-700.0
    )


def group_concs(
    n_free: float, log_w: np.ndarray, s: np.ndarray, mu: np.ndarray
) -> np.ndarray:
    """c_i = n_free * w_i * exp(s_i . mu) / (1 + sum_j w_j * exp(s_j . mu))
    for every variable state in a group, given the element chemical
    potentials `mu`."""
    exponents = log_w + s @ mu
    log_z = logsumexp(np.concatenate(([0.0], exponents)))
    return n_free * np.exp(exponents - log_z)


def content_and_hessian(
    group_data: GroupData, mu: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Total element content and its Jacobian d(content)/d(mu), summed
    over all groups, at element chemical potentials `mu`.

    Per group, the Jacobian is `n_free` times the per-site covariance
    of the stoichiometry vectors over the full state ensemble,
    including the empty state (stoichiometry zero, occupancy
    ``n_free - sum_i c_i``)."""
    K = len(mu)
    content = np.zeros(K)
    hessian = np.zeros((K, K))
    for n_free, log_w, s in group_data:
        c = group_concs(n_free, log_w, s, mu)
        content += s.T @ c
        mean_s = (s * c[:, None]).sum(axis=0) / n_free
        ds = s - mean_s
        hessian += (ds * c[:, None]).T @ ds
        hessian += max(n_free - c.sum(), 0.0) * np.outer(mean_s, mean_s)
    return content, hessian


def scaled_deviation(
    group_data: GroupData, mu: np.ndarray, remaining_vec: np.ndarray
) -> np.ndarray:
    """Per-element relative deviation of the content from its target at
    `mu`: ``|content / target - 1|`` for non-zero targets; for a zero
    target (a net-content balance over mixed-sign stoichiometries),
    ``|content|`` relative to the gross content ``sum_i |s_ie| c_i``,
    the natural scale of the balance. A zero target with zero gross
    content has no contributing occupancy at all and reads as exactly
    met (deviation 0)."""
    K = len(remaining_vec)
    content = np.zeros(K)
    gross = np.zeros(K)
    for n_free, log_w, s in group_data:
        c = group_concs(n_free, log_w, s, mu)
        content += s.T @ c
        gross += np.abs(s).T @ c
    deviation = np.empty(K)
    nonzero = remaining_vec != 0.0
    deviation[nonzero] = np.abs(content[nonzero] / remaining_vec[nonzero] - 1.0)
    zero = ~nonzero
    with np.errstate(invalid="ignore"):
        deviation[zero] = np.abs(content[zero]) / gross[zero]
    deviation[zero & (gross == 0.0)] = 0.0
    return deviation


def solve_along_direction(
    group_data: GroupData, mu: np.ndarray, direction: np.ndarray, target: float
) -> np.ndarray:
    """Solve ``direction @ content(mu + delta * direction) = target``
    for the scalar `delta` and return the updated `mu`. The projected
    content is non-decreasing in `delta` (its derivative is
    ``direction @ H @ direction >= 0`` with `H` positive
    semi-definite), so a sign change brackets a root. Bisection uses
    only signs, so the near-step profiles of saturating states cannot
    defeat it, but resolves `delta` only to an absolute floor; a Newton
    polish (below) recovers the precision a small net target needs. If no
    sign change exists within the expansion cap, the bracket edge is
    returned: the other directions must move first."""

    def excess_and_slope(delta: float) -> tuple[float, float]:
        content, hessian = content_and_hessian(group_data, mu + delta * direction)
        return float(direction @ content - target), float(
            direction @ hessian @ direction
        )

    f0, _ = excess_and_slope(0.0)
    if f0 == 0.0:
        return mu
    # Each failed probe becomes the near edge of the bracket, keeping
    # its width to the final step rather than the accumulated span.
    step = 1.0
    if f0 > 0:
        hi, lo = 0.0, -1.0
        while excess_and_slope(lo)[0] > 0:
            if lo < -1e7:
                return mu + lo * direction
            hi = lo
            step *= 2
            lo -= step
    else:
        lo, hi = 0.0, 1.0
        while excess_and_slope(hi)[0] < 0:
            if hi > 1e7:
                return mu + hi * direction
            lo = hi
            step *= 2
            hi += step
    delta = float(bisect(lambda d: excess_and_slope(d)[0], lo, hi, maxiter=200))

    # Bisection's absolute step floor (xtol ~ 2e-12) leaves the content
    # resolved to ~slope * 2e-12, far coarser than a small net target over
    # near-half-saturated sites needs (the answer is a cancellation of two
    # O(n_free) populations). One-dimensional Newton on the line --
    # g'(delta) = direction @ H @ direction, the exact slope -- contracts
    # quadratically from the bracketed point wherever the slope is
    # informative, and is simply skipped on the saturation plateaux where
    # the slope underflows and bisection is already at the float limit.
    g, slope = excess_and_slope(delta)
    best_delta, best_excess = delta, abs(g)
    for _ in range(50):
        if slope <= 0.0 or g == 0.0:
            break
        delta = delta - g / slope
        g, slope = excess_and_slope(delta)
        if abs(g) >= best_excess:
            break
        best_delta, best_excess = delta, abs(g)
    return mu + best_delta * direction


def bracketed_coordinate_solve(
    group_data: GroupData, remaining_vec: np.ndarray
) -> np.ndarray:
    """Scale-free fallback solve: cyclic exact line solves of the
    convex grand potential along each element's own chemical potential
    and along the collective all-ones direction, plus, each sweep, the
    regularised Newton direction ``(H + eps I)^-1 (target - content)``.
    No fixed set of directions spans every geometry's slow mode --
    plain coordinate descent ping-pongs when a shared species
    dominates several pools -- so the Newton line, which adapts to
    the local geometry and contracts quadratically wherever `H` is
    informative, carries the hard cases. Convergent for jointly
    feasible targets from any finite start; begins at the bounded,
    physically meaningful mu = 0 (the unconstrained populations).
    Returns the iterate with the smallest deviation; the caller
    verifies the targets and raises if they are unmet."""
    K = len(remaining_vec)
    mu = np.zeros(K)
    fixed_directions: list[tuple[np.ndarray, float]] = [
        (np.eye(K)[k], float(remaining_vec[k])) for k in range(K)
    ]
    if K > 1:
        fixed_directions.append((np.ones(K), float(remaining_vec.sum())))
    best = np.inf
    best_mu = mu.copy()
    stalled = 0
    for _ in range(500):
        for direction, target in fixed_directions:
            mu = solve_along_direction(group_data, mu, direction, target)
        content, hessian = content_and_hessian(group_data, mu)
        h_scale = float(np.trace(hessian)) / K
        if h_scale > 0:
            try:
                newton = np.linalg.solve(
                    hessian + 1e-12 * h_scale * np.eye(K), remaining_vec - content
                )
            except np.linalg.LinAlgError:
                # the regulariser underflows to zero against a Hessian at
                # denormal scale, leaving an exact zero mode unregularised
                newton = remaining_vec - content
        else:
            newton = remaining_vec - content
        norm = np.abs(newton).max()
        if norm > 0 and np.isfinite(norm):
            direction = newton / norm
            mu = solve_along_direction(
                group_data, mu, direction, float(direction @ remaining_vec)
            )
        deviation = scaled_deviation(group_data, mu, remaining_vec).max()
        if deviation >= 0.97 * best:
            stalled += 1
            if stalled >= 10:
                break
        else:
            stalled = 0
        if deviation < best:
            best, best_mu = deviation, mu.copy()
        if deviation <= 1e-10:
            break
    return best_mu


def solve_chemical_potentials(
    group_data: GroupData,
    elements: list[str],
    remaining_vec: np.ndarray,
    full_targets: list[float],
) -> np.ndarray:
    """Solve for one chemical potential per element such that the total
    content of element `elements[k]`, summed over `group_data`, equals
    ``remaining_vec[k]``. `full_targets` holds each element's full pool
    target (before subtracting fixed contributions), used only for
    diagnostics.

    The stationarity condition is solved as the per-element-scaled root
    problem ``content_X(mu) / remaining[X] - 1 = 0``, with the scaled
    Hessian as its analytic Jacobian (scipy's ``hybr``, Powell's hybrid
    method -- a damped Newton iteration), falling back to a log-residual
    stage of the same method and then to sign-based bracketing; the
    result of whichever stage produced it is verified against the
    targets.
    """
    K = len(elements)

    # feasibility: the most (least) of element X any group can supply is
    # n_free times the largest (most negative) stoichiometry among its
    # variable states.
    for k, elem in enumerate(elements):
        committed = full_targets[k] - remaining_vec[k]
        max_content = sum(
            n_free * max(0.0, s[:, k].max()) for n_free, _, s in group_data
        )
        if remaining_vec[k] > max_content:
            raise ElementPoolError(
                f"Element pool '{elem}': target {full_targets[k]:.3e} exceeds "
                f"the maximum achievable content "
                f"{committed + max_content:.3e} (committed "
                f"{committed:.3e} plus up to {max_content:.3e} from "
                "variable states)."
            )
        min_content = sum(
            n_free * min(0.0, s[:, k].min()) for n_free, _, s in group_data
        )
        if remaining_vec[k] < min_content:
            raise ElementPoolError(
                f"Element pool '{elem}': target {full_targets[k]:.3e} is "
                f"below the minimum achievable content "
                f"{committed + min_content:.3e} (committed {committed:.3e} "
                f"plus as little as {min_content:.3e} from variable "
                "states)."
            )

    def residual_and_jacobian(mu: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        content, hessian = content_and_hessian(group_data, mu)
        return content / remaining_vec - 1.0, hessian / remaining_vec[:, None]

    # Residuals measured in defects per cell span ~1e-30..1, and every
    # joint norm or step test a general-purpose solver applies to them
    # is dominated by the largest element, leaving dilute elements
    # unconverged at reported success (with trust-exact's absolute
    # gradient threshold this reached its extreme: doing nothing at all
    # counted as converged). Scaling each equation by its own target
    # keeps every residual O(1), however dilute the targets and however
    # many orders of magnitude separate them.
    #
    # The initial guess is one diagonal Newton step of the log-content
    # system ``log(content_X(mu)) = log(remaining[X])`` from mu = 0,
    # which inverts the dilute limit ``content_X ~ content_X(0) *
    # exp(sbar_X * mu_X)`` with characteristic stoichiometry
    # ``sbar_X = H_XX / content_X``. Started from mu = 0 itself, the
    # solver can sit on an underflow plateau where the Jacobian
    # vanishes.
    if (remaining_vec <= 0.0).any():
        # A zero net-content target cannot scale its own residual, and a
        # negative one has no log-space seed; the sign-based fallback
        # needs neither, so solve there directly.
        mu = bracketed_coordinate_solve(group_data, remaining_vec)
    else:
        c0, h0 = content_and_hessian(group_data, np.zeros(K))
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            x0 = np.log(remaining_vec / c0) * c0 / np.diag(h0)
        # 700 < log(largest double), so exp(x0) stays finite.
        x0 = np.clip(
            np.nan_to_num(x0, nan=700.0, posinf=700.0, neginf=-700.0),
            -700.0,
            700.0,
        )
        result = root(residual_and_jacobian, x0=x0, jac=True, method="hybr")
        mu = _sanitised_mu(result.x)
        deviation = scaled_deviation(group_data, mu, remaining_vec)

        if not (deviation.max() <= _element_pool_tolerance):
            # Second root-find stage, on log residuals: ln(content) is
            # nearly linear in mu throughout the dilute regime, and the
            # full Jacobian captures shared-species coupling that
            # defeats both the diagonal seed and coordinate sweeps (one
            # species carrying almost all of two pools' content). Run
            # second because the seeded scaled stage settles the common
            # regimes; ln(content) flattens at saturation, where neither
            # root stage is reliable and the bracketing sweep takes
            # over.
            def log_residual_and_jacobian(
                m: np.ndarray,
            ) -> tuple[np.ndarray, np.ndarray]:
                content, hessian = content_and_hessian(group_data, m)
                with np.errstate(divide="ignore", invalid="ignore"):
                    g = np.log(content / remaining_vec)
                    jacobian = hessian / content[:, None]
                return (
                    np.nan_to_num(g, nan=-700.0, posinf=700.0, neginf=-700.0),
                    np.nan_to_num(jacobian, nan=0.0, posinf=0.0, neginf=0.0),
                )

            log_result = root(
                log_residual_and_jacobian, x0=np.zeros(K), jac=True, method="hybr"
            )
            candidate = _sanitised_mu(log_result.x)
            candidate_deviation = scaled_deviation(
                group_data, candidate, remaining_vec
            )
            if candidate_deviation.max() < deviation.max():
                mu, deviation = candidate, candidate_deviation

        if not (deviation.max() <= _element_pool_tolerance):
            # The coordinate sweep needs no derivative -- each element
            # is solved by sign-based bracketing -- so plateaux where
            # the Jacobian underflows (saturated states, extreme
            # dilution) cannot stall it.
            mu = bracketed_coordinate_solve(group_data, remaining_vec)

    deviation = scaled_deviation(group_data, mu, remaining_vec)

    # The targets are verified independently of any solver's own
    # verdict; concentrations that silently miss a constraint must
    # never be returned. The comparison fails closed on NaN.
    if not (deviation.max() <= _element_pool_tolerance):
        worst = int(np.argmax(np.nan_to_num(deviation, nan=np.inf)))
        achieved, _ = content_and_hessian(group_data, mu)
        committed = full_targets[worst] - remaining_vec[worst]
        targets = ", ".join(
            f"{e}={r:.3e}" for e, r in zip(elements, remaining_vec, strict=True)
        )
        raise ElementPoolError(
            f"Element pool targets ({targets}) could not be satisfied: "
            f"element '{elements[worst]}' achieved "
            f"{committed + achieved[worst]:.6e} against target "
            f"{full_targets[worst]:.6e} (relative deviation "
            f"{deviation[worst]:.2e}). The targets may be jointly "
            "infeasible, or the solve may have failed to converge."
        )
    return mu
