import unittest
from unittest.mock import patch

import numpy as np
from scipy.optimize import OptimizeResult

from py_sc_fermi.element_pools import (
    ElementPoolError,
    bracketed_coordinate_solve,
    content_and_hessian,
    group_concs,
    scaled_deviation,
    solve_along_direction,
    solve_chemical_potentials,
)


def make_group(
    n_free: float,
    log_weights: list[float],
    stoich_rows: list[list[float]],
) -> tuple[float, np.ndarray, np.ndarray]:
    """A single exclusion group in the form the solver consumes:
    ``(n_free, log Boltzmann weights of shape (n,), stoichiometry rows of
    shape (n, K))``. Building these directly keeps the solver tests
    independent of `DefectSystem` and of any particular physics."""
    return (
        float(n_free),
        np.array(log_weights, dtype=float),
        np.array(stoich_rows, dtype=float),
    )


# A in pools X and Y, B in X only, C in Y only; weights 1 (log 0) so the
# sites sit at O(1) occupancy. Shared by the Jacobian and coupled-solve
# tests as a minimal two-element coupled system.
COUPLED_XY = [
    make_group(20.0, [0.0], [[1.0, 1.0]]),
    make_group(20.0, [0.0], [[1.0, 0.0]]),
    make_group(20.0, [0.0], [[0.0, 1.0]]),
]


class TestGroupConcs(unittest.TestCase):
    def test_matches_langmuir_formula(self):
        n_free, log_w, s = make_group(2.0, [0.3, -0.5], [[1.0], [2.0]])
        mu = np.array([0.4])
        exponents = log_w + s @ mu
        expected = n_free * np.exp(exponents) / (1.0 + np.exp(exponents).sum())
        np.testing.assert_allclose(group_concs(n_free, log_w, s, mu), expected)

    def test_dominant_weight_saturates_all_free_sites(self):
        c = group_concs(*make_group(5.0, [50.0], [[1.0]]), np.array([0.0]))
        np.testing.assert_allclose(c, [5.0])


class TestContentAndHessian(unittest.TestCase):
    def test_jacobian_matches_finite_differences_at_high_occupancy(self):
        """The analytic Jacobian of the element-content map must agree with
        finite differences at O(1) site occupancy, where the empty state's
        contribution to the per-site covariance is not negligible."""
        mu = np.array([-0.7, -1.6])
        _, jacobian = content_and_hessian(COUPLED_XY, mu)
        eps = 1e-6
        fd = np.zeros((2, 2))
        for k in range(2):
            d = np.zeros(2)
            d[k] = eps
            up, _ = content_and_hessian(COUPLED_XY, mu + d)
            down, _ = content_and_hessian(COUPLED_XY, mu - d)
            fd[:, k] = (up - down) / (2 * eps)
        np.testing.assert_allclose(jacobian, fd, rtol=1e-6)


class TestScaledDeviation(unittest.TestCase):
    def test_relative_to_nonzero_target(self):
        group_data = [make_group(2.0, [0.0], [[1.0]])]
        # content = 2 * 1 / (1 + 1) = 1.0 against a target of 2.0
        dev = scaled_deviation(group_data, np.array([0.0]), np.array([2.0]))
        np.testing.assert_allclose(dev, [0.5])

    def test_zero_target_scaled_by_gross_content(self):
        # one species adds the element (+1), one removes it (-1): the
        # net-content deviation is measured against the gross content
        group_data = [
            make_group(2.0, [0.0], [[1.0]]),
            make_group(2.0, [0.0], [[-1.0]]),
        ]
        mu = np.array([0.1])
        c_add = group_concs(*group_data[0], mu)[0]
        c_remove = group_concs(*group_data[1], mu)[0]
        expected = abs(c_add - c_remove) / (c_add + c_remove)
        dev = scaled_deviation(group_data, mu, np.array([0.0]))
        np.testing.assert_allclose(dev, [expected])

    def test_zero_target_with_no_contributing_state_reads_as_satisfied(self):
        # element Y (column 1) has no variable contribution anywhere, so a
        # zero target on it is met with gross content zero, not NaN
        group_data = [make_group(2.0, [0.0], [[1.0, 0.0]])]
        dev = scaled_deviation(
            group_data, np.array([0.0, 0.0]), np.array([0.5, 0.0])
        )
        self.assertEqual(dev[1], 0.0)


class TestSolveChemicalPotentials(unittest.TestCase):
    def test_single_element_hits_target(self):
        group_data = [make_group(1.0, [0.0], [[1.0]])]
        mu = solve_chemical_potentials(group_data, ["X"], np.array([0.3]), [0.3])
        content, _ = content_and_hessian(group_data, mu)
        np.testing.assert_allclose(content, [0.3], rtol=1e-10)

    def test_coupled_elements_hit_targets(self):
        mu = solve_chemical_potentials(
            COUPLED_XY, ["X", "Y"], np.array([8.0, 5.0]), [8.0, 5.0]
        )
        content, _ = content_and_hessian(COUPLED_XY, mu)
        np.testing.assert_allclose(content, [8.0, 5.0], rtol=1e-8)

    def test_target_beyond_capacity_raises(self):
        group_data = [make_group(1.0, [0.0], [[1.0]])]
        with self.assertRaises(ElementPoolError) as ctx:
            solve_chemical_potentials(group_data, ["X"], np.array([5.0]), [5.0])
        self.assertIn("X", str(ctx.exception))

    def test_guard_raises_when_solvers_do_not_converge(self):
        """The targets are verified independently of the solvers' verdicts:
        if both the root-find and the bracketing fallback return without
        reaching the target, an ElementPoolError names the unmet element
        rather than returning unconstrained concentrations."""
        group_data = [make_group(1.0, [np.log(2e-17)], [[1.0]])]
        def do_nothing(fun, x0, **kwargs):
            return OptimizeResult(x=np.zeros_like(x0), success=True, message="mocked")
        with (
            patch("py_sc_fermi.element_pools.root", do_nothing),
            patch(
                "py_sc_fermi.element_pools.bracketed_coordinate_solve",
                lambda group_data, remaining_vec: np.zeros(len(remaining_vec)),
            ),
        ):
            with self.assertRaises(ElementPoolError) as ctx:
                solve_chemical_potentials(
                    group_data, ["X"], np.array([4e-17]), [4e-17]
                )
        self.assertIn("'X'", str(ctx.exception))


class TestSolveAlongDirection(unittest.TestCase):
    def test_line_solve_hits_projected_target(self):
        group_data = [make_group(1.0, [0.0], [[1.0]])]
        mu = solve_along_direction(group_data, np.array([0.0]), np.array([1.0]), 0.3)
        content, _ = content_and_hessian(group_data, mu)
        np.testing.assert_allclose(content, [0.3], rtol=1e-10)


class TestBracketedCoordinateSolve(unittest.TestCase):
    """The derivative-free fallback in isolation: the entry point for
    zero/negative net-content targets and the recovery path when the
    Newton stages stall."""

    def test_solves_coupled_positive_targets(self):
        mu = bracketed_coordinate_solve(COUPLED_XY, np.array([8.0, 5.0]))
        content, _ = content_and_hessian(COUPLED_XY, mu)
        np.testing.assert_allclose(content, [8.0, 5.0], rtol=1e-8)

    def test_resolves_small_net_target_at_half_saturation(self):
        # O_i (+1) and V_O (-1) both half-occupied at mu = 0; a net target
        # of -1e-7 is a cancellation of two O(5) populations that the line
        # solve's Newton polish must resolve below the bisection floor
        group_data = [
            make_group(10.0, [0.0], [[1.0]]),
            make_group(10.0, [0.0], [[-1.0]]),
        ]
        mu = bracketed_coordinate_solve(group_data, np.array([-1e-7]))
        content, _ = content_and_hessian(group_data, mu)
        np.testing.assert_allclose(content, [-1e-7], rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
