import Mathlib.Data.Real.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.BigOperators.Ring
import Mathlib.Tactic.Ring

open BigOperators

/-!
# Network Synchronization and Spectral Sensitivity in Lean 4

This file formalizes mathematical foundations of synchronization and spectral 
outage sensitivity in coupled oscillator networks (Kuramoto models) representing 
high-voltage power transmission grids (e.g., the BSPTCL 24-bus network).

## Methodological Scope & Formal Verification Boundary
- **Formally Verified in Lean 4 (Machine-Checked, 0 custom axioms):**
  1. `two_le_of_ne`: Distinct bus indices `i ≠ j` structurally guarantee network size `2 ≤ n`,
     formally ruling out empty or trivial single-bus systems.
  2. `phase_sync_phase_difference_zero`: Complete phase synchrony (`θ i = θ j`) forces all
     pairwise phase differences to vanish identically, eliminating inter-bus power transfer.
  3. `power_damping_ratio_determines_frequency_sync`: In steady-state swing dynamics with
     vanishing coupling, uniform power-to-damping ratios (`P i / D i = P j / D j`) determine
     complete frequency synchronization (`ω i = ω j`), corresponding to the zero-discrepancy
     boundary `γ_c = 0` of the Dörfler–Bullo condition.
  4. `quadraticForm_eq_dotProduct`: Canonical bridge establishing that `quadraticForm M v`
     coincides exactly with Mathlib's standard `dotProduct v (Matrix.mulVec M v)`.
  5. `quadraticForm_rank_one`: The algebraic identity that for any incidence-like vector `e`,
     the quadratic form of the rank-1 dyadic matrix `(-K * e * eᵀ)` equals `-K * (vᵀ e)²`.
  6. `sum_mul_incidenceVector`: The inner product of any vector `v` with a line incidence
     vector `incidenceVector i j` evaluates exactly to `v i - v j`.
  7. `fiedler_outage_sensitivity`: The formal algebraic evaluation of the Rayleigh quotient
     perturbation kernel for line tripping: `quadraticForm (-K e_{ij} e_{ij}ᵀ) v = -K (v i - v j)²`.

- **Standard Matrix Analysis Boundary (Theorem 4 in Paper 2):**
  The first-order perturbation approximation `Δλ₂ ≈ v₂ᵀ ΔL v₂` for a simple eigenvalue
  under a symmetric matrix perturbation `ΔL = -K_{ij} e_{ij} e_{ij}ᵀ` is a standard classical
  result in matrix perturbation theory (e.g., Stewart & Sun, 1990). The algebraic evaluation
  of the resulting quadratic form is what is formally verified here without approximation.
-/

namespace DynamicalSystems

/-- Represent the phase angle state of a network of `n` coupled oscillators. -/
def PhaseVector (n : ℕ) := Fin n → ℝ

/-- Represent the frequency deviation (velocity) state of the network. -/
def FrequencyVector (n : ℕ) := Fin n → ℝ

/-- A coupling matrix representing transmission line susceptances or Laplacian perturbations. -/
def CouplingMatrix (n : ℕ) := Matrix (Fin n) (Fin n) ℝ

/-- Complete frequency synchronization (frequency locking) occurs when all frequency 
    deviations from the grid reference frequency are identical across all buses. -/
def IsFrequencySynchronized {n : ℕ} (ω : FrequencyVector n) : Prop :=
  ∀ i j : Fin n, ω i = ω j

/-- Complete phase synchronization (phase locking) occurs when all phase angles 
    across the grid nodes are identical. -/
def IsPhaseSynchronized {n : ℕ} (θ : PhaseVector n) : Prop :=
  ∀ i j : Fin n, θ i = θ j

/-- For any network containing distinct buses `i ≠ j : Fin n`, the network size
    is necessarily at least two (`2 ≤ n`). This formally excludes degenerate empty or
    single-bus networks from line outage analysis. -/
lemma two_le_of_ne {n : ℕ} {i j : Fin n} (h_ne : i ≠ j) : 2 ≤ n := by
  by_contra h
  have : i = j := by
    ext
    omega
  exact h_ne this

/-- When all phase angles are synchronized (`θ i = θ j` for all `i, j`),
    the pairwise phase difference between any pair of buses vanishes identically,
    eliminating all inter-bus sinusoidal coupling power transfer. -/
theorem phase_sync_phase_difference_zero {n : ℕ} (θ : PhaseVector n)
    (hPhase : IsPhaseSynchronized θ) (i j : Fin n) :
    θ i - θ j = 0 := by
  rw [hPhase i j, sub_self]

/-- In steady-state swing dynamics with balanced coupling, frequency deviations satisfy
    `P i - D i * ω i = 0`. Complete frequency synchronization (`ω i = ω j` for all `i, j`)
    holds when power-to-damping ratios are uniform across all grid buses (`P i / D i = P j / D j`),
    corresponding to the zero-discrepancy boundary `γ_c = 0` of the Dörfler-Bullo condition. -/
theorem power_damping_ratio_determines_frequency_sync {n : ℕ}
    (P D : Fin n → ℝ) (ω : FrequencyVector n)
    (hD : ∀ i, D i ≠ 0)
    (hEquil : ∀ i, P i - D i * ω i = 0)
    (hRatio : ∀ i j, P i / D i = P j / D j) :
    IsFrequencySynchronized ω := by
  intro i j
  have h_freq (k : Fin n) : ω k = P k / D k := by
    have hP : P k = D k * ω k := by
      have h := hEquil k
      exact sub_eq_zero.mp h
    rw [hP, mul_comm, mul_div_cancel_right₀ (ω k) (hD k)]
  rw [h_freq i, h_freq j]
  exact hRatio i j

/-- The incidence vector for a transmission line connecting bus `i` and bus `j`. -/
def incidenceVector {n : ℕ} (i j : Fin n) : Fin n → ℝ :=
  fun k => if k = i then 1 else if k = j then -1 else 0

/-- The quadratic form of a coupling matrix `M` with respect to a perturbation vector `v`. -/
def quadraticForm {n : ℕ} (M : CouplingMatrix n) (v : Fin n → ℝ) : ℝ :=
  ∑ r : Fin n, ∑ s : Fin n, v r * M r s * v s

/-- Bridge lemma: The standalone `quadraticForm` definition coincides with Mathlib's
    canonical matrix operations `dotProduct v (Matrix.mulVec M v)`. -/
lemma quadraticForm_eq_dotProduct {n : ℕ} (M : CouplingMatrix n) (v : Fin n → ℝ) :
    quadraticForm M v = dotProduct v (Matrix.mulVec M v) := by
  unfold quadraticForm dotProduct Matrix.mulVec dotProduct
  dsimp only
  apply Finset.sum_congr rfl
  intro i _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- Theorem: For any real vector `v` and any vector `e` representing line incidence,
    the quadratic form of the rank-1 perturbation matrix `(-K * e * eᵀ)` is exactly `-K * (vᵀ e)²`. -/
theorem quadraticForm_rank_one {n : ℕ} (K : ℝ) (e : Fin n → ℝ) (v : Fin n → ℝ) :
  quadraticForm (fun r s => -K * e r * e s) v = -K * (∑ k : Fin n, v k * e k) ^ 2 := by
  unfold quadraticForm
  have h1 : ∀ r s, v r * (-K * e r * e s) * v s = -K * (v r * e r) * (v s * e s) := by
    intro r s; ring
  simp_rw [h1]
  have h2 : ∀ r, (∑ s, -K * (v r * e r) * (v s * e s)) = -K * (v r * e r) * (∑ s, v s * e s) := by
    intro r
    rw [← Finset.mul_sum]
  simp_rw [h2]
  have h3 : (∑ r, -K * (v r * e r) * (∑ s, v s * e s)) = -K * (∑ s, v s * e s) * (∑ r, v r * e r) := by
    have h4 : ∀ r, -K * (v r * e r) * (∑ s, v s * e s) = (-K * (∑ s, v s * e s)) * (v r * e r) := by
      intro r; ring
    simp_rw [h4]
    rw [← Finset.mul_sum]
  rw [h3]
  ring

/-- Helper lemma: The inner product of any vector `v` with the incidence vector `incidenceVector i j`
    for distinct nodes `i ≠ j` evaluates exactly to `v i - v j`. -/
lemma sum_mul_incidenceVector {n : ℕ} (i j : Fin n) (h_ne : i ≠ j) (v : Fin n → ℝ) :
    (∑ k : Fin n, v k * incidenceVector i j k) = v i - v j := by
  have hi : incidenceVector i j i = 1 := by
    simp [incidenceVector]
  have hj : incidenceVector i j j = -1 := by
    simp [incidenceVector, h_ne.symm]
  rw [← Finset.add_sum_erase Finset.univ _ (Finset.mem_univ i)]
  have hj_mem : j ∈ Finset.univ.erase i := by
    rw [Finset.mem_erase]
    exact ⟨h_ne.symm, Finset.mem_univ j⟩
  rw [← Finset.add_sum_erase (Finset.univ.erase i) _ hj_mem]
  have h_zero : ∀ k ∈ (Finset.univ.erase i).erase j, v k * incidenceVector i j k = 0 := by
    intro k hk
    rw [Finset.mem_erase, Finset.mem_erase] at hk
    have hki : k ≠ i := hk.2.1
    have hkj : k ≠ j := hk.1
    simp [incidenceVector, hki, hkj]
  have h_sum_zero : ∑ k ∈ (Finset.univ.erase i).erase j, v k * incidenceVector i j k = 0 :=
    Finset.sum_eq_zero h_zero
  rw [h_sum_zero, add_zero, hi, hj]
  ring

/-- Theorem stating that for distinct nodes `i` and `j`, the projection of the Fiedler 
    eigenvector `v` onto the rank-1 outage perturbation matrix equals `-K * (v i - v j)²`. -/
theorem fiedler_outage_sensitivity {n : ℕ} (i j : Fin n) (h_ne : i ≠ j) (K : ℝ) (v : Fin n → ℝ) :
    quadraticForm (fun r s => -K * (incidenceVector i j r) * (incidenceVector i j s)) v = -K * (v i - v j) ^ 2 := by
  rw [quadraticForm_rank_one]
  rw [sum_mul_incidenceVector i j h_ne v]

end DynamicalSystems
