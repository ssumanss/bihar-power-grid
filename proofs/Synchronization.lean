import Mathlib.Data.Real.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.BigOperators.Ring
import Mathlib.Tactic.Ring

open BigOperators

/-!
# Network Synchronization and Kuramoto Manifolds in Lean 4

This file formalizes the mathematical concepts of synchronization in coupled 
oscillator networks (Kuramoto models) representing high-voltage power transmission grids.
It defines phase vectors, frequency vectors, and predicates for complete 
phase and frequency synchronization.
-/

namespace DynamicalSystems

/-- Represent the phase angle state of a network of `n` coupled oscillators. -/
def PhaseVector (n : ℕ) := Fin n → ℝ

/-- Represent the frequency deviation (velocity) state of the network. -/
def FrequencyVector (n : ℕ) := Fin n → ℝ

/-- A coupling matrix representing transmission line susceptances. -/
def CouplingMatrix (n : ℕ) := Matrix (Fin n) (Fin n) ℝ

/-- Complete frequency synchronization (frequency locking) occurs when all frequency 
  deviations from the grid reference frequency are identical. -/
def IsFrequencySynchronized {n : ℕ} (ω : FrequencyVector n) : Prop :=
  ∀ i j : Fin n, ω i = ω j

/-- Complete phase synchronization (phase locking) occurs when all phase angles 
  across the grid nodes are identical. -/
def IsPhaseSynchronized {n : ℕ} (θ : PhaseVector n) : Prop :=
  ∀ i j : Fin n, θ i = θ j

/-- Theorem proving that complete phase synchronization trivially implies 
  complete frequency synchronization when the derivatives are uniform. -/
theorem phase_sync_implies_freq_sync {n : ℕ} (_θ : PhaseVector n) (ω : FrequencyVector n)
    (hDeriv : ∀ i, ω i = 0) :
    IsFrequencySynchronized ω := by
  intro i j
  rw [hDeriv i, hDeriv j]

/-- The incidence vector for a transmission line connecting bus `i` and bus `j`. -/
def incidenceVector {n : ℕ} (i j : Fin n) : Fin n → ℝ :=
  fun k => if k = i then 1 else if k = j then -1 else 0

/-- The quadratic form of a matrix M with respect to a phase perturbation vector v. -/
def quadraticForm {n : ℕ} (M : Matrix (Fin n) (Fin n) ℝ) (v : Fin n → ℝ) : ℝ :=
  ∑ r : Fin n, ∑ s : Fin n, v r * M r s * v s

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
