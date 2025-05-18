# 📡 Microwave Amplifier

Welcome to **Microwave Amplifier** — a MATLAB project that helps you compute important metrics for RF transistor amplifiers like **Delta**, **Stability Factor (K)**, and **Maximum Transducer Gain (Gt)** using S-parameters.

All the heavy lifting is wrapped in clean, modular helper functions for clarity and reusability. Designed for RF enthusiasts, students, and engineers working in the microwave domain!

---

## 🚀 Features

- 📐 Converts polar to cartesian S-parameters
- 🔍 Calculates Delta (Δ) and Stability Factor (K)
- 🎯 Determines optimal source/load reflection coefficients (ρₛ and ρₗ)
- 📈 Computes maximum gain: Gs, Gl, Go, and Gt
- 💡 Modular function-based structure for clarity

---

## 🧪 Example Setup

```matlab
% Transistor AT-41411 @ 3GHz
S11_mag = 0.529;
S11_phase = 139.153; 
S21_mag = 2.896;
S21_phase = 55.280;
S12_mag = 0.098;
S12_phase = 60.012;
S22_mag = 0.335;
S22_phase = -27.192;

% Run the script and watch the magic happen!
