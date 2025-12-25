# EFD Curve Generator Toolbox

A MATLAB-based toolbox for generating and benchmarking 2D and 3D Elliptic Fourier Descriptor (EFD) curves. This project focuses on the mathematical generation of trajectories used in multi-robot formation control.

## 🚀 Features
- **2D/3D Generation:** Generate random closed-loop shapes using harmonic summation.
- **Dynamic 3D Animator:** Real-time morphing of 3D shapes with FPS and calculation time benchmarking.
- **Sticky Points:** Implementation of specific coordinate tracking (useful for multi-robot alignment).
- **Performance Metrics:** Benchmarks "Model Fitting Time" and "Point Generation Time" to evaluate suitability for real-time controllers.

## 🛠️ Requirements
- MATLAB (Tested on R2023b and later)

## 📖 Usage
1. Open MATLAB and navigate to the project folder.
2. Run `EFD_Curve_Generator_2D.m` for a basic 2D shape.
3. Run `Dynamic_EFD_3D_Sticky_Points.m` to see real-time 3D trajectory evolution.

## 🎓 Academic Context
This toolbox serves as the foundational geometry engine for my research in multi-robot formation control, specifically implementing and extending concepts from the work of **Prof. Mustafa Ünel** regarding holonomic and non-holonomic robot trajectories.
