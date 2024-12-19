# Geometric-Dynamic-Movement-Primitives

Here you can find the codes to play with Geometric Dynamic Motion Primitives (GDMP), featuring the transformation of the phase variable to the curvilinear abscissa of the demonstrated curve. All scripts were written in matlab code, including a Simulink file showing how to implement a Geometric DMP.

## The project
This study presents a novel parameterization approach for DMP that allows for the decoupling of the Transformation System from the timing law employed during the demonstration of a desired task. Consequently, it offers the flexibility to independently select the phase law for the Canonical System.

If you find the article useful, please cite:
```bibtex 
@misc{braglia2024phasefreedynamicmovementprimitives,
      title={Phase-free Dynamic Movement Primitives Applied to Kinesthetic Guidance in Robotic Co-manipulation Tasks}, 
      author={Giovanni Braglia and Davide Tebaldi and Luigi Biagiotti},
      year={2024},
      eprint={2401.08238},
      archivePrefix={arXiv},
      primaryClass={cs.RO},
      url={https://arxiv.org/abs/2401.08238}, 
}
``` 
## Authors
Giovanni Braglia, Davide Tebaldi and Luigi Biagiotti , from the Engineering Departement Enzo Ferrari, University of Modena and Reggio Emilia.
Check out our lab at: https://www.automatica.unimore.it/

## Important Resources

Here I leave attached the materials that helped us writing the codes.

- A. J. Ijspeert, J. Nakanishi, H. Hoffmann, P. Pastor and S. Schaal, "Dynamical Movement Primitives: Learning Attractor Models for Motor Behaviors," in Neural Computation, vol. 25, no. 2, pp. 328-373, Feb. 2013, doi: 10.1162/NECO_a_00393.
- Sidiropoulos, A., & Doulgeri, Z. (2021, May). A reversible dynamic movement primitive formulation. In 2021 IEEE International Conference on Robotics and Automation (ICRA) (pp. 3147-3153). IEEE.
- Andersson, J. A., Gillis, J., Horn, G., Rawlings, J. B., & Diehl, M. (2019). CasADi: a software framework for nonlinear optimization and optimal control. Mathematical Programming Computation, 11, 1-36.
- Biagiotti, L., & Melchiorri, C. (2008). Trajectory planning for automatic machines and robots. Springer Science & Business Media.

## Folders

- **Data/X_dim.mat** : an array of recorded Cartesian position used as an example for the codes;
- **Data/Q_dim.mat** : an array of recorded joint position used as an example for the codes; 
- **Codes/TrajParam_Workspace.m** : script that produces the symbolic functions used in the GDMP for Cartesian position;
- **Codes/TrajParam_Jointspace.m** : script that produces the symbolic functions used in the GDMP for joints position; 
- **Codes/OptimalPhase_Workspace.m** : script to optimize the phase variable of GDMP to minimize the task execution period and reach a constant feed-rate;
- **Codes/OptimalPhase_Jointspace.m** : script to optimize the phase variable of GDMP to minimize the task execution period and bound joints' velocity/acceleration;
- **Codes/SpatialSampling.m** : function that filters the input trajectory to retrieve a new one with equally distanced points;
- **Codes/(first_der.m,second_der.m)** : function to calculate analitically the first/second derivative of normalized Radial Basis Functions (RBF);
- **Codes/LagrangeApprox.m** : function to calculate analitically the weights coefficient to perform the regression of a desired signal as a sum of RBF;
- **Codes/Simulink_Parameters.m** : initializes the parameters for the Simulink files in the *Models* folder;
- **Models/(Geometric_DMP_R2020a.slx, Geometric_DMP_R2021a.slx, Geometric_DMP_R2022a.slx)** : Simulink models of a 1D GDMP block.

## Questions & Suggestions
All codes were tested on MATLAB 2022a. Due to internal dependencies in the codes, it is important to run first *Codes/TrajParam_Workspace.m*, successively *Codes/Workspace2Jointspace.m* and then *Codes/TrajParam_Jointspace.m*.
To test Simulink blocks, the related parameters are written in *Codes/Simulink_Paramters.m*.
For any doubt, question or suggestion, please feel free to email at:
giovanni.braglia@unimore.it
