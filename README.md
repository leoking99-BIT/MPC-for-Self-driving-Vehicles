# 无人驾驶车辆与模型预测控制(第2版)

无人驾驶车辆模型预测控制(第2版)随书仿真程序及扩展学习资料.

This repo holds the Simulink/CarSim codes for examples of Self-driving Vehicles and Model Predictive Contorl (2ed edition).


## Contents

1. Chapter-2: Vehicle model validation, including 4 examples
    * tire model validation
    * dynamic model derivation
    * kinematic model validation
    * dynamic model validation

2. Chapter-3: MPC for longitudinal control, including 4 examples
    * a simple example using MPC toolbox for speed tracking demo
    * MPC for speed tracking control with du as input
    * Extended: MPC for speed tracking control with u as input
    * Extended: MPC for Adaptive Cruise Control, involve multi-vehicle joint simulation

3. Chapter-4: MPC use kinematic model for path tracking, including 4 examples
    * Matlab code for given path tracking
    * Simulink/CarSim for given path tracking
    * Adaptive path fitting with 3rd-order Bezier curve
    * MPC for general path tracking

4. Chapter-5: MPC use dynamic model for Active Frontwheel Steering (AFS), including 2 examples
    * Tire Cornering Stiffness Estimation use RLS method
    * MPC for AFS 

5. Chapter-6: MPC for path tracking with local planning, including 1 examples
    * MPC for obstacle avoidance in local planning, then MPC for tracking

6. Chapter-7: MPC for handlig stability of high-speed self-driving Vehicles, including 3 examples
    * Different qp-solver compare test
    * MPC for high-speed self-driving vehicles, ignore handling stability
    * MPC for high-speed self-driving vehicles consider handling stability

7. Chapter-8:MPC for rollover prevention consider road terrain information, including 2 examples
    * Zero-Moment-Point based Rollover criteria validation
    * MPC for path tracking considering road terrain and rollover prevention

## Relatives links
IVRC: [Link for Inteligent Vehicle Research Center (IVRC) of BIT](https://github.com/bit-ivrc)

qpOASES: [Link for qpOASES](https://github.com/leoking99-BIT/qpOASES)

OSQP: [Link for OSQP](https://github.com/leoking99-BIT/osqp)

## Bug reports and support
Please report any issues via the [Github issue tracker](https://github.com/leoking99-BIT/Self-driving-Vehicles-and-Model-Predictive-Control/issues). All types of issues are welcome, including bug reports, documentation typos, feature requests and so on.

## Contact Me
Personal homepage: [https://leoking99-bit.github.io/](https://leoking99-bit.github.io/)

Email: leoking1025@bit.edu.cn

If you have any questions or suggestions, please don't hesitate to contact me. Any advices and comments would be highly appreciated. 