# 2020_ICRA_helix_RFT_swirl

This repository contains simulation codes related to the publication:

H. O. Caldag, S. Yesilyurt, Steering Control of Magnetic Helical Swimmers
in Swirling Flows due to Confinement, 2020 International Conference on
Robotics and Automation. https://doi.org/10.1109/ICRA40945.2020.9196521

Specifically, the repository contains the resistive force theory-based model described in this paper. It contains three pieces of Matlab scripts.

The model is for controlling the trajectories of magnetically actuated slender helices. These helices, when confined, exhibit oscillatory trajectories. The confinement effects here are represented with a swirling flow. Swimmer is controlled by tilting it up and down depending on its position relative to the desired position.

The model computes the resistance matrix for a slender helix and simulates its swimming in bulk fluid and applies control input depending on its instantaneous position.

The differential equations for the swimming kinematics are resolved using ode45 package of Matlab.

The package contains:

- helix_control_with_swirl.m : This is the main script for running the simulation of the helix under swirling flow.
- get_mobility_icra.m : A function that evaluates the components of the mobility matrix of a slender helix. The expressions are simplified versions of those presented in Man & Lauga (2013). The expressions used here are available in our publication as well.
- swirl_compute.m: This function computes the forces acting on the helix due to the swirling flow.

Refer to the codes for more details.

If you have any questions, write to hakanosmanc@gmail.com

REFERENCES

Y. Man, E. Lauga (2013). The wobbling-to-swimming transition of rotated helices. Physics of Fluids, https://doi.org/10.1063/1.4812637.
