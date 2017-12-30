# system_approximation
heuns.m is a program written in MatLab that's purpose is to numerically approximate the solutions to a system of first order homogenous differential equations of size N by N ( N being a variable number ). The equations are approximated by Heun's method. This functionality allows the program to also approximate the solution to a differential equation of N order if it is first broken down into its N by N matrix.

Heun's Method explained: http://calculuslab.deltacollege.edu/ODE/7-C-2/7-C-2-h.html

error_analysis.m is another MatLab program that illustrates how Heun's method is more accurate than Euler's method in approximating a system of two first order homogenous differential equations by computing the error of each value and comparing the two approximations accuracy.

Further development will allow the program to approximate nonhomogenous systems and nonlinear systems.
