# 5DoF_IK
An analytical solution for the Inverse Kinematics of a 5-DoF robotic arm with a prismatic joint, implemented in MATLAB.
The structure of the arm is RRPRR, where R denotes a revolute joint, and P denotes a prismatic joint. The DH parameters for the robot are:

![DH parameters](https://github.com/vighv/5DoF_IK/blob/master/DHparams.PNG)

The IK function takes in a 4x4 homogeneous transformation matrix T, and returns a 1x5 vector of joint variables. The length parameters (d<sub>i<\sub>, a<sub>i<\sub>) can be adapted to any other RRPRR arm with the same revolute joint ranges.
