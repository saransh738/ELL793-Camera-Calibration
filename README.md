# Camera-Calibration
## Problem Statement:
To find intrinsic and extrinsic camera calibration parameters of mobile phone's camera. Intrinsic calibration parameters include focal length, skew, radial distortion parameters, other distortion parameters, and the camera's optic center.
Extrinsic calibration parameters include Rotation and translation (scale) computation for every photo taken by the mobile phone.
We need to calibrate our mobile phone camera using a checkerboard object. 
Print two (or three) copies of the checkerboard image and stick them to two (or three) orthogonal planes (wall corner). 
Click a picture of the checkerboard from your phone. 
Camera calibration requires 2D and 3D correspondences. 
Create a dataset that contains XYZ coordinates of N points marked out on the wall checkerboard and also the XY coordinates of the corresponding points on the image. 
Now, write a program that estimates the 3 × 4 projection matrix P and then decompose it into the intrinsics and extrinsics.
