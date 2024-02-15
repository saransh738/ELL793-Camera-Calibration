# Camera-Calibration
## Problem Statement:
To find intrinsic and extrinsic camera calibration parameters of mobile phone's camera. 
* Intrinsic calibration parameters include focal length, skew, radial distortion parameters, other distortion parameters, and the camera's optic center.
* Extrinsic calibration parameters include Rotation and translation (scale) computation for every photo taken by the mobile phone.
  
We need to calibrate our mobile phone camera using a checkerboard object. 
Print two (or three) copies of the checkerboard image and stick them to two (or three) orthogonal planes (wall corner). 
Click a picture of the checkerboard from your phone. 
Camera calibration requires 2D and 3D correspondences. 
Create a dataset that contains XYZ coordinates of N points marked out on the wall checkerboard and also the XY coordinates of the corresponding points on the image. 
Now, write a program that estimates the 3 × 4 projection matrix P and then decompose it into the intrinsics and extrinsics.

## Approach:

* Normalize the data such that the centroid of 2D and 3D points are at the origin and the average Euclidean distance of 2D and 3D points from the origin is √2 and √3, respectively.
* Find the transformation matrices T and U that achieve this for 2D and 3D respectively, i.e., ˆx = Tx and Xˆ = UX where x and X are the unnormalized 2D and 3D points in homogeneous coordinates.
* Estimate the normalized projection matrix Pˆ using the DLT method. Denormalize the projection matrix Pˆ. (P = T−1PˆU).
* Decompose the projection matrix P = K[R| −RXo] into intrinsic matrix K, rotation matrix R, and the camera center Xo. K and R can be estimated using RQ decomposition.
* Verify that the projection matrix is correctly estimated by computing the RMSE between the 2D points marked by you and the estimated 2D projections of the marked 3D points. Visualize the points on the image and include them in the report. Also, mention why it is a good idea to normalize the points before performing DLT.
