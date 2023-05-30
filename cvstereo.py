#! /usr/bin/env python3

import sys
# Thanks to : https://www.tutorialspoint.com/how-to-create-a-depth-map-from-stereo-images-in-opencv-python
import cv2
img0 = cv2.imread(sys.argv[1] + ".ppm", 0)
imgL = cv2.resize(cv2.imread(sys.argv[1] + "-L.png", 0), (img0.shape[1], img0.shape[0]))
imgR = cv2.resize(cv2.imread(sys.argv[1] + "-R.png", 0), (img0.shape[1], img0.shape[0]))
stereo = cv2.StereoBM_create(numDisparities=16, blockSize=15)
depth = stereo.compute(imgL, imgR)
cv2.imwrite(sys.argv[1] + "-bumps.png", cv2.resize(depth, (img0.shape[1], img0.shape[0])))

