#! /usr/bin/env python3

import sys
# Thanks to : https://www.tutorialspoint.com/how-to-create-a-depth-map-from-stereo-images-in-opencv-python
import cv2
imgL = cv2.imread(sys.argv[1] + "-L.png", 0)
imgR = cv2.imread(sys.argv[1] + "-R.png", 0)
stereo = cv2.StereoBM_create(numDisparities=16, blockSize=15)
depth = stereo.compute(imgL, imgR)
cv2.imwrite(sys.argv[1] + "-bumps.png", depth)

