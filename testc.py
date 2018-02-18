#! /usr/bin/env python

import os
import sys
import subprocess

argv = sys.argv

if(len(argv) < 3):
  print "not much arguments."
  exit(- 1)

line = argv[2]
root, ext = os.path.splitext(line)
if(ext != ".ppm"):
  subprocess.call(["convert", line, "-compress", "none", line + ".ppm"])
  line = line + ".ppm"
  root, ext = os.path.splitext(line)

line2 = argv[3]
root2, ext2 = os.path.splitext(line2)
if(ext2 != ".ppm"):
  subprocess.call(["convert", argv[3], "-compress", "none", argv[3] + ".ppm"])
  line2 = line2 + ".ppm"

if(not os.path.exists(root + "-bump" + ext)):
  subprocess.call([argv[1], "bump", line, root + "-bump" + ext, root + ".obj"])
if(not os.path.exists(root + "-mask.obj.mtl")):
  subprocess.call(["cp", os.path.dirname(argv[1]) + "/material.mtl", root + ".obj.mtl"])
  f = open(root + ".obj.mtl", "a")
  f.write("map_Ka " + line + "\n")
  f.write("map_Kd " + line + "\n")
  f.close()
  subprocess.call([argv[1], "maskobj", line2, root + ".obj", root + "-mask.obj"])
  subprocess.call(["cp", root + ".obj.mtl", root + "-mask.obj.mtl"])

