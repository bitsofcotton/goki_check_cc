#! /usr/bin/env python

import os
import sys
import subprocess

argv = sys.argv

if(len(argv) < 4):
  print "not much arguments."
  exit(- 1)

root0, ext0 = os.path.splitext(argv[2])
root1, ext1 = os.path.splitext(argv[3])
subprocess.call([argv[1], "bump", argv[2], root0 + "-bump.ppm"])
subprocess.call(["convert", root0 + "-bump.ppm", "-blur", "8x8", "-compress", "none", root0 + "-bump-blur.ppm"])
subprocess.call([argv[1], "bump", argv[3], root1 + "-bump.ppm"])
subprocess.call(["convert", root1 + "-bump.ppm", "-blur", "8x8", "-compress", "none", root1 + "-bump-blur.ppm"])
subprocess.call([argv[1], "match",   argv[2], "match" + root0 + "-" + root1, argv[3], root0 + "-bump-blur.ppm", root1 + "-bump-blur.ppm"])
if(4 <= len(argv)):
  subprocess.call([argv[1], "match3d", argv[2], "match3d" + root0, argv[4], root0 + "-bump-blur.ppm"])
  subprocess.call([argv[1], "match3d", argv[3], "match3d" + root1, argv[4], root1 + "-bump-blur.ppm"])

