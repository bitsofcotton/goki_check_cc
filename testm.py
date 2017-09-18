#! /usr/bin/env python

import os
import sys
import glob
import subprocess

argv = sys.argv

if(len(argv) < 2):
  print "not much arguments."
  exit(- 1)

root0, ext0 = os.path.splitext(argv[2])
if(not os.path.exists(root0 + "-bump" + ext0)):
  subprocess.call([argv[1], "bump", argv[2], root0 + "-bump" + ext0])
if(not os.path.exists(root0 + "-bump-blur" + ext0)):
  subprocess.call(["convert", root0 + "-bump" + ext0, "-blur", "4x4", "-compress", "none", root0 + "-bump-blur" + ext0])
root1, ext1 = os.path.splitext(argv[3])
if(not os.path.exists(root1 + "-bump" + ext1)):
  subprocess.call([argv[1], "bump", argv[3], root1 + "-bump" + ext1])
if(not os.path.exists(root1 + "-bump-blur" + ext1)):
  subprocess.call(["convert", root1 + "-bump" + ext1, "-blur", "4x4", "-compress", "none", root1 + "-bump-blur" + ext1])
subprocess.call([argv[1], "match", root0 + "-bump-blur" + ext0, "match" + root0 + "-" + root1, root1 + "-bump-blur" + ext1, argv[3], argv[2]])
subprocess.call([argv[1], "match3d", root0 + "-bump-blur" + ext0, "match3d" + root0, argv[2], argv[4]])
subprocess.call([argv[1], "match3d", root1 + "-bump-blur" + ext1, "match3d" + root1, argv[3], argv[4]])

