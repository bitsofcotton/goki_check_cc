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
root1, ext1 = os.path.splitext(argv[3])
subprocess.call([argv[1], "match", argv[2], "match" + root0 + "-" + root1, argv[3]])
subprocess.call([argv[1], "match3d", argv[2], "match3d" + root0, argv[4]])
subprocess.call([argv[1], "match3d", argv[3], "match3d" + root1, argv[4]])

