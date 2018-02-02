#! /usr/bin/env python

import os
import sys
import subprocess

argv = sys.argv

if(len(argv) < 3):
  print "not much arguments."
  exit(- 1)

for idx in range(2, len(argv)):
  line = argv[idx]
  print line
  root, ext = os.path.splitext(line)
  if(ext != ".ppm"):
    subprocess.call(["convert", line, "-compress", "none", line + ".ppm"])
    line = line + ".ppm"
    root, ext = os.path.splitext(line)
  if(not os.path.exists(root + "-bump" + ext)):
    subprocess.call([argv[1], "bump2", line, root + "-bump" + ext])
    subprocess.call(["convert", "-resize", "128x", "-resize", "x128<", "-resize", "50%", "-gravity", "center", "-crop", "64x64+0+0", root + "-bump" + ext, root + "-bump-small" + ".jpeg"])

