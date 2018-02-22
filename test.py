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
    subprocess.call([argv[1], "bump", line, root + "-bump" + ext])
    subprocess.call(["convert", root + "-bump" + ext, "-blur", "16x16", "-compress", "none", root + "-bump-blur.ppm"])
    subprocess.call(["convert", line, root + "-bump-blur.ppm", "-alpha", "on", "-channel", "a", "-evaluate", "set", "30%", "-compose", "Multiply", "-composite", root + "-emph" + ext])
  if(not os.path.exists(root + ".obj")):
    if(not os.path.exists(root + "-mask.ppm")):
      subprocess.call(["convert", line, "-fill", "black", "-colorize", "100", "-compress", "none", root + "-mask.ppm"])
    subprocess.call([argv[1], "obj", root + "-bump-blur.ppm", root + "0.obj"])
    subprocess.call([argv[1], "maskobj", root + "-mask.ppm", root + "0.obj", root + ".obj", ".05", "2"])
  if(not os.path.exists(os.path.dirname(root) + "/material.mtl")):
    subprocess.call(["cp", os.path.dirname(argv[1]) + "/material.mtl", root + ".obj.mtl"])
    f = open(root + ".obj.mtl", "a")
    f.write("map_Ka " + line + "\n")
    f.write("map_Kd " + line + "\n\n")
    f.close()
  if(not os.path.exists(root + "-enl-75" + ext)):
    subprocess.call([argv[1], "collect", line, root + "-collect" + ext])
    subprocess.call([argv[1], "enlarge", line, root + "-enl" + ext])
    subprocess.call(["convert", root + "-enl" + ext, "-resize", "75%", root + "-enl-75" + ext])

