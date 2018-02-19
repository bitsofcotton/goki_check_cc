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
    subprocess.call(["convert", root + "-bump" + ext, "-blur", "8x8", "-compress", "none", root + "-bump-blur.ppm"])
  if(not os.path.exists(root + ".obj")):
    subprocess.call([argv[1], "obj", root + "-bump-blur.ppm", root + ".obj"])
  if(not os.path.exists(root + "-emph" + ext)):
    subprocess.call(["convert", root + "-bump" + ext, "-blur", "2x2", "-alpha", "on", "-channel", "a", "-evaluate", "set", "30%", root + "-bump-test.png"])
    subprocess.call(["convert", line, root + "-bump-test.png", "-compose", "Multiply", "-composite", root + "-emph" + ext])
  if(not os.path.exists(os.path.dirname(root) + "/material.mtl")):
    subprocess.call(["cp", os.path.dirname(argv[1]) + "/material.mtl", root + ".obj.mtl"])
    f = open(root + ".obj.mtl", "a")
    f.write("map_Ka " + line + "\n")
    f.write("map_Kd " + line + "\n")
    f.close()
  if(not os.path.exists(root + "-collect" + ext)):
    subprocess.call([argv[1], "collect", line, root + "-collect" + ext])
  if(not os.path.exists(root + "-enl-last" + ext)):
    subprocess.call([argv[1], "enlarge", line, root + "-enl" + ext])
    subprocess.call(["convert", line, "-matte", "-virtual-pixel", "Mirror", "-filter", "point", "+distort", "SRT", "45", "-compress", "none", root + "-d" + ext])
    subprocess.call([argv[1], "enlarge", root + "-d" + ext, root + "-enl-d" + ext])
    subprocess.call(["convert", root + "-enl-d" + ext, "-rotate", "-45", root + "-enl-d-margin" + ext])
    p = subprocess.Popen(["identify", "-format", "%wx%h", root + "-enl" + ext], stdout=subprocess.PIPE)
    szorig = p.communicate("\n")[0] + "+0+0"
    subprocess.call(["convert", root + "-enl-d-margin" + ext, "-gravity", "center", "-crop", szorig, root + "-enl-d-crop" + ext])
    #subprocess.call(["convert", root + "-enl" + ext, root + "-enl-d-crop" + ext, "-evaluate-sequence", "mean", "-resize", "75%", root + "-enl-last" + ext])
    subprocess.call(["convert", root + "-enl" + ext, root + "-enl-d-crop" + ext, "-evaluate-sequence", "mean", root + "-enl-last" + ext])

