#! /usr/bin/env python

import os
import sys
import subprocess
import numpy as np

argv = sys.argv

if(len(argv) < 4):
  print "not much argments."
  exit(- 1)

if(argv[2] == "match"):
  if(len(argv) < 5):
    print "no much argments."
    exit(- 1)
  root0, ext0 = os.path.splitext(argv[3])
  root1, ext1 = os.path.splitext(argv[4])
  if(ext0 != ".ppm"):
    subprocess.call(["convert", argv[3], "-compress", "none", root0 + ".ppm"])
  if(ext1 == ".obj"):
    subprocess.call([argv[1], "match3d", root0 + ".ppm", "match3dbase-" + root0 + "-" + root1 + "-", root1 + ".obj", root0 + "-bump.ppm", root0 + "-mask.ppm"])
  else:
    if(ext1 != ".ppm"):
      subprocess.call(["convert", argv[4], "-compress", "none", root1 + ".ppm"])
    if(len(argv) < 6):
      subprocess.call([argv[1], "match", root0 + ".ppm", "matchbase-" + root0 + "-" + root1 + "-", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm"])
    elif(os.path.splitext(argv[5])[1] == ".obj"):
      subprocess.call([argv[1], "match2dh3d", root0 + ".ppm", "match3dh2dbase-" + root0 + root1 + os.path.splitext(argv[5])[0] + "-", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm", argv[5]])
  exit(0)

pixels = 4
bhabit = ""
habit0 = ""
for line in argv[3:]:
  try:
    pixels = int(line)
    continue
  except:
    root, ext = os.path.splitext(line)
  if(ext != ".ppm"):
    subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
  if(argv[2] == "col"):
    subprocess.call([argv[1], "collect", root + ".ppm", root + "-collect.ppm"])
  elif(argv[2] == "enl"):
    subprocess.call([argv[1], "enlarge", root + ".ppm", root + "-enl.ppm"])
  elif(argv[2] == "bump"):
    subprocess.call([argv[1], "bump", root + ".ppm", root + "-bump.ppm"])
  elif(argv[2] == "emph"):
    subprocess.call(["convert", root + ".ppm", root + "-bump.ppm", "-alpha", "on", "-channel", "a", "-evaluate", "set", "30%", "-compose", "Multiply", "-composite", root + "-emph" + ext])
  elif(argv[2] == "mask"):
    subprocess.call(["convert", root + ".ppm", "-fill", "black", "-colorize", "100", "-compress", "none", root + "-mask.ppm"])
  elif(argv[2] == "obj"):
    subprocess.call([argv[1], "obj",     root + "-bump.ppm", root + "0.obj"])
    subprocess.call([argv[1], "maskobj", root + "-mask.ppm", root + "0.obj", root + ".obj", ".05", "2"])
  elif(argv[2] == "mtl"):
    subprocess.call(["cp", os.path.dirname(argv[1]) + "/material.mtl", root + ".obj.mtl"])
    subprocess.call(["cp", os.path.dirname(argv[1]) + "/material.mtl", root + "0.obj.mtl"])
    f = open(root + ".obj.mtl", "a")
    f.write("map_Ka " + line + "\n")
    f.write("map_Kd " + line + "\n\n")
    f.close()
    f = open(root + "0.obj.mtl", "a")
    f.write("map_Ka " + line + "\n")
    f.write("map_Kd " + line + "\n\n")
    f.close()
  elif(argv[2] == "tilt"):
    subprocess.call([argv[1], "tilt", root + ".ppm", root + "-tilt-base", root + "-bump.ppm"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-tilt-base-%d.ppm", "-r", "8", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
  elif(argv[2] == "jps"):
    subprocess.call([argv[1], "tilt2", root + ".ppm", root, root + "-bump.ppm"])
    subprocess.call(["montage", root + "-L.ppm", root + "-R.ppm", "-geometry", "100%x100%", root + "-stereo.jps"])
  elif(argv[2] == "64"):
    subprocess.call([argv[1], "bump2", root + ".ppm", root + "-bump64.ppm"])
    subprocess.call(["convert", "-resize", "128x", "-resize", "x128<", "-resize", "50%", "-gravity", "center", "-crop", "64x64+0+0", root + "-bump64.ppm", root + "-bump64e" + ".jpeg"])
  elif(argv[2] == "demosaic"):
    print pixels
    subprocess.call(["convert", line, "-resize", str(int(100. / pixels * 100) / 100.) + "%", "-compress", "none", root + "-0.ppm"])
    s0 = int(np.ceil(np.log(pixels) / np.log(2.)))
    for s in range(0, s0):
      subprocess.call([argv[1], "enlarge", root + "-" + str(s) + ".ppm", root + "-" + str(s + 1) + ".ppm"])
  elif(argv[2] == "habit"):
    if(len(bhabit) <= 0):
      bhabit = line
      habit0 = root + "-mask.ppm"
      subprocess.call(["cp", bhabit, bhabit + "-out.obj"])
    else:
      subprocess.call([argv[1], "habit", habit0, line + "-out.obj", bhabit + "-out.obj", line])
      bhabit = line

