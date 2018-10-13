#! /usr/bin/env python

import os
import sys
import subprocess
import glob
import numpy as np

argv = sys.argv

if(len(argv) < 4):
  print "not much argments."
  exit(- 1)

if(argv[2] == "match"):
  if(len(argv) < 4):
    print "no much argments."
    exit(- 1)
  root0, ext0 = os.path.splitext(argv[3])
  root1, ext1 = os.path.splitext(argv[4])
  if(ext0 != ".ppm"):
    subprocess.call(["convert", argv[3], "-compress", "none", root0 + ".ppm"])
  if(ext1 == ".obj" or ext1 == ".gltf"):
    subprocess.call([argv[1], "match", "4", "20", "8", "60", "15", root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], "match-" + root0 + "-" + argv[4]])
  elif(ext1 != ".ppm"):
    subprocess.call(["convert", argv[4], "-compress", "none", root1 + ".ppm"])
  subprocess.call([argv[1], "match", "4", "20", "8", "60", "15", root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", "match-" + root0 + "-" + root1])
  exit(0)

pixels = 4
bhabit = ""
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
    subprocess.call([argv[1], "enlarge", "1", root + ".ppm", root + "-enl-0.ppm"])
    subprocess.call(["convert", root + "-enl-0.ppm", "-sharpen", "4x4", root + "-enl.png"])
  elif(argv[2] == "bump"):
    subprocess.call([argv[1], "bump", root + ".ppm", root + "-bump.ppm"])
  elif(argv[2] == "pextend"):
    subprocess.call([argv[1], "pextend", "8", root + ".ppm", root + "-pextend.ppm"])
  elif(argv[2] == "emph"):
    subprocess.call(["convert", root + "-bump.ppm", root + ".ppm", "-compose", "hard-light", "-composite", root + "-emph" + ext])
  elif(argv[2] == "pnga"):
    subprocess.call(["convert", line, root + "-bump.ppm", "-channel-fx", "| gray=>alpha", root + "-alpha.png"])
  elif(argv[2] == "idet"):
    subprocess.call([argv[1], "idetect", root + ".ppm", root + "-idet.ppm"])
  elif(argv[2] == "mask0"):
    subprocess.call(["convert", root + ".ppm", "-fill", "black", "-colorize", "100", root + "-mask.png"])
  elif(argv[2] == "mask"):
    subprocess.call(["convert", root + "-mask.png", "-compress", "none", root + "-mask.ppm"])
  elif(argv[2] == "obj"):
    subprocess.call([argv[1], "obj", "0", "3", ".15", root + "-bump.ppm", root + ".obj"])
    subprocess.call([argv[1], "obj", "  .2", "3", ".15", root + "-bump.ppm", root + "-L.obj"])
    subprocess.call([argv[1], "obj", "- .2", "3", ".15", root + "-bump.ppm", root + "-R.obj"])
    subprocess.call([argv[1], "obj", "0", "3", ".15", root + "-bump.ppm", root + "-mask.ppm", root + "-mask.obj"])
    subprocess.call([argv[1], "obj", "stand", "3", ".2", "1", ".15", root + "-bump.ppm", root + "-mask.ppm", root + "-stand.obj"])
  elif(argv[2] == "scn"):
    subprocess.call(["xcrun", "scntool", "--convert", root + "-L.obj", "--format", "scn", "--output", root + "-L.scn"])
    subprocess.call(["xcrun", "scntool", "--convert", root + "-R.obj", "--format", "scn", "--output", root + "-R.scn"])
  elif(argv[2] == "mtl"):
    f = open(root + ".obj.mtl", "a")
    f.write("newmtl material0\n")
    f.write("Ka 1.000000 1.000000 1.000000\n")
    f.write("Kd 1.000000 1.000000 1.000000\n")
    f.write("Ks 0.000000 0.000000 0.000000\n")
    f.write("illum 1\n")
    f.write("map_Ka " + line + "\n")
    f.write("map_Kd " + line + "\n\n")
    f.close()
    subprocess.call(["cp", root + ".obj.mtl", root + "-L.obj.mtl"])
    subprocess.call(["cp", root + ".obj.mtl", root + "-R.obj.mtl"])
    subprocess.call(["cp", root + ".obj.mtl", root + "-mask.obj.mtl"])
    subprocess.call(["cp", root + ".obj.mtl", root + "-stand.obj.mtl"])
  elif(argv[2] == "jps"):
    subprocess.call([argv[1], "tilt", "1", "2", ".05", "0", root + ".ppm", root + "-bump.ppm", root + "-L.ppm"])
    subprocess.call([argv[1], "tilt", "0", "2", ".05", "0", root + ".ppm", root + "-bump.ppm", root + "-R.ppm"])
    subprocess.call(["montage", root + "-L.ppm", root + "-R.ppm", "-geometry", "100%x100%", root + "-stereo.jps"])
    subprocess.call(["montage", root + "-stereo.jps", "-geometry", "100%x100%", root + "-stereo.png"])
  elif(argv[2] == "tilt"):
    for s in range(0, 16):
      subprocess.call([argv[1], "tilt", str(s), "16", ".05", "0", root + ".ppm", root + "-bump.ppm", root + "-tilt-base-" + str(s) + ".ppm"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-tilt-base-%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
  elif(argv[2] == "btilt"):
    for s in range(0, 16):
      subprocess.call([argv[1], "tilt", "0", "2", str((s - 8) / 8. * .1), "0", root + ".ppm", root + "-bump.ppm", root + "-btilt-base-" + str(s) + ".ppm"])
      subprocess.call(["cp", root + "-btilt-base-" + str(s) + ".ppm", root + "-btilt-base-" + str(16 * 2 - s - 1) + ".ppm"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt-base-%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
  elif(argv[2] == "sbox"):
    for s in range(0, 8):
      subprocess.call([argv[1], "sbox", str(- s * 3), "1", root + ".ppm", root + ".obj", root + "-sbox-" + str(s) + ".ppm"])
  elif(argv[2] == "extend"):
    for tam in [.5 / 6., .75 / 6., 1. / 6.]:
      for s in range(0, 4):
        subprocess.call([argv[1], "tilt", str(s), "4", str(tam), "0", root + ".ppm", root + ".obj", root + "-tilt" + str(s) + ".ppm"])
        subprocess.call([argv[1], "bump", root + "-tilt" + str(s) + ".ppm", root + "-bumpext" + str(s) + ".ppm"])
        subprocess.call([argv[1], "obj", "0", "3", ".15", root + "-bumpext" + str(s) + ".ppm", root + "-bumpext" + str(s) + ".obj"])
      subprocess.call([argv[1], "habit", root + "-bumpext0.obj", root + "-bumpext2.obj", "2", "4", str(tam), root + "-bumpextA.obj"])
      subprocess.call([argv[1], "habit", root + "-bumpext1.obj", root + "-bumpext3.obj", "3", "4", str(tam), root + "-bumpextB.obj"])
      subprocess.call([argv[1], "habit", root + "-bumpextA.obj", root + "-bumpextB.obj", "0", "1", str(tam), root + "-bumpext.obj"])
      subprocess.call(["cp", root + "-bumpext.obj", root + ".obj"])
  elif(argv[2] == "demosaic"):
    print pixels
    s0  = int(np.ceil(np.log(pixels) / np.log(1.75)))
    subprocess.call(["convert", line, "-resize", str(int(10000 / float(pixels)) / 100.) + "%", "-compress", "none", root + "-0.ppm"])
    for s in range(0, s0):
      subprocess.call([argv[1], "enlarge", "1", root + "-" + str(s) + ".ppm", root + "-" + str(s + 1) + "-0.ppm"])
      subprocess.call(["convert", root + "-" + str(s + 1) + "-0.ppm", "-sharpen", "2x2", "-resize", "75%", "-compress", "none", root + "-" + str(s + 1) + ".ppm"])
  elif(argv[2] == "habit"):
    if(len(bhabit) <= 0):
      bhabit = line
      subprocess.call(["cp", bhabit, bhabit + "-out.obj"])
    else:
      subprocess.call([argv[1], "habit", bhabit + "-out.obj", line, line + "-out.obj"])
      bhabit = line
  elif(argv[2] == "prep"):
    subprocess.call(["convert", line, "-resize", str(pixels) + "@>", root + "-prep.png"])

