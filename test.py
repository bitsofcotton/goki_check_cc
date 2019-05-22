#! /usr/bin/env python

import os
import sys
import subprocess
import numpy as np

argv = sys.argv

if(len(argv) < 4):
  print "not much argments."
  exit(- 1)

if(argv[2] == "match" or argv[2] == "matcho"):
  if(len(argv) < 4):
    print "no much argments."
    exit(- 1)
  root0, ext0 = os.path.splitext(argv[3])
  root1, ext1 = os.path.splitext(argv[4])
  if(ext0 != ".ppm"):
    subprocess.call(["convert", argv[3], "-compress", "none", root0 + ".ppm"])
  if(ext1 == ".obj" or ext1 == ".gltf"):
    if(argv[2] == "match"):
      subprocess.call([argv[1], "match", "16", "40", "5", "5", root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], "match-" + root0 + "-" + argv[4]])
    else:
      for s in [0, .25, .5, .75, 1]:
        subprocess.call([argv[1], "matcho", argv[5], str(s), "5", "5", root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], "match-" + root0 + "-" + argv[4] + "-" + str(s)])
    exit(0)
  elif(ext1 != ".ppm"):
    subprocess.call(["convert", argv[4], "-compress", "none", root1 + ".ppm"])
  if(argv[2] == "match"):
    subprocess.call([argv[1], "match", "16", "40", "5", "5", root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm", "match-" + root0 + "-" + root1])
  else:
    for s in [0, .25, .5, .75, 1]:
      subprocess.call([argv[1], "matcho", argv[5], str(s), "5", "5", root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", "match-" + root0 + "-" + root1 + "-" + str(s)])
  exit(0)

if(argv[2] == "objtilt"):
  root0, ext0 = os.path.splitext(argv[3])
  if(ext0 != ".ppm"):
    subprocess.call(["convert", argv[3], "-compress", "none", root0 + ".ppm"])
  root = os.path.splitext(argv[4])[0]
  subprocess.call([argv[1], "drawr", root0 + ".ppm", argv[4], root + "-drawr-bump.ppm"])
  subprocess.call([argv[1], "drawm", root0 + ".ppm", argv[4], root + "-drawr-mask.ppm"])
  subprocess.call([argv[1], "obj", "0", "3", ".2", root + "-drawr-bump.ppm", root + "-drawr-mask.ppm", root + "-drawr.obj"])
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
    subprocess.call([argv[1], "enlarge", "1", root + ".ppm", root + "-enl.ppm"])
  elif(argv[2] == "enlp"):
    subprocess.call([argv[1], "enlarge", str(pixels), root + ".ppm", root + "-enl.ppm"])
  elif(argv[2] == "bump"):
    subprocess.call([argv[1], "bump", root + ".ppm", root + "-bump.ppm"])
  elif(argv[2] == "reshape"):
    subprocess.call([argv[1], "reshape", str(pixels), root + ".ppm", root + "-enl.ppm", root + "-enlout.ppm"])
  elif(argv[2] == "pextend"):
    subprocess.call([argv[1], "pextend", str(pixels), root + ".ppm", root + "-pextend.ppm"])
  elif(argv[2] == "emph"):
    subprocess.call(["convert", root + "-bump.ppm", root + ".ppm", "-compose", "hard-light", "-composite", root + "-emph" + ext])
  elif(argv[2] == "emphe"):
    subprocess.call(["convert", root + "-collect.ppm", root + ".ppm", "-compose", "hard-light", "-composite", root + "-emphe" + ext])
  elif(argv[2] == "pnga"):
    subprocess.call(["convert", line, root + "-bump.ppm", "-channel-fx", "| gray=>alpha", root + "-alpha.png"])
  elif(argv[2] == "idet"):
    subprocess.call([argv[1], "idetect", root + ".ppm", root + "-idet.ppm"])
  elif(argv[2] == "mask0"):
    subprocess.call(["convert", root + ".ppm", "-fill", "black", "-colorize", "100", root + "-mask.png"])
  elif(argv[2] == "mask"):
    subprocess.call(["convert", root + "-mask.png", "-compress", "none", root + "-mask.ppm"])
  elif(argv[2] == "obj"):
    subprocess.call([argv[1], "obj", "  0", "3", ".2", root + "-bump.ppm", root + ".obj"])
    subprocess.call([argv[1], "obj", " .2", "3", ".2", root + "-bump.ppm", root + "-L.obj"])
    subprocess.call([argv[1], "obj", "-.2", "3", ".2", root + "-bump.ppm", root + "-R.obj"])
    subprocess.call([argv[1], "obj", "  0", "3", ".2", root + "-bump.ppm", root + "-mask.ppm", root + "-mask.obj"])
    subprocess.call([argv[1], "obj", "stand", "3", ".2", "1", ".2", root + "-bump.ppm", root + "-mask.ppm", root + "-stand.obj"])
  elif(argv[2] == "scn"):
    subprocess.call(["xcrun", "scntool", "--convert", root + "-L.obj", "--format", "scn", "--output", root + "-L.scn"])
    subprocess.call(["xcrun", "scntool", "--convert", root + "-R.obj", "--format", "scn", "--output", root + "-R.scn"])
  elif(argv[2] == "mtl"):
    f = open(root + ".obj.mtl", "w")
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
    subprocess.call([argv[1], "tilt", "1", "4", ".005", "0", root + ".ppm", root + ".obj", root + "-L.ppm"])
    subprocess.call([argv[1], "tilt", "3", "4", ".005", "0", root + ".ppm", root + ".obj", root + "-R.ppm"])
    subprocess.call(["montage", root + "-R.ppm", root + "-L.ppm", "-geometry", "100%x100%", root + "-stereo.jps"])
    subprocess.call(["montage", root + "-stereo.jps", "-geometry", "100%x100%", root + "-stereo.png"])
  elif(argv[2] == "tilt"):
    for s in range(0, 16):
      subprocess.call([argv[1], "tilt", str(s), "16", ".05", "0", root + ".ppm", root + ".obj", root + "-tilt-base-" + str(s) + ".ppm"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-tilt-base-%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
  elif(argv[2] == "btilt"):
    for s in range(0, 32):
      subprocess.call([argv[1], "tilt", "1", "4", str((s - 16) / 16. * .1), "0", root + ".ppm", root + ".obj", root + "-btilt-base-" + str(s) + ".ppm"])
      subprocess.call(["cp", root + "-btilt-base-" + str(s) + ".ppm", root + "-btilt-base-" + str(32 * 2 - s - 1) + ".ppm"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt-base-%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-b.mp4"])
  elif(argv[2] == "flicker"):
    for s in range(0, 16):
      subprocess.call([argv[1], "obj", "0", "3", str(.15 / 16. * s), root + "-bump.ppm", root + "-flicker.obj"])
      subprocess.call([argv[1], "tilt", "1", "4", "0.25", "0", root + ".ppm", root + "-flicker.obj", root + "-flicker-base-" + str(s) + "-L.ppm"])
      subprocess.call([argv[1], "tilt", "3", "4", "0.25", "0", root + ".ppm", root + "-flicker.obj", root + "-flicker-base-" + str(s) + "-R.ppm"])
      subprocess.call(["montage", root + "-flicker-base-" + str(s) + "-R.ppm", root + "-flicker-base-" + str(s) + "-L.ppm", "-geometry", "100%x100%", root + "-flicker-base-" + str(s) + ".png"])
      subprocess.call(["cp", root + "-flicker-base-" + str(s) + ".png", root + "-flicker-base-" + str(16 * 2 - s - 1) + ".png"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-flicker-base-%d.png", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-ficker.mp4"])
  elif(argv[2] == "sbox"):
    for s in range(0, 4):
      subprocess.call([argv[1], "sbox", str(- s * 16), "1", root + ".ppm", root + ".obj", root + "-sbox-" + str(3 - s) + ".ppm"])
  elif(argv[2] == "extend"):
    for tami in range(1, 4):
      tam = tami / 360. * 8
      for s in range(0, 4):
        subprocess.call([argv[1], "tilt", str(s), "4", str(tam), "0", root + ".ppm", root + ".obj", root + "-tilt" + str(s) + ".ppm"])
        subprocess.call([argv[1], "bump", root + "-tilt" + str(s) + ".ppm", root + "-bumpext" + str(s) + ".ppm"])
        subprocess.call([argv[1], "obj", "0", "3", ".15", root + "-bumpext" + str(s) + ".ppm", root + "-bumpext" + str(s) + ".obj"])
      subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext0.obj", "0", "4", str(tam), root + "-bumpextA.obj"])
      subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext1.obj", "1", "4", str(tam), root + "-bumpextB.obj"])
      subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext2.obj", "2", "4", str(tam), root + "-bumpextC.obj"])
      subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext3.obj", "3", "4", str(tam), root + "-bumpextD.obj"])
      subprocess.call([argv[1], "habit", root + "-bumpextA.obj", root + "-bumpextC.obj", "0", "1", "0", root + "-bumpextAC.obj"])
      subprocess.call([argv[1], "habit", root + "-bumpextB.obj", root + "-bumpextD.obj", "0", "1", "0", root + "-bumpextBD.obj"])
      subprocess.call([argv[1], "habit", root + "-bumpextAC.obj", root + "-bumpextBD.obj", "0", "1", "0", root + "-bumpext.obj"])
      subprocess.call(["cp", root + "-bumpext.obj", root + ".obj"])
  elif(argv[2] == "demosaic"):
    print pixels
    subprocess.call(["convert", line, "-resize", str(int(10000 / float(pixels)) / 100.) + "%", "-compress", "none", root + "-demosaic0.ppm"])
    subprocess.call(["python", argv[0], argv[1], "enlp", str(int(np.ceil(np.log(pixels) / np.log(2)))), root + "-demosaic0.ppm"])
  elif(argv[2] == "habit"):
    if(len(bhabit) <= 0):
      bhabit = line
      subprocess.call(["cp", bhabit, bhabit + "-out.obj"])
    else:
      subprocess.call([argv[1], "habit", bhabit + "-out.obj", line, line + "-out.obj"])
      bhabit = line
  elif(argv[2] == "prep"):
    subprocess.call(["convert", line, "-resize", str(pixels) + "@>", root + "-prep.png"])

