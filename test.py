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

if(argv[2] == "match" or argv[2] == "pmatch"):
  if(len(argv) < 5):
    print "no much argments."
    exit(- 1)
  root0, ext0 = os.path.splitext(argv[3])
  root1, ext1 = os.path.splitext(argv[4])
  if(ext0 != ".ppm"):
    subprocess.call(["convert", argv[3], "-compress", "none", root0 + ".ppm"])
  if(ext1 == ".obj"):
    subprocess.call([argv[1], "match3d", root0 + ".ppm", "match3dbase-" + root0 + "-" + root1 + "-", argv[4], root0 + "-bump.ppm", root0 + "-mask.ppm"])
  elif(ext1 == ".gltf"):
    subprocess.call([argv[1], "match3dbone", root0 + ".ppm", "match3dbase-" + root0 + "-" + root1 + "-", argv[4], root0 + "-bump.ppm", root0 + "-mask.ppm"])
  else:
    if(ext1 != ".ppm"):
      subprocess.call(["convert", argv[4], "-compress", "none", root1 + ".ppm"])
    if(len(argv) < 6):
      subprocess.call([argv[1], argv[2], root0 + ".ppm", "matchbase-" + root0 + "-" + root1 + "-", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm"])
    elif(os.path.splitext(argv[5])[1] == ".obj"):
      subprocess.call([argv[1], "match2dh3d", root0 + ".ppm", "match3dh2dbase-" + root0 + root1 + os.path.splitext(argv[5])[0] + "-", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm", argv[5]])
    elif(os.path.splitext(argv[5])[1] == ".gltf"):
      subprocess.call([argv[1], "match2dh3dbone", root0 + ".ppm", "match3dh2dbase-" + root0 + root1 + os.path.splitext(argv[5])[0] + "-", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm", argv[5]])
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
    subprocess.call([argv[1], "enlarge", root + ".ppm", root + "-enl-0.ppm"])
    subprocess.call(["convert", root + "-enl-0.ppm", "-sharpen", "4x4", root + "-enl.png"])
  elif(argv[2] == "enl4"):
    subprocess.call([argv[1], "enlarge4", root + ".ppm", root + "-enl4-0.ppm"])
    subprocess.call(["convert", root + "-enl4-0.ppm", "-sharpen", "32x32", "-resize", "50%", root + "-enl4.png"])
  elif(argv[2] == "bump"):
    subprocess.call([argv[1], "bump", root + ".ppm", root + "-bump.ppm"])
  elif(argv[2] == "emph"):
    subprocess.call(["convert", root + "-bump.ppm", root + ".ppm", "-compose", "hard-light", "-composite", root + "-emph" + ext])
  elif(argv[2] == "mask0"):
    subprocess.call(["convert", root + ".ppm", "-fill", "black", "-colorize", "100", root + "-mask.png"])
  elif(argv[2] == "mask"):
    subprocess.call(["convert", root + "-mask.png", "-compress", "none", root + "-mask.ppm"])
  elif(argv[2] == "obj"):
    subprocess.call([argv[1], "obj",     root + "-bump.ppm", root + "0.obj"])
    subprocess.call([argv[1], "arobj",   root + "-bump.ppm", root + "1"])
    subprocess.call([argv[1], "maskobj", root + "-mask.ppm", root + "0.obj", root + ".obj", ".05", "2"])
  elif(argv[2] == "scn"):
    subprocess.call(["xcrun", "scntool", "--convert", root + "1-L.obj", "--format", "scn", "--output", root + "-L.scn"])
    subprocess.call(["xcrun", "scntool", "--convert", root + "1-R.obj", "--format", "scn", "--output", root + "-R.scn"])
  elif(argv[2] == "mtl"):
    subprocess.call(["cp", os.path.dirname(argv[1]) + "/material.mtl", root + ".obj.mtl"])
    f = open(root + ".obj.mtl", "a")
    f.write("map_Ka " + line + "\n")
    f.write("map_Kd " + line + "\n\n")
    f.close()
    subprocess.call(["cp", root + ".obj.mtl", root + "0.obj.mtl"])
    subprocess.call(["cp", root + ".obj.mtl", root + "1-L.obj.mtl"])
    subprocess.call(["cp", root + ".obj.mtl", root + "1-R.obj.mtl"])
  elif(argv[2] == "tilt"):
    subprocess.call([argv[1], "tilt", root + ".ppm", root + "-tilt-base", root + "0.obj"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-tilt-base-%d.ppm", "-r", "8", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
  elif(argv[2] == "btilt"):
    subprocess.call([argv[1], "tilt4", root + ".ppm", root + "-btilt-base", root + "0.obj"])
    files = []
    for s in range(0, 200):
      file = glob.glob(root + "-btilt-base-" + str(s) + ".ppm")
      if(len(file) < 1):
        break
      files.append(file[0])
    for s in range(0, len(files)):
      subprocess.call(["cp", files[len(files) - s - 1], root + "-btilt-base-" + str(s + len(files)) + ".ppm"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt-base-%d.ppm", "-r", "8", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
  elif(argv[2] == "btilt2"):
    subprocess.call([argv[1], "tilt5", root + ".ppm", root + "-btilt-base", root + "0.obj"])
    files = []
    for s in range(0, 200):
      file = glob.glob(root + "-btilt-base-" + str(s) + ".ppm")
      if(len(file) < 1):
        break
      files.append(file[0])
    for s in range(0, len(files)):
      subprocess.call(["cp", files[len(files) - s - 1], root + "-btilt-base-" + str(s + len(files)) + ".ppm"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt-base-%d.ppm", "-r", "8", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
  elif(argv[2] == "flicker"):
    subprocess.call([argv[1], "tiltp", root + ".ppm", root + "-tiltrot-base", root + "0.obj"])
    files = []
    for s in range(0, 200):
      file = glob.glob(root + "-tiltrot-base-" + str(s) + "-[LR].ppm")
      if(len(file) < 2):
        break
      files.append(root + "-" + str(s) + "-tr.png")
      subprocess.call(["montage", file[0], file[1], "-geometry", "100%x100%", files[- 1]])
    for s in range(0, len(files)):
      subprocess.call(["cp", files[len(files) - s - 1], root + "-" + str(s + len(files)) + "-tr.png"])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-%d-tr.png", "-r", "8", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-tr.mp4"])
  elif(argv[2] == "pextend"):
    subprocess.call([argv[1], "pextend", root + ".ppm", root + "-pextend.ppm"])
  elif(argv[2] == "extend"):
    for tam in [.5 / 6., .75 / 6., 1. / 6.]:
      subprocess.call([argv[1], "tilt3", root + ".ppm", root + "-tilt3", root + "0.obj", str(tam)])
      for s in range(0, 4):
        subprocess.call([argv[1], "bump", root + "-tilt3-" + str(s) + ".ppm", root + "-bumpext-" + str(s) + ".ppm"])
        subprocess.call([argv[1], "obj", root + "-bumpext-" + str(s) + ".ppm", root + "-bumpext-" + str(s) + ".obj"])
      subprocess.call([argv[1], "habit2", root + "-mask.ppm", root + "-bumpextA.obj", root + "-bumpext-0.obj", root + "0.obj", "0", "4", str(tam)])
      subprocess.call([argv[1], "habit2", root + "-mask.ppm", root + "-bumpextB.obj", root + "-bumpext-1.obj", root + "0.obj", "1", "4", str(tam)])
      subprocess.call([argv[1], "habit2", root + "-mask.ppm", root + "-bumpextC.obj", root + "-bumpext-2.obj", root + "0.obj", "2", "4", str(tam)])
      subprocess.call([argv[1], "habit2", root + "-mask.ppm", root + "-bumpextD.obj", root + "-bumpext-3.obj", root + "0.obj", "3", "4", str(tam)])
      subprocess.call([argv[1], "habit2", root + "-mask.ppm", root + "-bumpextAC.obj", root + "-bumpextA.obj-emph.obj", root + "-bumpextC.obj-emph.obj", "4", "4", str(tam)])
      subprocess.call([argv[1], "habit2", root + "-mask.ppm", root + "-bumpextBD.obj", root + "-bumpextB.obj-emph.obj", root + "-bumpextD.obj-emph.obj", "4", "4", str(tam)])
      subprocess.call([argv[1], "habit2", root + "-mask.ppm", root + "-bumpext.obj", root + "-bumpextAC.obj-emph.obj", root + "-bumpextBD.obj-emph.obj", "4", "4", str(tam)])
      subprocess.call([argv[1], "maskobj2", root + "-mask.ppm", root + "-bumpext.obj-emph.obj", root + "-bumpextmask.obj"])
      subprocess.call(["cp", root + "-bumpext.obj-emph.obj", root + "-" + str(tam) + ".obj"])
      subprocess.call(["cp", root + "-bumpext.obj-emph.obj", root + "0.obj"])
  elif(argv[2] == "jps"):
    subprocess.call([argv[1], "tilt2", root + ".ppm", root, root + "-bump.ppm"])
    subprocess.call(["montage", root + "-L.ppm", root + "-R.ppm", "-geometry", "100%x100%", root + "-stereo.jps"])
    subprocess.call(["montage", root + "-stereo.jps", "-geometry", "100%x100%", root + "-stereo.png"])
  elif(argv[2] == "pnga"):
    subprocess.call(["convert", line, root + "-bump.ppm", "-channel-fx", "| gray=>alpha", root + "-alpha.png"])
  elif(argv[2] == "64"):
    subprocess.call([argv[1], "bump2", root + ".ppm", root + "-bump64.ppm"])
    subprocess.call(["convert", "-resize", "128x", "-resize", "x128<", "-resize", "50%", "-gravity", "center", "-crop", "64x64+0+0", root + "-bump64.ppm", root + "-bump64e" + ".jpeg"])
  elif(argv[2] == "demosaic"):
    print pixels
    s0  = int(np.ceil(np.log(pixels) / np.log(8.)))
    subprocess.call(["convert", line, "-resize", str(int(10000 / float(pixels)) / 100.) + "%", root + "-0.png"])
    for s in range(0, s0):
      subprocess.call(["python", argv[0], argv[1], "enl4", root + "-" + str(s) + ".png"])
      subprocess.call(["mv", root + "-" + str(s) + "-enl4.png", root + "-" + str(s + 1) + ".png"])
  elif(argv[2] == "habit"):
    if(len(bhabit) <= 0):
      bhabit = line
      habit0 = root + "-mask.ppm"
      subprocess.call(["cp", bhabit, bhabit + "-out.obj"])
    else:
      subprocess.call([argv[1], "habit", habit0, line + "-out.obj", bhabit + "-out.obj", line])
      bhabit = line
  elif(argv[2] == "prep"):
    subprocess.call(["convert", line, "-resize", str(pixels) + "@>", root + "-prep.png"])

