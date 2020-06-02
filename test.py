#! /usr/bin/env python

import os
import sys
import subprocess
import linecache
import numpy as np

argv   = sys.argv
pixels = 4
bhabit = ""

if(len(argv) < 4):
  print "no much argments."
elif(argv[2] == "match" or argv[2] == "match0" or argv[2] == "matcho"):
  root0, ext0 = os.path.splitext(argv[3])
  root1, ext1 = os.path.splitext(argv[4])
  nemph       = 4
  if((argv[2] == "match" or argv[2] == "match0") and len(argv) > 5):
    pixels = int(argv[5])
  elif(argv[2] == "matcho" and len(argv) > 6):
    pixels = int(argv[6])
    if(len(argv) > 7):
      nemph = int(argv[7])
  if(ext1 == ".obj" and (argv[2] == "match" or argv[2] == "match0")):
    subprocess.call([argv[1], argv[2], "16", "40", str(pixels), str(pixels), root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], "match-" + root0 + "-" + argv[4]])
  elif(ext1 == ".obj" and argv[2] != "match"):
    subprocess.call([argv[1], argv[2], argv[5], str(nemph), str(pixels), str(pixels), root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], argv[5]])
  elif(argv[2] == "match" or argv[2] == "match0"):
    subprocess.call([argv[1], argv[2], "16", "40", str(pixels), str(pixels), root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm", "match-" + root0 + "-" + root1])
  else:
    subprocess.call([argv[1], argv[2], argv[5], str(nemph), str(pixels), str(pixels), root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", argv[5]])
    subprocess.call(["ffmpeg", "-loop", "1", "-i", argv[5] + "-%d-" + str(nemph) + ".ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", argv[5] + ".mp4"])
elif(argv[2] == "pmerge"):
  root0, ext0 = os.path.splitext(argv[3])
  root1, ext1 = os.path.splitext(argv[4])
  subprocess.call([argv[1], argv[2], "1", ".1", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-pose.txt"])
elif(argv[2] == "pcomp"):
  for line in argv[3:]:
    root, ext = os.path.splitext(line)
    if(ext != ".ppm"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    subprocess.call(["sh", "-c", argv[1] + " plist 1 .1 " + root + ".ppm > " + root + "-plist.txt"])
  caller = ["sh", "-c", argv[1] + " pcomp "]
  for line in argv[3:]:
    root, ext = os.path.splitext(line)
    caller[- 1] += " " + root + "-plist.txt"
  caller[- 1] += " > complement-list.txt"
  subprocess.call(caller)
  subprocess.call([argv[1], "pdraw", "complement-list.txt", "complement-pdraw.ppm", "600", "600", "1"])
elif(argv[2] == "pcomp0"):
  root0, ext0 = os.path.splitext(argv[3])
  if(ext0 != ".ppm"):
    subprocess.call(["convert", argv[3], "-compress", "none", root0 + ".ppm"])
  subprocess.call([argv[1], "pose", "1", ".1", root0 + "-pose.txt", root0 + ".ppm", root0 + "-bump.ppm", str(pixels), root0 + "-pose"])
  for line in argv[4:]:
    root, ext = os.path.splitext(line)
    if(ext != ".ppm"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    subprocess.call([argv[1], "pmerge", "1", ".1", root0 + "-bump.ppm", root + "-bump.ppm", root + "-pose-" + root0 + ".txt"])
  cache = linecache.getline(root0 + ".ppm", 2).split(" ")
  x     = int(cache[1])
  y     = int(cache[0])
  caller = ["sh", "-c", argv[1] + " pcomp " + root0 + "-pose.txt"]
  for line in argv[4:]:
    root, ext = os.path.splitext(line)
    caller[- 1] += " " + root + "-pose-" + root0 + ".txt"
  caller[- 1] += " > complement-list0.txt"
  subprocess.call(caller)
  subprocess.call([argv[1], "pdraw", "complement-list0.txt", "complement-pdraw.ppm", str(x), str(y), "0"])
  subprocess.call([argv[1], "poso", "1", ".1", "complement-list0.txt", root0 + ".ppm", root0 + "-bump.ppm", str(pixels), root0 + "-pose"])
  for s in range(0, pixels):
    subprocess.call(["cp", root + "-pose" + str(pixels - 1 - s) + ".ppm", root + "-pose" + str(pixels + s) + ".ppm"])
  subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-pose%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-pose.mp4"])
else:
  for line in argv[3:]:
    try:
      pixels = int(line)
      continue
    except:
      root, ext = os.path.splitext(line)
    if(ext != ".ppm"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    if(argv[2] == "sample"):
      cmd = ["python2", argv[0], argv[1], "col", line]
      subprocess.call(cmd)
      cmd[3] = "bump"
      subprocess.call(cmd)
      cmd[3] = "obj"
      subprocess.call(cmd)
      cmd[3] = "mtl"
      subprocess.call(cmd)
      cmd[3] = "jps"
      subprocess.call(cmd)
      cmd[3] = "tilt"
      subprocess.call(cmd)
      cmd[3] = "btilt"
      subprocess.call(cmd)
      cmd[3] = "btilt2"
      subprocess.call(cmd)
      cmd[3] = "flicker"
      subprocess.call(cmd)
      cmd[3] = "light"
      subprocess.call(cmd)
      subprocess.call(["mogrify", "-format", "png", "*.ppm"])
    elif(argv[2] == "col"):
      subprocess.call([argv[1], "collect", root + ".ppm", root + "-collect.ppm"])
    elif(argv[2] == "penetrate"):
      subprocess.call(["cp", root + ".ppm", root + "-penetrate-light.ppm"])
      for s in range(0, pixels):
        subprocess.call(["convert", root + "-penetrate-light.ppm", "-blur", "2x2", "-compress", "none", root + "-penetrate.ppm"])
        subprocess.call([argv[1], "light", "1", root + "-penetrate.ppm", root + "-penetrate-light.ppm"])
      subprocess.call(["convert", root + "-penetrate-light.ppm", root + "-penetrate.png"])
    elif(argv[2] == "enlp"):
      subprocess.call(["cp", root + ".ppm", root + "-enl.ppm"])
      for s in range(0, pixels):
        subprocess.call([argv[1], "enlarge", root + "-enl.ppm", root + "-enl0.ppm"])
        subprocess.call(["convert", root + "-enl0.ppm", "-resize", "75%", "-compress", "none", root + "-enl.ppm"])
      subprocess.call(["convert", root + "-enl.ppm", root + "-enl.png"])
    elif(argv[2] == "light"):
      subprocess.call([argv[1], "light", str(pixels), root + ".ppm", root + "-light.ppm"])
    elif(argv[2] == "bump"):
      subprocess.call([argv[1], "bump", str(pixels), root + ".ppm", root + "-bump.ppm"])
      subprocess.call(["convert", root + "-bump.ppm", root + ".ppm", "-compose", "hard-light", "-composite", root + "-emph" + ext])
      subprocess.call(["convert", line, root + "-bump.ppm", "-channel-fx", "| gray=>alpha", root + "-alpha.png"])
    elif(argv[2] == "pextend"):
      subprocess.call([argv[1], "pextend", str(pixels), root + ".ppm", root + "-pextend.ppm"])
    elif(argv[2] == "mask0"):
      subprocess.call(["convert", root + ".ppm", "-fill", "black", "-colorize", "100", root + "-mask.png"])
    elif(argv[2] == "mask"):
      subprocess.call(["convert", root + "-mask.png", "-compress", "none", root + "-mask.ppm"])
    elif(argv[2] == "obj"):
      subprocess.call([argv[1], "obj", str(pixels), "1", ".001", "0", root + "-bump.ppm", root + ".obj"])
      subprocess.call([argv[1], "obj", str(pixels), ".1", ".001", ".4", root + "-bump.ppm", root + "-mask.ppm", root + "-stand.obj"])
      subprocess.call(["xcrun", "scntool", "--convert", root + ".obj", "--format", "scn", "--output", root + ".scn"])
    elif(argv[2] == "jps"):
      subprocess.call([argv[1], "tilt", "1", "4", ".005", root + ".ppm", root + ".obj", root + "-L.ppm"])
      subprocess.call([argv[1], "tilt", "3", "4", ".005", root + ".ppm", root + ".obj", root + "-R.ppm"])
      subprocess.call(["montage", root + "-R.ppm", root + "-L.ppm", "-geometry", "100%x100%", root + "-stereo.jps"])
      subprocess.call(["montage", root + "-stereo.jps", "-geometry", "100%x100%", root + "-stereo.png"])
    elif(argv[2] == "tilt"):
      for s in range(0, 16):
        subprocess.call([argv[1], "tilt", str(s), "16", ".05", root + ".ppm", root + ".obj", root + "-tilt-base-" + str(s) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-tilt-base-%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
    elif(argv[2] == "btilt"):
      for s in range(0, 32):
        subprocess.call([argv[1], "tilt", "1", "4", str((s - 16) / 16. * .1), root + ".ppm", root + ".obj", root + "-btilt-base-" + str(s) + ".ppm"])
        subprocess.call(["cp", root + "-btilt-base-" + str(s) + ".ppm", root + "-btilt-base-" + str(32 * 2 - s - 1) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt-base-%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-b.mp4"])
    elif(argv[2] == "btilt2"):
      subprocess.call([argv[1], "tilt", "1", "4", ".0125", root + ".ppm", root + ".obj", root + "-btilt2-base-0.ppm"])
      subprocess.call([argv[1], "tilt", "3", "4", ".0125", root + ".ppm", root + ".obj", root + "-btilt2-base-1.ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt2-base-%d.ppm", "-r", "20", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-b2.mp4"])
    elif(argv[2] == "flicker"):
      for s in range(0, 16):
        subprocess.call([argv[1], "obj", "0", str(pixels), str(.15 / 16. * s), root + "-bump.ppm", root + "-flicker.obj"])
        subprocess.call([argv[1], "tilt", "1", "4", "0.25", root + ".ppm", root + "-flicker.obj", root + "-flicker-base-" + str(s) + "-L.ppm"])
        subprocess.call([argv[1], "tilt", "3", "4", "0.25", root + ".ppm", root + "-flicker.obj", root + "-flicker-base-" + str(s) + "-R.ppm"])
        subprocess.call(["montage", root + "-flicker-base-" + str(s) + "-R.ppm", root + "-flicker-base-" + str(s) + "-L.ppm", "-geometry", "100%x100%", root + "-flicker-base-" + str(s) + ".png"])
        subprocess.call(["cp", root + "-flicker-base-" + str(s) + ".png", root + "-flicker-base-" + str(16 * 2 - s - 1) + ".png"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-flicker-base-%d.png", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-ficker.mp4"])
    elif(argv[2] == "sbox"):
      for s in range(0, 4):
        subprocess.call([argv[1], "sbox", str(- s * 16), "1", root + ".ppm", root + ".obj", root + "-sbox-" + str(3 - s) + ".ppm"])
    elif(argv[2] == "demosaic"):
      subprocess.call(["convert", line, "-resize", str(int(10000 / float(pixels)) / 100.) + "%", "-compress", "none", root + "-demosaic0.ppm"])
      subprocess.call(["python2", argv[0], argv[1], "enlp", str(int(np.ceil(np.log(pixels) / np.log(1.5)))), root + "-demosaic0.ppm"])
    elif(argv[2] == "habit"):
      if(len(bhabit) <= 0):
        bhabit = line
        subprocess.call(["cp", bhabit, bhabit + "-out.obj"])
      else:
        subprocess.call([argv[1], "habit", bhabit + "-out.obj", line, line + "-out.obj"])
        bhabit = line
    elif(argv[2] == "extend"):
      for tami in range(1, 5):
        tam = tami / 360. * 15
        for s in range(0, 4):
          subprocess.call([argv[1], "tilt", str(s), "4", str(tam), root + ".ppm", root + ".obj", root + "-tilt" + str(s) + ".ppm"])
          subprocess.call([argv[1], "bump", "4", root + "-tilt" + str(s) + ".ppm", root + "-bumpext" + str(s) + ".ppm"])
          subprocess.call([argv[1], "obj", "1", "1", ".001", "0", root + "-bumpext" + str(s) + ".ppm", root + "-bumpext" + str(s) + ".obj"])
        subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext0.obj", "0", "4", str(tam), root + "-bumpextA.obj"])
        subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext1.obj", "1", "4", str(tam), root + "-bumpextB.obj"])
        subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext2.obj", "2", "4", str(tam), root + "-bumpextC.obj"])
        subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext3.obj", "3", "4", str(tam), root + "-bumpextD.obj"])
        subprocess.call([argv[1], "habit", root + "-bumpextA.obj", root + "-bumpextC.obj", "0", "1", "0", root + "-bumpextAC.obj"])
        subprocess.call([argv[1], "habit", root + "-bumpextB.obj", root + "-bumpextD.obj", "0", "1", "0", root + "-bumpextBD.obj"])
        subprocess.call([argv[1], "habit", root + "-bumpextAC.obj", root + "-bumpextBD.obj", "0", "1", "0", root + "-bumpext.obj"])
        subprocess.call(["cp", root + "-bumpext.obj", root + ".obj"])
    elif(argv[2] == "pose"):
      subprocess.call([argv[1], "pose", "1", ".1", root + "-pose.txt", root + ".ppm", root + "-bump.ppm", str(pixels), root + "-pose"])
    elif(argv[2] == "poso"):
      subprocess.call([argv[1], "poso", "1", ".1", root + "-pose.txt", root + ".ppm", root + "-bump.ppm", str(pixels), root + "-pose"])
      for s in range(0, pixels):
        subprocess.call(["cp", root + "-pose" + str(pixels - 1 - s) + ".ppm", root + "-pose" + str(pixels + s) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-pose%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-pose.mp4"])
    elif(argv[2] == "prep"):
      subprocess.call(["convert", line, "-resize", str(pixels) + "@>", root + "-prep.png"])

