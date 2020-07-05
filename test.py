#! /usr/bin/env python

import os
import sys
import subprocess

argv    = sys.argv
pixels  = 4
zratio  = .25
psi     = .05
bhabit  = ""

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
    subprocess.call([argv[1], argv[2], "16", "40", str(pixels), str(pixels), str(zratio), root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], "match-" + root0 + "-" + argv[4]])
  elif(ext1 == ".obj"):
    subprocess.call([argv[1], argv[2], argv[5], str(nemph), str(pixels), str(pixels), str(zratio), root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], argv[5]])
  elif(argv[2] == "match" or argv[2] == "match0"):
    subprocess.call([argv[1], argv[2], "16", "40", str(pixels), str(pixels), str(zratio), root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm", "match-" + root0 + "-" + root1])
  else:
    subprocess.call([argv[1], argv[2], argv[5], str(nemph), str(pixels), str(pixels), str(zratio), root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", argv[5]])
  if(argv[2] == "matcho"):
    subprocess.call(["ffmpeg", "-loop", "1", "-i", argv[5] + "-%d-" + str(nemph) + ".ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", argv[5] + ".mp4"])
elif(argv[2] == "pred"):
  cmd = [argv[1], argv[2], "pred.ppm"]
  for s in argv[3:]:
    r, e = os.path.splitext(s)
    cmd.append(r + ".ppm")
  subprocess.call(cmd)
elif(argv[2] == "pcopy" or argv[2] == "ppred"):
  roots = []
  exts  = []
  for s in argv[3:]:
    try:
      pixels = int(s)
      continue
    except:
      r, e = os.path.splitext(s)
      roots.append(r)
      exts.append(e)
  #cmd = [argv[1], argv[2], "1", ".175", str(zratio), str(pixels), "pose"]
  #cmd = [argv[1], argv[2], "4", ".05", str(zratio), str(pixels), "pose"]
  cmd = [argv[1], argv[2], "20", ".05", str(zratio), str(pixels), "pose"]
  if(argv[2] == "pcopy"):
    cmd.extend([roots[0] + ".ppm", roots[0] + "-bump.ppm", roots[1] + ".ppm", roots[1] + "-bump.ppm"])
  else:
    for s in roots:
      cmd.append(s + ".ppm")
      cmd.append(s + "-bump.ppm")
  subprocess.call(cmd)
  for s in range(0, pixels):
    subprocess.call(["cp", "pose" + str(pixels - 1 - s) + ".ppm", "pose" + str(pixels + s) + ".ppm"])
  subprocess.call(["ffmpeg", "-loop", "1", "-i", "pose%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", "pose.mp4"])
else:
  for line in argv[3:]:
    try:
      pixels = int(line)
      continue
    except:
      root, ext = os.path.splitext(line)
    if(ext != ".ppm"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    if(argv[2] == "collect" or argv[2] == "bump" or argv[2] == "enlarge" or argv[2] == "sharpen" or argv[2] == "pextend"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", str(pixels)])
    elif(argv[2] == "penetrate"):
      subprocess.call(["cp", root + ".ppm", root + "-penetrate-sharpen.ppm"])
      for s in range(0, pixels):
        subprocess.call(["convert", root + "-penetrate-sharpen.ppm", "-blur", "2x2", "-compress", "none", root + "-penetrate.ppm"])
        # this isn't enough, limit of this is enough but it has glitches.
        subprocess.call([argv[1], "sharpen", root + "-penetrate.ppm", root + "-penetrate-sharpen.ppm"])
      subprocess.call(["convert", root + "-penetrate-sharpen.ppm", root + "-penetrate.png"])
    elif(argv[2] == "obj"):
      subprocess.call([argv[1], "obj", str(pixels), "1",  str(zratio), "0", root + "-bump.ppm", root + ".obj"])
      subprocess.call([argv[1], "obj", str(pixels), ".1", str(zratio), ".4", root + "-bump.ppm", root + "-mask.ppm", root + "-stand.obj"])
      subprocess.call(["xcrun", "scntool", "--convert", root + ".obj", "--format", "scn", "--output", root + ".scn"])
    elif(argv[2] == "jps"):
      subprocess.call([argv[1], "tilt", "1", "4", str(psi), root + ".ppm", root + ".obj", root + "-L.ppm"])
      subprocess.call([argv[1], "tilt", "3", "4", str(psi), root + ".ppm", root + ".obj", root + "-R.ppm"])
      subprocess.call(["montage", root + "-R.ppm", root + "-L.ppm", "-geometry", "100%x100%", root + "-stereo.jps"])
      subprocess.call(["montage", root + "-stereo.jps", "-geometry", "100%x100%", root + "-stereo.png"])
    elif(argv[2] == "tilt"):
      for s in range(0, pixels):
        subprocess.call([argv[1], "tilt", str(s), str(pixels), str(psi), root + ".ppm", root + ".obj", root + "-tilt-base-" + str(s) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-tilt-base-%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
    elif(argv[2] == "btilt"):
      for s in range(0, pixels * 2):
        subprocess.call([argv[1], "tilt", "1", "4", str((s - pixels) / float(pixels) * psi), root + ".ppm", root + ".obj", root + "-btilt-base-" + str(s) + ".ppm"])
        subprocess.call(["cp", root + "-btilt-base-" + str(s) + ".ppm", root + "-btilt-base-" + str(pixels * 4 - s - 1) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt-base-%d.ppm", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-b.mp4"])
    elif(argv[2] == "btilt2"):
      subprocess.call([argv[1], "tilt", "1", "4", str(psi), root + ".ppm", root + ".obj", root + "-btilt2-base-0.ppm"])
      subprocess.call([argv[1], "tilt", "3", "4", str(psi), root + ".ppm", root + ".obj", root + "-btilt2-base-1.ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt2-base-%d.ppm", "-r", "20", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-b2.mp4"])
    elif(argv[2] == "flicker"):
      for s in range(0, pixels):
        subprocess.call([argv[1], "obj", "1", "1", str(s / float(pixels) * zratio), "0", root + "-bump.ppm", root + "-flicker.obj"])
        subprocess.call([argv[1], "tilt", "1", "4", str(psi), root + ".ppm", root + "-flicker.obj", root + "-flicker-base-" + str(s) + "-L.ppm"])
        subprocess.call([argv[1], "tilt", "3", "4", str(psi), root + ".ppm", root + "-flicker.obj", root + "-flicker-base-" + str(s) + "-R.ppm"])
        subprocess.call(["montage", root + "-flicker-base-" + str(s) + "-R.ppm", root + "-flicker-base-" + str(s) + "-L.ppm", "-geometry", "100%x100%", root + "-flicker-base-" + str(s) + ".png"])
        subprocess.call(["cp", root + "-flicker-base-" + str(s) + ".png", root + "-flicker-base-" + str(16 * 2 - s - 1) + ".png"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-flicker-base-%d.png", "-r", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-ficker.mp4"])
    elif(argv[2] == "sbox"):
      for s in range(0, pixels):
        subprocess.call([argv[1], argv[2], str(int(s - (pixels + 1) / 2.)), str(pixels), str(zratio), root + ".ppm", root + ".obj", root + "-sbox-" + str(pixels - s) + ".ppm"])
    elif(argv[2] == "sboxb2w"):
      subprocess.call([argv[1], "b2w", root + "-sbox-1.ppm", root + "-sbox-bw-1.ppm", "1"])
      for s in range(1, pixels):
        subprocess.call([argv[1], "b2wd", root + "-sbox-" + str(s + 1) + ".ppm", root + "-sbox-bw-" + str(s + 1) + ".ppm", root + "-sbox-" + str(s) + ".ppm"])
    elif(argv[2] == "demosaic"):
      subprocess.call(["convert", line, "-resize", str(int(10000. * pow(1.5, - pixels)) / 100.) + "%", "-compress", "none", root + "-demosaic0.ppm"])
      for s in range(0, pixels):
        subprocess.call([argv[1], "enlarge", root + "-demosaic0.ppm", root + "-demosaic-enl.ppm", "3"])
        subprocess.call(["convert", root + "-demosaic-enl.ppm", "-resize", str(150 / 3.) + "%", "-compress", "none", root + "-demosaic0.ppm"])
    elif(argv[2] == "extend"):
      for tami in range(1, pixels + 1):
        tam = tami / 360. * 60 / pixels
        for s in range(0, 4):
          subprocess.call([argv[1], "tilt", str(s), "4", str(tam), root + ".ppm", root + ".obj", root + "-tilt" + str(s) + ".ppm"])
          subprocess.call([argv[1], "bump", root + "-tilt" + str(s) + ".ppm", root + "-bumpext" + str(s) + ".ppm"])
          subprocess.call([argv[1], "obj", "1", "1", str(zratio), "0", root + "-bumpext" + str(s) + ".ppm", root + "-bumpext" + str(s) + ".obj"])
        subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext0.obj", "0", "4", str(tam), root + "-bumpextA.obj"])
        subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext1.obj", "1", "4", str(tam), root + "-bumpextB.obj"])
        subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext2.obj", "2", "4", str(tam), root + "-bumpextC.obj"])
        subprocess.call([argv[1], "habit", root + ".obj", root + "-bumpext3.obj", "3", "4", str(tam), root + "-bumpextD.obj"])
        subprocess.call([argv[1], "habit", root + "-bumpextA.obj", root + "-bumpextC.obj", "0", "1", "0", root + "-bumpextAC.obj"])
        subprocess.call([argv[1], "habit", root + "-bumpextB.obj", root + "-bumpextD.obj", "0", "1", "0", root + "-bumpextBD.obj"])
        subprocess.call([argv[1], "habit", root + "-bumpextAC.obj", root + "-bumpextBD.obj", "0", "1", "0", root + "-bumpext.obj"])
        subprocess.call(["cp", root + "-bumpext.obj", root + ".obj"])
    elif(argv[2] == "prep"):
      subprocess.call(["convert", line, "-resize", str(pixels) + "@>", root + "-prep.png"])
    elif(argv[2] == "prepsq"):
      subprocess.call(["convert", line, "-resize", str(pixels) + "x" + str(pixels) + "!", root + "-prepsq.png"])
    elif(argv[2] == "mask"):
      subprocess.call(["convert", root + "-mask.png", "-compress", "none", root + "-mask.ppm"])
    elif(argv[2] == "mask0"):
      subprocess.call(["convert", root + ".ppm", "-fill", "black", "-colorize", "100", root + "-mask.png"])

