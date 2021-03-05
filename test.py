#! /usr/bin/env python

import os
import sys
import subprocess

argv    = sys.argv
pixels  = 4
zratio  = .25
psi     = .1
rot     = 7
comp    = 4
tppred  = .0025
bhabit  = ""

if(len(argv) < 4):
  print "no much argments."
elif(argv[2] == "match" or argv[2] == "match0" or argv[2] == "matcho" or argv[2] == "rmatch" or argv[2] == "rmatch0"):
  root0, ext0 = os.path.splitext(argv[3])
  root1, ext1 = os.path.splitext(argv[4])
  nemph       = 4
  if((argv[2] == "match" or argv[2] == "match0") and len(argv) > 5):
    pixels = int(argv[5])
  elif((argv[2] == "matcho" or argv[2] == "rmatch" or argv[2] == "rmatch0") and len(argv) > 6):
    pixels = int(argv[6])
    if(len(argv) > 7):
      nemph = int(argv[7])
  if(ext1 == ".obj" and (argv[2] == "match" or argv[2] == "match0")):
    subprocess.call([argv[1], argv[2], "16", "40", str(pixels), str(pixels), str(zratio), root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], "match-" + root0 + "-" + argv[4]])
  elif(ext1 == ".obj" and (argv[2] == "rmatch" or argv[2] == "rmatch0")):
    subprocess.call([argv[1], argv[2], str(nemph), "40", str(pixels), str(pixels), str(zratio), root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], argv[5]])
  elif(ext1 == ".obj"):
    subprocess.call([argv[1], argv[2], argv[5], str(nemph), str(pixels), str(pixels), str(zratio), root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", argv[4], argv[5]])
  elif(argv[2] == "match" or argv[2] == "match0"):
    subprocess.call([argv[1], argv[2], "16", "40", str(pixels), str(pixels), str(zratio), root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm", "match-" + root0 + "-" + root1])
  elif(argv[2] == "rmatch" or argv[2] == "rmatch0"):
    subprocess.call([argv[1], argv[2], "16", "40", str(pixels), str(pixels), str(zratio), root0 + ".ppm", root0 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", root0 + "-mask.ppm", root1 + "-mask.ppm", argv[5]])
  else:
    subprocess.call([argv[1], argv[2], argv[5], str(nemph), str(pixels), str(pixels), str(zratio), root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", argv[5]])
  if(argv[2] == "matcho" or argv[2] == "rmatch" or argv[2] == "rmatch0"):
    subprocess.call(["ffmpeg", "-loop", "1", "-i", argv[5] + "-%d-" + str(nemph) + ".ppm", "-framerate", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", argv[5] + ".mp4"])
elif(argv[2] == "pred" or argv[2] == "predf" or argv[2] == "cat"):
  cmd = [argv[1], argv[2], argv[2] + ".ppm"]
  for s in argv[3:]:
    r, e = os.path.splitext(s)
    if(e != ".ppm"):
      subprocess.call(["convert", s, "-compress", "none", r + ".ppm"])
    cmd.append(r + ".ppm")
  subprocess.call(cmd)
elif(argv[2] == "ppred" or argv[2] == "ppredr"):
  roots = []
  for s in argv[3:]:
    try:
      pixels = int(s)
      continue
    except:
      r, e = os.path.splitext(s)
      roots.append(r)
  #cmd = [argv[1], argv[2], "1", str(tppred), str(zratio), str(pixels), str(comp), "pose"]
  #cmd = [argv[1], argv[2], "2", str(tppred), str(zratio), str(pixels), str(comp), "pose"]
  cmd = [argv[1], argv[2], "4", str(tppred), str(zratio), str(pixels), str(comp), "pose"]
  #cmd = [argv[1], argv[2], "20", str(tppred), str(zratio), str(pixels), str(comp), "pose"]
  for s in roots:
    cmd.append(s + ".ppm")
    if(argv[2] == "ppred"):
      cmd.append(s + "-bump.ppm")
    else:
      cmd.append(s + ".ppm")
  print " ".join(cmd)
  subprocess.call(cmd)
  subprocess.call(["ffmpeg", "-i", "pose%d.ppm", "-frame_drop_threshold", "0.0", "-framerate", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "pose.mp4"])
else:
  for line in argv[3:]:
    try:
      pixels = int(line)
      continue
    except:
      root, ext = os.path.splitext(line)
    if(ext != ".ppm" and argv[2] != "prep" and argv[2] != "prepsq"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    if(argv[2] == "bump"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + "0.ppm", str(pixels), str(rot)])
      subprocess.call(["convert", root + "-" + argv[2] + "0.ppm", "-blur", "32x32+16", "-compress", "none", root + "-" + argv[2] + ".ppm"])
    elif(argv[2] == "enlarge"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + "-global.ppm", str(pixels), str(rot)])
      subprocess.call(["convert", root + "-" + argv[2] + "-global.ppm", "-blur", str(int(pow(pixels, .5))), "-compress", "none", root + "-" + argv[2] + "-global-blur.ppm"])
      subprocess.call([argv[1], "sharpen", root + "-" + argv[2] + "-global-blur.ppm", root + "-" + argv[2] + "-global-blur-sharpen.ppm", str(int(pow(pixels, .5))), str(rot)])
      subprocess.call(["convert", root + "-" + argv[2] + "-global-blur-sharpen.ppm", "-equalize", root + "-" + argv[2] + "-gbse.png"])
      subprocess.call(["convert", root + ".ppm", "-resize", str(pixels * 100) + "%", "-equalize", root + "-" + argv[2] + "-local.png"])
      subprocess.call(["convert", root + "-" + argv[2] + "-gbse.png", root + "-" + argv[2] + "-local.png", "-average", "-equalize", root + "-" + argv[2] + ".png"])
    elif(argv[2] == "collect" or argv[2] == "sharpen" or argv[2] == "bump" or argv[2] == "enlarge" or argv[2] == "flarge" or argv[2] == "pextend"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", str(pixels), str(rot)])
    elif(argv[2] == "obj"):
      subprocess.call([argv[1], "obj", str(pixels), "1",  str(zratio), "0", root + "-bump.ppm", root + ".obj"])
      #subprocess.call([argv[1], "obj", str(pixels), ".1", str(zratio), ".4", root + "-bump.ppm", root + "-mask.ppm", root + "-stand.obj"])
      #subprocess.call(["xcrun", "scntool", "--convert", root + ".obj", "--format", "scn", "--output", root + ".scn"])
    elif(argv[2] == "jps"):
      subprocess.call([argv[1], "tilt", "1", "4", str(psi), root + ".ppm", root + ".obj", root + "-L.ppm"])
      subprocess.call([argv[1], "tilt", "3", "4", str(psi), root + ".ppm", root + ".obj", root + "-R.ppm"])
      subprocess.call(["montage", root + "-R.ppm", root + "-L.ppm", "-geometry", "100%x100%", root + "-stereo.jps"])
      subprocess.call(["montage", root + "-stereo.jps", "-geometry", "100%x100%", root + "-stereo.png"])
      subprocess.call(["cp", root + "-L.ppm", root + "-jps-0.ppm"])
      subprocess.call(["cp", root + "-R.ppm", root + "-jps-1.ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-jps-%d.ppm", "-framerate", "20", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-LR.mp4"])
    elif(argv[2] == "tilt"):
      for s in range(0, pixels):
        subprocess.call([argv[1], "tilt", str(s), str(pixels), str(psi), root + ".ppm", root + ".obj", root + "-tilt-base-" + str(s) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-tilt-base-%d.ppm", "-framerate", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + ".mp4"])
    elif(argv[2] == "btilt"):
      for s in range(0, pixels * 2):
        subprocess.call([argv[1], "tilt", "1", "4", str((s - pixels) / float(pixels) * psi), root + ".ppm", root + ".obj", root + "-btilt-base-" + str(s) + ".ppm"])
        subprocess.call(["cp", root + "-btilt-base-" + str(s) + ".ppm", root + "-btilt-base-" + str(pixels * 4 - s - 1) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt-base-%d.ppm", "-framerate", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-b.mp4"])
    elif(argv[2] == "flicker"):
      for s in range(0, pixels):
        subprocess.call([argv[1], "obj", "4", "1", str(s / float(pixels) * zratio), "0", root + "-bump.ppm", root + "-flicker.obj"])
        subprocess.call([argv[1], "tilt", "1", "4", str(psi), root + ".ppm", root + "-flicker.obj", root + "-flicker-base-" + str(s) + "-L.ppm"])
        subprocess.call([argv[1], "tilt", "3", "4", str(psi), root + ".ppm", root + "-flicker.obj", root + "-flicker-base-" + str(s) + "-R.ppm"])
        subprocess.call(["montage", root + "-flicker-base-" + str(s) + "-R.ppm", root + "-flicker-base-" + str(s) + "-L.ppm", "-geometry", "100%x100%", root + "-flicker-base-" + str(s) + ".png"])
        subprocess.call(["cp", root + "-flicker-base-" + str(s) + ".png", root + "-flicker-base-" + str(pixels * 2 - s - 1) + ".png"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-flicker-base-%d.png", "-framerate", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root + "-flicker.mp4"])
    elif(argv[2] == "sbox"):
      for s in range(0, pixels):
        subprocess.call([argv[1], argv[2], str(int(s - (pixels + 1) / 2.)), str(pixels * 2), str(zratio), root + ".ppm", root + ".obj", root + "-sbox-" + str(pixels - s) + ".ppm"])
    elif(argv[2] == "sboxb2w"):
      subprocess.call([argv[1], "b2w", root + "-sbox-1.ppm", root + "-sbox-bw-1.ppm", "1"])
      for s in range(1, pixels):
        subprocess.call([argv[1], "b2wd", root + "-sbox-" + str(s + 1) + ".ppm", root + "-sbox-bw-" + str(s + 1) + ".ppm", root + "-sbox-" + str(s) + ".ppm"])
    elif(argv[2] == "demosaic"):
      subprocess.call(["convert", line, "-resize", str(int(10000. / pixels) / 100.) + "%", "-compress", "none", root + "-demosaic.ppm"])
      subprocess.call(["python2", argv[0], argv[1], "enlarge", str(pixels), root + "-demosaic.ppm"])
    elif(argv[2] == "prep"):
      subprocess.call(["convert", line, "-resize", str(pixels) + "@>", root + "-prep.png"])
    elif(argv[2] == "prepsq"):
      subprocess.call(["convert", line, "-resize", str(pixels) + "x" + str(pixels) + "!", root + "-prepsq.png"])
    elif(argv[2] == "rot"):
      subprocess.call(["convert", line, "-rotate", "90", root + "-rot.png"])
    elif(argv[2] == "mask"):
      subprocess.call(["convert", root + "-mask.png", "-compress", "none", root + "-mask.ppm"])
    elif(argv[2] == "mask0"):
      subprocess.call(["convert", root + ".ppm", "-fill", "black", "-colorize", "100", root + "-mask.png"])
    elif(argv[2] == "nurie"):
      subprocess.call(["convert", root + ".ppm", "-modulate", "50", root + "-bump.ppm", "-compose", "softlight", "-composite", "-equalize", root + "-nurie.png"])
    elif(argv[2] == "illust"):
      subprocess.call([argv[1], "reshape", str(pixels), root + ".ppm", root + "-bump.ppm", root + "-illust.ppm", str(1. / pow(pixels, .5))])

