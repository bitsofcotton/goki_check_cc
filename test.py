#! /usr/bin/env python3

import os
import sys
import subprocess

argv   = sys.argv
pixels = 4
psi    = 1. / 60.
rot    = 2
gather = "5%"

if(len(argv) < 4):
  print("no much argments.")
elif(argv[1] == "senganext"):
  cmd = ["python3", argv[0], argv[2] , "diffraw", "8"]
  cmd.extend(argv[5:])
  subprocess.call(cmd)
  cmd = [argv[3]]
  for line in argv[5:]:
    root, ext = os.path.splitext(line)
    cmd.append(root + "-diffraw.ppm")
  subprocess.call(cmd)
  cmd[0] = argv[4]
  subprocess.call(cmd)
elif(argv[2] == "match"):
  root0, ext0 = os.path.splitext(argv[3])
  root1, ext1 = os.path.splitext(argv[4])
  vboxd = vboxs = nsub = nemph = 4
  if(len(argv) > 5):
    vboxd = int(argv[5])
  if(len(argv) > 6):
    vboxs = int(argv[6])
  if(len(argv) > 7):
    nsub = int(argv[7])
  if(len(argv) > 8):
    nemph = int(argv[8])
  subprocess.call([argv[1], argv[2], str(nsub), str(nemph), str(vboxd), str(vboxs), root0 + ".ppm", root1 + ".ppm", root0 + "-bump.ppm", root1 + "-bump.ppm", "match-" + root0 + "-" + root1])
  subprocess.call(["ffmpeg", "-loop", "1", "-i", "match-" + root0 + "-" + root1 + "-%d-" + str(nemph) + ".ppm", "-framerate", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "12", root0 + "-" + root1 + ".mp4"])
elif(argv[2] == "cat" or argv[2] == "catr" or argv[2] == "catb" or argv[2] == "catbr"):
  cmd = [argv[1], argv[2]]
  if(argv[2] == "cat" or argv[2] == "catb"):
    cmd[1] = "catr"
    cmd.append("dummy")
  elif(argv[2] == "catr" or argv[2] == "catbr"):
    cmd[1] = "cat"
    cmd.append("dummy")
  for s in argv[3:]:
    r, e = os.path.splitext(s)
    if(e != ".ppm"):
      subprocess.call(["convert", s, "-compress", "none", r + ".ppm"])
    if(argv[2] == "catb"):
      cmd.append(r + "-bump.ppm")
    elif(argv[2] == "catr"):
      cmd.append(r + "-represent.ppm")
    elif(argv[2] == "catbr"):
      cmd.append(r + "-bump-represent.ppm")
    else:
      cmd.append(r + ".ppm")
  subprocess.call(cmd)
elif(argv[2] == "seinsq" or argv[2] == "seinpdf"):
  files = []
  if(argv[2] == "seinsq"):
    exts = argv[4:]
  else:
    exts = ["pdf"]
  for root, dirs, filesw in os.walk(argv[3]):
    for f in filesw:
      r, e = os.path.splitext(f)
      for s in exts:
        if(e.lower() == "." + s.lower()):
          files.append(os.path.join(root, f))
          break
  ex0 = len(files)
  try:
    ex0 = min(ex0, int(argv[- 1]))
  except:
    pass
  ex = len(str(ex0))
  pixels = int(pow(float(ex0), 2. / 3.))
  if(argv[2] == "seinsq"):
    for t in range(0, ex0):
      subprocess.call(["convert", files[t], "-resize", str(pixels) + "x" + str(pixels) + "!", "-compress", "none", "sein-" + str(t).zfill(ex) + ".ppm"])
  else:
    for t in range(0, ex0):
      subprocess.call(["pdftopng", files[t], "seinpdf-" + str(t).zfill(ex)])
elif(argv[2] == "tilecat" or argv[2] == "tilecatb" or argv[2] == "tilecatr" or argv[2] == "tilecatbr"):
  pixels = int(argv[3])
  cmd = ["montage"]
  t   = 0
  for line in sys.stdin:
    if(len(line) <= 2):
      cmd.extend(["-tile", str(pixels) + "x" + str(pixels), "-geometry", "+0+0", "tilecat-" + str(t) + ".png"])
      subprocess.call(cmd)
      cmd = ["montage"]
      t  += 1
    elif(argv[2] == "tilecatb"):
      cmd.append(line[:- 1 - len("-bump.ppm")] + ".ppm")
    elif(argv[2] == "tilecatr"):
      cmd.append(line[:- 1 - len("-represent.ppm")] + ".ppm")
    elif(argv[2] == "tilecatbr"):
      cmd.append(line[:- 1 - len("-bump-represent.ppm")] + ".ppm")
    else:
      cmd.append(line[:- 1])
  cmd.extend(["-tile", str(pixels) + "x" + str(pixels), "-geometry", "+0+0", "tilecat-" + str(t) + ".png"])
  subprocess.call(cmd)
elif(argv[2] == "tile"):
  idx = 3
  try:
    pixels = int(argv[3])
    idx += 1
  except:
    pass
  for t in range(0, int((pow(len(argv) - idx, .5) + pixels) / pixels)):
    cmd = ["montage"]
    cmd.extend(argv[t * pixels * pixels + idx: min((t + 1) * pixels * pixels + idx, len(argv))])
    cmd.extend(["-tile", str(pixels) + "x" + str(pixels), "-geometry", "+0+0", "tile-" + str(t) + ".png"])
    subprocess.call(cmd)
elif(argv[2] == "apply"):
  sz = 1
  with open(argv[1] + ".txt") as f:
    sz = int(f.readline())
  list = []
  for line in argv[3:]:
    subprocess.call(["convert", line, "-resize", str(sz) + "x" + str(sz) + "!", "-compress", "none", line + "-genl.ppm"])
    list.append(line + "-genl.ppm")
  subprocess.call(["sh", "-c", argv[1] + " - " + " ".join(list) + " < " + argv[1] + ".txt"])
elif(argv[2] == "denlarge"):
  for line in argv[3:]:
    subprocess.call(["mogrify", "-resize", gather, "-compress", "none", line])
    subprocess.call([argv[1], "enlarge", line, line + "-enl.ppm", "2", "8"])
    subprocess.call(["python3", argv[0], argv[1], "penlarge", line + "-enl.ppm"])
elif(argv[2] == "getcontext"):
  cmd = [argv[1], "+"]
  cmd.extend(argv[3:])
  subprocess.call(["sh", "-c", " ".join(cmd) + " > " + argv[1] + ".txt"])
elif(argv[2] == "i2i"):
  idx = 3
  try:
    pixels = int(argv[3])
    idx += 1
  except:
    pass
  for line in argv[idx:]:
    root, ext = os.path.splitext(line)
    if(ext != ".ppm" and argv[2] != "prep" and argv[2] != "prepsq"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
  for linex in argv[idx:]:
    rootx, extx = os.path.splitext(linex)
    for liney in argv[idx:]:
      if(linex == liney): continue
      rooty, exty = os.path.splitext(liney)
      subprocess.call([argv[1], "recolor3", str(pixels), rooty + ".ppm", rootx + ".ppm", rooty + "-" + rootx + "-i2i0.ppm"])
      subprocess.call([argv[1], "recolor",  str(pixels), rootx + ".ppm", rooty + ".ppm", rooty + "-" + rootx + "-i2i1.ppm", "2.5"])
      subprocess.call([argv[1], "recolor3", str(pixels), rooty + "-" + rootx + "-i2i1.ppm", rooty + "-" + rootx + "-i2i0.ppm", rooty + "-" + rootx + "-i2i.ppm"])
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
    elif(argv[2] == "penlarge"):
      subprocess.call(["convert", line, "-resize", "200%", "-compress", "none", root + "-penl0.ppm"])
      subprocess.call([argv[1], "sharpen", root + "-penl0.ppm", root + "-penl1.ppm", str(pixels)])
      subprocess.call(["convert", root + "-penl1.ppm", "-resize", "50%", root + "-penl0.ppm", "-compose", "minus", "-composite", "-modulate", "50", "-negate", root + "-penl0.ppm", "-average", "-normalize", root + "-penl.png"])
    elif(argv[2] == "1to1enl"):
      subprocess.call(["cp", root + ".ppm", root + "-sharpen.ppm"])
      for t in range(0, pixels):
        subprocess.call(["cp", root + "-sharpen.ppm", root + "-sharpen0.ppm"])
        subprocess.call([argv[1], "sharpen", root + "-sharpen0.ppm", root + "-sharpen1.ppm", str(pixels)])
        subprocess.call([argv[1], "diffraw", root + "-sharpen1.ppm", root + "-sharpen2.ppm", str(pixels), str(rot)])
        subprocess.call(["convert", root + "-sharpen2.ppm", "-equalize", "-compress", "none", root + "-sharpen.ppm"])
    elif(argv[2] == "sharpen"):
      subprocess.call(["cp", root + ".ppm", root + "-sharpen.ppm"])
      for t in range(0, pixels):
        subprocess.call([argv[1], "sharpen", root + "-sharpen.ppm", root + "-sharpen-work.ppm", str(pixels)])
        subprocess.call(["convert", root + "-sharpen-work.ppm", "-resize", "50%", "-compress", "none", root + "-sharpen.ppm"])
      subprocess.call(["convert", root + ".ppm", root + "-sharpen.ppm", "-compose", "minus", "-composite", "-negate", "-normalize", "+contrast", root + "-sharpen-minus-normalize.png"])
    elif(argv[2] == "diffraw"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", str(pixels), str(pixels)])
    elif(argv[2] == "bump" or argv[2] == "represent" or argv[2] == "collect" or argv[2] == "flarge" or argv[2] == "blink" or argv[2] == "diffraw" or argv[2] == "enlarge"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", str(pixels), str(rot)])
    elif(argv[2] == "obj"):
      subprocess.call(["convert", root + "-bump0.ppm", "-resize", gather, "-normalize", "-compress", "none", root + "-bump1.ppm"])
      subprocess.call([argv[1], "obj", str(rot), root + "-bump1.ppm", root + ".obj"])
      with open(root + "-bump0.ppm") as f:
        f.readline()
        a = f.readline().split(" ")
        w = int(a[0])
        h = int(a[1])
      subprocess.call(["convert", root + "-bump1.ppm", "-resize", str(w) + "x" + str(h) + "!", "-compress", "none", root + "-bump.ppm"])
    elif(argv[2] == "jps"):
      subprocess.call([argv[1], "tilt", "1", "4", str(  psi), root + ".ppm", root + "-bump.ppm", root + "-L.ppm"])
      subprocess.call([argv[1], "tilt", "1", "4", str(- psi), root + ".ppm", root + "-bump.ppm", root + "-R.ppm"])
      subprocess.call(["convert", root + "-R.ppm", "-define", "trim:edges=east,west", "-trim", root + "-R.png"])
      subprocess.call(["convert", root + "-L.ppm", "-define", "trim:edges=east,west", "-trim", root + "-L.png"])
      subprocess.call(["montage", root + "-R.png", root + "-L.png", "-geometry", "100%x100%", root + "-stereo.jps"])
      subprocess.call(["montage", root + "-L.png", root + "-R.png", "-geometry", "100%x100%", root + "-stereoR.jps"])
      subprocess.call(["montage", root + "-stereo.jps", "-geometry", "100%x100%", root + "-stereo.png"])
      subprocess.call(["montage", root + "-stereoR.jps", "-geometry", "100%x100%", root + "-stereoR.png"])
      subprocess.call(["cp", root + "-L.ppm", root + "-jps-0.ppm"])
      subprocess.call(["cp", root + "-R.ppm", root + "-jps-1.ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-jps-%d.ppm", "-framerate", "20", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "60", root + "-LR.mp4"])
    elif(argv[2] == "btilt"):
      w = h = 0
      with open(root + "-bump.ppm") as f:
        f.readline()
        a = f.readline().split(" ")
        w = int(a[0])
        h = int(a[1])
      for s in range(0, pixels * 2):
        subprocess.call([argv[1], "tilt", "1", "4", str((s - pixels) / float(pixels) * psi), root + ".ppm", root + "-bump.ppm", root + "-btilt-base-" + str(s) + ".ppm"])
        subprocess.call(["mogrify", "-define", "trim:edges=east,west", "-trim", "-resize", str(w) + "x" + str(h) + "!", "-compress", "none", root + "-btilt-base-" + str(s) + ".ppm"])
        subprocess.call(["cp", root + "-btilt-base-" + str(s) + ".ppm", root + "-btilt-base-" + str(pixels * 4 - s - 1) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-btilt-base-%d.ppm", "-framerate", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "60", root + ".mp4"])
    elif(argv[2] == "sbox"):
      for s in range(0, pixels):
        subprocess.call([argv[1], argv[2], str(int(s - (pixels + 1) / 2.)), str(pixels * 2), "0", root + ".ppm", root + "-bump.ppm", root + "-sbox-" + str(pixels - s) + ".ppm"])
    elif(argv[2] == "sboxb2w"):
      subprocess.call([argv[1], "b2w", root + "-sbox-1.ppm", root + "-sbox-bw-1.ppm", "1"])
      for s in range(1, pixels):
        subprocess.call([argv[1], "b2wd", root + "-sbox-" + str(s + 1) + ".ppm", root + "-sbox-bw-" + str(s + 1) + ".ppm", root + "-sbox-" + str(s) + ".ppm"])
    elif(argv[2] == "prep"):
      subprocess.call(["convert", line, "-resize", str(pixels) + "x>", "-resize", "x" + str(pixels) + ">", root + "-prep.png"])
    elif(argv[2] == "prepsq"):
      subprocess.call(["convert", line, "-resize", str(pixels) + "x" + str(pixels) + "!", root + "-prepsq.png"])
    elif(argv[2] == "equalize"):
      subprocess.call(["convert", line, "-equalize", root + "-equalize.png"])
    elif(argv[2] == "nurie"):
      subprocess.call(["convert", root + ".ppm", "-modulate", "50", root + "-bump.ppm", "-compose", "softlight", "-composite", "-equalize", root + "-nurie.png"])
    elif(argv[2] == "illust"):
      subprocess.call([argv[1], "reshape", str(pixels), root + ".ppm", root + "-bump.ppm", root + "-illust.ppm", str(pow(float(pixels), - .5))])
    elif(argv[2] == "edge"):
      subprocess.call(["convert", root + "-collect.ppm", "-type", "GrayScale", "-negate", "-compress", "none", root + "-collect-negate.ppm"])
      subprocess.call(["convert", line, root + "-collect-negate.ppm", "-compose", "multiply", "-composite", root + "-edge.png"])
    elif(argv[2] == "gray"):
      subprocess.call(["convert", line, "-colorspace", "Gray", "-separate", "-average", root + "-gray.png"])

