#! /usr/bin/env python3

import os
import sys
import subprocess

argv   = sys.argv
pixels = 4
psi    = 1. / 3.
rot    = 0

if(len(argv) < 4):
  print("no much argments.")
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
  elif(argv[2] == "catr" or argv[2] == "catbr"):
    cmd[1] = "cat"
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
elif(argv[2] == "move"):
  src = argv[3].split("%")
  dst = argv[4].split("%")
  try:
    subprocess.call(["mv", src[0] + "0" + src[1], dst[0] + "0" + dst[1]])
  except:
    pass
  try:
    for t in range(1, int(argv[5])):
      subprocess.call(["mv", src[0] + str(t) + src[1], dst[0] + str(t) + dst[1]])
  except:
    pass
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
    pixels = argv[3]
    idx += 1
  except:
    pass
  cmd = ["montage"]
  cmd.extend(argv[idx:])
  cmd.extend(["-tile", pixels, "-geometry", "+0+0", "tile.png"])
  subprocess.call(cmd)
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
      subprocess.call([argv[1], "recolor2", str(pixels), rooty + "-" + rootx + "-i2i.ppm", rooty + "-" + rootx + "-i2i--02.ppm", "-.02"])
elif(argv[2] == "pextend"):
  for line in argv[3:]:
    root, ext = os.path.splitext(line)
    subprocess.call(["convert", line, "-compress", "none", root + "-pext.ppm"])
    w = h = 0
    with open(root + "-pext.ppm") as f:
      f.readline()
      a = f.readline().split(" ")
      w = int(a[0])
      h = int(a[1])
    subprocess.call(["mogrify", "-resize", str(w * w) + "x" + str(h) + "!", "-compress", "none", root + "-pext.ppm"])
    subprocess.call([argv[1], root + "-pext.ppm"])
    with open(root + "-pext.ppm") as f:
      f.readline()
      a = f.readline().split(" ")
      h = int(a[1])
    subprocess.call(["mogrify", "-resize", str(w) + "x" + str(h) + "!", "-format", "png", root + "-pext.ppm"])
elif(argv[2] == "pred"):
  w = h = 0
  for line in argv[4:]:
    root, ext = os.path.splitext(line)
    subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    with open(root + ".ppm") as f:
      f.readline()
      a = f.readline().split(" ")
      w0 = int(a[0])
      h0 = int(a[1])
      if(w == 0): w = w0
      if(h == 0): h = h0
      if(w != w0): sys.exit(- 1)
      if(h != h0): sys.exit(- 1)
  r  = pow((pow(len(argv[4:]) / int(argv[3]), 4.) / 3.) / w / h, .5)
  w *= r
  h *= r
  w  = int(w)
  h  = int(h)
  cmd = [argv[1]]
  for line in argv[4:]:
    subprocess.call(["mogrify", "-compress", "none", "-resize", str(w) + "x" + str(h) + "!", root + ".ppm"])
    cmd.append(root + ".ppm")
  subprocess.call(cmd)
else:
  for line in argv[3:]:
    try:
      pixels = int(line)
      continue
    except:
      root, ext = os.path.splitext(line)
    if(ext != ".ppm" and argv[2] != "prep" and argv[2] != "prepsq"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    if(argv[2] == "represent" or argv[2] == "collect" or argv[2] == "flarge" or argv[2] == "blink" or argv[2] == "enlarge" or argv[2] == "sharpen" or argv[2] == "limit" or argv[2] == "bit"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", str(pixels), str(rot)])
    elif(argv[2] == "bump"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + "0.ppm", str(pixels), str(rot)])
      subprocess.call(["python3", argv[0], argv[1], "cleansq", root + "-bump0.ppm"])
      subprocess.call(["convert", root + "-bump0-cleansq.png", "-format", "ppm", "-compress", "none", root + "-bump.ppm"])
    elif(argv[2] == "obj"):
      subprocess.call([argv[1], argv[2], root + "-bump.ppm", root + ".obj"])
    elif(argv[2] == "jps"):
      subprocess.call([argv[1], "tilt",  "1", "4", str(psi), root + ".ppm", root + "-bump.ppm", root + "-L.ppm"])
      subprocess.call([argv[1], "tilt", "-1", "4", str(psi), root + ".ppm", root + "-bump.ppm", root + "-R.ppm"])
      subprocess.call(["montage", root + "-R.ppm", root + "-L.ppm", "-geometry", "100%x100%", root + "-stereo.jps"])
      subprocess.call(["montage", root + "-L.ppm", root + "-R.ppm", "-geometry", "100%x100%", root + "-stereoR.jps"])
      subprocess.call(["montage", root + "-stereo.jps", "-geometry", "100%x100%", root + "-stereo.png"])
      subprocess.call(["montage", root + "-stereoR.jps", "-geometry", "100%x100%", root + "-stereoR.png"])
    elif(argv[2] == "tilt"):
      for s in range(0, pixels * 2):
        subprocess.call([argv[1], argv[2], "1", "4", str((s - pixels) / float(pixels) * psi), root + ".ppm", root + "-bump.ppm", root + "-tilt-" + str(s) + ".ppm"])
        subprocess.call(["cp", root + "-tilt-" + str(s) + ".ppm", root + "-tilt-" + str(pixels * 4 - s - 1) + ".ppm"])
      subprocess.call(["ffmpeg", "-loop", "1", "-i", root + "-tilt-%d.ppm", "-framerate", "6", "-an", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-vcodec", "libx264", "-pix_fmt", "yuv420p", "-t", "60", root + ".mp4"])
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
    elif(argv[2] == "cleans" or argv[2] == "cleansq"):
      w = h = 0
      with open(root + ".ppm") as f:
        f.readline()
        a = f.readline().split(" ")
        w = int(a[0])
        h = int(a[1])
      # Thanks:
      # R.B. https://www.imagemagick.org/Usage/filter/nicolas/#downsample
      # From googling, via https://qiita.com/yoya/items/b1590de289b623f18639 .
      if(argv[2][- 1] == "q"):
        subprocess.call(["convert", line, "-colorspace", "RGB", "-filter", "LanczosRadius", "-distort", "Resize", str(100. * pow(w * h, - .25)) + "%", "-colorspace", "sRGB", "-colorspace", "RGB", "+sigmoidal-contrast", "7.5", "-filter", "LanczosRadius", "-distort", "Resize", str(w) + "x" + str(h) + "!", "-sigmoidal-contrast", "7.5", "-colorspace", "sRGB", root + "-cleansq.png"])
      else:
        subprocess.call(["convert", line, "-colorspace", "RGB", "-filter", "LanczosRadius", "-distort", "Resize", str(100. * pow(w * h, - .25)) + "%", "-colorspace", "sRGB", root + "-cleans.png"])

