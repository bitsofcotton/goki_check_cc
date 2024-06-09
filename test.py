#! /usr/bin/env python3

import os
import sys
import subprocess

argv   = sys.argv
pixels = 4
psi    = 1. / 3.
rot    = 0

if(len(argv) < 4 and (argv[2] != "move")):
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
  curdir = os.path.basename(os.getcwd())
  s = 0
  while(True):
    try:
      b = subprocess.call(["mv", "predg-backward-" + str(s) + ".ppm", curdir + "-b" + str(s) + ".ppm"])
      f = subprocess.call(["mv", "predg-forward-" + str(s) + ".ppm", curdir + "-f" + str(s) + ".ppm"])
      if(b != 0 or f != 0):
        exit(0)
        break
    except:
      exit(0)
      break
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
elif(argv[2] == "predC"):
  upix = upixr = sz = sz2 = w = h = 0
  for argn in argv[3:]:
    with open(argn) as f:
      f.readline()
      a = f.readline().split(" ")
      w = int(a[0])
      h = int(a[1])
    upix   = pow(19683 / 8., .5)
    upixr  = upix / pow(w * h, .5)
    upixr /= 3.
    sz     = str(int(upixr * w)) + "x" + str(int(upixr * h))
    sz2    = str(int(upixr * w) * int(upixr * w)) + "x" + str(int(upixr * h) * int(upixr * h))
    subprocess.call(["convert", argn, "+sigmoidal-contrast", "7.5", "-filter", "LanczosRadius", "-distort", "Resize", sz, "-sigmoidal-contrast", "7.5", "-compress", "none", argn + "-root.ppm"])
    for x in range(0, int(upixr * w) ):
      for y in range(0, int(upixr * h) ):
        subprocess.call(["convert", argn, "+sigmoidal-contrast", "7.5", "-filter", "LanczosRadius", "-distort", "Resize", sz2, "-sigmoidal-contrast", "7.5", "-crop", sz + "+" + str(int(upixr * w) * x) + "+" + str(int(upixr * h) * y), "-compress", "none", argn + "-" + str(x) + "-" + str(y) + ".ppm"])
  cmd = ["predgmp"]
  for argn in argv[3:]:
    cmd.append(argn + "-root.ppm")
  subprocess.call(cmd)
  s = 0
  loop = True
  while(loop):
    t = 1
    while(True):
      try:
        b  = subprocess.call(["mogrify", "-despeckle", "-compress", "none", "predg-backward-" + str(s) + "-" + str(t) + ".ppm"])
        b2 = subprocess.call(["gokibinmp", "enlarge", "predg-backward-" + str(s) + "-" + str(t) + ".ppm", "predg-backward-root-" + str(s) + "-" + str(t) + ".ppm", "2", "12"])
        f  = subprocess.call(["mogrify", "-despeckle", "-compress", "none", "predg-forward-" + str(s) + "-" + str(t) + ".ppm"])
        f2 = subprocess.call(["gokibinmp", "enlarge", "predg-forward-" + str(s) + "-" + str(t) + ".ppm", "predg-forward-root-" + str(s) + "-" + str(t) + ".ppm", "2", "12"])
        if(b != 0 or f != 0):
          if(t == 1): loop = False
          break
        t += 1
      except:
        if(t == 1): loop = False
        break
    s += 1
  for x in range(0, int(upixr * w) ):
    for y in range(0, int(upixr * h) ):
      cmd = ["predgmp"]
      for argn in argv[3:]:
        cmd.append(argn + "-" + str(x) + "-" + str(y) + ".ppm")
      subprocess.call(cmd)
      s = 0
      loop = True
      while(loop):
        t = 1
        while(True):
          try:
            b  = subprocess.call(["mogrify", "-despeckle", "-compress", "none", "predg-backward-" + str(s) + "-" + str(t) + ".ppm"])
            b2 = subprocess.call(["gokibinmp", "enlarge", "predg-backward-" + str(s) + "-" + str(t) + ".ppm", "predg-backward-" + str(x) + "-" + str(y) + "-" + str(s) + "-" + str(t) + ".ppm", "2", "12"])
            f  = subprocess.call(["mogrify", "-despeckle", "-compress", "none", "predg-forward-" + str(s) + "-" + str(t) + ".ppm"])
            f2 = subprocess.call(["gokibinmp", "enlarge", "predg-forward-" + str(s) + "-" + str(t) + ".ppm", "predg-forward-" + str(x) + "-" + str(y) + "-" + str(s) + "-" + str(t) + ".ppm", "2", "12"])
            if(b != 0 or b2 != 0 or f != 0 or f2 != 0):
              if(t == 1): loop = False
              break
            t += 1
          except:
            if(t == 1): loop = False
            break
        s += 1
elif(argv[2] == "predCbond"):
  # XXX: LD_PRELOAD=... allocator python calls ksh causes imagemagick
  #      noerror exit.
  upix = upixr = sz = sz2 = w = h = 0
  for argn in argv[3:]:
    with open(argn) as f:
      f.readline()
      a = f.readline().split(" ")
      w = int(a[0])
      h = int(a[1])
    upix   = pow(19683 / 8., .5)
    upixr  = upix / pow(w * h, .5)
    upixr /= 3.
    sz     = str(int(upixr * w)) + "x" + str(int(upixr * h))
    sz2    = str(int(upixr * w) * int(upixr * w)) + "x" + str(int(upixr * h) * int(upixr * h))
  curdir = os.path.basename(os.getcwd())
  s = 0
  loop = True
  while(loop):
    t = 1
    while(True):
      try:
        cmd = ["montage"]
        for y in range(0, int(upixr * h)):
          for x in range(0, int(upixr * w)):
            cmd.append("predg-backward-" + str(x) + "-" + str(y) + "-" + str(s) + "-" + str(t) + ".ppm")
        cmd.extend(["-tile", str(int(upixr * w)) + "x" + str(int(upixr * h)), "-geometry", "+0+0", "predg-backward-extend-" + str(s) + "-" + str(t) + ".png"])
        bb = subprocess.call(cmd)
        b  = subprocess.call(["convert", "predg-backward-root-" + str(s) + "-" + str(t) + ".ppm", "-resize", sz2, "predg-backward-extend-" + str(s) + "-" + str(t) + ".png", "-resize", sz2, "-compose", "multiply", "-composite", curdir + "-b" + str(s) + "-" + str(t) + ".png"])
        cmd = ["montage"]
        for y in range(0, int(upixr * h)):
          for x in range(0, int(upixr * w)):
            cmd.append("predg-forward-" + str(x) + "-" + str(y) + "-" + str(s) + "-" + str(t) + ".ppm")
        cmd.extend(["-tile", str(int(upixr * w)) + "x" + str(int(upixr * h)), "-geometry", "+0+0", "predg-forward-extend-" + str(s) + "-" + str(t) + ".png"])
        ff = subprocess.call(cmd)
        f  = subprocess.call(["convert", "predg-forward-root-" + str(s) + "-" + str(t) + ".ppm", "-resize", sz2, "predg-forward-extend-" + str(s) + "-" + str(t) + ".png", "-resize", sz2, "-compose", "multiply", "-composite", curdir + "-f" + str(s) + "-" + str(t) + ".png"])
        if(b != 0 or bb != 0 or f != 0 or ff != 0):
          if(t == 1): loop = False
          break
        t += 1
      except:
        if(t == 1): loop = False
        break
    s += 1
else:
  for line in argv[3:]:
    try:
      pixels = int(line)
      continue
    except:
      root, ext = os.path.splitext(line)
    if(ext != ".ppm" and argv[2] != "prep" and argv[2] != "prepsq"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    if(argv[2] == "represent" or argv[2] == "collect" or argv[2] == "flarge" or argv[2] == "blink" or argv[2] == "enlarge" or argv[2] == "shrink" or argv[2] == "sharpen" or argv[2] == "limit" or argv[2] == "bit" or argv[2] == "rgb2xyz" or argv[2] == "xyz2rgb"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", str(pixels), str(rot)])
    elif(argv[2] == "bump"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + "0.ppm", str(pixels), str(rot)])
      subprocess.call(["python3", argv[0], argv[1], "cleanq", root + "-bump0.ppm"])
      subprocess.call(["convert", root + "-bump0-cleanq.png", "-format", "ppm", "-compress", "none", root + "-bump.ppm"])
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
    elif(argv[2] == "nurie"):
      subprocess.call(["convert", root + ".ppm", "-modulate", "50", root + "-bump.ppm", "-compose", "softlight", "-composite", "-equalize", root + "-nurie.png"])
    elif(argv[2] == "illust"):
      subprocess.call([argv[1], "reshape", str(pixels), root + ".ppm", root + "-bump.ppm", root + "-illust.ppm", str(pow(float(pixels), - .5))])
    elif(argv[2] == "gray"):
      subprocess.call(["convert", line, "-colorspace", "Gray", "-separate", "-average", root + "-gray.png"])
    elif(argv[2] == "cleanl" or argv[2] == "cleanq" or argv[2] == "cleanc" or argv[2] == "cleanC"):
      w = h = 0
      with open(root + ".ppm") as f:
        f.readline()
        a = f.readline().split(" ")
        w = int(a[0])
        h = int(a[1])
      # N.B.
      # With ddpmopt/README.md, we have 20kbit upper limit on function entropy.
      # In practical, 2 bit from MSB expecting input is enough.
      # However, we expect 8 bit color for each pixel.
      upix  = pow(19683 / 8., .5)
      upixr = upix / pow(w * w + h * h, .5)
      if(argv[2][- 1] == "C"): upixr /= pow(3. , .5)
      sz    = str(int(upixr * w) + 1) + "x" + str(int(upixr * h) + 1)
      # N.B.
      # refering en.wikipedia.org/wiki/Condorcet's-jury-theorem
      # n ~ 11 to get .95 from 2/3 probability.
      # This is enough if it's from probability based concerns.
      # However, in deterministic meaning, we might use cj == upix as well.
      cj     = int(pow(11, .5) + 1)
      # Thanks:
      # R.B. https://www.imagemagick.org/Usage/filter/nicolas/#downsample
      # From googling, via https://qiita.com/yoya/items/b1590de289b623f18639 .
      if(argv[2][- 1] == "q"):
        subprocess.call(["convert", line, "-despeckle", "-filter", "LanczosRadius", "-distort", "Resize", sz, "-despeckle", "+sigmoidal-contrast", "7.5", "-filter", "LanczosRadius", "-distort", "Resize", str(w) + "x" + str(h) + "!", "-sigmoidal-contrast", "7.5", root + "-cleanq.png"])
      elif(argv[2][- 1] == "l"):
        subprocess.call(["convert", line, "+sigmoidal-contrast", "7.5", "-filter", "LanczosRadius", "-distort", "Resize", str(int(cj * upixr * w + 1)) + "x" + str(int(cj * upixr * h + 1)), "-sigmoidal-contrast", "7.5", root + "-cleanl.png"])
      elif(argv[2][- 1] == "c" or argv[2][- 1] == "C"):
        subprocess.call(["convert", line, "-despeckle", "-filter", "LanczosRadius", "-distort", "Resize", sz, "-despeckle", "-compress", "none", root + "-" + argv[2] + ".ppm"])
        subprocess.call([argv[1], "enlarge", root + "-" + argv[2] + ".ppm", root + "-" + argv[2] + "e.ppm", "2", str(pixels)])
        subprocess.call(["convert", root + "-" + argv[2] + "e.ppm", "-equalize", "-despeckle", root + "-" + argv[2] + "e.png"])
      else:
        print("unknown command: ", argv[2])

