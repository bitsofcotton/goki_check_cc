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
  subprocess.call(["mv", "predg.ppm", curdir + ".ppm"])
  subprocess.call(["mv", "predgw.ppm", curdir + "w.ppm"])
  subprocess.call(["mv", "predgw4.ppm", curdir + "w4.ppm"])
elif(argv[2] == "pred" or argv[2] == "predg"):
  cmd = ["python3", argv[0], argv[4], "egeppred", argv[1], argv[5]]
  cmd.extend(argv[6:])
  subprocess.call(cmd)
  cmd = argv[0:2]
  if(argv[2][- 1] == "g"):
    cmd.append("wgpredg")
  else:
    cmd.append("wgpred")
  cmd.extend(argv[3:])
  subprocess.call(cmd)
elif(argv[2] == "egeppred"):
  list0 = argv[5:]
  loopb = 0
  while(True):
    score = []
    for t in range(0, len(list0) - 4):
      cmd = [argv[1], "t"]
      cmd.extend(list0[t:t + 4])
      score.append(abs(float(subprocess.check_output(cmd, encoding="utf-8").split(",")[0])) )
    score = sorted(score)[int(len(score) / 6.)]
    list = [list0[0], list0[1], list0[2] ]
    for t in range(3, len(list0)):
      cmd = [argv[1], "t"]
      cmd.extend(list[- 3:])
      cmd.append(list0[t])
      lscore = abs(float(subprocess.check_output(cmd, encoding="utf-8").split(",")[0]))
      if(score < lscore):
        list.append(list0[t])
    list.reverse()
    rlist = []
    for t in range(max(0, len(list) - 7), len(list)):
      cmd = [argv[1], "t"]
      cmd.extend(list[t:t + 4])
      lscore = abs(float(subprocess.check_output(cmd, encoding="utf-8").split(",")[0]))
      if(score < lscore):
        rlist.append(list[t + 3])
    list = list[:- 4]
    list.extend(rlist)
    list.reverse()
    if(len(list) <= len(argv[5:]) / 2 or len(list) == loopb): break
    print(len(list))
    list0 = list
    loopb = len(list)
  cmd = [argv[1], "c"]
  cmd.extend(list)
  subprocess.check_output(cmd)
  cmd = ["python3", argv[0], argv[3], "bit", str(int(argv[4]))]
  for idx in range(0, len(list)):
    cmd.append(list[idx] + "-c3.ppm")
  subprocess.check_output(cmd)
  cmd = [argv[1], "p"]
  for idx in range(0, len(list)):
    cmd.append(list[idx] + "-c3-bit.ppm")
  subprocess.check_output(cmd)
  subprocess.check_output(["python3", argv[0], argv[3], "move"])
  curdir = os.path.basename(os.getcwd())
  subprocess.call([argv[3], "bit", curdir + ".ppm", curdir + "-bit.ppm", str(- int(argv[4])), "1"])
elif(argv[2] == "wgpred" or argv[2] == "wgpredg"):
  if(argv[2][- 1] == "g"):
    pxs = int(len(argv[6:]) / int(argv[5]))
    ext = "pgm"
  else:
    pxs = int(len(argv[6:]) / float(int(argv[5])) / 3)
    ext = "ppm"
  for l in argv[6:]:
    subprocess.call(["convert", l, "-resize", str(pxs) + "@^", "-compress", "none", l + "-wg." + ext])
    subprocess.call([argv[1], "bit", l + "-wg." + ext, l + "-wg-bit." + ext, str(int(argv[5])), "0"])
  list0 = []
  for l in argv[6:]:
    list0.append(l + "-wg-bit." + ext)
  subprocess.call(["sh", "-c", argv[3] + " + " + " ".join(list0) + " > wgL.txt"])
  subprocess.call(["sh", "-c", argv[3] + " - " + " ".join(list0) + " < wgL.txt"])
  list4 = []
  listw = []
  for l in list0:
    list4.append(l + "-4.ppm")
    listw.append(l)
    listw.append(l + "-4.ppm")
  subprocess.call(["sh", "-c", argv[4] + " p " + " ".join(list4)])
  subprocess.call(["mv", "predg.ppm", "predgw4.ppm"])
  subprocess.call(["sh", "-c", argv[4] + " w " + " ".join(listw)])
  subprocess.call(["python3", argv[0], argv[2], "move"])
  curdir = os.path.basename(os.getcwd())
  subprocess.call([argv[1], "bit", curdir + "w.ppm", curdir + "w-bit.ppm", str(- int(argv[5])), "1"])
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
else:
  for line in argv[3:]:
    try:
      pixels = int(line)
      continue
    except:
      root, ext = os.path.splitext(line)
    if(ext != ".ppm" and argv[2] != "prep" and argv[2] != "prepsq"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
    if(argv[2] == "represent" or argv[2] == "flarge" or argv[2] == "blink" or argv[2] == "enlarge" or argv[2] == "shrink" or argv[2] == "sharpen" or argv[2] == "limit" or argv[2] == "bit" or argv[2] == "nbit" or argv[2] == "rgb2xyz" or argv[2] == "xyz2rgb"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", str(pixels), str(rot)])
    elif(argv[2] == "bump" or argv[2] == "blur" or argv[2] == "collect"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", "1", str(pixels)])
    elif(argv[2] == "obj"):
      subprocess.call([argv[1], argv[2], root + "-bump.ppm", root + ".obj"])
    elif(argv[2] == "jps"):
      subprocess.call([argv[1], "tilt",  "1", "4", str(psi / 2.), root + ".ppm", root + "-bump.ppm", root + "-L.ppm"])
      subprocess.call([argv[1], "tilt", "-1", "4", str(psi / 2.), root + ".ppm", root + "-bump.ppm", root + "-R.ppm"])
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
    elif(argv[2] == "sbit"):
      subprocess.call(["convert", line, "-gravity", "northwest", "-crop", "100%x100%+0+1", root + "-sbit.png"])

