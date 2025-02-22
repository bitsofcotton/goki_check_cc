#! /usr/bin/env python3

import os
import sys
import subprocess

argv   = sys.argv
pixels = 4
psi    = 1. / 3.
rot    = 0

def ifloat(x, offset = 0):
  try:
    return float(x)
  except:
    try:
      b = x.split("*")
      n = e = 0.
      tbl = {'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, \
        '8': 8, '9': 9, 'a': 10, 'b': 11, 'c': 12, 'd': 13, 'e': 14, 'f': 15, \
        ' ': -1, '\t': -1, '\n': -1, '-': 16}
      m = False
      for ff in b[0]:
        if(tbl[ff] < 0): continue
        if(tbl[ff] == 16):
          m = True
          continue
        n *= 16
        n += tbl[ff]
      if(m): n = - n
      m = False
      for ff in b[1][2:]:
        if(tbl[ff] < 0): continue
        if(tbl[ff] == 16):
          m = True
          continue
        e *= 16
        e += tbl[ff]
      if(m): e = - e
      return n * pow(2., e - offset)
    except:
      pass
  return 0.

if(len(argv) < 4 and (argv[2] != "move")):
  print("no much argments.")
elif(argv[2] == "move"):
  import numpy
  f = "{0:0" + str(int(numpy.ceil(numpy.log(float(int(len(argv[3:])))) / numpy.log(10.)))) + "d}"
  for idx in range(0, len(argv[3:])):
    subprocess.call(["mv", argv[idx + 3], argv[1] + f.format(idx) + os.path.splitext(argv[idx + 3])[1] ])
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
elif(argv[2] == "pred" or argv[2] == "pred+" or argv[2] == "pred++" or argv[2] == "predg" or argv[2] == "predg+" or argv[2] == "predg++" or argv[2] == "predq" or argv[2] == "predq+" or argv[2] == "predq++" or argv[2] == "predgq" or argv[2] == "predgq+" or argv[2] == "predgq++"):
  if(argv[2][- 1] == 'q' or argv[2][- 2] == 'q' or argv[2][- 3] == 'q'):
    mode = 1
  else:
    mode = 0
  # N.B. Condorcet jury doesn't work on our sample run.
  if(argv[2][- 2] == '+'):
    # N.B. upper bound for #f(x,y,z) when completely separatable {x, y, z}
    #      however, {x,y,z,w} is too large to exist without internal relations.
    #      we count them as a operation, so we take output each result pixels.
    lsz  = pow(pow(3., pow(3., 3.)), .5)
    # N.B. however, with predictions, only 1/3 will be warranted to get
    #      results as in use.
    lsz /= 3.
  elif(argv[2][- 1] == '+'):
    # N.B. f(x, y) == 0, #f pure function case.
    lsz  = 19683
  else:
    # N.B. context saturation case.
    lsz  = 768
  # N.B. we treat each bit condition to them.
  lsz /= 8.
  if(argv[2][- 1] == 'g' or argv[2][- 2] == 'g' or argv[2][- 3] == 'g' or argv[2][- 4] == 'g'):
    bits = max(1, min(8, int(len(argv[5:]) / 4 / 4)))
    pxs  = int(len(argv[5:]) / bits / 4)
    ext  = "pgm"
  else:
    bits = max(1, min(8, int(len(argv[5:]) / 3 / 4 / 4)))
    pxs  = int(len(argv[5:]) / bits / 3 / 4)
    ext  = "ppm"
    lsz /= 3.
  if(mode == 1):
    lsz  = pow(lsz, .5)
  lsz = int(lsz)
  for f in argv[5:]:
    if(mode == 1):
      subprocess.call(["convert", f, "-resize", str(lsz) + "x" + str(lsz) + "!", "-compress", "none", f + "-wgL." + ext])
    else:
      subprocess.call(["convert", f, "-resize", str(lsz) + "@", "-compress", "none", f + "-wgL." + ext])
  list0 = []
  for f in argv[5:]:
    list0.append(f + "-wgL." + ext)
  list  = []
  loopb = 0
  while(True):
    score = []
    for t in range(0, len(list0) - 4):
      cmd = [argv[4], "t"]
      cmd.extend(list0[t:t + 4])
      score.append(abs(ifloat(subprocess.check_output(cmd, encoding="utf-8").split(",")[0])) )
    # XXX: magic number, ideally, is [1] but too slow to run.
    score = sorted(score)[int(len(score) / 5.)]
    list = [list0[0], list0[1], list0[2] ]
    for t in range(3, len(list0)):
      cmd = [argv[4], "t"]
      cmd.extend(list[- 3:])
      cmd.append(list0[t])
      lscore = abs(ifloat(subprocess.check_output(cmd, encoding="utf-8").split(",")[0]))
      if(score < lscore):
        list.append(list0[t])
    list.reverse()
    rlist = []
    for t in range(max(0, len(list) - 5), len(list)):
      cmd = [argv[4], "t"]
      cmd.extend(list[t:t + 4])
      lscore = abs(ifloat(subprocess.check_output(cmd, encoding="utf-8").split(",")[0]))
      if(score < lscore):
        rlist.append(list[t + 3])
    list = list[:- 4]
    list.extend(rlist)
    list.reverse()
    if(len(list) <= len(argv[5:]) / 2 or len(list) == loopb): break
    print(len(list))
    list0 = list
    loopb = len(list)
  cmd = [argv[4], "c"]
  cmd.extend(list)
  subprocess.check_output(cmd)
  list1 = []
  for f in list:
    list1.append(f + "-c3.ppm")
  list2 = []
  list3 = []
  for f in list1:
    subprocess.call(["convert", f, "-resize", str(pxs) + "@^", "-compress", "none", f + "-wg." + ext])
    subprocess.call([argv[1], "bit", f + "-wg." + ext, f + "-wg-bit." + ext, str(bits), "0"])
    subprocess.call([argv[1], "bit", f, f + "-wgL-bit." + ext, str(bits), "0"])
    list2.append(f + "-wg-bit." + ext)
    list3.append(f + "-wgL-bit." + ext)
  cmd = [argv[4], "p"]
  cmd.extend(list3)
  subprocess.check_output(cmd)
  curdir = os.path.basename(os.getcwd())
  subprocess.call([argv[1], "bit", "predg.ppm", curdir + "-bit.ppm", str(- bits), "1"])
  subprocess.call([argv[1], "nbit", "predg.ppm", curdir + "-nbit.ppm", str(- bits), "1"])
  # XXX: white space, delimiter, should use Popen with pipe.
  subprocess.call(["sh", "-c", argv[3] + " + " + " ".join(list2) + " > wgL.txt"])
  subprocess.call(["sh", "-c", argv[3] + " - " + " ".join(list2) + " < wgL.txt"])
  list4 = [argv[4], "p"]
  listl = [argv[4], "w"]
  listr = [argv[4], "w"]
  for l in list2:
    list4.append(l + "-4.ppm")
  for idx in range(0, len(list3)):
    listl.append(list3[idx])
    listl.append(list2[idx] + "-4.ppm")
  for idx in range(0, len(list3) - 1):
    listr.append(list3[idx])
    listr.append(list2[idx + 1] + "-4.ppm")
  listr.append(list3[- 1])
  listr.append("predgw4.ppm")
  subprocess.call(list4)
  subprocess.call(["mv", "predg.ppm", "predgw4.ppm"])
  subprocess.call(listr)
  subprocess.call([argv[1], "bit", "predgw.ppm", curdir + "w4-bit.ppm", str(- bits), "1"])
  subprocess.call([argv[1], "nbit", "predgw.ppm", curdir + "w4-nbit.ppm", str(- bits), "1"])
  subprocess.call(listl)
  subprocess.call([argv[1], "bit", "predgw.ppm", curdir + "w-bit.ppm", str(- bits), "1"])
  subprocess.call([argv[1], "nbit", "predgw.ppm", curdir + "w-nbit.ppm", str(- bits), "1"])
elif(argv[2] == "crossarg"):
  step = int(argv[3])
  spl  = 0
  for idx in range(4, len(argv)):
    if(argv[idx] == ","):
      spl = idx
      break
  if(spl == 0):
    print("ERR")
  else:
    arg1 = argv[4:spl]
    arg2 = argv[spl + 1:]
    sw   = 1
    for s in range(0, min(len(arg1), len(arg2)) * 2):
      if(s % step == 0):
        sw = - sw
      if(0 < sw):
        print(arg1[int(s / step / 2) * step + (s % step) - min(len(arg1), len(arg2))])
      else:
        print(arg2[int(s / step / 2) * step + (s % step) - min(len(arg1), len(arg2))])
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
    if(argv[2] == "represent" or argv[2] == "flarge" or argv[2] == "blink" or argv[2] == "enlarge" or argv[2] == "shrink" or argv[2] == "sharpen" or argv[2] == "limit" or argv[2] == "bit" or argv[2] == "nbit" or argv[2] == "slide" or argv[2] == "nslide" or argv[2] == "rgb2xyz" or argv[2] == "xyz2rgb"):
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

