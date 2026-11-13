#! /usr/bin/env python3

import os
import sys
import math
import subprocess

argv   = sys.argv
pixels = 4
psi    = 1. / 3.
rot    = 0

if(len(argv) < 3):
  print("no much argments.")
elif(argv[2] == "repeat"):
  # rep == 11 for 2/3 -> .9?
  rep = 3
  n   = int(len(argv[5:]) / 21)
  l   = int(len(argv[5:]) / n)
  # N.B. we crush to last as ind2vd limit which internally *= 2.
  for s in range(0, n):
    cmd = argv[3:5]
    # N.B. shorten input for large output, we cannot match entire meaning
    #      better than average, last subtracted series.
    if(s == 0): cmd.extend(argv[5:][- l:])
    else: cmd.extend(argv[5:- s][::-1][::s + 1][::-1][- l:])
    for t in range(0, rep):
      subprocess.call(cmd)
      subprocess.call(["mv", argv[1], argv[1] + "-" + str(s) + "-" + str(t) + ".ppm"])
elif(argv[2] == "move"):
  t = 0
  for f in argv[3:]:
    root, ext = os.path.splitext(f)
    subprocess.call(["mv", f, argv[1] + str(t).zfill(len(str(len(argv[3:]))) ) + ext])
    t += 1
elif(argv[2] == "seinp" or argv[2] == "seinq" or argv[2] == "seinP" or argv[2] == "seinQ"):
  files = []
  for root, dirs, filesw in os.walk(argv[3]):
    for f in filesw:
      r, e = os.path.splitext(f)
      for s in argv[4:]:
        if(e.lower() == "." + s.lower()):
          files.append(os.path.join(root, f))
          break
  if(argv[2][len("sein")] == 'p' or argv[2][len("sein")] == 'q'):
    sz  = len(files) / 3
    ext = ".ppm"
  else:
    sz  = len(files)
    ext = ".pgm"
  if(argv[2][len("sein")] == 'p' or argv[2][len("sein")] == 'P'):
    rsz = str(int(sz)) + "@"
  else:
    sz  = pow(sz, .5)
    rsz = str(int(sz)) + "x" + str(int(sz)) + "!"
  ex = len(str(len(files)))
  for idx in range(0, len(files)):
    subprocess.call(["convert", files[idx], "-compress", "none", "-resize", rsz, "sein-" + str(idx).zfill(ex) + ext])
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
elif(argv[2] == "tilecat"):
  pixels = int(argv[3])
  cmd = ["montage"]
  t   = 0
  for line in sys.stdin:
    if(len(line) <= 2):
      cmd.extend(["-tile", str(pixels) + "x" + str(pixels), "-geometry", "+0+0", "tilecat-" + str(t) + ".png"])
      subprocess.call(cmd)
      cmd = ["montage"]
      t  += 1
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
elif(argv[2] == "check"):
  # Thanks to claude.ai 2026/11/14
  import numpy as np
  from PIL import Image
  from scipy.stats import rankdata
  B = 4
  rng = np.random.default_rng(0)
  def prep(a):
    H, W = a.shape
    return a[ : H // B * B, : W // B * B]
  def rank01(a):
    a = np.asarray(a, float).ravel()
    return (rankdata(a, method="average") - 0.5) / a.size
  def resid(y):
    return y - y.mean()
  def bsum(a):
    return a.reshape(H, W).reshape(hb, B, wb, B).sum((1, 3))
  def stats(xv, full=False):
    rX = resid(rank01(xv))
    nX = np.linalg.norm(rX)
    c  = bsum(rX * rP)
    sxx = bsum(rX * rX)
    u  = c / (nX * nP)
    ok = (sxx * syy) > 0
    rb = np.zeros_like(c)
    rb[ok] = c[ok] / np.sqrt(sxx[ok] * syy[ok])
    v  = np.array([u.sum(), u[u > 0].sum(), u[u < 0].sum(),
      (rb[ok] ** 2).mean() if ok.any() else 0.0])
    return (v, c, sxx, u) if full else v
  def phase_surrogate(x, rng):
    F = np.fft.fft2(x - x.mean())
    ph = np.angle(np.fft.fft2(rng.standard_normal(x.shape)))
    return np.real(np.fft.ifft2(np.abs(F) * np.exp(1j * ph)))
  def pvals(obs, null, two_sided=False):
    m = null.mean()
    up = (1 + (null >= obs).sum()) / (1 + len(null))
    if two_sided:
      return up, (1 + (np.abs(null - m) >= abs(obs - m)).sum()) / (1 + len(null))
    return up
  x = prep(np.asarray(Image.open(argv[1]).convert("L"), dtype=float))
  p = prep(np.asarray(Image.open(argv[3]).convert("L"), dtype=float))
  H,  W  = x.shape
  hb, wb = H // B, W // B
  k   = 0
  rP  = resid(rank01(p))
  nP  = np.linalg.norm(rP)
  syy = bsum(rP * rP)
  obs, c, sxx, u = stats(x.ravel(), full=True)
  null = np.array([stats(phase_surrogate(x, rng).ravel()) for _ in range(20000)])
  m, s = null[:, k].mean(), null[:, k].std()
  out = {}
  out["rho"] = dict(obs=float(obs[k]), null_mean=float(m), null_sd=float(s),
    z=float((obs[k] - m) / s) if s > 0 else float("nan"))
  out["rho"]["p_one_sided"], out["rho"]["p_two_sided"] = [float(v) for v in
    pvals(obs[k], null[:, k], two_sided=True)]
  idx   = rng.integers(0, c.size, (2000, c.size))
  boots = c.ravel()[idx].sum(1) / np.sqrt(sxx.ravel()[idx].sum(1) * syy.ravel()[idx].sum(1))
  lo, hi = np.percentile(boots, [2.5, 97.5])
  out["ci95"] = (float(lo), float(hi))
  print(out)
elif(argv[2] == "i2i"):
  idx = 3
  try:
    pixels = int(argv[3])
    idx += 1
  except:
    pass
  for line in argv[idx:]:
    root, ext = os.path.splitext(line)
    if(ext != ".ppm"):
      subprocess.call(["convert", line, "-compress", "none", root + ".ppm"])
  for linex in argv[idx:]:
    rootx, extx = os.path.splitext(linex)
    for liney in argv[idx:]:
      if(linex == liney): continue
      rooty, exty = os.path.splitext(liney)
      subprocess.call([argv[1], "recolor3", str(pixels), rooty + ".ppm", rootx + ".ppm", rooty + "-" + rootx + "-i2i0.ppm"])
      subprocess.call([argv[1], "recolor",  str(pixels), rootx + ".ppm", rooty + ".ppm", rooty + "-" + rootx + "-i2i1.ppm", ".00000001"])
      subprocess.call([argv[1], "recolor3", str(pixels), rooty + "-" + rootx + "-i2i1.ppm", rooty + "-" + rootx + "-i2i0.ppm", rooty + "-" + rootx + "-i2i.ppm"])
      subprocess.call([argv[1], "recolor2", str(pixels), rooty + "-" + rootx + "-i2i.ppm", rooty + "-" + rootx + "-i2i-emph.ppm", "-.000000005"])
else:
  for line in argv[3:]:
    try:
      pixels = int(line)
      continue
    except:
      root, ext = os.path.splitext(line)
    if(ext != ".ppm" and ext != ".pgm" and argv[2] != "prep" and argv[2] != "prepsq"):
      subprocess.call(["convert", line, "-compress", "none", "-depth", "16", root + ".ppm"])
    if(argv[2] == "enlarge" or argv[2] == "shrink" or argv[2] == "blink" or argv[2] == "limit" or argv[2] == "bit" or argv[2] == "nbit" or argv[2] == "bump" or argv[2] == "normalize"):
      subprocess.call([argv[1], argv[2], root + ".ppm", root + "-" + argv[2] + ".ppm", str(pixels), str(rot)])
    elif(argv[2] == "collect" or argv[2] == "fcollect" or argv[2] == "sharpen" or argv[2] == "blur" or argv[2] == "abs"):
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
      subprocess.call(["convert", line, "-resize", str(pixels) + "x>", "-resize", "x" + str(pixels) + ">", "-depth", "16", root + "-prep.png"])
    elif(argv[2] == "prepsq"):
      subprocess.call(["convert", line, "-resize", str(pixels) + "x" + str(pixels) + "!", "-depth", "16", root + "-prepsq.png"])
    elif(argv[2] == "nurie"):
      subprocess.call(["convert", root + ".ppm", "-modulate", "50", root + "-bump.ppm", "-compose", "softlight", "-composite", "-equalize", root + "-nurie.png"])
    elif(argv[2] == "illust"):
      subprocess.call([argv[1], "reshape", str(pixels), root + ".ppm", root + "-bump.ppm", root + "-illust.ppm", str(pow(float(pixels), - .5))])

