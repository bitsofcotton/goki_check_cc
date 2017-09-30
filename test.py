#! /usr/bin/env python

import os
import sys
import glob
import subprocess

argv = sys.argv

if(len(argv) < 2):
  print "not much arguments."
  exit(- 1)

for idx in range(2, len(argv)):
  line = argv[idx]
  print line
  root, ext = os.path.splitext(line)
  if(not os.path.exists(root + "-enl1" + ext)):
    subprocess.call([argv[1], "enlargeds", line, root + "-enl0" + ext])
    subprocess.call(["convert", root + "-enl0" + ext, "-resize", "75%", "-unsharp", "10x4+1+0", "-compress", "none", root + "-enl1" + ext])
  if(not os.path.exists(root + "-collect" + ext)):
    subprocess.call([argv[1], "collect", line, root + "-collect" + ext])
  if(not os.path.exists(root + "-bump" + ext)):
    subprocess.call([argv[1], "bump", line, root + "-bump" + ext])
  if(not os.path.exists(root + "-bump-blur" + ext)):
    subprocess.call([argv[1], "bumpobj", line, root + "-bump-blur" + ext, root + ".obj"])
  if(not os.path.exists(root + "-tilt-1" + ext) and not os.path.exists(root + ".mp4")):
    subprocess.call([argv[1], "tilt", line, root + "-tilt", root + "-bump-blur" + ext])
  if(not os.path.exists(root + ".mp4")):
    os.system("ffmpeg -r 25 -loop 1 -i " + root + "-tilt-%d.ppm" + " -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -c:v libx264 -preset veryfast -pix_fmt yuv420p -t 00:00:15 " + root + ".mp4 && rm -f " + root + "-tilt-*.ppm")

