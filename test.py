#! /usr/bin/env python

import os
import sys
import subprocess

argv = sys.argv

if(len(argv) < 3):
  print "not much arguments."
  exit(- 1)

for idx in range(2, len(argv)):
  line = argv[idx]
  print line
  root, ext = os.path.splitext(line)
  if(ext != ".ppm"):
    subprocess.call(["convert", line, "-compress", "none", line + ".ppm"])
    line = line + ".ppm"
    root, ext = os.path.splitext(line)
  if(not os.path.exists(root + "-enl-last" + ext)):
    subprocess.call([argv[1], "enlarge", line, root + "-enl" + ext])
    subprocess.call(["convert", line, "-matte", "-virtual-pixel", "Mirror", "-filter", "point", "+distort", "SRT", "45", "-compress", "none", root + "-d" + ext])
    subprocess.call([argv[1], "enlarge", root + "-d" + ext, root + "-enl-d" + ext])
    subprocess.call(["convert", root + "-enl-d" + ext, "-rotate", "-45", root + "-enl-d-margin" + ext])
    p = subprocess.Popen(["identify", "-format", "%wx%h", root + "-enl" + ext], stdout=subprocess.PIPE)
    szorig = p.communicate("\n")[0] + "+0+0"
    subprocess.call(["convert", root + "-enl-d-margin" + ext, "-gravity", "center", "-crop", szorig, root + "-enl-d-crop" + ext])
    subprocess.call(["convert", root + "-enl" + ext, root + "-enl-d-crop" + ext, "-evaluate-sequence", "mean", "-resize", "75%", root + "-enl-last" + ext])
  if(not os.path.exists(root + "-collect" + ext)):
    subprocess.call([argv[1], "collect", line, root + "-collect" + ext])
  if(not os.path.exists(root + "-bump-blur" + ext)):
    subprocess.call([argv[1], "bump", line, root + "-bump" + ext, root + "-bump-blur" + ext, root + ".obj"])
  if(not os.path.exists(root + "-emph" + ext)):
    subprocess.call(["convert", root + "-bump" + ext, "-blur", "2x2", "-alpha", "on", "-channel", "a", "-evaluate", "set", "30%", root + "-bump-test.png"])
    subprocess.call(["convert", line, root + "-bump-test.png", "-compose", "Multiply", "-composite", root + "-emph" + ext])
  if(not os.path.exists(root + "-tilt-1" + ext) and not os.path.exists(root + ".mp4")):
    subprocess.call([argv[1], "tilt", line, root + "-tilt", root + "-bump" + ext])
  if(not os.path.exists(root + ".mp4")):
    os.system("ffmpeg -r 16 -loop 1 -i " + root + "-tilt-%d.ppm" + " -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -c:v libx264 -preset veryfast -pix_fmt yuv420p -t 00:00:8 " + root + ".mp4 && rm -f " + root + "-tilt-*.ppm")

