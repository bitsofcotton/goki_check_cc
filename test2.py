import sys
import numpy as np

block = 44100

def light(d, r):
  dr = []
  di = []
  for uu in d:
    dr.append(uu.real)
    di.append(uu.imag)
  dr = np.fft.fft(dr)
  di = np.fft.fft(di)
  for uu in range(1, len(dr)):
    dr[uu] /= pow(np.exp(np.pi * 1.j * uu / len(dr)) - 1., r)
    di[uu] /= pow(np.exp(np.pi * 1.j * uu / len(di)) - 1., r)
  dr = np.fft.ifft(dr)
  di = np.fft.ifft(di)
  for uu in range(0, len(d)):
    d[uu] = dr[uu].real + 1.j * di[uu].real
  return d

def diff(d, r):
  for uu in range(0, len(d)):
    d[uu] *= pow(- 2.j * np.pi * uu / len(d), r)
  return d

data = []
d0   = []
for line in sys.stdin:
  if(line[0] == ';'):
    print line,
    continue
  l = [s for s in line.split(" ") if s != '' and s != '\r\n']
  if(len(data) < len(l)):
    for s in range(0, len(l) - len(data)):
      data.append([])
  d0.append(l[0])
  for s in range(1, len(l)):
    data[s - 1].append(float(l[s]))

for v in range(0, int(sys.argv[2])):
  out = []
  MM  = 0.
  for u in data:
    if(len(u) == 0):
      continue
    out.append([])
    for s in range(0, len(u) / block + 1):
      d = np.fft.fft(u[s * block:min((s + 1) * block, len(u))])
      if(sys.argv[1] == "ld"):
        d = diff(light(d, 1), 1)
      elif(sys.argv[1] == "l"):
        d = light(d, 1)
      elif(sys.argv[1] == "L"):
        d = light(u[s * block:min((s + 1) * block, len(u))], 1)
      elif(sys.argv[1] == "d"):
        d = diff(d, 1)
      if(sys.argv[1] != "L"):
        d = np.fft.ifft(d)
      for v in d:
        out[- 1].append(v.real)
    MM = max(MM, sorted(out[- 1])[len(out[- 1]) * 7 / 8] * 2.)
  data = out
  for u in data:
    for v in u:
      v = max(- 1, min(1, v / MM))
for u in range(0, len(data[0])):
  print d0[u], " ",
  for s in range(0, len(data)):
    if(u < len(data[s])):
      print max(- 1, min(1, data[s][u] / MM)), " ",
    else:
      print "0 ",
  print

