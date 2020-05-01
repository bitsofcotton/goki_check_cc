import sys
import numpy as np

block = 44100

def light(d, r):
  dr = []
  di = []
  for u in d:
    dr.append(u.real)
    di.append(u.imag)
  dr = np.fft.fft(dr)
  di = np.fft.fft(di)
  for u in range(1, len(dr)):
    dr[u] /= pow(np.exp(np.pi * 1.j * u / len(dr)) - 1., r)
    di[u] /= pow(np.exp(np.pi * 1.j * u / len(di)) - 1., r)
  dr = np.fft.ifft(dr)
  di = np.fft.ifft(di)
  for u in range(0, len(d)):
    d[u] = dr[u].real + 1.j * di[u].real
  return d

def diff(d, r):
  for u in range(0, len(d)):
    d[u] *= pow(- 2.j * np.pi * u / len(d), r)
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
    d2 = []
    for s in range(0, len(u) / block + 1):
      d2.append(np.fft.fft(u[s * block:min((s + 1) * block, len(u))]))
    for s in range(0, len(d2[0])):
      d = []
      for v in range(0, len(d2)):
        if(s < len(d2[v])):
          d.append(d2[v][s])
      if(sys.argv[1] == "ld"):
        d = diff(light(d, 1), 1)
      elif(sys.argv[1] == "l"):
        d = light(d, 1)
      elif(sys.argv[1] == "d"):
        d = diff(d, 1)
      for v in range(0, len(d)):
        d2[v][s] = d[v]
    for s in range(0, len(u) / block + 1):
      d = np.fft.ifft(d2[s])
      for v in d:
        out[- 1].append(float(v.real))
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

