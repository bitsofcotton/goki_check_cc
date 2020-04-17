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

data0 = []
for line in sys.stdin:
  data0.append(float(line.split(",")[0]))

out = []
MM  = 0.
for s in range(0, len(data0) / block + 1):
  d = np.fft.fft(data0[s * block:min((s + 1) * block, len(data0))])
  if(sys.argv[1] == "ld"):
    for u in range(0, int(sys.argv[2])):
      d = diff(light(d, 1), 1)
  elif(sys.argv[1] == "l"):
    d = light(d, int(sys.argv[2]))
  elif(sys.argv[1] == "d"):
    d = diff(d, int(sys.argv[2]))
  d = np.fft.ifft(d)
  for u in d:
    out.append(u.real)
    MM = max(MM, abs(u.real))
for s in out:
  print s / MM, "\r"

