import numpy as np

ACASE48 = 'KIAM'
SPIDER = 'SPIDER'
IDUM = 3
MESHR = 128
MESHZ = 257

folder = 'Files'
filename = 'eqdsk_equilx'
with open("%s/%s" % (folder, filename), 'w') as file:
  file.write("%8s%8s%36d%4d%4d" % (ACASE48, SPIDER, IDUM, MESHR, MESHZ))
