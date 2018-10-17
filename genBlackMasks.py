import numpy as np
from skimage import io
import os
import sys

folder = sys.argv[1]

print folder
root = os.path.split(folder)[0]
print root

for f in os.listdir(folder):
    print f
    if '.jpg' in f:
        img = io.imread(os.path.join(folder, f))
        zeros = np.zeros(img.shape)
        maskname = os.path.join(root, "Masks", f.split('.')[0] + '_mask.jpg')
        io.imsave(maskname, zeros)


