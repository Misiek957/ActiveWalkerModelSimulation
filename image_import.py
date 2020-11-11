#%%

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from PIL import Image
import cv2
# %matplotlib inline

im = Image.open("cc_fin.bmp")
im2 = mpimg.imread('cc_fin.bmp')
im3 = cv2.imread('cc_fin.bmp')
print(im3.shape)
#%%
b = im3.copy()
# set green and red channels to 0
b[:, :, 1] = 0
b[:, :, 2] = 0
# plt.imshow(b)
cv2.imshow('B-RGB', b)
#%%
g = im3.copy()
# set blue and red channels to 0
g[:, :, 0] = 0
g[:, :, 2] = 0
# plt.imshow(g)
cv2.imshow('G-RGB', g)
#%%
r = im3.copy()
# set blue and green channels to 0
r[:, :, 0] = 0
r[:, :, 1] = 0
# plt.imshow(r)
cv2.imshow('R-RGB', r)