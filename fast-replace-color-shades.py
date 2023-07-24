import numpy as np
import argparse
from tqdm import tqdm
from PIL import Image
from colorama import init
init()
from colorama import Style, Fore, Back
from traceback_with_variables import activate_by_import
import os
import pyprodigy.main as pdg
#from numba import njit

pdg.check("header", "Image Color Swap")

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input file")
parser.add_argument("output", help="Location to save the output file")
parser.add_argument("colormap", help="Location of colormap : a csv file with 6 columns corresponding to input RGB and output RGB values")
parser.add_argument("delta", help="The color delta")

args = parser.parse_args()

input_ = args.input
file_name = input_.split("/")[-1]
output = args.output
colormap_loc = args.colormap
colour_delta = float(args.delta)

# Load Image
pdg.check("pass", "Image Loaded")
image = np.asarray(Image.open(input_))

# Load colormap : format = dictionary with input as keys and output as values
pdg.check("pass", "Colour map loaded")
colormap_raw = np.genfromtxt(colormap_loc, delimiter=",", dtype=int)
colormap = {}

for colors in colormap_raw:
    in_col = tuple(colors[:3])
    out_col = colors[3:]

    colormap[in_col] = [out_col]

print("COLORMAP = ", str(colormap))


new_image = np.zeros(shape=image.shape)
rows = image.shape[0]
cols = image.shape[1]

print("Shape of the image is : ", rows, " x ", cols)

#@njit
def in_delta(in_pixel, out_pixel, cutoff):
    in_pixel_ = np.array(in_pixel)
    out_pixel_ = np.array(out_pixel)

    delta = np.average(np.absolute(out_pixel_-in_pixel_))
    if delta < cutoff:
        return True, out_pixel
    else:
        return False, None

keys = colormap.keys()

# Transform Image

for key in tqdm(keys, colour="blue", leave=False):
    transfer_image = np.zeros(shape=image.shape) + np.array([[list(key)+[255]]])
    delta_array = image - transfer_image
    delta_array = np.average(delta_array[:,:,:3].reshape(-1,3), axis=1).reshape(rows, cols)
    change_indices = np.argwhere(delta_array > colour_delta)
    
    colour = np.array([[list(key) + [255]]])

    for index in tqdm(change_indices, leave=False):
        #print(index)
        image[tuple(index)] = colour
    
pdg.check("pass", "Image transformation complete")
#new_image[row_idx][col_idx] = new_pixel

# Save Image
out_image = Image.fromarray(image.astype("uint8"), "RGBA")
out_image.save(os.path.join(output, file_name), "PNG")

pdg.check("pass", "Image saved.")
pdg.check("finish")
