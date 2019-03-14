import cv2 as cv
import numpy as np
import tkinter as tk
from PIL import ImageTk, Image
import sys

imgCV = cv.imread(sys.argv[1])

print(imgCV.shape)

root = tk.Tk()
geometry = "%dx%d+0+0"%(imgCV.shape[0], imgCV.shape[1])
root.geometry()

def leftclick(event):
    global click_nr
    global skip
    if skip:
        skip = False
    else:
        sys.stdout.write('\n')
        click_nr += 1

    sys.stdout.write("{},{},{}".format(click_nr, event.x, event.y))
    sys.stdout.flush()
    

def rightclick(event):
    global skip
    skip = True
    sys.stdout.write('\r---                   \r')
    sys.stdout.flush()

# import image
click_nr = 0
skip = False

img = ImageTk.PhotoImage(Image.open(sys.argv[1]))
panel = tk.Label(root, image = img)
panel.bind("<Button-1>", leftclick)
panel.bind("<Button-3>", rightclick)
#panel.pack(side = "bottom", fill = "both", expand = "no")
panel.pack(fill = "both", expand = 1)

root.mainloop()
