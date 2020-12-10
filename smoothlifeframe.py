# smoothlifecanvas.py
#
# An object that extends the tkinter Canvas widget.

import os
import sys
import tkinter
from tkinter import ttk
import numpy as np

DEFAULT_CANVAS_WIDTH = 400
DEFAULT_CANVAS_HEIGHT = 400
DEFAULT_BUTTON_FRAME_WIDTH = 100
DEFAULT_BUTTON_FRAME_HEIGHT = 400

class SmoothLifeFrame(ttk.Frame):

    def __init__(self, x0=0, xn=1, y0=0, yn=1, numx=None, numy=None, r0=1, r1=3, b=1, master=None, **kwargs):
        """
        SmoothLifeFrame.__init__()
        
        Constructor fot the SmoothLifeFrame class, which extends the ttk frame.
        
        Inputs:
        
            x0 - The lower left corner of the canvas in Cartesian coordinates.
            xn - The lower right corner of the canvas in Cartesian coordinates.
        
        Outputs:
        
            A new SmoothLifeFrame instance.
        """
        
        super().__init__(master=master, **kwargs)
        
        # Set the members of the class. Most of these settings are for the canvas.
        self.x0 = x0
        self.xn = xn
        self.y0 = y0
        self.yn = yn
        self.numx = (DEFAULT_CANVAS_WIDTH // 2 + 1) if numx is None else numx
        self.numy = (DEFAULT_CANVAS_HEIGHT // 2 + 1) if numy is None else numy
        self.hx = (self.xn - self.x0) / (self.numx - 1)
        self.hy = (self.yn - self.y0) / (self.numy - 1)
        
        # Set the inner and outer radii of the annulus.
        self.r0 = r0
        self.r1 = r1
        
        # Set the antialising width to be the smaller of half the inner radius
        # or 1/4 of the thickness of the annulus.
        self.b = 0
        
        # Set the birth and death intervals.
        self.bint = (0.0, 0.6)
        self.dint = (0.4, 1.0)
        
        # Create the canvas object.
        self._create_canvas()
        self._create_lines()
        
        # Create the button panel on the side.
        self._create_button_panel()
        
        # Pack this frame into the root frame.
        self.grid(row=0, column=0,
                  sticky=(tkinter.N, tkinter.S, tkinter.E, tkinter.W))

    #-------------------------------------------------------------------------
    # Private methods for creating components of the GUI and drawings within
    # them.
    #-------------------------------------------------------------------------

    def _create_canvas(self):
        self.canvas_frame = ttk.Frame(master=self)
        self.canvas = tkinter.Canvas(master=self.canvas_frame,
                                     width=DEFAULT_CANVAS_WIDTH,
                                     height=DEFAULT_CANVAS_HEIGHT,
                                     background='white')

        self.canvas.grid(row=0, column=0,
                         sticky=(tkinter.N, tkinter.S, tkinter.E, tkinter.W))
        self.canvas_frame.grid(row=0, column=0)

    def _create_button_panel(self):
        # Create the panel itself.
        self.button_frame = ttk.Frame(master=self,
                                      width=DEFAULT_BUTTON_FRAME_WIDTH,
                                      height=DEFAULT_BUTTON_FRAME_HEIGHT)

        # Create the buttons to go into the grame.
        self.start_stop_button = ttk.Button(master=self.button_frame,
                                            text='Start')
        self.step_button = ttk.Button(master=self.button_frame,
                                      text='Step')
        self.clear_button = ttk.Button(master=self.button_frame,
                                       text='Clear')
        self.randomize_button = ttk.Button(master=self.button_frame,
                                           text='Randomize')

        self.start_stop_button.grid(row=0, column=0)
        self.step_button.grid(row=1, column=0)
        self.clear_button.grid(row=2, column=0)
        self.randomize_button.grid(row=3, column=0)
        self.button_frame.grid(row=0, column=1, sticky=(tkinter.E,))

    def _create_lines(self):
        width = float(self.canvas.cget('width'))
        height = float(self.canvas.cget('width'))
        
        print('width, height = {0}, {1}'.format(width, height))
        
        # Create the vertical lines.
        py0 = 0
        pyn = height
        
        # The "slope" for the x transformation.
        m = width / (self.xn - self.x0)
        
        for i in range(1, self.numx - 1):
            px = m * i * self.hx
            
            self.canvas.create_line(px, py0, px, pyn, fill='#000000', width=1, tag='vline')

        # Create the vertical lines.
        px0 = 0
        pxn = width

        # The "slope" for the y transformation.
        m = -height / (self.yn - self.y0)
        
        for j in range(1, self.numy - 1):
            py = m * j * self.hy + height
            
            self.canvas.create_line(px0, py, pxn, py, fill='#000000', width=1, tag='hline')

    def _create_cells(self):
        pass

    #-------------------------------------------------------------------------
    # Class methods to add and remove items from the canvas.
    #-------------------------------------------------------------------------
    
    def clear_all(self):
        self.delete('vline')
        self.delete('hline')
        self.delete('cell')
        
    def clear_cells(self):
        self.canvas.delete('cell')
    
    #-------------------------------------------------------------------------
    # Class methods for event bindings.
    #-------------------------------------------------------------------------
    
    def on_resize(self, event):
        pass

    def on_left_click(self, event):
        pass

    def on_center_click(self, event):
        pass

    def on_right_click(self, event):
        pass

    def on_key_press(self, event):
        pass

    def on_button_press(self, event):
        pass

    #-------------------------------------------------------------------------
    # Class methods for updating the grid at each step.
    #-------------------------------------------------------------------------

    def step(self):
        pass

    def update(self):
        pass

if __name__ == '__main__':
    x0 = 0
    xn = 10
    y0 = 0
    yn = 10
    numx = 80
    numy = 80
    width = 800
    height = 800
    r0 = 1
    r1 = 2

    master = tkinter.Tk()
    cnvs = SmoothLifeFrame(master=master, x0=x0, xn=xn, y0=y0, yn=yn,
                           numx=numx, numy=numy, r0=r0, r1=r1)

    tkinter.mainloop()