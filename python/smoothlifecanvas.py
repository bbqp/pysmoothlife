"""
smoothlifecanvas.py

An extension of the tkinter canvas widget, but for SmoothLife
"""

import tkinter as tk


class SmoothLifeCanvas(tk.Canvas):
	def __init__(self, x0, xn, y0, yn, numx, numy, r0, rn, aa, b0, bn, d0, dn, master=None, cnf={}, **kw):
		super().__init__(master, cnf, **kw)
        
        self.linetag = 'line'
        self.celltag = 'cell'
        
        self.x0 = x0
        self.xn = xn
        self.numx = numx
        self.hx = (self.xn - self.x0) / (self.numx - 1)
        
        self.y0 = y0
        self.yn = yn
        self.numy = numy
        self.hy = (self.yn - self.y0) / (self.numy - 1)
        
        # Set the radii for the annulus.
        self.r0 = r0
        self.rn = rn
        self.aa = aa
        
        # Set the birth and death intervals.
        self.b0 = b0
        self.bn = bn
        self.d0 = d0
        self.dn = dn
        
        # Add the lines and pack the canvas.
        self.add_lines()
        self.pack()
        
        # Create arrays to keep track of the current and next iterates.
        self.state = np.zeros(shape=(self.numx - 1, self.numy - 1), dtype=np.float32)

    def init_canvas(self, seedval=0):
        self.add_lines()
    
        np.random.seed(seedval)
        self.state = np.random.rand()

	def add_lines(self):
        # Create the x-lines.
        py0 = 0
        pyn = self.height
        for i in range(self.numx):
            px = ((i * self.hx) / (self.xn - self.x0)) * self.width
            self.create_line(px, py0, px, pyn, tags=self.linetag)

        # Create the x-lines.
        px0 = 0
        pxn = self.width
        for j in range(self.numy):
            py = ((j * self.hy) / (self.yn - self.y0)) * self.height
            self.create_line(px0, py, pxn, py, tags=self.linetag)
    
    def clear(self):
        self.delete(self.celltag)

    def clearall(self):
        self.delete('all')

    def update(self):
        self.clear()

if __name__ == '__main__':
	# init tk
	root = tkinter.Tk()

	# create canvas
	myCanvas = tkinter.Canvas(root, bg="white", height=300, width=300)

	# draw arcs
	coord = 10, 10, 300, 300
	arc = myCanvas.create_arc(coord, start=0, extent=150, fill="red")
	arv2 = myCanvas.create_arc(coord, start=150, extent=215, fill="green")

	# add to window and show
	myCanvas.pack()
	root.mainloop()