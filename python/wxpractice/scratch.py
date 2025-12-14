class SmoothLifeOpenGLCanvas(wx.lib.agw.speedmeter.BufferedWindow):
    def __init__(self, parent, size, style=wx.NO_FULL_REPAINT_ON_RESIZE, bufferedstyle=wx.lib.agw.speedmeter.SM_BUFFERED_DC):
        wx.lib.agw.speedmeter.BufferedWindow.__init__(self, parent, size=size, style=style, bufferedstyle=bufferedstyle)

        self.delta = 100
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.update, self.timer)
        
        self.Bind(wx.EVT_PAINT, self.OnPaint)

        # # Disable background erasing (flicker-licious)
        # def disable_event(*pargs,**kwargs):
            # pass # the sauce, please
        # self.Bind(wx.EVT_ERASE_BACKGROUND, disable_event) 
 
        self.width = size[0]
        self.height = size[1]
        
        self.nrows = self.height
        self.ncols = self.width
        
        self.dpx = self.width // self.nrows
        self.dpy = self.height // self.ncols
        
        self.state = np.random.rand(self.nrows, self.ncols)
        
        self.colormap = mpl.colormaps['viridis']
        
        self._create_color_list()
        self._create_rect_list()
        self._create_pen_list()
        self._create_brush_list()

    def start(self, e):
        self.timer.Start(milliseconds=self.delta)

    def step(self, e):
        self.timer.Start(milliseconds=self.delta, oneShot=wx.TIMER_ONE_SHOT)
        
    def stop(self, e):
        self.timer.Stop()
        
    def clear(self, e):
        self.state[:,:] = 0
        self._update_color_list()
        self._update_pen_list()
        self._update_brush_list()
        
        self.Refresh()

    def populate(self, e):
        self.state[:,:] = np.random.randint(0, 2, size=(self.nrows, self.ncols))
        self._update_color_list()
        self._update_pen_list()
        self._update_brush_list()
        
        self.Refresh()

    def update(self, e):
        self.state[:,:] = np.random.randint(0, 2, size=(self.nrows, self.ncols))
        self._update_color_list()
        self._update_pen_list()
        self._update_brush_list()
        
        self.UpdateDrawing()
        self.Refresh()

    def _create_color_list(self):
        self.color_list = []

        for j in range(self.ncols):
            for i in range(self.nrows):
                color = [int(c * 255) for c in self.colormap(self.state[i, j])]
                self.color_list.append(color)

    def _create_rect_list(self):
        self.rect_list = []
        
        for j in range(self.ncols):
            y0 = j * self.dpy
            for i in range(self.nrows):
                x0 = i * self.dpx
                self.rect_list.append((x0, y0, self.dpx, self.dpy))
        
    def _create_pen_list(self):
        self.pen_list = []
        
        for k in range(len(self.color_list)):
            self.pen_list.append(wx.ThePenList.FindOrCreatePen(wx.Colour(*self.color_list[k])))
        
    def _create_brush_list(self):
        self.brush_list = []
        
        for k in range(len(self.color_list)):
            self.brush_list.append(wx.TheBrushList.FindOrCreateBrush(wx.Colour(*self.color_list[k])))
            
    def _update_color_list(self):
        for j in range(self.ncols):
            for i in range(self.nrows):
                color = [int(c * 255) for c in self.colormap(self.state[i, j])]
                self.color_list[j * self.nrows + i] = color
        
    def _update_pen_list(self):
        for k in range(len(self.color_list)):
            self.pen_list[k] = wx.ThePenList.FindOrCreatePen(wx.Colour(*self.color_list[k]))

    def _update_brush_list(self):
        for k in range(len(self.color_list)):
            self.brush_list[k] = wx.TheBrushList.FindOrCreateBrush(wx.Colour(*self.color_list[k]))

    def OnSize(self, e):
        pass

    def OnPaint(self, e):
        self.draw_grid()
        
    def draw_grid(self):
        dc = wx.BufferedPaintDC(self)  # This is WAY faster than paintDC.
        dc.Clear()
        
        self.Draw(dc)
        
    def Draw(self, dc):
        dc.DrawRectangleList(self.rect_list, self.pen_list, self.brush_list)