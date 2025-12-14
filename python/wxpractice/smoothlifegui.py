import wx
import numpy as np
import matplotlib as mpl
import smoothlifestate as sls


class SmoothLifeCanvas(wx.Panel):
    def __init__(self, parent, size, name='SmoothLife Canvas', style=wx.NO_FULL_REPAINT_ON_RESIZE):
        wx.Panel.__init__(self, parent, size=size, name=name, style=style)

        self.delta = 100
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.update, self.timer)
        
        self.Bind(wx.EVT_PAINT, self.OnPaint)

        # Disable background erasing (flicker-licious)
        def disable_event(*pargs,**kwargs):
            pass # the sauce, please
        self.Bind(wx.EVT_ERASE_BACKGROUND, disable_event) 
 
        self.width = size[0]
        self.height = size[1]
        
        self.nrows = self.height // 10
        self.ncols = self.width // 10
        
        self.dpx = self.width // self.ncols
        self.dpy = self.height // self.nrows
        
        self.state = sls.SmoothlifeState(numx=self.ncols+1, numy=self.nrows+1)
        self.state.set_random()
        #self.state.set_value(1)
        self.state_data = self.state.get_state_data_copy()
        print(self.state_data)
        
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
        self.state.set_value(0)
        self.state_data = self.state.get_state_data_copy()
        print(self.state_data)
        self._update_color_list()
        self._update_pen_list()
        self._update_brush_list()
        
        self.Refresh()

    def populate(self, e):
        self.state.set_random()
        self.state_data = self.state.get_state_data_copy()
        print(self.state_data)
        self._update_color_list()
        self._update_pen_list()
        self._update_brush_list()
        
        self.Refresh()

    def update(self, e):
        self.state.update()
        self.state_data = self.state.get_state_data_copy()
        print(self.state_data)
        
        self._update_color_list()
        self._update_pen_list()
        self._update_brush_list()
        
        self.Refresh()

    def _create_color_list(self):
        self.color_list = []

        for j in range(self.ncols):
            for i in range(self.nrows):
                color = [int(c * 255) for c in self.colormap(self.state_data[i, j])]
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
                color = [int(c * 255) for c in self.colormap(self.state_data[i, j])]
                self.color_list[j * self.nrows + i] = color
        
    def _update_pen_list(self):
        for k in range(len(self.color_list)):
            self.pen_list[k] = wx.ThePenList.FindOrCreatePen(wx.Colour(*self.color_list[k]))

    def _update_brush_list(self):
        for k in range(len(self.color_list)):
            self.brush_list[k] = wx.TheBrushList.FindOrCreateBrush(wx.Colour(*self.color_list[k]))


    def OnPaint(self, e):
        self.draw_grid()
        
    def draw_grid(self):
        # dc = wx.PaintDC(self)
        dc = wx.BufferedPaintDC(self)  # This is WAY faster than paintDC.
        dc.Clear()
        
        dc.DrawRectangleList(self.rect_list, self.pen_list, self.brush_list)


class SmoothLifeFrame(wx.Frame):
    def __init__(self, parent=None, title='SmoothLife'):
        wx.Frame.__init__(self, parent, title=title)

        # A pointer to the game state.
        #self.state = GameOfLifeState(nrows, ncols)

        # Create a sizer for laying out the panels in the frame with a 10-pixel gap.
        self.grid_rows = 2
        self.grid_cols = 1
        self.grid_gap_px = 10
        #self.panel_sizer = wx.GridSizer(self.grid_rows, self.grid_cols, self.grid_gap_px)
        self.panel_sizer = wx.BoxSizer(orient=wx.VERTICAL)

        self._setup_canvas_panel()
        self._setup_button_panel()

        self.panel_sizer.Add(self.canvas)
        self.panel_sizer.Add(self.button_panel)
        self.SetSizer(self.panel_sizer)
        self.Fit()

    def _setup_canvas_panel(self):
        # Canvas for drawing.
        self.canvas = SmoothLifeCanvas(self, size=wx.Size(800, 800))

    def _setup_button_panel(self):
        # Create a sizer for the botton layout in the bottom panel.
        self.button_panel = wx.Panel(self, name="Button Panel")

        self.populate_button = wx.Button(self.button_panel, id=wx.ID_ANY, label="Populate", name="Populate Button")
        self.clear_button = wx.Button(self.button_panel, id=wx.ID_ANY, label="Clear", name="Clear Button")
        self.start_button = wx.Button(self.button_panel, id=wx.ID_ANY, label="Start", name="Start Button")
        self.step_button = wx.Button(self.button_panel, id=wx.ID_ANY, label="Step", name="Step Button")
        self.stop_button = wx.Button(self.button_panel, id=wx.ID_ANY, label="Stop", name="Stop Button")
        self.exit_button = wx.Button(self.button_panel, id=wx.ID_ANY, label="Exit", name="Exit Button")

        self.populate_button.Bind(wx.EVT_BUTTON, self.canvas.populate)
        self.clear_button.Bind(wx.EVT_BUTTON, self.canvas.clear)
        self.start_button.Bind(wx.EVT_BUTTON, self.canvas.start)
        self.step_button.Bind(wx.EVT_BUTTON, self.canvas.step)
        self.stop_button.Bind(wx.EVT_BUTTON, self.canvas.stop)
        self.exit_button.Bind(wx.EVT_BUTTON, self.on_exit)

        self.button_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.button_sizer.Add(self.populate_button, 1, wx.ALL)
        self.button_sizer.Add(self.clear_button, 1, wx.ALL)
        self.button_sizer.Add(self.start_button, 1, wx.ALL)
        self.button_sizer.Add(self.step_button, 1, wx.ALL)
        self.button_sizer.Add(self.stop_button, 1, wx.ALL)
        self.button_sizer.Add(self.exit_button, 1, wx.ALL)
        
        self.button_panel.SetSizer(self.button_sizer)
        
    def on_exit(self, e):
        self.Close(True)

if __name__ == '__main__':
    app = wx.App(False)
    frame = SmoothLifeFrame()
    frame.Show(True)
    app.MainLoop()
    # app = wx.App(False)
    # frame = SmoothLifeFrame()
    # frame.Show(True)
    # app.MainLoop()














