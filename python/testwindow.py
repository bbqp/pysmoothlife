import moderngl_window as mglw
import numpy as np


class Test(mglw.WindowConfig):
    gl_version = (4, 6)
    window_size = (800, 800)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Do initialization here
        #self.prog = self.ctx.program(...)
        #self.vao = self.ctx.vertex_array(...)
        self.texture = self.ctx.texture(self.wnd.size, 4)

        self._init_vertex_data()

    def on_render(self, time: float, frametime: float):
        # This method is called every frame
        #self.vao.render()
        pass

    def _init_vertex_data(self):
        N = 21
        x = np.linspace(0, Test.window_size[0], N)
        y = np.linspace(0, Test.window_size[1], N)
        c = np.random.uniform(size=N)

        self.vertex_data = np.dstack([x, y, c])

if __name__ == '__main__':
    # Blocking call entering rendering/event loop
    Test.run()
