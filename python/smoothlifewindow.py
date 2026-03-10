import moderngl_window as mglw

class Test(mglw.WindowConfig):
    gl_version = (4, 6)
    window_size = (800, 600)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Do initialization here
        #self.prog = self.ctx.program(...)
        #self.vao = self.ctx.vertex_array(...)
        self.texture = self.ctx.texture(self.wnd.size, 4)

    def on_render(self, time: float, frametime: float):
        # This method is called every frame
        #self.vao.render()
        pass

if __name__ == '__main__':
    # Blocking call entering rendering/event loop
    Test.run()
