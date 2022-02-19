"""Animation tests."""

from manim import *
from manim_physics import *
from utils import *
import numpy as np
import math

class LogoIntro(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi=55 * DEGREES,
            theta=0
        )

        def magnetic_field_func(p, i, location=ORIGIN):
            p = p - location
            r = np.linalg.norm(p)
            if False:
                return p
            else:
                return np.array([i * (1/r**2) * p[1], -i * (1/r**2) * p[0], 0])

        fil = Line(10 * IN, 10 * OUT, stroke_width=6)
        stream_lines = StreamLines(
            func=lambda p: magnetic_field_func(p, 4),
            virtual_time=10,
            # color=BLUE,
            y_range=[
                math.floor(-config["frame_width"] / 2),
                math.ceil(config["frame_width"] / 2)
            ],
        )
        self.add(stream_lines, fil)
        stream_lines.start_animation(
            flow_speed=1.5,
            time_width=0.5
        )
        self.wait(2)

        t = ValueTracker()

        def theta_func(t):
            return t * (450 * DEGREES) / 3

        def zoom_func(t):
            return 1 - (t * 0.9 / 3)
            
        # self.camera.theta_tracker.add_updater(
        #     lambda tr: tr.set_value(theta_func(t.get_value()))
        # )
        self.camera.zoom_tracker.add_updater(
            lambda tr: tr.set_value(zoom_func(t.get_value()))
        )
        self.add(self.camera.zoom_tracker)

        self.play(
            LaggedStart(
            t.animate.set_value(3),
            self.camera._frame_center.animate.move_to(15 * DOWN)
            ),
            run_time=3
        )
        self.wait()

class Test(Scene):
    def construct(self):
        # self.camera.background_color = "#FFFFFF"
        im = [[255 * math.sqrt(2) / np.linalg.norm([x, y])
               for x in np.linspace(-100, 100, 400)] for y in np.linspace(-50, 50, 200)]

        temp_image = ImageMobject(np.uint8(im))
        temp_image.height = 7
        self.add(temp_image)


        

