"""
OpenING manim source code, used to produce physics related videos.
See: https://www.polytech-reseau.org/opening/
"""

from manim import *
from manim_physics import *
from utils import *
import numpy as np

YELLOW_G = "#FCBA03"

class Anim14_1(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=45 * DEGREES
        )
        axes1 = ThreeDAxes(
            x_length=config.frame_height + 1.5
        ).set_stroke(width=0.7).scale(0.45).shift(3 * DR)
        axes2 = ThreeDAxes().set_stroke(width=0.7).scale(0.45)
        axes3 = ThreeDAxes(
            y_length=config.frame_height + 1.5
        ).set_stroke(width=0.7).scale(0.45).shift(3 * UL)

        x_label1 = axes1.get_x_axis_label('x').scale(0.7)
        y_label1 = axes1.get_y_axis_label('y').scale(0.7)
        z_label1 = axes1.get_z_axis_label('z', rotation=0).scale(0.7)

        x_label2 = axes2.get_x_axis_label('x').scale(0.7)
        y_label2 = axes2.get_y_axis_label('y').scale(0.7)
        z_label2 = axes2.get_z_axis_label('z', rotation=0).scale(0.7)

        x_label3 = axes3.get_x_axis_label('x').scale(0.7)
        y_label3 = axes3.get_y_axis_label('y').scale(0.7).shift(1.7 * UP)
        z_label3 = axes3.get_z_axis_label('z', rotation=0).scale(0.7)

        labels1 = VGroup(x_label1, y_label1, z_label1)
        labels2 = VGroup(x_label2, y_label2, z_label2)
        labels3 = VGroup(x_label3, y_label3, z_label3)

        line11 = DashedLine(
            axes1.c2p(0, 2.5, 0),
            axes1.c2p(2.5, 2.5, 0),
            stroke_width=3
        )
        line12 = DashedLine(
            axes1.c2p(2.5, 0, 0),
            axes1.c2p(2.5, 2.5, 0),
            stroke_width=3
        )
        line13 = DashedLine(
            axes1.c2p(2.5, 2.5, 0),
            axes1.c2p(2.5, 2.5, 2),
            stroke_width=3
        )
        line14 = DashedLine(
            axes1.c2p(2.5, 2.5, 2),
            axes1.c2p(0, 0, 2),
            stroke_width=3
        )
        m1 = Tex('M').scale(0.7).next_to(axes1.c2p(2.5, 2.5, 2), LEFT, buff=SMALL_BUFF)
        xm1 = MathTex("x_M").scale(0.5).next_to(line11, UP)
        ym1 = MathTex("y_M").scale(0.5).next_to(line12)
        zm1 = MathTex("z_M").scale(0.5).next_to(line13, LEFT, buff=SMALL_BUFF)

        line21 = DashedLine(
            axes2.c2p(0, 0, 0),
            axes2.c2p(2.5, 3, 0),
            stroke_width=3
        )
        line22 = DashedLine(
            axes2.c2p(2.5, 3, 0),
            axes2.c2p(2.5, 3, 3),
            stroke_width=3
        )
        line23 = Line(
            axes2.c2p(0, 0, 0),
            axes2.c2p(2.5, 3, 3),
            stroke_width=3
        )
        angle2 = Angle(
            Line(axes2.c2p(0, 0, 0), axes2.c2p(1, 0, 0)),
            line21
        )

        m2 = Tex('M').scale(0.7).next_to(axes2.c2p(2.5, 3, 3), LEFT, buff=SMALL_BUFF)
        thetam2 = MathTex(r"\theta_M").scale(0.5).move_to(0.4 * UR + 0.1 * RIGHT)
        rm2 = MathTex("r_M").scale(0.5).move_to(0.9 * UR)
        zm2 = MathTex("z_M").scale(0.5).next_to(line22, LEFT, buff=SMALL_BUFF)

        line31 = DashedLine(
            axes3.c2p(0, 0, 0),
            axes3.c2p(2.5, 3.5, 0),
            stroke_width=3
        )
        line32 = DashedLine(
            axes3.c2p(2.5, 3.5, 0),
            axes3.c2p(2.5, 3.5, 3),
            stroke_width=3
        )
        line33 = Line(
            axes3.c2p(0, 0, 0),
            axes3.c2p(2.5, 3.5, 3),
            stroke_width=3
        )

        angle31 = Angle(
            Line(axes3.c2p(0, 0, 0), axes3.c2p(1, 0, 0)),
            line31
        )
        angle32 = Arc(
            radius=0.8,
            angle=cartesian_to_spherical([2.5, 3.5, 3])[1]
        ).rotate(
            PI / 2,
            axis=DOWN,
            about_point=UP
        ).move_to(axes3.c2p(0, 0.72, 1.8))

        m3 = Tex('M').scale(0.7).move_to(axes3.c2p(2.7, 3.5, 3.5))
        phim3 = MathTex(r"\varphi_M").scale(0.5).move_to(axes3.c2p(1.8, 1, 0))
        thetam3 = MathTex(r"\theta_M").scale(0.5).move_to(axes3.c2p(0, 1, 3))
        rm3 = MathTex("r_M").scale(0.5).move_to(axes3.c2p(2.5, 2.5, 1.4))

        self.add_fixed_orientation_mobjects(*labels1, xm1, ym1, zm1, m1)
        self.play(
            LaggedStart(
                AnimationGroup(
                    FadeIn(axes1),
                    FadeIn(labels1)
                ),
                AnimationGroup(
                    *map(Create, [line11, line12, line13, line14]),
                    *map(FadeIn, [xm1, ym1, zm1, m1])
                ),
                lag_ratio=0.6
            )
        )
        self.add_fixed_orientation_mobjects(*labels2, thetam2, rm2, zm2, m2)
        self.play(
            LaggedStart(
                AnimationGroup(
                    FadeIn(axes2),
                    FadeIn(labels2)
                ),
                AnimationGroup(
                    *map(Create, [line21, line22, line23, angle2]),
                    *map(FadeIn, [thetam2, rm2, zm2, m2])
                ),
                lag_ratio=0.6
            )
        )
        self.add_fixed_orientation_mobjects(*labels3, phim3, thetam3, rm3, m3)
        self.play(
            LaggedStart(
                AnimationGroup(
                    FadeIn(axes3),
                    FadeIn(labels3)
                ),
                AnimationGroup(
                    *map(Create, [line31, line32, line33, angle31, angle32]),
                    *map(FadeIn, [phim3, thetam3, rm3, m3])
                ),
                lag_ratio=0.6
            )
        )
        self.wait()

        cartesian_t = Tex("Cartésien").move_to(4.3 * LEFT + 2 * DOWN)
        cylindric_t = Tex("Cylindrique").move_to(2 * DOWN)
        spheric_t = Tex("Sphérique").move_to(4.5 * RIGHT + 2 * DOWN)

        self.add_fixed_in_frame_mobjects(cartesian_t)
        self.play(
            FadeIn(cartesian_t),
            run_time=0.8
        )
        self.add_fixed_in_frame_mobjects(cylindric_t)
        self.play(
            FadeIn(cylindric_t),
            run_time=0.8
        )
        self.add_fixed_in_frame_mobjects(spheric_t)
        self.play(
            FadeIn(spheric_t),
            run_time=0.8
        )

        self.wait()

class Anim14_2(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=35 * DEGREES
        )
        axes = ThreeDAxes().set_stroke(width=0.6)
        x_label = axes.get_x_axis_label('x')
        y_label = axes.get_y_axis_label('y').shift(2 * UP)
        z_label = axes.get_z_axis_label('z',buff=LARGE_BUFF, rotation=0)

        e_x = Arrow(ORIGIN, 1.2 * RIGHT, buff=0)
        e_y = Arrow(ORIGIN, 1.2 * UP, buff=0)
        e_z = Arrow(ORIGIN, 1.2 * OUT, buff=0)
        e_x_t = MathTex(r"\overrightarrow{e_x}").scale(0.8).next_to(e_x.get_tip(), DOWN + 0.2 * LEFT)
        e_y_t = MathTex(r"\overrightarrow{e_y}").scale(0.8).next_to(e_y.get_tip(), 1.5 * OUT)
        e_z_t = MathTex(r"\overrightarrow{e_z}").scale(0.8).next_to(e_z.get_tip(), OUT)


        labels = VGroup(x_label, y_label, z_label)
        cartesian_base = VGroup(e_x, e_y, e_z)
        cartesian_base_t = VGroup(e_x_t, e_y_t, e_z_t)
        self.add_fixed_orientation_mobjects(*labels, *cartesian_base_t)
        self.play(
            Create(axes),
            Write(labels),
            Create(cartesian_base),
            Write(cartesian_base_t)
        )
        self.wait()

        m_dot = Dot(shade_in_3d=True)
        line_x = DashedLine(ORIGIN, m_dot)
        line_x_t = MathTex("x_M").scale(0.8).next_to(line_x, 1.3 * UP)

        def line_x_updater(line):
            line.become(
                DashedLine(
                    [0, m_dot.get_center()[1], 0],
                    [m_dot.get_center()[0], m_dot.get_center()[1], 0]
                )
            )
        
        def line_x_t_updater(t):
            t.next_to(line_x, 1.3 * UP)
        
        self.add_fixed_orientation_mobjects(line_x_t)
        self.play(
            Create(m_dot),
            Create(line_x),
            FadeIn(line_x_t)
        )
        line_x.add_updater(line_x_updater)
        line_x_t.add_updater(line_x_t_updater)
        self.play(
            m_dot.animate.move_to(3 * RIGHT)
        )
        self.wait(0.7)

        line_y = DashedLine([m_dot.get_center()[0], 0, 0], m_dot)
        line_y_t = MathTex("y_M").scale(0.8).next_to(line_y, 1.3 * RIGHT)

        def line_y_updater(line):
            line.become(
                DashedLine(
                    [m_dot.get_center()[0], 0, 0],
                    [m_dot.get_center()[0], m_dot.get_center()[1], 0]
                )
            )

        def line_y_t_updater(t):
            t.next_to(line_y, 1.3 * RIGHT)
        
        self.add(line_y)
        self.add_fixed_orientation_mobjects(line_y_t)
        self.play(
            FadeIn(line_y_t)
        )
        line_y.add_updater(line_y_updater)
        line_y_t.add_updater(line_y_t_updater)

        self.play(
            m_dot.animate.shift(2.5 * UP)
        )
        self.wait(0.7)

        line_z = DashedLine(m_dot, m_dot)
        line_z_t = MathTex("z_M").scale(0.8).next_to(line_z, 1.3 * UP)

        def line_z_updater(line):
            line.become(
                DashedLine([m_dot.get_center()[0], m_dot.get_center()[1], 0], m_dot.get_center())
            )

        def line_z_t_updater(t):
            t.next_to(line_z, 1.3 * UP)
        
        self.add(line_z)
        self.add_fixed_orientation_mobjects(line_z_t)
        self.play(
            FadeIn(line_z_t)
        )
        line_z.add_updater(line_z_updater)
        line_z_t.add_updater(line_z_t_updater)

        self.play(
            m_dot.animate.shift(1.1 * OUT)
        )
        m_t = Tex('M').scale(0.8).next_to(m_dot, UP)
        self.add_fixed_orientation_mobjects(m_t)
        self.play(
            FadeIn(m_t)
        )
        self.wait()

        m_def = MathTex("M(x_M, y_M, z_M)").scale(0.9).to_corner(UL)
        self.add_fixed_in_frame_mobjects(m_def)
        self.play(
            TransformFromCopy(Dot().move_to(0.42 * RIGHT + 0.34 * DOWN), m_def)
        )
        self.wait()
        
        self.play(
            Indicate(
                VGroup(*cartesian_base, *cartesian_base_t),
                scale_factor=1.3,
                color=YELLOW_G,
                run_time=0.85
            )
        )
        self.wait()

        om_line = DashedLine(ORIGIN, m_dot.get_center(), color=BLUE_D)
        om_t = MathTex(
            r"OM=\sqrt{{x_M}^2+{y_M}^2+{z_M}^2}",
            color=BLUE_D
        ).scale(0.8).next_to(m_def, DOWN, aligned_edge=LEFT)

        self.move_camera(
            phi=60 * DEGREES,
            theta=20 * DEGREES,
            added_anims=[Create(om_line)]
        )
        self.add_fixed_in_frame_mobjects(om_t)
        self.play(
            FadeIn(om_t)
        )
        self.wait(2)
        line_x.remove_updater(line_x_updater)
        line_y.remove_updater(line_y_updater)
        line_z.remove_updater(line_z_updater)
        self.move_camera(
            zoom=2,
            frame_center=2 * UP,
            added_anims=[
                *map(FadeOut, [m_def, om_line, om_t, line_x, line_x_t, line_y, line_y_t, line_z, line_z_t]),
                m_t.animate.scale(1.2).shift(0.2 * OUT + 0.7 * DOWN),
                m_dot.animate.scale(0.8)
            ]
        )
        self.wait()

        line_dx = Line(
            m_dot.get_center(),
            m_dot.get_center() + 0.7 * RIGHT
        )
        dx_t = MathTex(r"\mathrm{d}x").next_to(line_dx, DOWN, buff=SMALL_BUFF)

        surf_dx_dy = Rectangle(
            height=0,
            width=0.7,
            fill_opacity=0.5,
            fill_color=BLUE_D
        ).move_to(m_dot.get_center() + 0.35 * RIGHT)

        dy_t = MathTex(
            r"\mathrm{d}y"
        ).next_to(line_dx, DOWN, buff=SMALL_BUFF).shift(0.6 * RIGHT + 0.7 * UP)

        y_value = ValueTracker(0)

        def surf_updater(surf):
            surf.become(
                Rectangle(
                    height=y_value.get_value(),
                    width=0.7,
                    fill_opacity=surf.get_fill_opacity(),
                    fill_color=surf.get_fill_color()
                ).move_to(m_dot.get_center() + 0.35 * RIGHT + y_value.get_value() / 2 * UP)
            )

        self.add_fixed_orientation_mobjects(dx_t)

        self.play(
            LaggedStart(
                Create(line_dx),
                FadeIn(dx_t),
                lag_ratio=0.3
            )
        )
        self.wait(0.5)
        self.add(surf_dx_dy)
        surf_dx_dy.add_updater(surf_updater)
        self.add_fixed_orientation_mobjects(dy_t)
        self.play(
            LaggedStart(
                y_value.animate.set_value(0.7),
                FadeIn(dy_t),
                lag_ratio=0.3
            )
        )
        self.wait(0.5)

        vol_dxyz = Prism(
            dimensions=[0.7, 0.7, 0.7],
            fill_color=BLUE_D,
            stroke_width=DEFAULT_STROKE_WIDTH
        ).move_to(m_dot.get_center() + 0.35 * UR)

        dz_t = MathTex(
            r"\mathrm{d}z"
        ).next_to(surf_dx_dy, UP, buff=SMALL_BUFF).shift(0.5 * OUT + 0.7 * RIGHT + 0.1 * UP)

        z_value = ValueTracker()

        def vol_updater(vol):
            vol.become(
                Prism(
                    dimensions=[0.7, 0.7, z_value.get_value()],
                    fill_color=BLUE_D,
                    stroke_width=DEFAULT_STROKE_WIDTH
                ).move_to(m_dot.get_center() + 0.35 * UR + z_value.get_value() / 2 * OUT)
            )

        self.add(vol_dxyz)
        self.remove(surf_dx_dy)
        self.remove(line_dx)
        vol_dxyz.add_updater(vol_updater)
        self.add_fixed_orientation_mobjects(dz_t)

        self.play(
            LaggedStart(
                z_value.animate.set_value(0.35),
                FadeIn(dz_t),
                lag_ratio=0.3
            )
        )
        self.wait()

        vol_dxyz.remove_updater(vol_updater)
        self.play(
            LaggedStart(
                AnimationGroup(
                    FadeOut(m_dot),
                    FadeOut(m_t)
                ),
                VGroup(vol_dxyz, dx_t, dy_t, dz_t).animate.scale(0.1),
                lag_ratio=0.4
            )
        )
        self.wait()

        self.play(
            self.camera._frame_center.animate(run_time=2).move_to(
                vol_dxyz.get_center()
            ),
            self.camera.zoom_tracker.animate(run_time=3).set_value(15),
        )

        self.wait()

        self.play(
            LaggedStart(
                self.camera.zoom_tracker.animate(run_time=1.5).set_value(2),
                VGroup(vol_dxyz, dx_t, dy_t, dz_t).animate(run_time=1.5).scale(10),
                lag_ratio=0.3
            )
        )

        self.wait()

        dv_t = MathTex(r"\mathrm{d}^3V").move_to(3 * RIGHT + DOWN)
        dtau_t = MathTex(r"\mathrm{d}\tau").next_to(dv_t, DOWN, aligned_edge=RIGHT)

        self.add_fixed_in_frame_mobjects(dv_t)

        self.play(
            Write(dv_t)
        )
        self.wait()
        self.add_fixed_in_frame_mobjects(dtau_t)
        self.play(
            Write(dtau_t)
        )
        self.wait()
        self.play(
            FadeOut(dtau_t)
        )
        self.wait()
        self.play(
            Indicate(
                dv_t[0][1],
                color=YELLOW_G
            )
        )
        self.wait()
        
        d_product_t = MathTex(r"=\mathrm{d}x\mathrm{d}y\mathrm{d}z").next_to(dv_t, aligned_edge=DOWN)
        self.add_fixed_in_frame_mobjects(d_product_t)
        self.play(
            FadeIn(d_product_t, lag_ratio=0.1)
        )
        self.wait()

class Anim14_3(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=35 * DEGREES
        )
        axes = ThreeDAxes().set_stroke(width=0.6)
        x_label = axes.get_x_axis_label('x')
        y_label = axes.get_y_axis_label('y').shift(2 * UP)
        z_label = axes.get_z_axis_label('z',buff=LARGE_BUFF, rotation=0)

        labels = VGroup(x_label, y_label, z_label)
        self.add_fixed_orientation_mobjects(*labels)
        self.play(
            Create(axes),
            Write(labels)
        )
        self.wait()

        m_dot = Dot(shade_in_3d=True)
        line_r = DashedLine(ORIGIN, m_dot.get_center())

        def line_r_updater(line):
            line.become(
                DashedLine(
                    ORIGIN,
                    [m_dot.get_center()[0], m_dot.get_center()[1], 0]
                )
            )
        
        
        self.play(
            Create(m_dot),
            Create(line_r)
        )
        line_r.add_updater(line_r_updater)
        self.play(
            m_dot.animate.move_to(2.5 * UR)
        )
        line_r_t = MathTex("r_M").scale(0.8).move_to(m_dot.get_center() + 0.7 * DOWN)
        self.add_fixed_orientation_mobjects(line_r_t)
        self.play(
            FadeIn(line_r_t)
        )
        self.wait(0.7)

        angle = Angle(
            Line(ORIGIN, RIGHT),
            line_r,
            radius=0.8
        )
        angle_t = MathTex(r"\theta_M").scale(0.8).move_to(1.3 * RIGHT + 0.5 * UP)
        self.add_fixed_orientation_mobjects(angle_t)
        self.play(
            Create(angle),
            FadeIn(angle_t)
        )
        def angle_updater(angle):
            angle.become(
                Angle(
                    Line(ORIGIN, RIGHT),
                    line_r,
                    radius=angle.radius
                )
            )
        
        def r_m_updater(t):
            t.move_to(m_dot.get_center() + 0.7 * DOWN)
        
        angle.add_updater(angle_updater)
        line_r_t.add_updater(r_m_updater)
        self.wait()

        self.play(
            MoveAlongPath(
                m_dot,
                Arc(
                    radius=np.linalg.norm(m_dot.get_center()),
                    start_angle=PI / 4,
                    angle=PI / 8
                )
            ),
        )
        self.wait()
        self.play(
            MoveAlongPath(
                m_dot,
                Arc(
                    radius=np.linalg.norm(m_dot.get_center()),
                    start_angle=PI / 4 + PI / 8,
                    angle=-PI / 8
                )
            ),

        )
        self.wait()

        line_r_t.clear_updaters()

        line_z = DashedLine(m_dot, m_dot)
        line_z_t = MathTex("z_M").scale(0.8).next_to(line_z, 1.3 * UP)

        def line_z_updater(line):
            line.become(
                DashedLine([m_dot.get_center()[0], m_dot.get_center()[1], 0], m_dot.get_center())
            )

        def line_z_t_updater(t):
            t.next_to(line_z, 1.3 * UP)
        
        self.add(line_z)
        self.add_fixed_orientation_mobjects(line_z_t)
        self.play(
            FadeIn(line_z_t)
        )
        line_z.add_updater(line_z_updater)
        line_z_t.add_updater(line_z_t_updater)

        self.play(
            m_dot.animate.shift(1.1 * OUT)
        )
        m_t = Tex('M').scale(0.8).next_to(m_dot, UP)
        self.add_fixed_orientation_mobjects(m_t)
        self.play(
            FadeIn(m_t)
        )
        self.wait()

        m_def = MathTex(R"M(r_M, \theta_M, z_M)").scale(0.9).to_corner(UL)
        self.add_fixed_in_frame_mobjects(m_def)
        self.play(
            TransformFromCopy(Dot().move_to(0.72 * RIGHT + 0.24 * DOWN), m_def)
        )
        self.wait()

        om_line = DashedLine(ORIGIN, m_dot.get_center(), color=BLUE_D)
        om_t = MathTex(
            R"OM=\sqrt{{r_M}^2+{z_M}^2}",
            color=BLUE_D
        ).scale(0.8).next_to(m_def, DOWN, aligned_edge=LEFT)

        self.move_camera(
            phi=60 * DEGREES,
            theta=20 * DEGREES,
            added_anims=[Create(om_line)]
        )
        self.add_fixed_in_frame_mobjects(om_t)
        self.play(
            FadeIn(om_t)
        )
        self.wait(2)
        line_r.remove_updater(line_r_updater)
        line_z.remove_updater(line_z_updater)
        angle.remove_updater(angle_updater)
        self.move_camera(
            zoom=2,
            frame_center=2 * UP,
            added_anims=[
                *map(FadeOut, [m_def, om_line, om_t]),
                m_t.animate.scale(1.2).shift(0.2 * OUT + 0.7 * DOWN),
                m_dot.animate.scale(0.8)
            ]
        )
        self.wait(2)

        e_r, e_theta, e_z = create_cylindric_base(m_dot.get_center(), length=0.8)
        e_r_t = MathTex(r"\overrightarrow{e_r}").scale(0.8).next_to(e_r.get_tip(), 0.2 * RIGHT)
        e_theta_t = MathTex(r"\overrightarrow{e_\theta}").scale(0.8).next_to(e_theta.get_tip(), 1.5 * OUT)
        e_z_t = MathTex(r"\overrightarrow{e_z}").scale(0.8).next_to(e_z.get_tip(), OUT)
        cylindric_base = VGroup(e_r, e_theta, e_z)
        cylindric_base_t = VGroup(e_r_t, e_theta_t, e_z_t)
        self.add_fixed_orientation_mobjects(*cylindric_base_t)
        self.play(
            *map(Create, cylindric_base),
            *map(Create, cylindric_base_t)
        )
        self.wait(2)

        self.play(
            *map(FadeOut, [*cylindric_base, *cylindric_base_t, line_r_t, line_z_t])
        )
        self.wait()

        line_dr = Line(
            m_dot.get_center(),
            cylindric_to_cartesian(
                cartesian_to_cylindric(m_dot.get_center()) + 0.7 * RIGHT
            )
        )
        dr_t = MathTex(r"\mathrm{d}r").move_to(line_dr.get_center() + 0.2 * RIGHT)

        self.add_fixed_orientation_mobjects(dr_t)

        self.play(
            LaggedStart(
                Create(line_dr),
                FadeIn(dr_t),
                lag_ratio=0.3
            )
        )

        self.wait(0.7)

        drdtheta_surf = Surface(
            lambda u, v: ORIGIN,
            u_range=[0, 0],
            v_range=[0, 0]
        )

        angle_dtheta = Arc(
            start_angle=PI / 4,
            angle=PI / 14
        )
        theta_vt = ValueTracker(PI / 4)

        line_dtheta = DashedLine(
            ORIGIN,
            cylindric_to_cartesian([angle_dtheta.radius, theta_vt.get_value(), 0])
        )
        line_dtheta.add_updater(lambda l: l.become(
            DashedLine(
                ORIGIN,
                cylindric_to_cartesian([angle_dtheta.radius, theta_vt.get_value(), 0])
            )
        ))
        dtheta_t = MathTex(r"\mathrm{d}\theta").move_to(UP + 0.8 * RIGHT)


        def surf_updater(surf):
            surf.become(
                Surface(
                    func=lambda theta, r: np.array([
                        r * np.cos(theta),
                        r * np.sin(theta),
                        m_dot.get_center()[2]
                    ]),
                    u_range=[PI / 4, theta_vt.get_value()],
                    v_range=[
                        np.linalg.norm(m_dot.get_center()[:2]),
                        np.linalg.norm(m_dot.get_center()[:2]) + 0.7
                    ],
                    resolution=(1, 1),
                    checkerboard_colors=[BLUE_D],
                    fill_opacity=0.5,
                    stroke_width=4
                )
            )
        self.add(drdtheta_surf, line_dtheta)
        drdtheta_surf.add_updater(surf_updater)
        self.add_fixed_orientation_mobjects(dtheta_t)
        self.play(
            theta_vt.animate.increment_value(PI / 14),
            Create(angle_dtheta),
            FadeIn(dtheta_t)
        )
        
        self.wait()

        self.move_camera(
            phi=0,
            theta=0,
            zoom=1,
            frame_center=ORIGIN,
            added_anims=[
                FadeOut(z_label, run_time=1.5),
                dr_t.animate(run_time=1.5).shift(0.1 * RIGHT).scale(0.7),
                m_t.animate(run_time=1.5).shift(0.1 * RIGHT).scale(0.7),
                dtheta_t.animate(run_time=1.5).shift(0.05 * DL).scale(0.7)
            ],
            run_time=1.5
        )
        self.wait()

        arc_dtheta = Arc(
            radius=np.linalg.norm(m_dot.get_center()[:2]),
            start_angle=PI / 4,
            angle=PI / 14
        ).shift(m_dot.get_center()[2] * OUT)


        self.add(arc_dtheta)
        self.play(
            Indicate(arc_dtheta, scale_factor=2, color=YELLOW_G, run_time=1.5)
        )
        self.play(
            FadeOut(arc_dtheta)
        )
        self.wait()

        circle_ex = Arc(
            radius=np.linalg.norm(m_dot.get_center()[:2]),
            start_angle=PI / 4,
            angle=TAU,
            color=RED
        ).shift(m_dot.get_center()[2] * OUT)


        circ_t = MathTex(r"P=2\pi\times r", color=RED).to_corner(UL)
        new_angle_t = MathTex(r"\mathrm{d}\theta").next_to(circ_t[0][:2], aligned_edge=DOWN)

        self.play(
            Create(circle_ex)
        )
        self.wait()

        self.add_fixed_in_frame_mobjects(circ_t)
        self.play(
            FadeIn(circ_t, lag_ratio=0.1)
        )
        self.wait()
        self.play(
            Circumscribe(circ_t.get_part_by_tex(r"2\pi"))
        )
        self.wait()
        arc_dtheta.set_color(RED)
        self.add(arc_dtheta)
        self.add_fixed_in_frame_mobjects(new_angle_t)
        self.play(
            Uncreate(circle_ex),
            FadeOut(circ_t[0][2:4]),
            FadeIn(new_angle_t),
            run_time=1.8
        )
        self.wait()

        self.move_camera(
            phi=60 * DEGREES,
            theta=20 * DEGREES,
            zoom=2,
            frame_center=2 * UP,
            added_anims=[
                FadeIn(z_label, run_time=1.5),
                dr_t.animate(run_time=1.5).shift(0.1 * LEFT).scale(1 / 0.7),
                m_t.animate(run_time=1.5).shift(0.1 * LEFT).scale(1 / 0.7),
                dtheta_t.animate(run_time=1.5).shift(0.05 * UR).scale(1 / 0.7)
            ],
            run_time=1.5
        )
        self.wait(0.5)
        rdtheta_t = MathTex(r"r\mathrm{d}\theta").move_to(dr_t.get_center() + LEFT + 0.3 * DOWN)
        self.add_fixed_orientation_mobjects(rdtheta_t)
        self.play(
            FadeOut(arc_dtheta),
            FadeOut(circ_t[0][:2]),
            FadeOut(circ_t[0][4:]),
            FadeOut(new_angle_t),
            FadeIn(rdtheta_t)
        )
        self.wait()

        z_vt = ValueTracker(0)

        surfz1 = always_redraw(lambda:
            Surface(
                func=lambda theta, z: np.array([
                    np.linalg.norm(m_dot.get_center()[:2]) * np.cos(theta),
                    np.linalg.norm(m_dot.get_center()[:2]) * np.sin(theta),
                    z
                ]),
                u_range=[PI / 4, theta_vt.get_value()],
                v_range=[
                    m_dot.get_center()[2],
                    m_dot.get_center()[2] + z_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )
        surfz2 = always_redraw(lambda:
            Surface(
                func=lambda theta, z: np.array([
                    (np.linalg.norm(m_dot.get_center()[:2]) + 0.7) * np.cos(theta),
                    (np.linalg.norm(m_dot.get_center()[:2]) + 0.7) * np.sin(theta),
                    z
                ]),
                u_range=[PI / 4, theta_vt.get_value()],
                v_range=[
                    m_dot.get_center()[2],
                    m_dot.get_center()[2] + z_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )
        surfz3 = always_redraw(lambda:
            Surface(
                func=lambda r, z: np.array([
                    r * np.cos(PI / 4),
                    r * np.sin(PI / 4),
                    z
                ]),
                u_range=[
                    np.linalg.norm(m_dot.get_center()[:2]),
                    np.linalg.norm(m_dot.get_center()[:2]) + 0.7
                ],
                v_range=[
                    m_dot.get_center()[2],
                    m_dot.get_center()[2] + z_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )
        surfz4 = always_redraw(lambda:
            Surface(
                func=lambda r, z: np.array([
                    r * np.cos(theta_vt.get_value()),
                    r * np.sin(theta_vt.get_value()),
                    z
                ]),
                u_range=[
                    np.linalg.norm(m_dot.get_center()[:2]),
                    np.linalg.norm(m_dot.get_center()[:2]) + 0.7
                ],
                v_range=[
                    m_dot.get_center()[2],
                    m_dot.get_center()[2] + z_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )
        surfz5 = always_redraw(lambda:
            Surface(
                func=lambda theta, r: np.array([
                    r * np.cos(theta),
                    r * np.sin(theta),
                    m_dot.get_center()[2] + z_vt.get_value()
                ]),
                u_range=[PI / 4, theta_vt.get_value()],
                v_range=[
                    np.linalg.norm(m_dot.get_center()[:2]),
                    np.linalg.norm(m_dot.get_center()[:2]) + 0.7
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )

        surfaces = VGroup(surfz1, surfz2, surfz3, surfz4, surfz5)

        z_t = MathTex(R"\mathrm{d}z").move_to(dr_t.get_center() + 0.3 * OUT + 0.1 * UP)

        self.add(*surfaces)
        self.play(
            z_vt.animate.set_value(0.5)
        )
        self.add_fixed_orientation_mobjects(z_t)
        self.play(
            FadeIn(z_t)
        )
        self.wait()

        drdtheta_surf.clear_updaters()
        self.play(
            Indicate(drdtheta_surf, color=YELLOW_G),
            Indicate(dv_t[0][4:10], color=YELLOW_G),
            run_time=2
        )
        self.wait()

        dv_t = MathTex(
            R"\mathrm{d}^3V=\mathrm{d}r.r\mathrm{d}\theta.\mathrm{d}z"
        ).move_to(3 * RIGHT + 1.5 * DOWN)

        self.add_fixed_in_frame_mobjects(dv_t)
        self.play(
            FadeIn(dv_t, lag_ratio=0.1)
        )
        self.wait()


class Anim14_4(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=35 * DEGREES
        )
        axes = ThreeDAxes().set_stroke(width=0.6)
        x_label = axes.get_x_axis_label('x')
        y_label = axes.get_y_axis_label('y').shift(2 * UP)
        z_label = axes.get_z_axis_label('z',buff=LARGE_BUFF, rotation=0)

        labels = VGroup(x_label, y_label, z_label)
        self.add_fixed_orientation_mobjects(*labels)
        self.play(
            Create(axes),
            Write(labels)
        )
        self.wait()
        
        theta_vt = ValueTracker()
        r_vt = ValueTracker(3)
        phi_vt = ValueTracker(70 * DEGREES)

        m_dot = always_redraw(lambda: Dot(
            spherical_to_cartesian([r_vt.get_value(), phi_vt.get_value(), theta_vt.get_value()]),
            shade_in_3d=True
        ))
        m_t = Tex('M').next_to(m_dot, UP)
        def m_t_updater(t):
            t.next_to(m_dot, UP)
            self.add_fixed_orientation_mobjects(t)
        line_r = always_redraw(lambda: Line(ORIGIN, m_dot.get_center()))

        angle_theta = always_redraw(lambda: Arc(
            radius=0.8,
            angle=theta_vt.get_value()
        ).rotate(
            PI / 2,
            axis=DOWN,
            about_point=UP
        ).rotate(
            PI / 2 - phi_vt.get_value(),
            axis=IN,
            about_point=IN
        ))

        self.add_fixed_orientation_mobjects(m_t)
        self.play(
            *map(FadeIn, [m_dot, line_r, angle_theta, m_t])
        )
        m_t.add_updater(m_t_updater)
        self.wait(0.3)
        self.play(
            theta_vt.animate.set_value(PI / 2 - PI / 7)
        )
        theta_t = MathTex(R"\theta_M").move_to(0.5 * UP + OUT).scale(0.8)
        self.add_fixed_orientation_mobjects(theta_t)
        self.play(
            FadeIn(theta_t)
        )
        self.wait()
        self.play(
            r_vt.animate.increment_value(1)
        )
        self.wait(0.2)
        r_t = MathTex("r_M").move_to(0.9 * UP + 0.2 * OUT).scale(0.8)
        self.add_fixed_orientation_mobjects(r_t)
        self.play(
            r_vt.animate.increment_value(-1),
            FadeIn(r_t)
        )
        self.wait()

        line_z = always_redraw(lambda: DashedLine(
            [*m_dot.get_center()[:2], 0],
            m_dot.get_center()
        ))
        line_xy = always_redraw(lambda: DashedLine(
            ORIGIN,
            [*m_dot.get_center()[:2], 0]
        ))
        angle_phi = always_redraw(lambda: Angle(
            Line(ORIGIN, RIGHT),
            line_xy,
            radius=0.8
        ))
        phi_t = MathTex(R"\varphi_M").move_to(UR).scale(0.8)
        self.play(
            FadeIn(line_z),
            FadeIn(line_xy)
        )
        self.add_fixed_orientation_mobjects(phi_t)
        self.play(
            FadeIn(angle_phi),
            FadeIn(phi_t)
        )
        self.wait()
        self.play(
            phi_vt.animate.set_value(120 * DEGREES)
        )
        self.wait(0.2)
        self.play(
            phi_vt.animate.set_value(70 * DEGREES)
        )
        self.wait()

        m_def = MathTex(R"M(r_M, \theta_M, \varphi_M)").scale(0.9).to_corner(UL)
        self.add_fixed_in_frame_mobjects(m_def)
        self.play(
            TransformFromCopy(Dot().move_to(1.72 * RIGHT + 0.43 * UP).scale(0.2), m_def)
        )
        self.wait()

        om_t = MathTex(
            R"OM=r_m",
            color=BLUE_D
        ).next_to(m_def, DOWN, aligned_edge=LEFT)
        self.add_fixed_in_frame_mobjects(om_t)
        self.play(
            FadeIn(om_t)
        )
        self.wait()

        self.play(
            self.camera.zoom_tracker.animate.set_value(2),
            self.camera._frame_center.animate.move_to(2 * UP + OUT),
            theta_vt.animate.set_value(PI / 4),
            r_t.animate.shift(0.6 * OUT),
            run_time=1.5
        )
        self.wait()

        e_r, e_theta, e_phi = create_spheric_base(m_dot.get_center())
        e_r_t = MathTex(R"\overrightarrow{e_r}").scale(0.8).next_to(e_r.get_tip(), 0.4 * OUT)
        e_theta_t = MathTex(R"\overrightarrow{e_\theta}").scale(0.8).next_to(e_theta.get_tip(), 0.35 * UP)
        e_phi_t = MathTex(R"\overrightarrow{e_\varphi}").scale(0.8).next_to(e_phi.get_tip(), OUT)

        spheric_base = VGroup(e_r, e_theta, e_phi)
        spheric_base_t = VGroup(e_r_t, e_theta_t, e_phi_t)
        self.add_fixed_orientation_mobjects(*spheric_base_t)
        self.play(
            *map(Create, spheric_base),
            *map(Create, spheric_base_t)
        )
        self.wait(2)
        self.play(
            *map(FadeOut, [*spheric_base, *spheric_base_t, om_t, m_def])
        )
        self.play(
            theta_vt.animate.set_value(PI / 2 - PI / 7),
            self.camera.theta_tracker.animate.set_value(20 * DEGREES),
            self.camera.phi_tracker.animate.set_value(80 * DEGREES),
            r_t.animate.shift(0.55 * IN + 0.2 * UP)
        )
        self.wait()

        dtheta_vt = ValueTracker(theta_vt.get_value())
        dphi_vt = ValueTracker(phi_vt.get_value())
        dr_vt = ValueTracker(r_vt.get_value())
        surf_dthetadphi = always_redraw(lambda: Surface(
            func=lambda theta, phi: np.array([
                r_vt.get_value() * np.sin(theta) * np.cos(phi),
                r_vt.get_value() * np.sin(theta) * np.sin(phi),
                r_vt.get_value() * np.cos(theta)
            ]),
            u_range=[theta_vt.get_value(), dtheta_vt.get_value()],
            v_range=[
                phi_vt.get_value(),
                dphi_vt.get_value()
            ],
            resolution=(1, 1),
            checkerboard_colors=[BLUE_D],
            fill_opacity=0.5,
            stroke_width=4
        ))

        angle_dtheta = always_redraw(lambda: Arc(
            radius=0.6,
            start_angle=theta_vt.get_value(),
            angle=dtheta_vt.get_value() - theta_vt.get_value()
        ).rotate(
            PI / 2,
            axis=DOWN,
            about_point=UP
        ).rotate(
            PI / 2 - phi_vt.get_value(),
            axis=IN,
            about_point=IN
        ))
        dtheta_t = MathTex(R"\mathrm{d}\theta").scale(0.7).move_to(angle_dtheta).shift(0.2 * UP)


        self.add(surf_dthetadphi, angle_dtheta)
        self.add_fixed_orientation_mobjects(dtheta_t)
        self.play(
            dtheta_vt.animate.increment_value(PI / 14),
            FadeIn(dtheta_t)
        )
        self.wait()
        m_t.clear_updaters()
        self.move_camera(
            phi=70 * DEGREES,
            theta=60 * DEGREES,
            added_anims=[
                m_t.animate.shift(0.4 * OUT),
                r_t.animate.shift(0.3 * RIGHT + 0.2 * OUT),
                FadeOut(angle_dtheta),
                FadeOut(dtheta_t)
            ]
        )
        self.wait(0.5)

        angle_dphi = always_redraw(lambda: Arc(
            radius=0.6,
            start_angle=phi_vt.get_value(),
            angle=dphi_vt.get_value() - phi_vt.get_value()
        ))
        dphi_t = MathTex(R"\mathrm{d}\varphi").scale(0.6).move_to(angle_dphi).shift(0.3 * UP + 0.1 * LEFT)
        self.add_fixed_orientation_mobjects(dphi_t)
        self.add(angle_dphi)
        self.play(
            dphi_vt.animate.increment_value(PI / 12),
            FadeIn(dphi_t)
        )
        self.wait()
        self.play(
            Indicate(line_xy, color=YELLOW_G)
        )
        self.wait()
        self.move_camera(
            phi=65 * DEGREES,
            theta=15 * DEGREES,
            added_anims=[
                FadeOut(angle_dphi, dphi_t),
                m_t.animate.shift(0.6 * DOWN + 0.2 * IN),
                r_t.animate.shift(0.2 * IN)
            ]
        )
        self.wait()

        surfr1 = always_redraw(lambda:
            Surface(
                func=lambda theta, phi: np.array([
                    dr_vt.get_value() * np.sin(theta) * np.cos(phi),
                    dr_vt.get_value() * np.sin(theta) * np.sin(phi),
                    dr_vt.get_value() * np.cos(theta)
                ]),
                u_range=[theta_vt.get_value(), dtheta_vt.get_value()],
                v_range=[
                    phi_vt.get_value(),
                    dphi_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )
        surfr2 = always_redraw(lambda:
            Surface(
                func=lambda r, theta: np.array([
                    r * np.sin(theta) * np.cos(phi_vt.get_value()),
                    r * np.sin(theta) * np.sin(phi_vt.get_value()),
                    r * np.cos(theta)
                ]),
                u_range=[r_vt.get_value(), dr_vt.get_value()],
                v_range=[
                    theta_vt.get_value(),
                    dtheta_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )
        surfr3 = always_redraw(lambda:
            Surface(
                func=lambda r, theta: np.array([
                    r * np.sin(theta) * np.cos(dphi_vt.get_value()),
                    r * np.sin(theta) * np.sin(dphi_vt.get_value()),
                    r * np.cos(theta)
                ]),
                u_range=[r_vt.get_value(), dr_vt.get_value()],
                v_range=[
                    theta_vt.get_value(),
                    dtheta_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )
        surfr4 = always_redraw(lambda:
            Surface(
                func=lambda r, phi: np.array([
                    r * np.sin(theta_vt.get_value()) * np.cos(phi),
                    r * np.sin(theta_vt.get_value()) * np.sin(phi),
                    r * np.cos(theta_vt.get_value())
                ]),
                u_range=[r_vt.get_value(), dr_vt.get_value()],
                v_range=[
                    phi_vt.get_value(),
                    dphi_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )
        surfr5 = always_redraw(lambda:
            Surface(
                func=lambda r, phi: np.array([
                    r * np.sin(dtheta_vt.get_value()) * np.cos(phi),
                    r * np.sin(dtheta_vt.get_value()) * np.sin(phi),
                    r * np.cos(dtheta_vt.get_value())
                ]),
                u_range=[r_vt.get_value(), dr_vt.get_value()],
                v_range=[
                    phi_vt.get_value(),
                    dphi_vt.get_value()
                ],
                resolution=(1, 1),
                checkerboard_colors=[BLUE_D],
                fill_opacity=0.5,
                stroke_width=4
            )
        )

        surfaces = VGroup(surfr1, surfr2, surfr3, surfr4, surfr5)
        self.add(*surfaces)
        self.play(
            dr_vt.animate.increment_value(0.25)
        )
        self.wait()

        dv_t = MathTex(
            R"\mathrm{d}^3V=\mathrm{d}r.r\mathrm{d}\theta.r\sin{\theta}\mathrm{d}\varphi"
        ).move_to(3 * RIGHT + 2 * UP)

        self.add_fixed_in_frame_mobjects(dv_t)
        self.play(
            FadeIn(dv_t, lag_ratio=0.1)
        )
        self.wait()

        self.play(
            Indicate(
                Line(
                    m_dot.get_center(),
                    spherical_to_cartesian(cartesian_to_spherical(m_dot.get_center()) + 0.25 * RIGHT),
                    stroke_width=4
                ),
                color=YELLOW_G,
                scale_factor=1.8
            ),
            Indicate(dv_t[0][4:6], color=YELLOW_G, scale_factor=1.8),
            run_time=2
        )
        self.play(
            Indicate(
               Surface(
                    func=lambda theta, phi: np.array([
                        r_vt.get_value() * np.sin(theta) * np.cos(phi),
                        r_vt.get_value() * np.sin(theta) * np.sin(phi),
                        r_vt.get_value() * np.cos(theta)
                    ]),
                    u_range=[theta_vt.get_value(), dtheta_vt.get_value()],
                    v_range=[phi_vt.get_value(), phi_vt.get_value()],
                    resolution=(1, 1),
                    checkerboard_colors=[BLUE_D],
                    fill_opacity=0.5,
                    stroke_width=6
                ),
                color=YELLOW_G,
                scale_factor=1.8
            ),
            Indicate(dv_t[0][7:10], color=YELLOW_G, scale_factor=1.8),
            run_time=2
        )
        self.play(
            Indicate(
                Surface(
                    func=lambda theta, phi: np.array([
                        r_vt.get_value() * np.sin(theta) * np.cos(phi),
                        r_vt.get_value() * np.sin(theta) * np.sin(phi),
                        r_vt.get_value() * np.cos(theta)
                    ]),
                    u_range=[theta_vt.get_value(), theta_vt.get_value()],
                    v_range=[phi_vt.get_value(), dphi_vt.get_value()],
                    resolution=(1, 1),
                    checkerboard_colors=[BLUE_D],
                    fill_opacity=0.5,
                    stroke_width=6
                ),
                color=YELLOW_G,
                scale_factor=1.8
            ),
            Indicate(dv_t[0][11:], color=YELLOW_G, scale_factor=1.8),
            run_time=2
        )
        self.wait()




        # self.interactive_embed()
