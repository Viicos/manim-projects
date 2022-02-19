"""
Behind the Gauss theorem.
Currently pending.
"""

from manim import *
from manim_physics import *
import numpy as np
import math
import os

YELLOW_G = "#FCBA03"

class Introduction(Scene):
    def construct(self):
        esint_template = TexTemplate()
        esint_template.add_to_preamble(r'\usepackage{esint}')

        treatise_picture = ImageMobject(
            os.path.join('assets', 'A_Treatise_on_Electricity_and_Magnetism_Volume_1_003.jpg')
        )
        maxwell_picture = ImageMobject(
            os.path.join('assets', 'PSM_V17_D008_James_Clerk_Maxwell.jpg')
        )
        treatise_picture.width = 4
        treatise_picture.move_to(3.5 * LEFT)
        maxwell_picture.width = 4
        maxwell_picture.next_to(
            treatise_picture,
            direction=RIGHT,
            buff=2.5 * RIGHT,
            aligned_edge=UP
        )

        james_name = Tex(
            r"James Clerk Maxwell\\(1832\textendash 1879)"
        ).next_to(maxwell_picture, DOWN).align_to(treatise_picture, DOWN)
        self.play(
            FadeIn(maxwell_picture),
            FadeIn(treatise_picture),
            FadeIn(james_name),
            run_time=2
        )
        self.wait()

        eqs = MathTex(
            r"&\vec{\nabla}\cdot\vec{E}=\frac{\rho}{\varepsilon_0}\\",
            r"&\vec{\nabla}\cdot\vec{B}=0\\",
            r"&\vec{\nabla}\times\vec{E}=-\frac{\partial\vec{B}}{\partial t}\\",
            r"&\vec{\nabla}\times\vec{B}=\mu_0\vec{j} + \mu_0\varepsilon_0\frac{\partial\vec{E}}{\partial t}"
        ).move_to(maxwell_picture).shift(DOWN).scale(0.8)
        brace = Brace(eqs, LEFT)
        self.play(
            FadeOut(maxwell_picture),
            FadeOut(james_name)
        )
        self.play(
            FadeIn(eqs, lag_ratio=0.1, run_time=2.5),
            FadeIn(brace, run_time=2)
        )
        self.wait()

        self.play(
            Indicate(
                eqs.get_part_by_tex(
                    r"\vec{\nabla}\cdot\vec{E}=\frac{\rho}{\varepsilon_0}"
                ),
                scale_factor=1.3,
                color=YELLOW_G,
                run_time=1.5
            )
        )
        self.wait()

        self.play(
            *map(FadeOut, [treatise_picture, eqs, brace])
        )
        self.wait()
        int_eq = MathTex(
            r"{{\oiint_S}}{{ \overrightarrow{E} }}.\overrightarrow{\mathrm{d^2S}} = \frac{1}{ {{ \varepsilon_0 }} }{{ \iiint_V }}{{ \rho }}\mathrm{d^3V}",
            tex_template=esint_template
        ).scale(1.3)
        self.play(
            FadeIn(int_eq, run_time=1.5)
        )
        self.wait()

        self.play(
            LaggedStartMap(
                lambda m: Indicate(m, color=YELLOW_G),
                [int_eq[0], int_eq[1], int_eq[3], int_eq[5], int_eq[6]],
                lag_ratio=1
            ),
            run_time=4
        )
        self.wait()
        self.play(
            int_eq.animate(run_time=1.5, rate_func=there_and_back, lag_ratio=0.04).set_color(YELLOW_G)
        )
        self.wait(0.3)

        behind_t1 = Tex(
            "Derrière le théorème",
            font_size=1.3 * DEFAULT_FONT_SIZE
        ).shift(UP)
        behind_t2 = Tex(
            "de Gauss",
            font_size = 3.25 * DEFAULT_FONT_SIZE
        ).next_to(behind_t1, DOWN, aligned_edge=RIGHT)

        charge1 = Charge(point=2.65 * RIGHT + 0.72 * DOWN)
        charge2 = Charge(magnitude=-1, point=3.2 * LEFT + 0.2 * UP)
        electric_field = ElectricField(charge1, charge2)


        self.play(
            FadeOut(int_eq)
        )
        self.add_foreground_mobjects(behind_t1, behind_t2, charge1, charge2)
        self.play(
            LaggedStart(
                AnimationGroup(
                    FadeIn(behind_t1, shift=DOWN),
                    FadeIn(behind_t2, shift=RIGHT),
                    run_time=2
                ),
                AnimationGroup(
                    FadeIn(charge1),
                    FadeIn(charge2),
                    FadeIn(electric_field)
                ),
                lag_ratio=0.5
            )
        )
        def update_charge1(c, dt):
            c.shift(0.05 * dt * UR)

        def update_charge2(c, dt):
            c.shift(0.05 * dt * DL)
        
        def update_efield(efield):
            efield.become(ElectricField(charge1, charge2))
        charge1.add_updater(update_charge1)
        charge2.add_updater(update_charge2)
        electric_field.add_updater(update_efield)
        self.wait(4)

class ConceptsAbordes(MovingCameraScene):
    def construct(self):
        line_v = Line(3.5 * UP, 3.5 * DOWN, stroke_width=1.5 * DEFAULT_STROKE_WIDTH)
        line_h = Line(7 * LEFT, 7 * RIGHT, stroke_width=1.5 * DEFAULT_STROKE_WIDTH)
        text_1 = Tex(
            r"Champs scalaires\\et vectoriels",
            font_size=0.8 * DEFAULT_FONT_SIZE
        ).to_corner(UL).shift(1.4 * RIGHT + 0.4 * UP)
        text_2 = Tex(
            r"Systèmes de coordonnées",
            font_size=0.8 * DEFAULT_FONT_SIZE
        ).to_corner(UR).shift(0.7 * LEFT + 0.4 * UP)
        text_3 = Tex(
            r"Volumes élémentaires",
            font_size=0.8 * DEFAULT_FONT_SIZE
        ).to_corner(DL).shift(3 * UP + 0.9 * RIGHT)
        text_4 = Tex(
            r"Divergence",
            font_size=0.8 * DEFAULT_FONT_SIZE
        ).to_corner(DR).shift(2.9 * UP + 2.1 * LEFT)
        self.play(
            GrowFromCenter(line_v),
            GrowFromCenter(line_h)
        )
        self.play(
            *map(lambda t: FadeIn(t, lag_ratio=0.1), [text_1, text_2, text_3, text_4]),
            run_time=0.8
        )
        self.wait()

        self.play(
            self.camera.frame.animate.scale(0.5).move_to(
                config["frame_width"] / 3.9 * LEFT + config["frame_height"] / 3.9 * UP
            )
        )
        self.wait(3)


class Champs(MovingCameraScene):
    def construct(self):
        self.camera.frame.save_state()
        axes = Axes(
            x_length=14,
            y_length=8,
            tips=False
        )
        plane = NumberPlane(
            axis_config={"include_numbers": True},
            background_line_style={
                "stroke_opacity": 0.4
            }
        )
        arrow_m = Arrow(
            ORIGIN,
            RIGHT + 1.5 * UP,
            buff=0,
            color=RED
        )
        point_m = Dot(RIGHT + 1.5 * UP)
        m_t = Tex('M').next_to(arrow_m.get_tip(), RIGHT).shift(0.5 * UP)
        self.play(
            Create(plane),
            Create(axes, run_time=0.8)
        )
        self.wait()
        self.play(
            LaggedStart(
                GrowArrow(arrow_m),
                Create(point_m),
                Write(m_t),
                lag_ratio=0.5
            )
        )
        self.wait(0.2)
        arrow_func = Arrow(
            ORIGIN,
            0.5 * RIGHT,
            buff=0,
            max_stroke_width_to_length_ratio=10
        ).next_to(m_t)
        ruler = SVGMobject(
            os.path.join("assets", "ruler.svg"),
            width=1.3  
        ).next_to(arrow_func)
        self.play(
            LaggedStart(
                self.camera.frame.animate.move_to(2 * RIGHT + 1.5 * UP).scale(0.6),
                GrowArrow(arrow_func),
                Create(ruler, lag_ratio=0.3),
                lag_ratio=0.5
            )
        )
        self.wait()

        number = Tex("12,3").next_to(arrow_func)
        vector2d = Matrix([[2], [1]]).next_to(arrow_func)
        self.play(
            FadeOut(ruler, shift=DOWN),
            FadeIn(number, shift=DOWN)
        )
        self.wait()
        self.play(
            FadeOut(number, shift=DOWN),
            FadeIn(vector2d, shift=DOWN)
        )
        self.wait(2)

        def scalar_field_func(p):
            return 30 / np.linalg.norm(p)
        
        def m_value_updater(m):
            m.next_to(point_m)
            m.set_value(scalar_field_func(point_m.get_center()))

        m_value = DecimalNumber(
            scalar_field_func(point_m.get_center()),
            unit='^\\circ C',
            num_decimal_places=1,
            include_background_rectangle=False
        ).next_to(point_m)
        m_value.add_updater(m_value_updater)
        self.play(
            Restore(self.camera.frame),
            *map(FadeOut, [vector2d, arrow_func, arrow_m, m_t]),
            FadeIn(m_value)
        )
        self.wait()
        self.play(
            Succession(
                ApplyMethod(point_m.move_to, 2 * DR, run_time=1.3),
                ApplyMethod(point_m.move_to, 2 * LEFT + 0.5 * UP, run_time=1.3),
                ApplyMethod(point_m.move_to, 4 * LEFT + 2 * DOWN, run_time=1.3),
                lag_ratio=1.3
            )
        )
        self.wait()



        
        

        




