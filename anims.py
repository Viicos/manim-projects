"""Animations from 2019-2020."""

from manim import *
from utils import *
import random
import numpy as np
import math

class Test(Scene):
    def construct(self):
        esint_template = TexTemplate()
        esint_template.add_to_preamble(r'\usepackage{esint}')
        flux_el_t1 = MathTex(
            r"{{\delta\Phi=}} {{\overrightarrow{E} . \overrightarrow{\mathrm{d^2S}}}}",
            tex_template=esint_template,
            tex_to_color_map={
                r"\overrightarrow{E}": "#FCBA03"
            }
        )
        flux_el_t1_int = MathTex(
            r"\iint_S {{\delta\Phi=}} \iint_S {{\overrightarrow{E} . \overrightarrow{\mathrm{d^2S}}}}",
            tex_template=esint_template,
            tex_to_color_map={
                r"\overrightarrow{E}": "#FCBA03"
            }
        )
        self.add(flux_el_t1)
        self.play(
            TransformMatchingTex(flux_el_t1, flux_el_t1_int)
        )
        self.wait()


class DensiteGlissante(ThreeDScene):
    @staticmethod
    def get_random_coords(x_range, distance=0.4):
        x = random.uniform(-distance, distance)
        y = random.uniform(*x_range)
        z = random.uniform(-distance, distance)

        return (x, y, z)
    
    @staticmethod
    def shift_sphere(sphere, dt):
        sphere.move_to(sphere.get_center() + (0, dt, 0))

    def construct(self):
        
        random.seed()

        self.set_camera_orientation(
            phi=80 * DEGREES,
            theta=-100 * DEGREES
        )
        axes = ThreeDAxes(y_range=(-8, 8, 1))

        self.add(axes)

        spheres_coords = [self.get_random_coords((-8, -1)) for _ in range(40)] + [self.get_random_coords((1, 8)) for _ in range(40)]

        np.save('spheres_coords', spheres_coords)

        spheres_el = [
            Dot3D(
                resolution=(6, 6),
                color=RED,
                fill_opacity=0.6
            ).move_to(coords)
            for coords in spheres_coords
        ]

        sphere_vol = Sphere(
            radius=0.5,
            fill_color=BLUE,
            fill_opacity=0.5
        ).move_to(axes.y_range[1] * DOWN)

        self.add(*spheres_el, sphere_vol)

        self.wait()

        sphere_vol.add_updater(self.shift_sphere)

        self.move_camera(
            theta=-260 * DEGREES,
            run_time=16,
            rate_func=lambda t: smooth(t, inflection=5.0)
        )

        self.wait()

class DensiteGlissanteGraph(Scene):

    @staticmethod
    def num_points_in_sphere(x, radius, sphere_coords):
        cp = x * UP
        return sum(map(lambda c: math.dist(cp, c) <= radius, sphere_coords))

    def construct(self):
        axes = Axes(
            x_range=[0, 8, 8],
            y_range=[0, 8, 6],
            tips=False,
            x_length=2.7,
            y_length=2.3,
            axis_config={"stroke_width": 4}
        )
        axes_label = axes.get_axis_labels(y_label=r"\rho(x)")
        self.play(
            Create(axes),
            Create(axes_label)
        )
        sphere_coords = np.load("I:\manim-cm-projects\spheres_coords.npy")
        
        x_array = np.arange(0, 8.5, 0.5)
        y_array = [self.num_points_in_sphere((x - 4) * 2, math.sqrt(2) - 0.5, sphere_coords) for x in x_array]

        pfit = np.polyfit(x_array, y_array, 10)
        pfunc = axes.plot(
            np.poly1d(pfit) / 2 - 0.7,
            color=BLUE
        )
        self.play(
            Create(pfunc),
            run_time=16,
            rate_func = lambda t: smooth(t, inflection=5.0)
        )
        self.wait()

class DensiteGlissanteZoom(ThreeDScene):
    @staticmethod
    def get_random_coords(x_range, distance=0.4):
        x = random.uniform(*x_range)
        y = random.uniform(-distance, distance)
        z = random.uniform(-distance, distance)

        return (x, y, z)
    
    @staticmethod
    def color_updater(sphere, sphere_dv):
        ct_dv = sphere_dv.get_center()
        ct = sphere.get_center()
        if math.dist(ct_dv, ct) <= sphere_dv.radius - sphere.radius:
            new_sphere = Dot3D(
                point=ct,
                radius=sphere.radius,
                color=GREEN_E,
                fill_opacity=0.6
            )
        else:
            new_sphere = Dot3D(
                point=ct,
                radius=sphere.radius,
                color=RED,
                fill_opacity=0.6
            )
        sphere.become(new_sphere)
    
    @staticmethod
    def num_points_in_sphere(sphere, sphere_coords):
        return sum(map(lambda c: math.dist(sphere.get_center(), c) <= sphere.radius - 0.25, sphere_coords))
    
    def construct(self):
        random.seed()

        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=100 * DEGREES
        )

        axes = ThreeDAxes().set_stroke(width=4)

        sphere = Sphere(
            center=3 * LEFT,
            radius=1.5,
            fill_opacity=0.4
        )

        spheres_coords = [self.get_random_coords([-5, 5]) for _ in range(14)]

        spheres_el = [
            Dot3D(
                point=coords,
                radius=0.25,
                resolution=(8, 8),
                color=RED,
                fill_opacity=0.6
            ).set_color(RED)\
                .add_updater(lambda s: self.color_updater(s, sphere))
            for coords in spheres_coords
        ]

        point_m = Dot(3 * LEFT)
        point_m_t = Tex('M').move_to(3 * LEFT + DOWN)
        density_t = MathTex(r"\rho(M)=\frac{\delta^{3}q(M)}{\mathrm{d^{3}V}}=").move_to(2 * UP + RIGHT)
        density_t_frac = MathTex(r"\frac{\quad e^{+}}{\mathrm{d^{3}V}}").next_to(density_t)
        
        q_vt = ValueTracker(self.num_points_in_sphere(sphere, spheres_coords))
        q_value = DecimalNumber(
            number=self.num_points_in_sphere(sphere, spheres_coords),
            num_decimal_places=0,
            font_size=44
        ).next_to(density_t, aligned_edge=UP).shift(0.12 * DR)

        q_vt.add_updater(
            lambda m: m.set_value(self.num_points_in_sphere(sphere, spheres_coords))
        )

        def q_value_updater(m):
            m.set_value(self.num_points_in_sphere(sphere, spheres_coords))
            self.add_fixed_in_frame_mobjects(m)

        self.add_fixed_orientation_mobjects(point_m_t)
        self.add_fixed_in_frame_mobjects(density_t, density_t_frac, q_value)


        self.play(
            Create(axes.x_axis),
            Create(sphere),
            Write(point_m),
            Write(point_m_t),
            Write(density_t),
            Write(density_t_frac),
            Write(q_value),
            *map(FadeIn, spheres_el)	
        )
        q_value.add_updater(q_value_updater)
        self.wait()

        self.play(
            (sphere + point_m + point_m_t).animate.shift(6 * RIGHT),
            run_time=8,
            rate_func=lambda t: smooth(t, inflection=5.0)
        )
        self.wait()

class DistributionCarre(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=40 * DEGREES
        )

        axes = ThreeDAxes().set_stroke(width=0.5)

        sphere = Sphere(
            radius=2.5,
            resolution=(18, 18),
            fill_opacity=0.4
        ).set_color(BLUE_E)

        ct_p = Dot().scale(1.2)
        ct_p_t = Tex('O').shift(0.3 * RIGHT + 0.15 * DOWN)
        self.add_fixed_orientation_mobjects(ct_p_t)

        self.play(
            *map(Create, [axes, sphere, ct_p]),
            Write(ct_p_t)
        )
        self.wait()

        self.move_camera(
            phi=90 * DEGREES,
            theta=90 * DEGREES,
            run_time=1.7
        )
        self.wait()

        # Partie 2D :
        graph_axes = Axes(
            x_range=[0, 1, 0.5],
            y_range=[0, 1, 0.5],
            tips=False,
            x_length=2.7,
            y_length=2.3,
            axis_config={"stroke_width": 4}
        ).move_to((3.5, 2, 0))

        axis_labels = graph_axes.get_axis_labels(x_label="r", y_label=r"\rho(r)")

        circles = [
            Circle(
                radius=r,
                color=RED_C,
                stroke_opacity=(r**4 / 40)
            )
            for r in np.arange(0, 2.5, 1 / 70)
        ]

        self.add_fixed_in_frame_mobjects(*circles)
        self.play(
            FadeOut(sphere),
            FadeOut(ct_p_t),
            *map(FadeIn, circles)
        )
        self.wait()

        self.add_fixed_in_frame_mobjects(graph_axes, *axis_labels)
        self.play(
            Create(graph_axes),
            Write(axis_labels)
        )
        self.wait()

        r_func = graph_axes.plot(
            lambda r: r**2,
            x_range=[0, 1],
            color=BLUE
        )
        f_label = graph_axes.get_graph_label(
            graph=r_func,
            label=r"\rho_{0}\left(\frac{r}{R}\right)^{2}",
            color=BLUE
        )

        self.add_fixed_in_frame_mobjects(r_func, f_label)
        self.play(
            Create(r_func),
            Write(f_label)
        )
        self.wait()

        self.wait()
        self.play(
            *map(FadeOut, circles + [graph_axes, axis_labels, r_func, f_label]),
            FadeIn(sphere),
            FadeIn(ct_p_t)
        )
        self.wait()
        self.move_camera(
            phi=70 * DEGREES,
            theta=50 * DEGREES,
            run_time=1.4
        )
        self.wait()

        m_coords = np.array([1, 3.5, 2.5])
        m = Dot(m_coords).scale(1.2)
        m_t = Tex('M').move_to((1.4, 3.5, 2.5))
        self.add_fixed_orientation_mobjects(m_t)

        self.play(
            Create(m),
            Write(m_t)
        )
        self.wait()

        spheric_base = create_spheric_base(m_coords)

        e_r_t = MathTex(r"\overrightarrow{e_r}")\
            .move_to(spherical_to_cartesian(cartesian_to_spherical(spheric_base[0].get_tip().get_center()) + 0.7 * RIGHT))
        e_theta_t = MathTex(r"\overrightarrow{e_\theta}")\
            .move_to(spherical_to_cartesian(cartesian_to_spherical(spheric_base[1].get_tip().get_center()) + 0.1 * OUT))
        e_phi_t = MathTex(r"\overrightarrow{e_\phi}")\
            .move_to(spherical_to_cartesian(cartesian_to_spherical(spheric_base[2].get_tip().get_center()) + 0.1 * UP))

        spheric_base_t = VGroup(e_r_t, e_theta_t, e_phi_t)
        self.add_fixed_orientation_mobjects(*spheric_base_t)

        self.play(
            *map(Create, spheric_base),
            *map(Write, spheric_base_t)
        )
        self.wait()

        pi1_t = MathTex(r"\pi_1 = \left\{M; \overrightarrow{e_r}; \overrightarrow{e_\theta}\right\}")\
            .to_corner(UR).scale(1.2)
        
        self.play(
            *map(lambda m: m.animate.set_color('#FCBA03'), [spheric_base[0], spheric_base[1], e_r_t, e_theta_t, m, m_t])
        )

        self.wait()
        self.add_fixed_in_frame_mobjects(pi1_t)
        self.play(
            Write(pi1_t, run_time=0.75)
        )
        self.wait()
        self.play(
            *map(lambda m: m.animate.set_color(WHITE), [spheric_base[0], spheric_base[1], e_r_t, e_theta_t, m, m_t])
        )
        self.wait()

        v1 = spheric_base[0].get_end() - m_coords
        v2 = spheric_base[1].get_end() - m_coords
        v3 = spheric_base[2].get_end() - m_coords

        pi1 = Surface(
            lambda u, v: parametric_plane_func(u, v, v1, v2, m_coords),
            u_range=[-6, 1],
            v_range=[-3, 3],
            resolution=(16, 16),
            checkerboard_colors=[RED_D, RED_E],
            fill_opacity=0.5
        )

        self.play(Create(pi1))
        self.wait(2)

        self.play(
            sphere.animate.fade(0.85),
            pi1.animate.fade(0.45)
        )
        self.wait()

        p_coords = np.array([1.7, 0, 0.7])
        p = Dot3D(
            point=p_coords,
            radius=0.07,
            color=WHITE,
            fill_opacity=0.5,
            stroke_width=0
        )
        p_t = Tex('P').next_to(p, RIGHT)

        line1 = DashedLine(p, ORIGIN)
        line1_t = MathTex('a').next_to(line1, OUT, buff=SMALL_BUFF)

        p2_coords = get_symmetric_point(p_coords, np.cross(v1, v2))
        p2 = Dot3D(
            point=p2_coords,
            radius=0.07,
            color=WHITE,
            fill_opacity=0.5,
            stroke_width=0
        )
        p2_t = Tex("P'").next_to(p2, LEFT)

        line2 = DashedLine(ORIGIN, p2_coords)
        line2_t = MathTex('a').next_to(line2, OUT, buff=SMALL_BUFF)
        self.add_fixed_orientation_mobjects(p_t)

        self.play(
            FadeIn(p),
            Write(p_t)
        )
        self.wait()
        self.play(
            Create(line1)
        )
        self.begin_ambient_camera_rotation(rate=0.08)
        self.add_fixed_orientation_mobjects(line1_t)
        self.play(
            Write(line1_t)
        )
        self.wait(2.5)
        self.add_fixed_orientation_mobjects(p2_t)
        self.play(
            FadeIn(p2),
            Write(p2_t)
        )
        self.add_fixed_orientation_mobjects(line2_t)
        self.play(
            Create(line2),
            Write(line2_t)
        )
        self.wait(2)

        dens_t = MathTex(r"\rho(P)=\rho(P')=\rho_0\left(\frac{a}{R}\right)^{2}").next_to(pi1_t, DOWN).scale(0.8)
        self.add_fixed_in_frame_mobjects(dens_t)
        self.play(
            Write(dens_t)
        )
        self.wait()
        self.stop_ambient_camera_rotation()
        self.wait()

        champ_t1 = MathTex(r"\overrightarrow{E(M)}=E_r(M)\overrightarrow{e_r}+E_\theta(M)\overrightarrow{e_\theta}")\
            .next_to(pi1_t, DOWN).shift(0.7 * LEFT).scale(0.9)
        
        e_r_theta = Arrow(
            start=m_coords,
            end=m_coords + spheric_base[0].get_vector() + spheric_base[1].get_vector(),
            buff=0,
            color="#FCBA03"
        )

        proj_e_r = DashedLine(e_r_theta.get_end(), spheric_base[0].get_end(), color="#FCBA03")
        proj_e_theta = DashedLine(e_r_theta.get_end(), spheric_base[1].get_end(), color="#FCBA03")

        self.add_fixed_in_frame_mobjects(champ_t1)
        self.play(
            FadeOut(dens_t, run_time=0.5),
            FadeIn(champ_t1, run_time=0.5),
            Create(e_r_theta),
            Create(proj_e_r),
            Create(proj_e_theta)
        )
        self.wait(3)

        self.play(
            *map(FadeOut, [champ_t1, p_t, p, p2_t, p2, line1, line1_t, line2, line2_t, pi1, e_r_theta, proj_e_r, proj_e_theta]),
            sphere.animate.set_opacity(0.4)
        )
        self.move_camera(
            theta=50 * DEGREES
        )

        pi2 = Surface(
            lambda u, v: parametric_plane_func(u, v, v1, v3, m_coords),
            u_range=[-6, 1],
            v_range=[-3, 3],
            resolution=(16, 16),
            checkerboard_colors=[RED_D, RED_E],
            fill_opacity=0.5
        )

        pi2_t = MathTex(r"\pi_2 = \left\{M; \overrightarrow{e_r}; \overrightarrow{e_\phi}\right\}")\
            .next_to(pi1_t, DOWN, buff=LARGE_BUFF + 1.5).scale(1.1)
        
        self.play(
            *map(lambda m: m.animate.set_color('#FCBA03'), [spheric_base[0], spheric_base[2], e_r_t, e_phi_t, m, m_t])
        )
        self.add_fixed_in_frame_mobjects(pi2_t)
        self.play(
            Write(pi2_t, run_time=0.75)
        )
        self.wait()
        self.play(
            *map(lambda m: m.animate.set_color(WHITE), [spheric_base[0], spheric_base[2], e_r_t, e_phi_t, m, m_t])
        )
        self.wait()
        self.camera.light_source_start_point = 9 * UP + 7 * LEFT + 10 * IN
        self.play(
            Create(pi2)
        )
        self.wait()
        self.move_camera(
            theta=-30 * DEGREES,
            run_time=6
        )
        self.wait(2)
        self.play(
            *map(FadeOut, [axes, sphere, ct_p, ct_p_t, m, m_t, *spheric_base, *spheric_base_t, pi2])
        )
        self.wait()
        self.play(
            pi1_t.animate.move_to(3 * LEFT + 1.5 * UP),
            pi2_t.animate.move_to(3 * RIGHT + 1.5 * UP)
        )
        self.play(
            pi1_t.animate.scale(1.2),
            pi2_t.animate.scale(1.2)
        )
        self.wait()

        arrow1 = Vector(1.25 * DOWN).next_to(pi1_t, DOWN)
        arrow2 = Vector(1.25 * DOWN).next_to(pi2_t, DOWN)
        
        champ_t1 = MathTex(r"\overrightarrow{E(M)}=E_r(M)\overrightarrow{e_r} + E_\theta(M)\overrightarrow{e_\theta}")\
            .next_to(arrow1, DOWN).shift(0.6 * LEFT)
        champ_t2 = MathTex(r"\overrightarrow{E(M)}=E_r(M)\overrightarrow{e_r} + E_\phi(M)\overrightarrow{e_\phi}")\
            .next_to(arrow2, DOWN).shift(0.6 * RIGHT)
        champ_f_t = MathTex(r"\overrightarrow{E(M)}=E(M)\overrightarrow{e_r}").scale(1.2)
        self.add_fixed_in_frame_mobjects(arrow1, arrow2, champ_t1, champ_t2)

        self.play(
            *map(Write, [arrow1, arrow2, champ_t1, champ_t2])
        )
        self.wait()
        champ_t1.scale_in_place
        self.play(
            *map(
                lambda m: m.animate.set_color("#FCBA03"),
                [pi1_t[0][6:10], pi2_t[0][6:10], champ_t1[0][14:18], champ_t2[0][14:18]]
            )
        )
        self.play(
            *map(
                lambda m: m.animate(rate_func=there_and_back).scale(1.5),
                [pi1_t[0][6:10], pi2_t[0][6:10], champ_t1[0][14:18], champ_t2[0][14:18]]
            )
        )
        self.play(
            *map(
                lambda m: m.animate.set_color(WHITE),
                [pi1_t[0][6:10], pi2_t[0][6:10], champ_t1[0][14:18], champ_t2[0][14:18]]
            )
        )
        self.wait()
        self.add_fixed_in_frame_mobjects(champ_f_t)
        self.play(
            *map(FadeOut, [arrow1, arrow2, champ_t1, champ_t2, pi1_t, pi2_t]),
            FadeIn(champ_f_t)
        )
        self.wait(2)

class ExplainFlux(ThreeDScene):
    def construct(self):
        esint_template = TexTemplate()
        esint_template.add_to_preamble(r'\usepackage{esint}')

        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=-30 * DEGREES
        )

        axes = ThreeDAxes(z_range=[-3, 3, 1]).set_stroke(width=0.5)
        x_label = axes.get_x_axis_label(label=r"\overrightarrow{e_x}")
        y_label = axes.get_y_axis_label(label=r"\overrightarrow{e_y}").shift(2 * UP)
        z_label = axes.get_z_axis_label(label=r"\overrightarrow{e_z}", rotation=0, buff=LARGE_BUFF)
        cartesian_base_t = VGroup(x_label, y_label, z_label)
        self.add_fixed_orientation_mobjects(x_label, y_label, z_label)

        self.play(
            Create(axes, run_time=0.8),
            *map(Write, cartesian_base_t)
        )

        surface = Surface(
            lambda u, v: parametric_plane_func(u, v, RIGHT, OUT, ORIGIN),
            u_range=[-2, 2],
            v_range=[-1.5, 1.5],
            resolution=(1, 1),
            checkerboard_colors=[GREY],
            fill_opacity=0.1,
            stroke_width=4
        )
        ds_surface = Surface(
            lambda u, v: parametric_plane_func(u, v, RIGHT, OUT, ORIGIN),
            u_range=[-0.25, 0.25],
            v_range=[-0.25, 0.25],
            resolution=(1, 1),
            fill_opacity=0,
            stroke_width=3
        ).shift(0.2 * OUT)

        ds = Arrow(0.2 * OUT, 1.8 * UP + 0.2 * OUT, buff=0)
        ds_t1 = MathTex(r"\overrightarrow{\mathrm{d^2S}}=\mathrm{d^2S}\overrightarrow{e_y}").to_corner(UR)

        self.play(
            Write(surface),
            run_time=0.8
        )
        self.move_camera(
            theta=30 * DEGREES,
            run_time=4
        )

        ds_t2 = MathTex(r"\overrightarrow{\mathrm{d^2S}}").next_to(ds.get_end(), buff=MED_LARGE_BUFF)
        self.add_fixed_in_frame_mobjects(ds_t1)
        self.add_fixed_orientation_mobjects(ds_t2)
        self.play(
            Write(ds_surface),
            Write(ds_t1),
            Write(ds_t2),
            Create(ds)
        )
        self.wait()

        ex_field_c = [np.array([x, -0.25, z]) for x in (1.2, 0, -1.2) for z in (0.85, 0, -0.85)]
        ex_field = [Arrow(c, c + 2 * UP, color="#FCBA03") for c in ex_field_c]
        flux_t = MathTex(
            r"\Phi(\overrightarrow{E}) =",
            r"\iint_S \overrightarrow{E}\overrightarrow{\mathrm{d^2S}}",
            tex_template=esint_template,
            # tex_to_color_map={
            #     "\overrightarrow{E}": "#FCBA03"
            # }
        ).next_to(ds_t1, DOWN, aligned_edge=RIGHT)
        flux_t[0][3:5].set_color("#FCBA03")
        flux_t[1][3:5].set_color("#FCBA03")

        self.play(
            *map(Create, ex_field)
        )
        self.wait()
        self.add_fixed_in_frame_mobjects(flux_t)
        self.play(
            Write(flux_t, run_time=0.8)
        )
        self.wait()
        ds_grp = VGroup(ds_t2, ds_surface)

        self.play(
            *map(FadeOut, ex_field),
            ds_grp.animate.shift(OUT + 1.6 * RIGHT),
            ds.animate.shift(OUT + 1.6 * RIGHT).set_opacity(0.5),
            FadeOut(flux_t)
        )
        self.wait()

        flux_el_t1 = MathTex(
            r"{{\delta\Phi =}} {{\overrightarrow{E} \overrightarrow{\mathrm{d^2S}}}}",
            tex_to_color_map={
                r"\overrightarrow{E}": "#FCBA03"
            }
        ).next_to(ds_t1, DOWN, aligned_edge=RIGHT)
        flux_el_t2 = MathTex(
            r"=\Vert\overrightarrow{E}\Vert\,\Vert\overrightarrow{\mathrm{d^2S}}\Vert\cos(\theta)\,\mathrm{u.a.}",
            tex_to_color_map={
                r"\Vert\overrightarrow{E}\Vert\,\Vert\overrightarrow{\mathrm{d^2S}}\Vert": "#FCBA03",
                r"\cos(\theta)": GREEN_C
            }
        ).next_to(flux_el_t1, DOWN, aligned_edge=RIGHT)
        flux_el_t3 = MathTex(
            r"=\qquad\times\qquad\mathrm{u.a.}"
        ).next_to(flux_el_t2, DOWN, aligned_edge=RIGHT)

        ex_arrow = Arrow(ds.get_start(), ds.get_end() + 0.5 * UP, color="#FCBA03", buff=0)
        self.add_fixed_in_frame_mobjects(flux_el_t1, flux_el_t2)
        self.play(
            FadeIn(ex_arrow),
            LaggedStartMap(FadeIn, [flux_el_t1, flux_el_t2], lag_ratio=1)
        )
        self.wait()

        e_vt = ValueTracker(1)
        e_value = DecimalNumber(
            1,
            num_decimal_places=1,
            include_sign=True,
            color="#FCBA03"
        ).next_to(flux_el_t3[0][0], RIGHT, buff=SMALL_BUFF)
        e_value.add_updater(lambda d: d.set_value(e_vt.get_value()))

        cos_vt = ValueTracker()
        cos_value = DecimalNumber(
            1,
            num_decimal_places=1,
            include_sign=True,
            color=GREEN_C
        ).next_to(flux_el_t3[0][1], RIGHT, buff=SMALL_BUFF)
        cos_value.add_updater(lambda d: d.set_value(np.cos(cos_vt.get_value())))

        self.add_fixed_in_frame_mobjects(flux_el_t3, e_value, cos_value)
        self.play(
            *map(FadeIn, [flux_el_t3, e_value, cos_value])
        )
        e_value.add_updater(lambda d: self.add_fixed_in_frame_mobjects(d))
        cos_value.add_updater(lambda d: self.add_fixed_in_frame_mobjects(d))

        self.wait()
        self.play(
            e_vt.animate.set_value(2.5),
            ex_arrow.animate.put_start_and_end_on(ex_arrow.get_start(), ex_arrow.get_end() + UP)
        )
        self.wait()
        self.play(
            e_vt.animate.set_value(1),
            ex_arrow.animate.put_start_and_end_on(ex_arrow.get_start(), ex_arrow.get_end() - UP)
        )
        cos_t = Integer(color=GREEN_C).move_to(0.8 * UP + 0.5 * LEFT).scale(1.2)
        cos_t.add_updater(lambda d: d.set_value(cos_vt.get_value() / DEGREES))
        deg_t = MathTex(
            r"\circ",
            color=GREEN_C
        ).next_to(cos_t, RIGHT).scale(1.2).shift(0.1 * UP)
        deg_t.add_updater(lambda d: d.next_to(cos_t, RIGHT).shift(0.1 * UP))

        self.add_fixed_in_frame_mobjects(cos_t, deg_t)
        self.play(
            FadeIn(cos_t),
            FadeIn(deg_t),
            run_time=0.5
        )
        cos_t.add_updater(lambda d: self.add_fixed_in_frame_mobjects(d))

        self.play(
            cos_vt.animate.set_value(45 * DEGREES),
            Rotate(ex_arrow, 45 * DEGREES, axis=RIGHT, about_point=ex_arrow.get_start()),
            run_time=1.5
        )
        self.wait()
        self.play(
            cos_vt.animate.set_value(90 * DEGREES),
            Rotate(ex_arrow, 45 * DEGREES, axis=RIGHT, about_point=ex_arrow.get_start()),
            run_time=1.5
        )
        self.wait()
        self.play(
            cos_vt.animate.set_value(180 * DEGREES),
            Rotate(ex_arrow, 90 * DEGREES, axis=RIGHT, about_point=ex_arrow.get_start()),
            run_time=1.5
        )
        self.wait(2)
        cos_value.clear_updaters()
        e_value.clear_updaters()
        deg_t.clear_updaters()
        cos_t.clear_updaters()
        self.play(
            *map(FadeOut, [cos_t, deg_t, cos_value, e_value, flux_el_t2, flux_el_t3]),
            Rotate(ex_arrow, -180 * DEGREES, axis=RIGHT, about_point=ex_arrow.get_start())
        )
        self.wait()

        ds_grp.add(ds, ex_arrow)

        for _ in range(10):
            self.play(ds_grp.animate.shift(0.15 * LEFT), run_time=0.3)
        for _ in range(12):
            self.play(ds_grp.animate.shift(0.15 * LEFT), run_time=0.065)
        self.play(ds_grp.animate.shift(0.4 * IN), run_time=0.2)
        for _ in range(13):
            self.play(ds_grp.animate.shift(0.15 * RIGHT), run_time=0.065)
        self.play(ds_grp.animate.set_opacity(0), run_time=0.2)

        self.wait()

        self.play(
            *map(Create, ex_field),
            FadeOut(flux_el_t1),
            FadeIn(flux_t)
        )

        flux_vt = ValueTracker(10)

        flux_value = DecimalNumber(
            10,
            num_decimal_places=1,
            include_sign=True
        )
        flux_value.next_to(flux_t[0], RIGHT).shift(0.08*DOWN)
        flux_value.add_updater(lambda d: d.set_value(flux_vt.get_value()))
        ua_t = MathTex(r"\mathrm{u.a.}").next_to(flux_value)
        self.add_fixed_in_frame_mobjects(flux_value, ua_t)
        self.play(
            FadeIn(flux_value),
            FadeOut(flux_t[1]),
            FadeIn(ua_t)
        )
        flux_value.add_updater(lambda d: self.add_fixed_in_frame_mobjects(d))

        self.wait()

        self.play(
            flux_vt.animate.set_value(15),
            *[arr.animate.put_start_and_end_on(arr.get_start(), arr.get_end() + UP) for arr in ex_field],
            run_time=2
        )
        self.wait()
        self.play(
            flux_vt.animate.set_value(10),
            *[arr.animate.put_start_and_end_on(arr.get_start(), arr.get_end() - UP) for arr in ex_field],
            run_time=2
        )
        self.wait()
        self.play(
            flux_vt.animate.set_value(8.5),
            *[Rotate(arr, 20 * DEGREES, axis=RIGHT, about_point=arr.get_start()) for arr in ex_field],
            run_time=2
        )
        self.wait()
        self.play(
            flux_vt.animate.set_value(0),
            *[Rotate(arr, 70 * DEGREES, axis=RIGHT, about_point=arr.get_start()) for arr in ex_field],
            run_time=2
        )
        self.wait()
        self.play(
            flux_vt.animate.set_value(-10),
            *[Rotate(arr, 90 * DEGREES, axis=RIGHT, about_point=arr.get_start()) for arr in ex_field],
            run_time=2
        )
        self.wait(2)

class Induction41(ThreeDScene):
    def construct(self):
        esint_template = TexTemplate()
        esint_template.add_to_preamble(r"\usepackage{physics}")
        esint_template.add_to_preamble(r"\usepackage{esint}")

        self.set_camera_orientation(
			phi =70 * DEGREES,
			theta=30 * DEGREES
		)
        field_t = MathTex(
            r"\overrightarrow{B}=B_0\cos\left(\frac{\pi r}{2R}\right)\cos(\omega t)\overrightarrow{e_z}",
            tex_to_color_map={
                r"\overrightarrow{B}": "#FCBA03"
            }
        )
        self.add_fixed_in_frame_mobjects(field_t)
        self.play(
            FadeIn(field_t),
            run_time=1.2
        )
        self.wait()
        self.play(
            field_t.animate.to_corner(UR).scale(0.8)
        )
        self.wait()

        axes = ThreeDAxes()
        cylinder = Cylinder(
            radius=2,
            height=8,
            direction=Y_AXIS,
            show_ends=False,
            fill_opacity=0.7
        )
        y_label = axes.get_y_axis_label(label=r"\overrightarrow{e_z}").shift(2 * UP)
        self.add_fixed_orientation_mobjects(y_label)
        self.play(
            FadeIn(axes),
            FadeIn(cylinder),
            Write(y_label)
        )
        self.wait()

        b_arrows = [
            Arrow(
                start=np.roll(cylindric_to_cartesian((r, theta, z)), 2),
                end=np.roll(cylindric_to_cartesian((r, theta, z + 2 * np.cos((PI * r) / (2 * cylinder.radius)))), 2),
                buff=0,
                color="#FCBA03",
                stroke_opacity=0.5
            )
            for r in np.linspace(0.5, 1.5, 2)
            for theta in np.linspace(0, TAU, 6)
            for z in np.linspace(-3, 3, 3)
        ]
        for arrow in b_arrows:
            arrow.original_length = arrow.get_length()
        self.play(
            *map(Create, b_arrows)
        )
        self.wait()

        t_vt = ValueTracker()
        t_value = DecimalNumber(num_decimal_places=1, unit="s").next_to(field_t, DOWN, aligned_edge=RIGHT)
        t_t = MathTex("t=").next_to(t_value, LEFT, aligned_edge=DOWN)
        t_value.add_updater(lambda t: t.set_value(t_vt.get_value()))
        self.add_fixed_in_frame_mobjects(t_value, t_t)
        self.play(
            *map(lambda a: a.animate.set_opacity(1), b_arrows),
            cylinder.animate.set_opacity(0.35),
            FadeIn(t_value),
            FadeIn(t_t)
        )
        t_value.add_updater(lambda t: self.add_fixed_in_frame_mobjects(t))

        self.wait()


        for arrow in b_arrows:
            arrow.add_updater(
                lambda a: a.become(
                    Arrow(
                        a.get_start(),
                        a.get_start() + np.cos(2 * t_vt.get_value()) * a.original_length * UP,
                        buff=0,
                        color="#FCBA03"
                    )
                )
            )
        
        # self.play(
        #     t_vt.animate.set_value(10),
        #     run_time=4,
        #     rate_func=linear
        # )
        self.wait()
        self.move_camera(
			theta=200 * DEGREES,
			phi=60 * DEGREES,
			run_time=4 * PI,
			rate_func=lambda t: smooth(t, inflection = 5.0),
            added_anims=[t_vt.animate(rate_func=linear, run_time=4 * PI).set_value(4 * PI)]
		)
        self.wait(2)

        circle = Circle(
            radius=1.3,
            stroke_width=4,
            fill_opacity=0.45,
            fill_color="#FCBA03",
            stroke_color=GREEN
        ).rotate(90 * DEGREES, axis=X_AXIS).shift(2.8 * DOWN)

        gamma_t = MathTex(r"\Gamma", color=GREEN).move_to(3.5 * LEFT + 3 * DOWN)

        ds = Arrow(2.8 * DOWN + 0.2 * OUT, 1.3 * DOWN + 0.2 * OUT, buff=0)
        ds_t = MathTex(r"\overrightarrow{\mathrm{d^2S}}").next_to(ds.get_end() + 0.6 * OUT)

        eq_t = MathTex(
            r"e=\oint_{\Gamma}\overrightarrow{E}.\overrightarrow{\mathrm{dl}}=", r"-\dv{}{t}", r"\iint_S \overrightarrow{B}.\overrightarrow{\mathrm{d^2S}}",
            tex_template=esint_template,
            # tex_to_color_map={
            #     r"\Gamma": GREEN
            # }
        ).scale(0.8).to_corner(UL)
        eq_t[0][3].set_color(GREEN)
        self.add_fixed_in_frame_mobjects(eq_t)
        self.add_fixed_orientation_mobjects(gamma_t, ds_t)
        self.move_camera(
            zoom=0.9,
            theta=220 * DEGREES,
            frame_center=UP,
            added_anims=[FadeIn(eq_t), Create(circle), Write(gamma_t), Create(ds), Write(ds_t)]
        )
        self.wait()

        self.play(
            eq_t[2].animate(rate_func=there_and_back, run_time=1.8).set_color("#FCBA03").scale(1.3),
            circle.animate(rate_func=there_and_back, run_time=1.8).scale(1.4)
        )
        self.wait()
        eq_t_flux = MathTex(
            r"-\dv{\Phi_S(t)}{t}"
        ).scale(0.8).next_to(eq_t[0], RIGHT)
        self.add_fixed_in_frame_mobjects(eq_t_flux)
        self.play(
            FadeOut(eq_t[1:3]),
            FadeIn(eq_t_flux),
            # Transform(eq_t[1:3], eq_t_flux)
        )
        self.wait()

        cos_axes = Axes(
            x_range=[0, 4 * PI, PI],
            y_range=[-2, 2, 2],
            tips=False,
            x_length=2.5,
            y_length=1.6,
            axis_config={"stroke_width": 4}
        ).next_to(eq_t, DOWN, aligned_edge=LEFT, buff=LARGE_BUFF)

        cos_axis_labels = cos_axes.get_axis_labels(
            x_label=MathTex("t").scale(0.8),
            y_label=MathTex(r"\Phi(t)", color="#FCBA03").scale(0.8)
        )
        self.add_fixed_in_frame_mobjects(cos_axes, cos_axis_labels)
        self.play(
            Create(cos_axes),
            Write(cos_axis_labels)
        )
        self.wait()
        cos_func = cos_axes.plot(
            lambda t: np.cos(2 * t),
            x_range=[0, 4 * PI],
            color="#FCBA03"
        )
        self.add_fixed_in_frame_mobjects(cos_func)
        self.begin_3dillusion_camera_rotation(rate=-0.05)
        self.play(
            t_vt.animate.set_value(8 * PI),
            Create(cos_func),
            rate_func=linear,
            run_time=4 * PI
        )
        self.wait()

        cos_func_l = cos_axes.plot(
            lambda t: 2 / 1.3 * np.cos(2 * t),
            x_range=[0, 4 * PI],
            color="#FCBA03"
        )
        # self.add_fixed_in_frame_mobjects(cos_func_l)
        self.play(
            circle.animate.scale(2 / 1.3),
            Transform(cos_func, cos_func_l),
            run_time=2.5
        )
        self.wait()
        self.play(
            circle.animate.scale(1.5),
            run_time=2
        )
        self.wait()
        self.play(
            circle.animate.scale(1 / 1.5),
            run_time=1.5
        )
        d_arrow = Vector(DOWN).next_to(cos_axes, DOWN)
        dt_t = MathTex(r"-\dv{}{t}").scale(0.8).next_to(d_arrow, LEFT)

        sin_axes = Axes(
            x_range=[0, 4 * PI, PI],
            y_range=[-2, 2, 2],
            tips=False,
            x_length=2.5,
            y_length=1.6,
            axis_config={"stroke_width": 4}
        ).next_to(cos_axes, DOWN, aligned_edge=LEFT, buff=1.6)

        sin_axis_labels = sin_axes.get_axis_labels(
            x_label=MathTex("t").scale(0.8),
            y_label=MathTex("e", color=RED).scale(0.8)
        )
        self.add_fixed_in_frame_mobjects(d_arrow, dt_t, sin_axes, sin_axis_labels)
        self.play(
            Create(d_arrow),
            FadeIn(dt_t),
            Create(sin_axes),
            Write(sin_axis_labels),
            FadeOut(cos_func)
        )
        self.wait()

        sin_func_l = sin_axes.plot(
            lambda t: 2 / 1.3 * np.sin(2 * t),
            x_range=[0, 4 * PI],
            color=RED
        )

        e_length = 1.2
        e_arrows = [
            Arrow(
                np.roll(cylindric_to_cartesian((2, theta, circle.get_center()[1])), 2),
                np.roll(cylindric_to_cartesian((math.sqrt(4 + e_length**2), theta + np.arctan2(e_length, 4), circle.get_center()[1])), 2),
                buff=0,
                color=RED
            )
            for theta in np.linspace(0, TAU, 10)
        ]

        for arrow in e_arrows:
            arrow.unit_vector = arrow.get_unit_vector()
            arrow.add_updater(
                lambda a: a.become(
                    Arrow(
                        a.get_start(),
                        a.get_start() + np.sin(2 * t_vt.get_value()) * e_length * a.unit_vector,
                        buff=0,
                        color=RED
                    )
                )
            )
        
        self.play(
            *map(Create, e_arrows)
        )

        self.add_fixed_in_frame_mobjects(sin_func_l)
        self.play(
            t_vt.animate.set_value(12 * PI),
            Create(cos_func),
            Create(sin_func_l),
            rate_func=linear,
            run_time=4 * PI
        )
        self.wait()

class ExplainCurl(ThreeDScene):
    def construct(self):
        self.set_camera_orientation(
            phi=70 * DEGREES,
            theta=30 * DEGREES
        )

        axes = ThreeDAxes(z_range=[-3, 3, 1]).set_stroke(width=0.5)
        x_label = axes.get_x_axis_label(label=r"\overrightarrow{e_x}")
        y_label = axes.get_y_axis_label(label=r"\overrightarrow{e_y}").shift(2 * UP)
        z_label = axes.get_z_axis_label(label=r"\overrightarrow{e_z}", rotation=0, buff=LARGE_BUFF)
        cartesian_base_t = VGroup(x_label, y_label, z_label)

        fil = Line(10 * IN, 10 * OUT, stroke_width=6)
        i_t = MathTex("I=").to_corner(UL)
        amp_vt = ValueTracker(2)
        amp_value = DecimalNumber(
            number=2,
            num_decimal_places=1,
            include_sign=True,
            # include_background_rectangle=True,
            unit="A"
        ).next_to(i_t, RIGHT)
        amp_value.add_updater(lambda d: d.set_value(amp_vt.get_value()))

        self.add_fixed_orientation_mobjects(x_label, y_label, z_label)
        self.play(
            Create(axes),
            *map(Write, cartesian_base_t),
            run_time=0.8
        )
        self.wait()
        self.add_fixed_in_frame_mobjects(i_t, amp_value)

        self.play(
            Create(fil, run_time=3),
            FadeIn(amp_value),
            FadeIn(i_t)
        )
        amp_value.add_updater(lambda d: self.add_fixed_in_frame_mobjects(d))

        def shift_arrow(arrow, dt):
            if arrow.get_start()[2] >= 4:
                arrow.move_to(-7 * Z_AXIS)
            else:
                arrow.move_to(arrow.get_center() + dt * Z_AXIS)
        
        i_arrows = [
            Arrow(
                start=(0, 0, z),
                end=(0, 0, z + 2),
                buff=0,
                color=RED
            )
            for z in [-5, -2, 1, 4]
        ]
        self.play(*map(Create, i_arrows))
        for arrow in i_arrows:
            arrow.add_updater(shift_arrow)
        self.wait(3)

        def magnetic_field_func(p, i):
            r = np.linalg.norm(p)
            scaled_p = p / r
            if r < 0.4:
                return ORIGIN
            else:
                return np.array([i * (1/r**2) * scaled_p[1], -i * (1/r**2) * scaled_p[0], 0])

        magnetic_field = ArrowVectorField(
            func=lambda p: magnetic_field_func(p, amp_vt.get_value()),
            y_range=[
                math.floor(-config["frame_width"] / 2),
                math.ceil(config["frame_width"] / 2)
            ],
            # length_func=lambda norm: 0.5 * sigmoid(norm)
        )
        self.play(
            Create(magnetic_field)
        )
        self.wait(2)

        self.move_camera(
            phi=0,
            theta=0,
            run_time=1.7
        )
        for arrow in i_arrows:
            arrow.clear_updaters()
        self.wait()

        circle = Circle(radius=0.28, color=WHITE)
        circle.add(Dot(radius=0.02))
        curl_norm_t = MathTex(
            r"\Vert\overrightarrow{rot}\,\overrightarrow{B}\Vert="
        )
        curl_norm_t.move_to(UP).add_background_rectangle()
        curl_value = DecimalNumber(
            1,
            num_decimal_places=1
        )
        curl_value.next_to(curl_norm_t, RIGHT).shift(0.08*DOWN)

        self.add_fixed_in_frame_mobjects(curl_norm_t, curl_value)

        self.play(
            Create(circle),
            FadeIn(curl_norm_t),
            FadeIn(curl_value)
        )
        curl_value.add_updater(lambda d: d.set_value(amp_vt.get_value() / 2))
        curl_value.add_updater(lambda d: self.add_fixed_in_frame_mobjects(d))
        self.wait()

        magnetic_field_l = ArrowVectorField(
            func=lambda p: magnetic_field_func(p, 4 * amp_vt.get_value()),
            y_range=[
                math.floor(-config["frame_width"] / 2),
                math.ceil(config["frame_width"] / 2)
            ],
            # length_func=lambda norm: 0.5 * sigmoid(norm)
        )

        self.play(
            amp_vt.animate.set_value(4),
            magnetic_field.animate.become(magnetic_field_l),
            run_time=3
        )
        self.wait()

        curl_value.clear_updaters()
        self.play(
            FadeOut(curl_norm_t),
            FadeOut(curl_value)
        )

        for arrow in i_arrows:
            arrow.add_updater(shift_arrow)
        self.move_camera(
            phi=70 * DEGREES,
            theta=30 * DEGREES,
            zoom=1.2
        )
        self.play(
            *map(FadeOut, i_arrows),
            FadeOut(fil)
        )
        self.wait()

        magnetic_field_i = ArrowVectorField(
            func=lambda p: magnetic_field_func(p, -4 * amp_vt.get_value()),
            y_range=[
                math.floor(-config["frame_width"] / 2),
                math.ceil(config["frame_width"] / 2)
            ],
            # length_func=lambda norm: 0.5 * sigmoid(norm)
        )

        curl_vector = Arrow(ORIGIN, 2 * IN, buff=0)
        curl_vector_t = MathTex(
            r"\overrightarrow{rot}\,\overrightarrow{B}"
        ).next_to(curl_vector.get_end(), OUT, buff=-0.6)
        self.add_fixed_orientation_mobjects(curl_vector_t)
        self.play(
            Create(curl_vector),
            Write(curl_vector_t)
        )
        self.wait()
        self.play(
            amp_vt.animate.set_value(-4),
            curl_vector.animate.become(Arrow(ORIGIN, 2 * OUT, buff=0)),
            magnetic_field.animate.become(magnetic_field_i),
            curl_vector_t.animate.shift(5.2 * OUT),
            run_time=3
        )
        self.wait()

        self.play(
            FadeOut(curl_vector),
            FadeOut(curl_vector_t)
        )

        div_t = MathTex(
            r"div\,\overrightarrow{B}=0"
        ).move_to(UP).add_background_rectangle()

        self.move_camera(
            phi=0,
            theta=0,
            run_time=1.7
        )

        self.add_fixed_in_frame_mobjects(div_t)
        self.play(
            FadeIn(div_t)
        )
        self.wait()
        self.play(
            circle.animate.move_to(3 * UP),
            div_t.animate.shift(2 * RIGHT)
        )

        arrows_div = [
            Vector(
                0.7 * LEFT,
                color=RED,
                stroke_width=6
            ).move_to(circle.get_center() + 0.8 * RIGHT + offset * UP)
            for offset in [-0.3, 0, 0.3]
        ]

        self.play(
            *map(Create, arrows_div),
            run_time=0.5
        )
        self.play(
            *[arrow.animate.shift(0.5 * LEFT) for arrow in arrows_div],
            run_time=0.5
        )
        self.play(
            *map(FadeOut, arrows_div),
            run_time=0.5

        )
        for arrow in arrows_div:
            arrow.shift(0.5 * LEFT)
        self.play(
            *map(Create, arrows_div),
            run_time=0.5
        )
        self.play(
            *[arrow.animate.shift(0.5 * LEFT) for arrow in arrows_div],
            run_time=0.5
        )
        self.play(
            *map(FadeOut, arrows_div),
            run_time=0.5
        )
        self.wait()

        self.move_camera(
            phi=70 * DEGREES,
            theta=30 * DEGREES,
            frame_center=2 * UP,
            added_anims=[
                FadeOut(circle),
                FadeOut(div_t)
            ]
        )
        self.wait()

        prism = Prism(
            dimensions=[2.5, 5, 2],
            fill_color=BLUE_E,
            stroke_width=2,
            fill_opacity=0.5
        ).move_to(1.75 * UP)

        prism_t = MathTex(
            r"\mathscr{P}=\{S_1, \dots, S_6\}"
        ).to_corner(UR)

        self.begin_ambient_camera_rotation()

        self.add_fixed_in_frame_mobjects(prism_t)
        self.play(
            Write(prism_t),
            FadeOut(i_t)
        )
        self.wait()
        surfaces_t_loc = [[-1, 2, 0.7], [0, -1, 0], [1, 2, 0.7], [0, 5, 1], [0, 3, 2], [0, 3, -2]]
        surfaces_t = [MathTex(f"S_{i}").move_to(p) for i, p in enumerate(surfaces_t_loc, start=1)]

        for i, t in zip([2, 5, 3, 4, 1, 0], surfaces_t):
            self.add_fixed_orientation_mobjects(t)
            self.play(
                Create(prism[i]),
                Write(t),
                run_time=0.7
            )
            if i == 4:
                self.wait(1)
        self.stop_ambient_camera_rotation()
        
        flux_t = MathTex(
            r"\Phi_{\mathscr{P}}(\overrightarrow{B})=\Phi_{S_1}+\Phi_{S_2}+\Phi_{S_3}+\Phi_{S_4}+\Phi_{S_5}+\Phi_{S_6}"
        ).to_corner(UL)
        self.add_fixed_in_frame_mobjects(flux_t)
        self.move_camera(
            theta=90 * DEGREES,
            run_time=2,
            added_anims=[
                FadeIn(flux_t, run_time=0.7)
            ]
        )
        self.wait()

