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
