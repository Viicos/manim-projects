"""Animations for 2021."""

import math
from utils import *
from manim import *
from manim.utils.space_ops import cartesian_to_spherical, spherical_to_cartesian


def CylindricToCartesian(point):  # point = (r, theta, z)
	r, theta, z = point
	return [r * np.cos(theta), r * np.sin(theta), z]

class Test2(ThreeDScene):
	def construct(self):
		self.set_camera_orientation(
			phi=70 * DEGREES,
			theta=50 * DEGREES,
			zoom=1.2
		)
		axes = ThreeDAxes()
		m_coords = np.array([0, 3.5, 1.7])
		theta_angle = ArcBetweenPoints(
			OUT,
			m_coords / np.linalg.norm(m_coords)
		).rotate(90 * DEGREES, axis=OUT - m_coords / np.linalg.norm(m_coords), about_point=m_coords / np.linalg.norm(m_coords))
		theta_angle = Arc(radius=0.8).rotate(PI / 2, axis=DOWN, about_point=UP)
		self.set_camera_orientation(
			phi=90 * DEGREES,
			theta=0 * DEGREES
		)

		self.add(axes, theta_angle, Dot(UP))
		# self.play(
		# 	theta_angle.animate.rotate(90 * DEGREES, axis=DOWN, about_point=m_coords / np.linalg.norm(m_coords))
		# )
		# self.play(
		# 	theta_angle.animate.rotate(PI / 4, axis=OUT, about_point=m_coords / np.linalg.norm(m_coords))
		# )
		# self.wait()


class GaussSquared(ThreeDScene):
	def construct(self):
		esint_template = TexTemplate()
		esint_template.add_to_preamble(r'\usepackage{esint}')

		self.set_camera_orientation(
			phi=70 * DEGREES,
			theta=50 * DEGREES
		)

		charged_sphere = Sphere(radius=2.5, resolution=(18, 18)).set_fill(BLUE_E, opacity=0.4)
		axes = ThreeDAxes().set_stroke(width=0.5)
		center_point = Dot3D(radius=0.06).scale(1.2)
		center_point_t = Tex('O').shift(0.3 * RIGHT + 0.15 * DOWN)
		
		self.add_fixed_orientation_mobjects(center_point_t)
		self.play(Create(charged_sphere), Create(axes), Create(center_point), Write(center_point_t, run_time=0.7))
		self.wait()

		# r0_t = Tex('$R(S) = R_0$').to_corner(UL)
		# self.add_fixed_in_frame_mobjects(r0_t)
		# self.play(Write(r0_t))

		point_m_coords = (0, 3, 2.5)
		point_m = Dot3D(point_m_coords, radius=0.06).scale(1.2)
		point_m_t = Tex('M').move_to((0.4, 3, 2.5)).scale(0.8)
		self.add_fixed_orientation_mobjects(point_m_t)
		self.play(Create(point_m), Create(point_m_t))
	
		self.wait()

		# er_vector = Arrow3D(
		# 	start=point_m_coords,
		# 	end=spherical_to_cartesian(*(cartesian_to_spherical(point_m_coords) + 1.3 * RIGHT)),
		# 	thickness=0.01,
		# 	height=0.2,
		# 	base_radius=0.05
		# )
		print(cartesian_to_spherical(point_m_coords), 1.3 * RIGHT, cartesian_to_spherical(point_m_coords) + 1.3 * RIGHT)
		er_vector = Arrow(
			start=point_m_coords,
			end=spherical_to_cartesian((cartesian_to_spherical(point_m_coords) + 1.3 * RIGHT)),
			buff=0,
			max_tip_length_to_length_ratio=0.5
		)


		e_r_t = Tex('$\\overrightarrow{e_r}$')\
                    .move_to(er_vector.end)\
					.shift(0.4 * IN + 0.3 * LEFT)\
                    .rotate(90 * DEGREES, axis=RIGHT)\
                    .rotate(90 * DEGREES - self.camera.get_theta(), axis=IN)\
                    .flip(IN)\
                    .scale(0.8)

		self.play(Create(er_vector), Write(e_r_t))

		self.wait()

		# e_vector = Arrow3D(
		# 	start=point_m_coords,
		# 	end=spherical_to_cartesian(*(cartesian_to_spherical(point_m_coords) + 1.8 * RIGHT))
		# ).set_color('#fcba03')
		e_vector = Arrow(
			start=point_m_coords,
			end=spherical_to_cartesian((cartesian_to_spherical(point_m_coords) + 1.8 * RIGHT)),
			buff=0,
			max_tip_length_to_length_ratio=0.5
		).set_color('#fcba03')
		e_t = Tex('$\\overrightarrow{E(M)} = E(R)\\overrightarrow{e_r}$')\
                    .move_to(e_vector.end)\
                    .shift(0.4 * IN + 2.2 * LEFT)\
                    .rotate(90 * DEGREES, axis=RIGHT)\
                    .rotate(90 * DEGREES - self.camera.get_theta(), axis=IN)\
                    .flip(IN)\
					.scale(0.7)\
					.set_color('#fcba03')
		
		# self.add_fixed_orientation_mobjects(e_r_t, e_t)
		self.play(Create(e_vector), Write(e_t))

		self.wait()

		self.play(
			ApplyMethod(
				e_vector.put_start_and_end_on,
				e_vector.get_start(),
				spherical_to_cartesian((cartesian_to_spherical(point_m_coords) + 2.2 * RIGHT))
			)
		)

		self.wait()

		self.play(
			ApplyMethod(
				e_vector.put_start_and_end_on,
				e_vector.get_start(),
				spherical_to_cartesian((cartesian_to_spherical(point_m_coords) + 0.8 *  LEFT))
			)
		)

		self.wait()

		self.play(
			*map(FadeOut, [e_vector, e_t, er_vector, e_r_t])
		)

		gauss_sphere = Sphere(radius=math.dist(ORIGIN, point_m_coords), resolution=(18, 18)).set_fill(RED_E, opacity=0.5)

		self.move_camera(
			zoom=0.89,
			frame_center=(2, 0, 0),
			added_anims=[Create(gauss_sphere)]
		)
		self.wait()

		ds_surface = Sphere(
			radius=math.dist(ORIGIN, point_m_coords) + 0.02,
			u_range=(0.75, 0.9), v_range=(2.25, 2.35), resolution=(1,1),
			fill_opacity=0.5, checkerboard_colors=[GREY], stroke_width=3
		)
		# ds_vector = Arrow3D(
		# 	start=ds_surface.get_center(),
		# 	end=spherical_to_cartesian(*(cartesian_to_spherical(ds_surface.get_center()) + 1.3 * RIGHT)),
		# 	thickness=0.01,
		# 	height=0.2,
		# 	base_radius=0.05
		# )
		ds_vector = Arrow(
			start=ds_surface.get_center(),
			end=spherical_to_cartesian((cartesian_to_spherical(ds_surface.get_center()) + 1.3 * RIGHT)),
			buff=0,
			max_tip_length_to_length_ratio=0.5
		)
		ds_t = Tex('$\\overrightarrow{\\mathrm{dS}}$')\
                    .move_to(ds_vector.end)\
                    .shift(0.2 * IN + 0.6 * LEFT)\
                    .scale(0.8)
		
		self.add_fixed_orientation_mobjects(ds_t)
		self.play(
			Create(ds_surface),
			Create(ds_vector),
			Write(ds_t)
		)
		self.begin_ambient_camera_rotation(rate=0.06)
		self.wait(2)

		gauss_t = MathTex(
			r'\oiint_S\overrightarrow{E}.\overrightarrow{\mathrm{dS}}=\frac{Q_{int}}{\epsilon_{0}}',
			tex_template=esint_template,
			tex_to_color_map={
				'_S': RED_B,
				'\overrightarrow{E}': '#fcba03' 
			}
		).to_edge(LEFT).shift(2 * UP)

		self.add_fixed_in_frame_mobjects(gauss_t)

		self.play(Write(gauss_t))

		self.wait()

		q_t = Tex(r'$\frac{Q_{int}}{\epsilon_{0}} = $').next_to(gauss_t, direction=DOWN, aligned_edge=LEFT)
		q_vt = ValueTracker(10)
		q_value = DecimalNumber(
			number=10,
			num_decimal_places=1,
			include_sign=True,
			background_stroke_opacity=0
		).next_to(q_t, direction=RIGHT)

		def q_value_updater(m):
			m.set_value(q_vt.get_value())
			self.add_fixed_in_frame_mobjects(m)

		self.add_fixed_in_frame_mobjects(q_t, q_value)
		self.play(Write(q_t), Write(q_value))
		q_value.add_updater(q_value_updater)
		
		self.wait()

		sphere_radius = math.dist(ORIGIN, point_m_coords)

		e_field = [
			Arrow(
				start=spherical_to_cartesian((sphere_radius, theta, phi)),
				end=spherical_to_cartesian((sphere_radius + 1.2, theta, phi)),
				buff=0,
				max_tip_length_to_length_ratio=0.5
			).set_color('#fcba03')
			for theta in np.linspace(PI / 8, PI - PI / 8, num=4)
			for phi in np.linspace(0, TAU, num=8)
		]

		self.play(*map(Create, e_field))

		self.wait()

		self.play(
			ApplyMethod(q_vt.set_value, 20),
			*map(lambda a: ApplyMethod(
				a.put_start_and_end_on,
				a.get_start(),
				spherical_to_cartesian((cartesian_to_spherical(a.get_end()) + 0.8 * RIGHT))
				),
				e_field
			)
		)

		self.wait(3)

		self.play(
			ApplyMethod(q_vt.set_value, -10),
			*map(lambda a: ApplyMethod(
				a.put_start_and_end_on,
				a.get_start(),
				spherical_to_cartesian((cartesian_to_spherical(a.get_start()) + LEFT))
				),
				e_field
			)
		)

		self.wait(3)

		self.play(
			ApplyMethod(q_vt.set_value, 20),
			*map(lambda a: ApplyMethod(
				a.put_start_and_end_on,
				a.get_start(),
				spherical_to_cartesian((cartesian_to_spherical(a.get_end()) + 3 * RIGHT))
                ),
				e_field
			)
		)

		self.wait(2)

		self.play(
			*map(FadeOut, [point_m, point_m_t, ds_surface, ds_vector, ds_t])
		)

		self.move_camera(
			zoom=0.70,
			added_anims=[
				gauss_sphere.animate(run_time=5).scale(1.5),
				*map(lambda a: ApplyMethod(
					a.put_start_and_end_on,
					spherical_to_cartesian((cartesian_to_spherical(a.get_start()) + 1.9 * RIGHT)),
					spherical_to_cartesian((cartesian_to_spherical(a.get_end()) + RIGHT)),
					run_time=5
					),
					e_field
				)
			],
			run_time=5
		)
		self.wait(2)

class AnaDonnees(Scene):
	def construct(self):
		point_list = [
			(1, 2, 0), (1.3, 1.5, 0), (1.8, 1, 0), (2, 1.8, 0),
			(-1, -2, 0), (-1.5, -1, 0), (-1.8, -2, 0)
		]
		a = list(map(sum, zip(*point_list)))
		print(a)
		
		axes = Axes(
			x_length=config.frame_width,
			y_length=config.frame_height,
			tips=False
		).add_coordinates()

		plane = NumberPlane(
			x_length=config.frame_width,
			y_length=config.frame_height
		)

		dots = map(lambda coords: Dot(plane.c2p(*coords), color=RED), point_list)
		
		self.add(plane)
		self.play(*map(Create, dots))
		self.wait()

class Hall(ThreeDScene):
	def construct(self):
		self.set_camera_orientation(
			phi=60 * DEGREES,
			theta=-60 * DEGREES
		)

		axes = ThreeDAxes().set_stroke(width=0.5)
		z_axis_label = axes.get_z_axis_label("z")
		x_y_axis_label = axes.get_axis_labels()

		fil = Prism().set_fill(RED_E, opacity=0.5)

		self.play(*map(Create, [axes, z_axis_label, x_y_axis_label, fil]))
		# self.wait(10)

class Invariance(ThreeDScene):
	def construct(self):
		self.set_camera_orientation(
			phi=70 * DEGREES,
			theta=20 * DEGREES,
			zoom=1.2
		)
		
		axes = ThreeDAxes().set_stroke(width=0.5)
		sphere = Sphere(
			radius=2.5,
			resolution=(18, 18),
			fill_opacity=0.4
		).set_color(BLUE_E)

		circles = [
			VGroup(*[
			Circle(
				radius=r,
				color=RED_C,
				stroke_opacity=(r / (math.sqrt(sphere.radius**2 - d**2)))**4
			).rotate(PI / 2, axis=UP).shift(d * RIGHT)
			for r in np.linspace(0, math.sqrt(sphere.radius**2 - d**2), 80)
			])
			for d in np.linspace(sphere.radius, -sphere.radius, 60)
		]
		circles_theta = circles[30].copy()

		distribution_t = MathTex(
			r"\forall P \in \mathscr{S},\,\rho(P)=\rho_0\left(\frac{r}{R}\right)^2"
		).to_corner(UL).scale(0.8)
		self.add_fixed_in_frame_mobjects(distribution_t)

		self.play(
			Create(axes),
			Create(sphere),
			Write(distribution_t),
			run_time=0.7
		)
		self.wait()

		self.move_camera(
			theta=50 * DEGREES,
			added_anims=[LaggedStartMap(
				FadeIn, circles, rate_func=there_and_back, run_time=4.5, lag_ratio=0.5
				)],
			run_time=4.5,
			rate_func=lambda t: smooth(t, inflection=5.0)
		)
		self.wait(2)


		m_coords = np.array([0, 3.5, 1.7])
		theta_vt = ValueTracker(cartesian_to_spherical(m_coords)[2])

		# Updaters:
		def m_theta_updater(m, vt):
			new_coords = cartesian_to_spherical(m.get_center())
			new_coords[2] = vt.get_value()
			m.move_to(spherical_to_cartesian(new_coords))

		
		def m_t_updater(t, dot):
			t.next_to(dot, UP)
		
		def line_updater(l, dot):
			l.become(DashedLine(ORIGIN, dot.get_center()))
		
		def theta_updater(a, vt):
			a.become(Arc(
					radius=0.8,
                    angle=vt.get_value()
				).rotate(
					PI / 2,
					axis=DOWN,
					about_point=UP
				)
			)

		m = Dot3D(m_coords, radius=0.06).scale(1.2)
		m_t = Tex('M').next_to(m, UP)
		line = DashedLine(ORIGIN, m_coords)
		theta_angle = Arc(
			radius=0.8,
			angle=theta_vt.get_value()
		).rotate(
			PI / 2,
			axis=DOWN,
			about_point=UP
		)
		theta_t = MathTex(r"\theta_M").next_to(theta_angle, IN, buff=SMALL_BUFF)
		m.add_updater(lambda m: m_theta_updater(m, theta_vt))
		m_t.add_updater(lambda t: m_t_updater(t, m))
		line.add_updater(lambda l: line_updater(l, m))
		theta_angle.add_updater(lambda a: theta_updater(a, theta_vt))
		self.add_fixed_orientation_mobjects(m_t, theta_t)


		self.play(
			FadeIn(circles_theta),
			Create(m),
			Create(theta_angle),
			Create(line),
			Write(m_t),
			Write(theta_t)
		)
		self.wait(2)

		self.move_camera(
			theta=0,
			phi=90 * DEGREES,
			added_anims=[FadeOut(sphere)]
		)
		self.wait()
		self.play(
			theta_vt.animate.set_value(3 * PI / 4),
			rate_func=lambda t: smooth(t, inflection=5.0),
			run_time=3
		)
		self.wait(0.7)
		self.play(
			theta_vt.animate.set_value(PI / 3),
			rate_func=lambda t: smooth(t, inflection=5.0),
			run_time=3
		)
		self.wait()

		self.move_camera(
			phi=70 * DEGREES,
			theta=50 * DEGREES,
			added_anims=[FadeIn(sphere)]
		)
		for mobject in (theta_angle, line, m):
			mobject.clear_updaters()
		self.play(
			*map(FadeOut, [circles_theta, theta_angle, line, theta_t])
		)

		circles_phi = VGroup(*[
			Circle(
				radius=r,
				color=RED_C,
				stroke_opacity=(r / (math.sqrt(sphere.radius**2 - 1.5**2)))**4
			).shift(1.5 * OUT)
			for r in np.linspace(0, math.sqrt(sphere.radius**2 - 1.5**2), 80)
			])
		
		self.play(
			Create(circles_phi)
		)
		self.wait()

		def phi_updater(a, vt):
			a.become(Arc(
					radius=0.8,
                    angle=vt.get_value()
				)
			)
		
		def m_phi_updater(m, vt):
			new_coords = cartesian_to_spherical(m.get_center())
			new_coords[1] = vt.get_value()
			m.move_to(spherical_to_cartesian(new_coords))

		phi_vt = ValueTracker(PI / 2)
		phi_angle = Arc(
			radius=0.8,
			angle=phi_vt.get_value()
		)

		phi_t = MathTex(r"\phi_M").next_to(phi_angle, RIGHT, buff=SMALL_BUFF)
		self.add_fixed_orientation_mobjects(phi_t)

		self.move_camera(
			phi=0,
			theta=0,
			added_anims=[
				FadeOut(sphere),
				*map(FadeIn, [phi_t, line, phi_angle])
			]
		)
		m.add_updater(lambda m: m_phi_updater(m, phi_vt))
		line.add_updater(lambda l: line_updater(l, m))
		phi_angle.add_updater(lambda a: phi_updater(a, phi_vt))
		self.wait()

		self.play(
			phi_vt.animate.set_value(5 * PI / 6),
			rate_func=lambda t: smooth(t, inflection=5.0),
			run_time=3
		)
		self.wait(0.7)
		self.play(
			phi_vt.animate.set_value(PI / 4),
			rate_func=lambda t: smooth(t, inflection=5.0),
			run_time=3
		)
		
		self.wait(3)