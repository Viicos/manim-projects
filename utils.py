"""Utility methods."""

from manim import *
import numpy as np
import math


def cylindric_to_cartesian(point):  # point = (r, theta, z)
    r, theta, z = point
    return np.array([r * np.cos(theta), r * np.sin(theta), z])

def cartesian_to_cylindric(point):
    x, y, z = point
    return np.array([math.sqrt(x**2 + y**2), math.atan2(y, x) , z])

def create_spheric_base(point, length=1):
    _, theta, phi = cartesian_to_spherical(point)

    e_r = Arrow(
        start=point,
        end=spherical_to_cartesian(cartesian_to_spherical(point) + length * RIGHT),
        buff=0
    )
    e_theta = Arrow(
        start=point,
        end=np.array(point) + length * np.array([
            np.cos(theta) * np.cos(phi),
            np.sin(theta) * np.cos(phi),
            -np.sin(phi)
            ]),
        buff=0
    )
    e_phi = Arrow(
        start=point,
        end=np.array(point) + length * np.array([-np.sin(theta), np.cos(theta), 0]),
        buff=0
    )
    return (e_r, e_theta, e_phi)

def create_cylindric_base(point, length=1):
    _, theta, _ = cartesian_to_cylindric(point)

    e_r = Arrow(
        start=point,
        end=cylindric_to_cartesian(cartesian_to_cylindric(point) + length * RIGHT),
        buff=0
    )
    e_theta = Arrow(
        start=point,
        end=np.array(point) + length * np.array([
                -np.sin(theta),
                np.cos(theta),
                0
            ]),
        buff=0
    )
    e_z = Arrow(
        start=point,
        end=np.array(point) + length * OUT,
        buff=0
    )
    return (e_r, e_theta, e_z)

def parametric_plane_func(u, v, v1, v2, point):
    return v1 * u + v2 * v + point

def parametric_plane_to_cartesian(v1, v2, point):
    """
    Get the cartesian equation of a plane from its parametric equation:
    - original equation: `u * v1 + v * v2 + point`
    - target: `ax + by + cz + d = 0`
    """
    normal_vector = np.cross(v1, v2)
    d = -sum(np.multiply(normal_vector, point))
    return (*normal_vector, d)

def get_intersection_of_point_normal_to_plane(point, plane_coords):
    """
    Get the intersection of the point normal to the plane.
    See https://math.stackexchange.com/questions/1215760/intersection-of-point-normal-to-plane
    """
    return point - sum(np.multiply(point, plane_coords)) / sum(np.square(plane_coords)) * plane_coords

def get_symmetric_point(point, plane_coords):
    normal_vector = point - get_intersection_of_point_normal_to_plane(point, plane_coords)

    return point - 2 * abs(sum(np.multiply(point, plane_coords))) / np.linalg.norm(plane_coords) * (normal_vector / np.linalg.norm(normal_vector))


class UnchargedParticle(Circle):
    """A neutral particle.

    Parameters
    ----------
    point : Union[:class:`list`, :class:`numpy.ndarray`], optional
        The location of the dot.
    radius : Optional[:class:`float`]
        The radius of the dot.
    fill_color : :class:`~.Colors`, optional
        The fill color of the particle.
    """

    def __init__(
        self,
        point=ORIGIN,
        radius=0.12,
        fill_color=LIGHT_GREY
    ):
        super().__init__(
            arc_center=point,
            radius=radius,
            stroke_width=1.3,
            stroke_color=WHITE,
            fill_color=fill_color,
            fill_opacity=0.8
        )
        
class Proton(Circle):
    """A positive particle.

    Parameters
    ----------
    point : Union[:class:`list`, :class:`numpy.ndarray`], optional
        The location of the dot.
    radius : Optional[:class:`float`]
        The radius of the dot.
    fill_color : :class:`~.Colors`, optional
        The fill color of the particle.
    """

    def __init__(
        self,
        point=ORIGIN,
        radius=0.12,
        fill_color=RED
    ):
        super().__init__(
            arc_center=point,
            radius=radius,
            stroke_width=1.3,
            stroke_color=WHITE,
            fill_color=fill_color,
            fill_opacity=0.8
        )

class Ruler(Rectangle):
    """A ruler consisting of a rectangle and dashes.
    
    """
    def __init__(
        self,
        width: float = 2.0,
        num_dashes: int = 3,
        **kwargs
        ):
            super().__init__(width=width, height=width / 2, **kwargs)
            dashes = VGroup(*[
                Line(ORIGIN, (width / 8) * DOWN)
                for _ in range(num_dashes)
            ]).arrange(buff=width / num_dashes).shift((width / 5.5) * UP)
            self.add(dashes)
