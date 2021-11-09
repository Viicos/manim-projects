from manim import *
import numpy as np
import math

def create_spheric_base(point):
    _, theta, phi = cartesian_to_spherical(point)

    e_r = Arrow(
        start=point,
        end=spherical_to_cartesian(cartesian_to_spherical(point) + 1.6 * RIGHT),
        buff=0
    )
    e_theta = Arrow(
        start=point,
        end=np.array(point) + 1.3 * np.array([
            np.cos(theta) * np.cos(phi),
            np.sin(theta) * np.cos(phi),
            -np.sin(phi)
            ]),
        buff=0
    )
    e_phi = Arrow(
        start=point,
        end=np.array(point) + 1.3 * np.array([-np.sin(theta), np.cos(theta), 0]),
        buff=0
    )
    return (e_r, e_theta, e_phi)

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
