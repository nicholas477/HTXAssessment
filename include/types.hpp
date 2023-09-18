#pragma once

#include <array>
#include <cmath>
#include <sstream>
#include <vector>

struct vec
{
	union
	{
		struct
		{
			float x, y, z;
		};
		float xyz[3];
	};

	vec(float _x = 0.f, float _y = 0.f, float _z = 0.f)
	    : x(_x)
	    , y(_y)
	    , z(_z)
	{
	}

	float dist(const vec& other) const
	{
		vec vec_btwn;
		vec_btwn.x = x - other.x;
		vec_btwn.y = y - other.y;
		vec_btwn.z = z - other.z;
		return vec_btwn.length();
	}

	float length() const
	{
		return std::sqrt(squared_length());
	}

	float squared_length() const
	{
		return x * x + y * y + z * z;
	}

	float dot(const vec& other) const
	{
		return (x * other.x) + (y * other.y) + (z * other.z);
	}

	vec cross(const vec& other) const
	{
		vec out;
		out.xyz[0] = xyz[1] * other.xyz[2] - xyz[2] * other.xyz[1];
		out.xyz[1] = -(xyz[0] * other.xyz[2] - xyz[2] * other.xyz[0]);
		out.xyz[2] = xyz[0] * other.xyz[1] - xyz[1] * other.xyz[0];
		return out;
	}

	std::string to_string() const
	{
		std::stringstream out_stream;
		out_stream << "(";
		out_stream << "x: " << x;
		out_stream << ", y: " << y;
		out_stream << ", z: " << z;
		out_stream << ")";
		return out_stream.str();
	}

	vec operator/(float rhs) const
	{
		vec out;
		out.x = x / rhs;
		out.y = y / rhs;
		out.z = z / rhs;
		return out;
	}

	vec operator*(float rhs) const
	{
		vec out;
		out.x = x * rhs;
		out.y = y * rhs;
		out.z = z * rhs;
		return out;
	}

	vec normalize() const
	{
		return (*this) / length();
	}

	vec operator-(const vec& rhs) const
	{
		vec out;
		out.x = x - rhs.x;
		out.y = y - rhs.y;
		out.z = z - rhs.z;
		return out;
	}

	vec operator+(const vec& rhs) const
	{
		vec out;
		out.x = x + rhs.x;
		out.y = y + rhs.y;
		out.z = z + rhs.z;
		return out;
	}
};

static vec operator*(float lhs, const vec& rhs)
{
	vec out;
	out.x = lhs * rhs.x;
	out.y = lhs * rhs.y;
	out.z = lhs * rhs.z;
	return out;
}

// A plane that bisects 3d space infinitely.
// Imagine this as a normal vector based on the origin, multiplied by w. The tip of that vector is where space is
// bisected in half along a plane orthogonal to that vector
struct plane
{
	float x, y, z, w;

	plane(const vec& base, const vec& normal)
	{
		x = normal.x;
		y = normal.y;
		z = normal.z;
		w = base.dot(normal);
	}
};

// Returns the distance between a point and a plane
// Copied from unreal engine source code
static float point_plane_distance(const vec& point, const vec& plane_base, const vec& plane_normal)
{
	return (point - plane_base).dot(plane_normal);
}

// Projects a point onto a plane
// Copied from unreal engine source code
static vec point_plane_project(const vec& point, const vec& plane_base, const vec& plane_normal)
{
	//Find the distance of X from the plane
	//Add the distance back along the normal from the point
	return point - (point_plane_distance(point, plane_base, plane_normal) * plane_normal);
}

struct intersection_result
{
	int num_intersections = 0;
	std::array<vec, 2> intersections;
};

// Gives you the number of times a line passes through a sphere and where it intersects
// This can either be 0 times, 1 times (tangent), or 2 times
static intersection_result line_sphere_intersections(const vec& line_start, const vec& line_end, const vec& sphere_p, float sphere_r)
{
	intersection_result result;

	// Calculate the direction vector of the line
	vec line_direction = line_end - line_start;

	// Calculate the vector from the sphere center to the line start point
	vec sphere_to_line_start = line_start - sphere_p;

	// Calculate the coefficients for the quadratic equation
	float a = line_direction.dot(line_direction);
	float b = 2.0 * line_direction.dot(sphere_to_line_start);
	float c = sphere_to_line_start.dot(sphere_to_line_start) - (sphere_r * sphere_r);

	// Calculate the discriminant
	float discriminant = b * b - 4 * a * c;

	if (discriminant > 0)
	{
		// Calculate the two possible t-values (parameter along the line)
		float t1 = (-b + sqrt(discriminant)) / (2 * a);
		float t2 = (-b - sqrt(discriminant)) / (2 * a);

		// Calculate the intersection points

		// The t value could be past the beginning or the end of the lines, so make sure they're in range
		if (t1 >= 0.f && t1 <= 1.f)
		{
			result.intersections[result.num_intersections] = line_start + (line_direction * t1);
			result.num_intersections++;
		}

		if (t2 >= 0.f && t2 <= 1.f)
		{
			result.intersections[result.num_intersections] = line_start + (line_direction * t2);
			result.num_intersections++;
		}
	}
	else if (discriminant == 0)
	{
		// One intersection point (tangent)
		result.num_intersections = 1;

		float t = -b / (2 * a);

		if (t >= 0.f && t <= 1.f)
		{
			result.intersections[result.num_intersections] = line_start + (line_direction * t);
			result.num_intersections++;
		}
	}
	else
	{
		// No intersection
		result.num_intersections = 0;
	}

	return result;
}

// Solves cubic equations in the form: x^3 + a*x^2 + b*x + c = 0
// Only returns the real roots

// http://math.ivanovo.ac.ru/dalgebra/Khashin/poly/index.html
// return 3: 3 real roots vec[0], vec[1], vec[2]
// return 1: 1 real root vec[0]. Does not return the other two complex roots
static std::vector<double> SolveP3(double a, double b, double c)
{
	constexpr double eps = 1e-14;
	const double two_pi  = 6.28318530718;

	double a2 = a * a;
	double q  = (a2 - 3 * b) / 9;
	double r  = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
	// equation y^3 - 3q*y + r/2 = 0 where x = y-a/3
	if (fabs(q) < eps)
	{
		if (fabs(r) < eps)
		{ // three identical roots
			return {-a / 3};
		}
		return {std::cbrt(-r / 2)};
	}
	double r2 = r * r;
	double q3 = q * q * q;
	double A, B;
	if (r2 <= (q3 + eps))
	{
		double t = r / sqrt(q3);
		if (t < -1)
			t = -1;
		if (t > 1)
			t = 1;
		t = acos(t);
		a /= 3;
		q = -2 * sqrt(q);
		return {q * cos(t / 3) - a, q * cos((t + two_pi) / 3) - a, q * cos((t - two_pi) / 3) - a};
	}
	else
	{
		A = -std::cbrt(fabs(r) + sqrt(r2 - q3));
		if (r < 0)
			A = -A;
		B = (A == 0 ? 0 : B = q / A);

		a /= 3;
		if (fabs(0.5 * sqrt(3.) * (A - B)) < eps)
		{
			return {(A + B) - a, -0.5 * (A + B) - a};
		}
		return {(A + B) - a};
	}
}

// Returns the distance between a point (px, py) and a parabola y = ax^2 + bx + c
static double parabola_point_distance(double a, double b, double c, double px, double py, double x)
{
	const double y = a * x * x + b * x + c;
	return sqrt((px - x) * (px - x) + (py - y) * (py - y));
}

// Finds the minimum distance between a parabola and a point (px,py)
// https://stackoverflow.com/questions/9800324/how-to-find-the-distance-between-a-point-and-a-parabola-in-code
// Finds the real roots of the derivative of the squared distance function, and then returns the minimum distance
//
// TODO: add in code for rejecting distances with x < 0 since our parabola doesn't extend pass x < 0
static double parabola_minimum_distance(double a, double b, double c, double px, double py)
{
	// Find the x^3, x^2, x coefficients
	// Full equation for the derivative of distance squared is:
	// 2bc - 2p(x) - 2bp(y) + 2x + 2b^2x + 4acx - 4axp(y) + 6abx^2 + 4a^2x^3
	const double x_cube_coefficient = 4 * a * a;
	double x_square_coefficient     = 6 * a * b;
	double x_coefficient            = (4 * a * c) - (4 * a * py) + (2 * b * b) + 2;
	double constants                = 2 * b * c - 2 * b * py - 2 * px;

	// Change it to the form x^3 + ax^2 + bx + c by dividing by the x cubed coefficient
	// We do this so we can feed it into the cubic equation solver
	x_square_coefficient /= x_cube_coefficient;
	x_coefficient /= x_cube_coefficient;
	constants /= x_cube_coefficient;

	// Now we have the roots of the derivative, so plug it back into the function and find the minimum distance
	std::vector<double> roots = SolveP3(x_square_coefficient, x_coefficient, constants);
	double min_distance       = parabola_point_distance(a, b, c, px, py, roots[0]);
	for (int i = 1; i < roots.size(); ++i)
	{
		min_distance = std::min(min_distance, parabola_point_distance(a, b, c, px, py, roots[i]));
	}

	return min_distance;
}

// Gets the radius of a sphere slice. The sphere slice is defined by the distance from the center of the sphere
static double get_sphere_slice_radius(double dist_from_center, double sphere_r)
{
	const double cap_h = sphere_r - std::abs(dist_from_center);
	return std::sqrt(cap_h * ((2 * sphere_r) - cap_h));
}