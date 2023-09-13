#pragma once

#include <array>
#include <cmath>
#include <sstream>

struct vec
{
	float x, y, z;

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
		return std::sqrt(x * x + y * y + z * z);
	}

	float dot(const vec& other) const
	{
		return (x * other.x) + (y * other.y) + (z * other.z);
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

// 4x4 matrix
struct mat
{
	std::array<float, 16> items;
};

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

// Calculates the intersection between a parabola and a circle. This is in 2D
static intersection_result parabola_circle_intersections(float a, float b, float c, const vec& center, float radius)
{
	intersection_result result;

	// Calculate the discriminant of the quadratic equation
	float discriminant = b * b - 4 * a * c;

	if (discriminant >= 0)
	{
		// Calculate the two possible x-values where the parabola intersects the circle
		float x1 = (-b + sqrt(discriminant)) / (2 * a);
		float x2 = (-b - sqrt(discriminant)) / (2 * a);

		// Calculate the corresponding y-values
		float y1 = a * x1 * x1 + b * x1 + c;
		float y2 = a * x2 * x2 + b * x2 + c;

		// Check if the intersection points are within the circle's boundary
		float distance1 = sqrt(pow(x1 - center.x, 2) + pow(y1 - center.y, 2));
		float distance2 = sqrt(pow(x2 - center.x, 2) + pow(y2 - center.y, 2));

		if (distance1 <= radius)
		{
			result.intersections[result.num_intersections] = vec {x1, y1, 0.f};
			result.num_intersections++;
		}
		if (distance2 <= radius)
		{
			result.intersections[result.num_intersections] = vec {x2, y2, 0.f};
			result.num_intersections++;
		}
	}

	return result;
}