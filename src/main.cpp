#include "types.hpp"
#include <array>
#include <assert.h>
#include <iostream>
#include <memory>

// Question 1. How do you determine if a 3D point lies inside a 3D sphere?
// A: If the distance from the center of the sphere to the point is less than
// the radius, the point is inside the sphere.
bool point_inside_sphere(const vec& p, const vec& sphere_p, float sphere_r)
{
	return p.dist(sphere_p) < sphere_r;
}

// Question 2: How do you determine if the path of a projectile passes within 1
// foot of a location? A: Center a sphere on the location, check if the path
// crosses the sphere twice.

// If the line is straight, a line passes WITHIN a foot of a location if it intersects a sphere
// centered at the location with the radius of 1ft twice.
bool line_passes_through_location(const vec& lineStart, const vec& lineEnd, const vec& sphereCenter)
{
	// Special case where the line starts or ends inside of the location
	if (lineStart.dist(sphereCenter) < 1.f || lineEnd.dist(sphereCenter))
	{
		return true;
	}

	// Otherwise check that the line intersects twice.
	return line_sphere_intersections(lineStart, lineEnd, sphereCenter, 1.f).num_intersections == 2;
}



// Question 3: The following diagram contains 5 squares in a 2x2 configuration
// (four 1x1 + one 2x2): Implement a function in either C++ that takes a
// parameter of N, where N is the size of the NxN matrix of any size greater
// than 0 and returns the total number of squares in the matrix.

// Returns the number of squares that can fit inside of a square NxN matrix
// without overlapping. This is this integer sequence: https://oeis.org/A222548
int squares(int N)
{
	int sum = 0;
	for (int i = 1; i <= N; ++i)
	{
		sum += (N / i) * (N / i);
	}
	return sum;
}

// Returns the number of squares that can fit inside of a square NxN matrix WITH
// overlapping.
int squares_overlapping(int N) { return (N * (N + 1) * ((2 * N) + 1)) / 6; }

// Question 4: Design the abstract class for a game particle system in C++. The
// system should support creating multiple particle groups in the game world.
// Show only the interface to the system, without implementation. Assume the
// existence of any math libraries or other systems you need for the interface.
//      a. Each particle group has a type id that defines (via a lookup table
//      not shown here) the number of particles, the lifetime of each particle,
//      the velocities of the particles, etc. 
//		b. Given a type id, particle groups can be created.
//		c. Particle groups are attached to models, placed at a location, or exist only in screen space.
//		d. Particle groups can be scaled and rotated.
//		e. Particle groups can be ended early.

class transform_interface
{
public:
	virtual mat get_world_transform() = 0;
	virtual mat get_local_transform() = 0;
};

class particle_system
{
public:
	struct particle_group&
	get_particle_group(const struct particle_group_type_id& type_id);
	struct particle_group*
	create_particle_group(const struct particle_group_type_id& type_id);
};

struct particle_group final: public transform_interface
{
	void set_parent(transform_interface* parent);

	virtual mat get_local_transform() override;
	virtual mat get_world_transform() override;
};

int main()
{
	std::cout << (point_inside_sphere(vec {0.5, 0, 0}, vec {0, .1, 0}, 0.6) ? "true" : "false")
	          << "\n";
	return 0;
}