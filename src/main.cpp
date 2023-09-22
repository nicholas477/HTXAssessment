#include "types.hpp"
#include <array>
#include <assert.h>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include "gtest/gtest.h"

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
	if (lineStart.dist(sphereCenter) < 1.f || lineEnd.dist(sphereCenter) < 1.f)
	{
		return true;
	}

	// Otherwise check that the line intersects twice.
	return line_sphere_intersections(lineStart, lineEnd, sphereCenter, 1.f).num_intersections == 2;
}

// If the line is a parabola, the parabola passes WITHIN a foot of a location if the minimum distance
// between the location and the parabola is less than a foot
static bool parabola_passes_location(const vec& starting_position, const vec& launch_velocity, const vec& sphere_p, float sphere_r, const vec& gravity)
{
	// All 3d parabolic arcs are really 2d, so just find the two vectors that define the 2d plane of the arc
	// TODO: reject straight line parabolas (parabolas where the X and Y basis vectors are parallel)
	const vec arc_up_vec           = gravity.normalize();
	const vec arc_forward_velocity = launch_velocity - (arc_up_vec * (arc_up_vec.dot(launch_velocity)));
	const vec arc_forward_vec      = arc_forward_velocity.normalize();
	const vec arc_left_vec         = arc_up_vec.cross(arc_forward_vec);

	// The arc up vector will define h and arc_forward_vector will define t
	// The parabola is defined by h(t) = at^2 + bt + c, and it starts at starting_position
	double a      = gravity.length();
	double b      = arc_up_vec.dot(launch_velocity);
	const float c = 0.f; // Since C is just a vertical offset set it to zero

	// Convert the parabola to a function of x as opposed to t
	// x = t * velocity, so just divide by velocity
	a /= (arc_forward_velocity.squared_length());
	b /= arc_forward_velocity.length();

	// Check the distance from the parabola plane to the sphere, and if its greater than the radius we can early out
	const float dist_from_sphere = point_plane_distance(sphere_p, starting_position, arc_left_vec);
	if (abs(dist_from_sphere) >= sphere_r)
	{
		return false;
	}

	// Now slice the sphere along the parabola plane to get a circle
	const vec sphere_plane_intersect = point_plane_project(sphere_p, starting_position, arc_left_vec);
	const vec circle_2dpos           = vec(arc_forward_vec.dot(sphere_plane_intersect - starting_position), arc_up_vec.dot(sphere_plane_intersect - starting_position), 0.f);
	const float circle_r             = get_sphere_slice_radius(dist_from_sphere, sphere_r);

	// Find the minimum distance between the circle center and the parabola
	const float minimum_distance = parabola_minimum_distance(a, b, c, circle_2dpos.x, circle_2dpos.y);
	return minimum_distance < circle_r;
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

// Assume this class has support for moving/scaling/rotating using matricies, vectors, etc
class transform
{
public:
	void set_parent(transform* parent); // Attaches this component to another component. Maintains local transform. Pass in nullptr to detach.
	transform* get_parent() const;
};

// Just to get this to compile
using FName = std::string;

namespace particles
{
	// Particle prefab type. Multiple emitters can be of the same type.
	// For example, if you had a fire particle emitter, it might be named "fire_01"
	struct type
	{
		FName id;

		struct description
		{
			// Num particles, velocities, lifetimes, etc here
		};

		bool operator==(const type& other)
		{
			return id == other.id;
		}
	};

	// Each emitter spawns an instance of the provided particle type
	class emitter
	{
	public:
		transform transform;

		bool create(FName type);
		void release();

		bool is_alive() const;
		bool is_screenspace() const;
		float get_lifetime() const;
	};

	class manager
	{
	public:
		// Holds references to all of the emitters in the scene for rendering
		std::unordered_set<emitter*> emitters;

		// Lookup table containing all of the registered particle types
		std::unordered_set<type> ParticleTypes;

		// Register a new particle type for use with particle emitters.
		bool add_type(const type& new_type);

		static manager& get();
	};
}

// Question 5: Define a list of unit tests to validate an implementation of the particle system interface from the previous question.

TEST(ParticleSystemTest, SpawnParticleEmitters)
{
	particles::type new_particle_type;
	new_particle_type.id = "particles/fire";

	EXPECT_EQ(particles::manager::get().add_type(new_particle_type), true);

	particles::emitter emitter;
	emitter.create("new_particle_type.id");
	EXPECT_EQ(emitter.is_alive(), true);
	EXPECT_EQ(emitter.get_lifetime() > 0.f, true);
}


// TODO: test transform
TEST(ParticleSystemTest, TestParticleGroupTransform)
{
}

// Question 6: Write a boid system using any engine you prefer. You can use any libraries or tech you want.
// The only stipulation is that the system needs to simulate a flying insect swarm instead of the normal
// birds/fish. The deliverable on this is the code as well as a screenshot (or video) that you think best
// exemplifies the behavior.

// Answer: https://github.com/nicholas477/InsectBoids/tree/main



// Question 7: Consider a parkour based movement system. What would be the MVP for such a system? What
// features would you add immediately after MVP, and what features could you postpone?

int main()
{

	return 0;
}