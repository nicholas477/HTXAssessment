#pragma once

#include <cmath>

struct vec
{
    float x,y,z;

    float dist(const vec& other) const
    {
        vec vec_btwn;
        vec_btwn.x = x - other.x;
        vec_btwn.y = y - other.y;
        vec_btwn.z = z - other.z;
        return std::sqrt(vec_btwn.x * vec_btwn.x + vec_btwn.y * vec_btwn.y + vec_btwn.z * vec_btwn.z);
    }

    float dot(const vec& other) const
    {
        return (x * other.x) + (y * other.y) + (z * other.z);
    }
};

// 4x4 matrix
struct mat
{
    std::array<float, 16> items;
};

static vec operator-(const vec& lhs, const vec& rhs)
{
    vec out;
    out.x = lhs.x - rhs.x;
    out.y = lhs.y - rhs.y;
    out.z = lhs.z - rhs.z;
    return out;
}