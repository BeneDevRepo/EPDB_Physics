#pragma once

#include "BMath/vector.h"

#include <vector>


struct PointMass;

extern std::vector<PointMass> points;


struct PointMass {
    float mass;
    bmath::vec<2, double> prevPos;
    bmath::vec<2, double> pos;
    bmath::vec<2, double> vel{}; // recalculated every frame
    bmath::vec<2, double> force{};

    inline PointMass():
        mass(1.f), prevPos(), pos()  { }
    inline PointMass(const float px, const float py, const float mass = 1.f):
        mass(mass), prevPos({px, py}), pos({px, py}) { }
    inline PointMass(const bmath::vec<2, double>& pos, const float mass = 1.f):
        mass(mass), prevPos(pos), pos(pos) { }

    inline void calcVel(const float dt) { // calculates velocity based on prevPos, pos and dt, and resets force
        vel = (pos - prevPos) / dt;
        force = {};
    }
    inline void applyForce(const bmath::vec<2, double>& force) {
        this->force += force;
    }
    inline void update(const float dt) { // moves particle based on vel and dt
        // integrate acceleration:
        const bmath::vec<2, double> accel = force / mass;
        vel += accel * dt;
    
        // integrate velocity:
        prevPos =  pos;
        pos += vel * dt;
    }

	// inline size_t id() const {
	// 	return this - points.data();
	// }
};