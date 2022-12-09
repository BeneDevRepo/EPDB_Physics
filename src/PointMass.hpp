#pragma once

#include "BMath/vector.h"


struct PointMass {
    float mass;
    bmath::vec2 prevPos;
    bmath::vec2 pos;
    bmath::vec2 vel{}; // recalculated every frame
    bmath::vec2 force{};

    inline PointMass():
        prevPos(), pos(), mass(1.f) { }
    inline PointMass(const float px, const float py, const float mass = 1.f):
        prevPos({px, py}), pos({px, py}), mass(mass) { }
    inline PointMass(const bmath::vec2& pos, const float mass = 1.f):
        prevPos(pos), pos(pos), mass(mass) { }

    inline void calcVel(const float dt) { // calculates velocity based on prevPos, pos and dt, and resets force
        vel = (pos - prevPos) / dt;
        force = {};
    }
    void applyForce(const bmath::vec2& force) {
        this->force += force;
    }
    inline void update(const float dt) { // moves particle based on vel and dt
        // integrate acceleration:
        const bmath::vec2 accel = force / mass;
        vel += accel * dt;
    
        // integrate velocity:
        prevPos =  pos;
        pos += vel * dt;
    }
};