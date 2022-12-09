#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cmath>
// #include <thread>
// #include <mutex>

#include "BWindow/GDIWindow.h"
#include "BMath/vector.h"

#include "PointMass.hpp"

#include "windows.h"


struct Constraint {
	virtual ~Constraint() = default;
	virtual inline void solve(const float dt) = 0;
};

struct Triangle {
	// bmath::vec2 p[3]{};
	PointMass* p[3];
	inline const bmath::vec2& p0() const { return p[0]->pos; }
	inline const bmath::vec2& p1() const { return p[1]->pos; }
	inline const bmath::vec2& p2() const { return p[2]->pos; }
	inline bmath::vec2& p0() { return p[0]->pos; }
	inline bmath::vec2& p1() { return p[1]->pos; }
	inline bmath::vec2& p2() { return p[2]->pos; }
};

// float calcArea(const Triangle& tri) {
float calcArea(const bmath::vec2& p0, const bmath::vec2& p1, const bmath::vec2& p2) {
	bmath::vec2 tmp(p1 - p0);
	std::swap(tmp[0], tmp[1]);
	tmp[1] *= -1;
	const bool positive = tmp.dot(p2 - p0) < 0;
	return (p1 - p0).cross(p2 - p0) * .5f;// * (positive ?: -1);
}

bmath::vec2 dirLinePoint(const bmath::vec2& l0, const bmath::vec2& l1, const bmath::vec2& p) {
	const bmath::vec2 ldir = l1 - l0; // direction l0 -> l1
	const bmath::vec2 lP = l0 + ldir * ldir.dot(p - l0); // point on line, that is closest to p
	const bmath::vec2 dir = p - (lP);
	return dir * calcArea(l0, l1, p) * 2 / dir.mag(); // rescale to match cross product magnitude
}

struct TriangleConstraint : Constraint {
private:
	Triangle target;
	const float area;
public:
	inline TriangleConstraint(Triangle target): target(target), area(calcArea(target.p0(), target.p1(), target.p2())) { }
	inline TriangleConstraint(Triangle target, const float area): target(target), area(area) { }
	virtual inline void solve(const float dt) {
		// const float invMass[3]{ 1.f, 1.f, 1.f };
		const float invMass[3] {
			1.f / target.p[0]->mass,
			1.f / target.p[1]->mass,
			1.f / target.p[2]->mass
		};

		const float currentArea = calcArea(target.p0(), target.p1(), target.p2());

		bmath::vec2 gradient0 = dirLinePoint(target.p1(), target.p2(), target.p0());
		bmath::vec2 gradient1 = dirLinePoint(target.p2(), target.p0(), target.p1());
		bmath::vec2 gradient2 = dirLinePoint(target.p0(), target.p1(), target.p2());

		// gradient0 *= currentArea * 2 / gradient0.mag();
		// gradient1 *= currentArea * 2 / gradient1.mag();
		// gradient2 *= currentArea * 2 / gradient2.mag();

		float w = 0.f;
		w += invMass[0] * gradient0.magSquared();
		w += invMass[1] * gradient1.magSquared();
		w += invMass[2] * gradient2.magSquared();

		const float compliance = 5.f;
		const float c = currentArea - area; // constraint value
		const float alpha = compliance / (dt * dt);
		const float lambda = -(2 * c) / (w + alpha);
		// target.p0() += gradient0 * s;
		std::cout << "Error: " << c << "  after: ";

		target.p[0]->pos += gradient0 * lambda * invMass[0];
		target.p[1]->pos += gradient1 * lambda * invMass[1];
		target.p[2]->pos += gradient2 * lambda * invMass[2];

		// target.p[0].prevPos += gradient0 * lambda * invMass[0];
		// target.p[1].prevPos += gradient1 * lambda * invMass[1];
		// target.p[2].prevPos += gradient2 * lambda * invMass[2];

		std::cout << (calcArea(target.p0(), target.p1(), target.p2()) - area) << "\n";
	}
};


struct Spring {
private:
	// bmath::vec2& p0;
	// bmath::vec2& p1;
	PointMass& p0;
	PointMass& p1;
	const float restLength;
	const float strength;
public:
	inline Spring(PointMass& p0, PointMass& p1, const float strength):
		p0(p0), p1(p1), restLength((p1.pos - p0.pos).mag()), strength(strength) {}
	// inline Spring(bmath::vec2& p0, bmath::vec2& p1, const float strength, const float restLength):
	// 	p0(p0), p1(p1), restLength(restLength), strength(strength) {}
	inline void update() {
		// const float invMass[2] { 1.f / m0, 1.f / m1};
		// const float invMass[2] { 1.f, 1.f};
		const bmath::vec2 p0p1 = p1.pos - p0.pos;
		const bmath::vec2 dir01 = p0p1 / p0p1.mag();
		const float currentLength = p0p1.mag();
		const float deviation = currentLength - restLength;
		float f = deviation * strength;
		p0.applyForce(f * dir01);
		p1.applyForce(-f * dir01);
	}
};


int main() {
	GDIWindow win(800, 800);

	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);

	LARGE_INTEGER lastTime;
	QueryPerformanceCounter(&lastTime);

	PointMass p0({3.f, 1.5f}, 1.f / 0.f);
	PointMass p1({6.f, 3.5f}, .2f);
	PointMass p2({2.f, 5.5f}, .2f);
	PointMass p3({4.f, 7.5f}, .2f);

	Triangle t0{ &p0, &p1, &p2 };
	Triangle t1{ &p2, &p1, &p3 };
	TriangleConstraint c0(t0);
	TriangleConstraint c1(t1);

	Spring s0(p0, p1, 50.f);
	Spring s1(p1, p2, 50.f);
	Spring s2(p2, p0, 50.f);

	Spring s3(p1, p3, 50.f);
	Spring s4(p2, p3, 50.f);

	// t0.p0()[0] += .3;

	float prevDt = .016f;
	float dt;
	while(!win.shouldClose()) {

		LARGE_INTEGER currentTime;
		QueryPerformanceCounter(&currentTime);

		// dt = (currentTime.QuadPart - lastTime.QuadPart) * 1. / freq.QuadPart;
		dt = 1.f / 120.f;
		lastTime = currentTime;

		constexpr size_t NUM_SUBSTEPS = 5;

		static float t = 0;
		t += dt / NUM_SUBSTEPS;
		t = fmodf(t, 1.f);


		for(size_t i = 0; i < NUM_SUBSTEPS; i++) {
			t += dt / NUM_SUBSTEPS;
			t = fmodf(t, 1.f);

			// t0.p0()[0] = 2.f + fabsf(2.f * t - 1.f);


			p0.calcVel(prevDt);
			p1.calcVel(prevDt);
			p2.calcVel(prevDt);
			p3.calcVel(prevDt);

			// std::cout << "dt: " << dt << "s\n";

			
			p0.applyForce({0.f, 9.81f});
			p1.applyForce({0.f, 9.81f});
			p2.applyForce({0.f, 9.81f});
			p3.applyForce({0.f, 9.81f});


			s0.update();
			s1.update();
			s2.update();
			s3.update();
			s4.update();


			p0.update(dt / NUM_SUBSTEPS);
			p1.update(dt / NUM_SUBSTEPS);
			p2.update(dt / NUM_SUBSTEPS);
			p3.update(dt / NUM_SUBSTEPS);


			c0.solve(dt / NUM_SUBSTEPS);
			c1.solve(dt / NUM_SUBSTEPS);

			prevDt = dt / NUM_SUBSTEPS;
		}


		win.graphics.clear(0xFFFFFFFF);

		// for(size_t i = 0; i < 3; i++)
		// 	win.graphics.line(
		// 		t0.p[i]->pos[0] * 100, t0.p[i]->pos[1] * 100,
		// 		t0.p[(i + 1) % 3]->pos[0] * 100, t0.p[(i + 1) % 3]->pos[1] * 100,
		// 		0xFF000000);
		win.graphics.line(p0.pos[0] * 100, p0.pos[1] * 100, p1.pos[0] * 100, p1.pos[1] * 100, 0xFF000000);
		win.graphics.line(p1.pos[0] * 100, p1.pos[1] * 100, p2.pos[0] * 100, p2.pos[1] * 100, 0xFF000000);
		win.graphics.line(p2.pos[0] * 100, p2.pos[1] * 100, p0.pos[0] * 100, p0.pos[1] * 100, 0xFF000000);
		win.graphics.line(p1.pos[0] * 100, p1.pos[1] * 100, p3.pos[0] * 100, p3.pos[1] * 100, 0xFF000000);
		win.graphics.line(p2.pos[0] * 100, p2.pos[1] * 100, p3.pos[0] * 100, p3.pos[1] * 100, 0xFF000000);

		for(size_t i = 0; i < 3; i++)
			win.graphics.fillCircle(t0.p[i]->pos[0] * 100, t0.p[i]->pos[1] * 100, 3, 0xFFFF0000);

		win.updateScreen();



		win.pollMsg();
	}

	return 0;
}