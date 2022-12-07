#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cmath>
// #include <thread>
// #include <mutex>

#include "BWindow/GDIWindow.h"
#include "BMath/vector.h"

#include "windows.h"


struct Constraint {
	virtual ~Constraint() = default;
	virtual inline void solve(const float dt) = 0;
};

struct Triangle {
	bmath::vec2 p[3]{};
	inline const bmath::vec2& p0() const { return p[0]; }
	inline const bmath::vec2& p1() const { return p[1]; }
	inline const bmath::vec2& p2() const { return p[2]; }
	inline bmath::vec2& p0() { return p[0]; }
	inline bmath::vec2& p1() { return p[1]; }
	inline bmath::vec2& p2() { return p[2]; }
};

float calcArea(const Triangle& tri) {
	return (tri.p1() - tri.p0()).cross(tri.p2() - tri.p0()) * .5f;
}

bmath::vec2 dirLinePoint(const bmath::vec2& l0, const bmath::vec2& l1, const bmath::vec2& p) {
	const bmath::vec2 ldir = l1 - l0;
	return p - (l0 + ldir * ldir.dot(p - l0));
}

struct TriangleConstraint : Constraint {
private:
	Triangle &target;
	const float area;
public:
	inline TriangleConstraint(Triangle& target): target(target), area(calcArea(target)) { }
	inline TriangleConstraint(Triangle& target, const float area): target(target), area(area) { }
	virtual inline void solve(const float dt) {
		// const float invMass[3]{ 1.f / 3.f, 1.f / 3.f, 1.f / 3.f };
		// const float invMass[3]{ 1.f, 1.f, 1.f };
		const float invMass[3]{ 0.f, 1.f, 1.f };

		const bmath::vec2 gradient0 = dirLinePoint(target.p0(), target.p1(), target.p2());
		const bmath::vec2 gradient1 = dirLinePoint(target.p1(), target.p2(), target.p0());
		const bmath::vec2 gradient2 = dirLinePoint(target.p2(), target.p0(), target.p1());

		float w = 0.f;
		w += invMass[0] * gradient0.magSquared();
		w += invMass[1] * gradient1.magSquared();
		w += invMass[2] * gradient2.magSquared();

		const float compliance = 0.f;
		const float c = (calcArea(target) - area); // constraint value
		const float alpha = compliance / (dt * dt);
		const float lambda = -(3 * c) / (w + alpha);
		// target.p0() += gradient0 * s;
		std::cout << "Error: " << c << "  after: ";

		target.p0() += gradient0 * lambda * invMass[0];
		target.p1() += gradient1 * lambda * invMass[1];
		target.p2() += gradient2 * lambda * invMass[2];

		std::cout << (calcArea(target) - area) << "\n";
	}
};


int main() {
	GDIWindow win(800, 800);

	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq); // TODO

	LARGE_INTEGER lastTime;
	QueryPerformanceCounter(&lastTime);

	Triangle t0 {
		bmath::vec2{ 1.f, 3.5f },
		bmath::vec2{ 5.f, 1.5f },
		bmath::vec2{ 5.f, 5.5f }
	};

	TriangleConstraint c0(t0);

	// t0.p0()[0] += .3;

	// float density = 1000.f;
	while(!win.shouldClose()) {
		LARGE_INTEGER currentTime;
		QueryPerformanceCounter(&currentTime);

		const double dt = (currentTime.QuadPart - lastTime.QuadPart) * 1. / freq.QuadPart;
		lastTime = currentTime;

		// std::cout << "dt: " << dt << "s\n";

		static float t = 0;
		t += dt;
		t = fmodf(t, 1.f);
		// while(t > 1)
		// 	t -= 1;

		// t0.p0()[0] = 2.f + 2.f * t;
		t0.p0()[0] = 2.f + fabsf(2.f * t - 1.f);

		c0.solve(dt);


		win.graphics.clear(0xFFFFFFFF);

		for(size_t i = 0; i < 3; i++)
			win.graphics.line(
				t0.p[i][0] * 100, t0.p[i][1] * 100,
				t0.p[(i + 1) % 3][0] * 100, t0.p[(i + 1) % 3][1] * 100,
				0xFF000000);

		for(size_t i = 0; i < 3; i++)
			win.graphics.fillCircle(t0.p[i][0] * 100, t0.p[i][1] * 100, 3, 0xFFFF0000);

		win.updateScreen();



		win.pollMsg();
	}

	return 0;
}