#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <vector>
// #include <thread>
// #include <mutex>

#include "BWindow/GDIWindow.h"
#include "BMath/vector.h"

#include "PointMass.hpp"
#include "Spring.hpp"
#include "LineConstraint.hpp"
#include "TriangleConstraint.hpp"

#include "windows.h"


GDIWindow* gwin;

std::vector<PointMass> points;
std::vector<LineConstraint> lines;
std::vector<Spring> springs;
std::vector<TriangleConstraint> triangles;



int main() {
	GDIWindow win(800, 800);
	gwin = &win;

	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);

	LARGE_INTEGER lastTime;
	QueryPerformanceCounter(&lastTime);


	bmath::vec<2, double> origin; // pixel coordinates of physics origin
	double scale = 100; // pixels per meter




	points.push_back(PointMass({3.f, 1.5f}, 1. / 0.));
	size_t base = 1;

	const double mass = 1.f;
	// points.push_back(PointMass({3.f, 2.5f}, 1. / 0.));
	points.push_back(PointMass({3.f, 2.5f}, mass));
	points.push_back(PointMass({6.f, 3.5f}, mass));
	points.push_back(PointMass({2.f, 5.5f}, mass));

	points.push_back(PointMass({4.f, 7.5f}, mass));


	triangles.push_back(TriangleConstraint(base + 0, base + 1, base + 2));
	triangles.push_back(TriangleConstraint(base + 2, base + 1, base + 3));

	springs.push_back(Spring(0, 1, 550.f));
	// lines.push_back(LineConstraint(0, 1));

	const double strength = 5.f;
	// springs.push_back(Spring(base + 0, base + 1, strength));
	// springs.push_back(Spring(base + 1, base + 2, strength));
	// springs.push_back(Spring(base + 2, base + 0, strength));

	// springs.push_back(Spring(base + 1, base + 3, strength));
	// springs.push_back(Spring(base + 2, base + 3, strength));

	const double restLength = 3.f;
	springs.push_back(Spring(base + 0, base + 1, strength, restLength));
	// springs.push_back(Spring(base + 1, base + 2, strength, restLength));
	springs.push_back(Spring(base + 2, base + 0, strength, restLength));

	springs.push_back(Spring(base + 1, base + 3, strength, restLength));
	springs.push_back(Spring(base + 2, base + 3, strength, restLength));


	int pmouseX, pmouseY;
	int mouseX = win.win.mouseX, mouseY = win.win.mouseY;

	int pscroll;
	int scroll = win.win.scroll;


	constexpr size_t NUM_SUBSTEPS = 100;

	
	while(!win.shouldClose()) {
		pmouseX = mouseX;
		pmouseY = mouseY;
		mouseX = win.win.mouseX;
		mouseY = win.win.mouseY;

		pscroll = scroll;
		scroll = win.win.scroll;

		if(GetAsyncKeyState(VK_MBUTTON) & 0x8000) {
			origin[0] += mouseX - pmouseX;
			origin[1] += mouseY - pmouseY;
		}
		if((scroll != pscroll) && GetAsyncKeyState(VK_CONTROL) & 0x8000) {
			const double diff = scroll - pscroll;
			const double scaleFac = (diff > 0) ? (1.1) : (.9);
			scale *= scaleFac;
			if(scale < 5) scale = 5;
			// if(scroll != pscroll)
			// std::cout << "Scroll\n";
		}
		if(GetAsyncKeyState(VK_LBUTTON) & 0x8000) {
			points[0].pos[0] = (std::max<double>(0., std::min<double>(mouseX, win.width))  - origin[0]) / scale;
			points[0].pos[1] = (std::max<double>(0., std::min<double>(mouseY, win.height)) - origin[1]) / scale;
			// points[0].prevPos[0] = points[0].pos[0];
			// points[0].prevPos[1] = points[0].pos[1];
		}
		points[0].vel *= 0;
		points[0].force *= 0;

		win.graphics.clear(0xFFFFFFFF);

		LARGE_INTEGER currentTime;
		QueryPerformanceCounter(&currentTime);

		// dt = (currentTime.QuadPart - lastTime.QuadPart) * 1. / freq.QuadPart;
		double dt = 1.f / 120.f;
		lastTime = currentTime;

		static double t = .5;

		for(size_t i = 0; i < NUM_SUBSTEPS; i++) {
			t += dt * .2f / NUM_SUBSTEPS;
			t = fmodf(t, 1.f);

			// points[0].pos[0] = 2.f + fabsf(2.f * t - 1.f);
			// points[0].pos[0] = 2.f + (sinf(2.f * 3.1415926535f * t) + .5f) * 2.f;


			// for(PointMass& p : points)
			// 	p.calcVel(prevDt);

			// std::cout << "dt: " << dt << "s\n";
			

			// accumulate forces:
			constexpr float gravity = 9.81f * .5f;
			for(PointMass& p : points)
				if(p.mass != 1.f / 0.f) // dont accelerate points with infinite mass
					p.applyForce({0.f, gravity * p.mass});

			for(PointMass& p : points)
				if(p.mass < 1.f / 0.f) // dont accelerate points with infinite mass
					p.applyForce(-p.vel * .15f);

			for(Spring& s : springs)
				s.update();


			// apply velocity:
			for(PointMass& p : points)
				p.update(dt / NUM_SUBSTEPS);


			// constraints:
			for(TriangleConstraint& t : triangles)
				t.solve(dt / NUM_SUBSTEPS);
			for(LineConstraint& l : lines)
				l.solve(dt / NUM_SUBSTEPS);

			// std::cout << "Correcting p0: " << points[0].pos[0] << ", " << points[0].pos[1];
			for(PointMass& p : points) {
				// p.pos[0] = std::max<double>(0., std::min<double>(p.pos[0], win.width / 100.));
				// p.pos[1] = std::max<double>(0., std::min<double>(p.pos[1], win.height / 100.));

				p.pos[1] = std::min<double>(p.pos[1], 8.);
			}
			// std::cout << "  After: " << points[0].pos[0] << ", " << points[0].pos[1] << "\n";
			


			// calculate velocity:
			for(PointMass& p : points)
				p.calcVel(dt / NUM_SUBSTEPS);

			// prevDt = dt / NUM_SUBSTEPS;
		}


		// for(size_t i = 0; i < 3; i++)
		// 	win.graphics.line(
		// 		t0.p[i]->pos[0] * 100, t0.p[i]->pos[1] * 100,
		// 		t0.p[(i + 1) % 3]->pos[0] * 100, t0.p[(i + 1) % 3]->pos[1] * 100,
		// 		0xFF000000);
		// win.graphics.line(p0.pos[0] * 100, p0.pos[1] * 100, p1.pos[0] * 100, p1.pos[1] * 100, 0xFF000000);
		// win.graphics.line(p1.pos[0] * 100, p1.pos[1] * 100, p2.pos[0] * 100, p2.pos[1] * 100, 0xFF000000);
		// win.graphics.line(p2.pos[0] * 100, p2.pos[1] * 100, p0.pos[0] * 100, p0.pos[1] * 100, 0xFF000000);
		// win.graphics.line(p1.pos[0] * 100, p1.pos[1] * 100, p3.pos[0] * 100, p3.pos[1] * 100, 0xFF000000);
		// win.graphics.line(p2.pos[0] * 100, p2.pos[1] * 100, p3.pos[0] * 100, p3.pos[1] * 100, 0xFF000000);
	

		for(const Spring& s : springs)
			win.graphics.line(
				origin[0] + points[s.p0].pos[0] * scale,
				origin[1] + points[s.p0].pos[1] * scale,
				origin[0] + points[s.p1].pos[0] * scale,
				origin[1] + points[s.p1].pos[1] * scale,
				0xFF000000);

		for(const PointMass& p : points)
			win.graphics.fillCircle(
				origin[0] + p.pos[0] * scale,
				origin[1] + p.pos[1] * scale,
				3,
				0xFFFF0000);

		// horizontal origin:
		win.graphics.line(
			0,
			origin[1],
			win.width,
			origin[1],
			0xFFFF0000);

		// vertical origin:
		win.graphics.line(
			origin[0],
			0,
			origin[0],
			win.height,
			0xFF00FF00);

		win.graphics.line(
			0,
			origin[1] + 8 * scale,
			win.width,
			origin[1] + 8 * scale,
			0xFF0000FF);

		win.updateScreen();



		win.pollMsg();
	}

	return 0;
}