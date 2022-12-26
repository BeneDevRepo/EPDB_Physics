#pragma once

#include "Constraint.hpp"
#include "PointMass.hpp"

#include "BWindow/GDIWindow.h"

extern GDIWindow* gwin;


inline double calcArea(const bmath::vec<2, double>& p0, const bmath::vec<2, double>& p1, const bmath::vec<2, double>& p2) {
	return (p1 - p0).cross(p2 - p0) * .5f;// * (positive ?: -1);
}

// inline bool positiveArea(const bmath::vec<2, double>& p0, const bmath::vec<2, double>& p1, const bmath::vec<2, double>& p2) {
// 	bmath::vec<2, double> tmp(p1 - p0);
// 	std::swap(tmp[0], tmp[1]);
// 	tmp[1] *= -1;
// 	const bool positive = tmp.dot(p2 - p0) < 0;
// 	return positive;
// }

inline bmath::vec<2, double> dirLinePoint(const bmath::vec<2, double>& l0, const bmath::vec<2, double>& l1, const bmath::vec<2, double>& p) {
	const bmath::vec<2, double> ldir = l1 - l0; // direction l0 -> l1
	const bmath::vec<2, double> l0p = p - l0;
	const float t = ldir.dot(l0p) / ldir.mag();
	const bmath::vec<2, double> lP = l0 + ldir / ldir.mag() * t; // point on line, that is closest to p
	const bmath::vec<2, double> dir = p - lP;
	// return dir / dir.mag() * calcArea(l0, l1, p) * 2; // rescale to match cross product magnitude
	return dir / dir.mag() * calcArea(l0, l1, p) * 2;// * (positiveArea(l0, l1, p) ?: -1); // rescale to match cross product magnitude
}


struct TriangleConstraint : Constraint {
private:
	const size_t p0, p1, p2;
	const double area;
public:
	// inline TriangleConstraint(Triangle target): target(target), area(calcArea(target.p0(), target.p1(), target.p2())) { }
	// inline TriangleConstraint(Triangle target, const double area): target(target), area(area) { }
	inline TriangleConstraint(const size_t p0, const size_t p1, const size_t p2):
		p0(p0), p1(p1), p2(p2),
		area(calcArea(points[p0].pos, points[p1].pos, points[p2].pos)) { }
	inline TriangleConstraint(const size_t p0, const size_t p1, const size_t p2, const double area):
		p0(p0), p1(p1), p2(p2),
		area(area) { }
	virtual inline void solve(const double dt) {
		// const double invMass[3]{ 1.f, 1.f, 1.f };
		const double invMass[3] {
			1.f / points[p0].mass,
			1.f / points[p1].mass,
			1.f / points[p2].mass
		};

		const double currentArea = calcArea(points[p0].pos, points[p1].pos, points[p2].pos);

		bmath::vec<2, double> gradient0 = dirLinePoint(points[p1].pos, points[p2].pos, points[p0].pos);
		bmath::vec<2, double> gradient1 = dirLinePoint(points[p2].pos, points[p0].pos, points[p1].pos);
		bmath::vec<2, double> gradient2 = dirLinePoint(points[p0].pos, points[p1].pos, points[p2].pos);

		// gwin->graphics.line(
		// 	points[p0].pos[0] * 100,
		// 	points[p0].pos[1] * 100,
		// 	points[p0].pos[0] * 100 + gradient0[0] * 10,
		// 	points[p0].pos[1] * 100 + gradient0[1] * 10,
		// 	0xFF00FF00
		// );
		// gwin->graphics.line(
		// 	points[p1].pos[0] * 100,
		// 	points[p1].pos[1] * 100,
		// 	points[p1].pos[0] * 100 + gradient1[0] * 10,
		// 	points[p1].pos[1] * 100 + gradient1[1] * 10,
		// 	0xFF00FF00
		// );
		// gwin->graphics.line(
		// 	points[p2].pos[0] * 100,
		// 	points[p2].pos[1] * 100,
		// 	points[p2].pos[0] * 100 + gradient2[0] * 10,
		// 	points[p2].pos[1] * 100 + gradient2[1] * 10,
		// 	0xFF00FF00
		// );

		// gradient0 *= currentArea * 2 / gradient0.mag();
		// gradient1 *= currentArea * 2 / gradient1.mag();
		// gradient2 *= currentArea * 2 / gradient2.mag();

		double w = 0.f;
		w += invMass[0] * gradient0.magSquared();
		w += invMass[1] * gradient1.magSquared();
		w += invMass[2] * gradient2.magSquared();

		const double compliance = 0.f;
		const double c = currentArea - area; // constraint value
		const double alpha = compliance / (dt * dt);
		const double lambda = -(2 * c) / (w + alpha);
		// target.p0() += gradient0 * s;

		// std::cout << "Error: " << c << "  after: ";
		// std::cout << "Area: " << currentArea << "  Error after: ";

		points[p0].pos += lambda * gradient0 * invMass[0];
		points[p1].pos += lambda * gradient1 * invMass[1];
		points[p2].pos += lambda * gradient2 * invMass[2];

		// std::cout << (calcArea(points[p0].pos, points[p1].pos, points[p2].pos) - area) << "\n";
	}
};

