#pragma once

#include "BMath/vector.h"

#include "PointMass.hpp"


struct Spring {
// private:
	// bmath::vec<2, double>& p0;
	// bmath::vec<2, double>& p1;
	// PointMass& p0;
	// PointMass& p1;
	size_t p0, p1;
	double restLength;
	double strength;
public:
	// inline Spring(PointMass& p0, PointMass& p1, const double strength):
	// 	p0(p0), p1(p1), restLength((p1.pos - p0.pos).mag()), strength(strength) {}
	inline Spring(const size_t p0, const size_t p1, const double strength):
		p0(p0), p1(p1), restLength((points[p1].pos - points[p0].pos).mag()), strength(strength) {}
	inline Spring(const size_t p0, const size_t p1, const double strength, const double restLength):
		p0(p0), p1(p1), restLength(restLength), strength(strength) {}
	// inline Spring(bmath::vec<2, double>& p0, bmath::vec<2, double>& p1, const double strength, const double restLength):
	// 	p0(p0), p1(p1), restLength(restLength), strength(strength) {}
	inline void update() {
		// const double invMass[2] { 1.f / m0, 1.f / m1};
		// const double invMass[2] { 1.f, 1.f};
		const bmath::vec<2, double> p0p1 = points[p1].pos - points[p0].pos;
		const bmath::vec<2, double> dir01 = p0p1 / p0p1.mag();
		const double currentLength = p0p1.mag();
		const double deviation = currentLength - restLength;
		double f = deviation * strength;
		points[p0].applyForce(f * dir01);
		points[p1].applyForce(-f * dir01);
	}
};