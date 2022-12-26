#pragma once

#include "Constraint.hpp"
#include "PointMass.hpp"


struct LineConstraint : Constraint {
private:
	const size_t p0, p1;
	const double length;
public:
	inline LineConstraint(const size_t p0, const size_t p1):
		p0(p0), p1(p1),
		length((points[p0].pos - points[p1].pos).mag()) { }
	inline LineConstraint(const size_t p0, const size_t p1, const double length):
		p0(p0), p1(p1),
		length(length) { }
	virtual inline void solve(const double dt) {
		// const double invMass[3]{ 1.f, 1.f, 1.f };
		const double invMass[2] {
			1.f / points[p0].mass,
			1.f / points[p1].mass,
			// 1.f / points[p2].mass
		};

		const double currentLength = (points[p0].pos - points[p1].pos).mag();

		bmath::vec<2, double> gradient0 = (points[p0].pos - points[p1].pos) / currentLength;
		bmath::vec<2, double> gradient1 = -gradient0;


		// gradient0 *= currentArea * 2 / gradient0.mag();
		// gradient1 *= currentArea * 2 / gradient1.mag();
		// gradient2 *= currentArea * 2 / gradient2.mag();

		double w = 0.f;
		w += invMass[0] * gradient0.magSquared();
		w += invMass[1] * gradient1.magSquared();
		// w += invMass[2] * gradient2.magSquared();

		const double compliance = 0.f;
		const double c = currentLength - length; // constraint value
		const double alpha = compliance / (dt * dt);
		const double lambda = -c / (w + alpha);
		// target.p0() += gradient0 * s;
		// std::cout << "Error: " << c << "  after: ";

		points[p0].pos += lambda * gradient0 * invMass[0];
		points[p1].pos += lambda * gradient1 * invMass[1];
		// points[p2].pos += lambda * gradient2 * invMass[2];

		// target.p[0].prevPos += gradient0 * lambda * invMass[0];
		// target.p[1].prevPos += gradient1 * lambda * invMass[1];
		// target.p[2].prevPos += gradient2 * lambda * invMass[2];

		// std::cout << (calcArea(points[p0].pos, points[p1].pos, points[p2].pos) - area) << "\n";
	}
};

