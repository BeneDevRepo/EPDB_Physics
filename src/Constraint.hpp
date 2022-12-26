#pragma once


struct Constraint {
	virtual ~Constraint() = default;
	virtual inline void solve(const double dt) = 0;
};