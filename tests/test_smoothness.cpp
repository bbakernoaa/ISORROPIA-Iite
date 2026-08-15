#include <gtest/gtest.h>
#include "Isorropia/Solver.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

TEST(IsorropiaSmoothnessTest, Case1DifferentiabilitySweep) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    input.temp = 298.15;
    input.rh = 0.80; // standard Relative Humidity
    input.iprob = 0; // Forward solver

    double sulfate = 1e-6;
    input.w[1] = sulfate; // Sulfate

    double h = 1e-5;
    std::vector<double> sulrats;
    std::vector<double> derivatives;

    // Sweep across sulrat from 0.80 to 2.20
    for (double sulrat = 0.80; sulrat <= 2.20; sulrat += 0.002) {
        input.w[2] = (sulrat - h) * sulfate; // Ammonia (left side)
        Isorropia::State state_left;
        solver.solve(input, state_left);

        input.w[2] = (sulrat + h) * sulfate; // Ammonia (right side)
        Isorropia::State state_right;
        solver.solve(input, state_right);

        double dWater_dAmmonia = (state_right.water - state_left.water) / (2.0 * h * sulfate);

        EXPECT_FALSE(std::isnan(dWater_dAmmonia)) << "NaN derivative at sulrat = " << sulrat;
        EXPECT_FALSE(std::isinf(dWater_dAmmonia)) << "Inf derivative at sulrat = " << sulrat;
        
        sulrats.push_back(sulrat);
        derivatives.push_back(dWater_dAmmonia);
    }

    // Verify derivative properties:
    // 1. Boundedness
    for (size_t i = 0; i < derivatives.size(); ++i) {
        EXPECT_LT(std::abs(derivatives[i]), 1e4) << "Derivative blew up at sulrat = " << sulrats[i] << " with value " << derivatives[i];
    }

    // 2. Smoothness/Differentiability (no sudden large jumps in derivative)
    for (size_t i = 1; i < derivatives.size(); ++i) {
        double diff = std::abs(derivatives[i] - derivatives[i-1]);
        EXPECT_LT(diff, 500.0) << "Discontinuity/jump in derivative found between sulrat = " 
                               << sulrats[i-1] << " and " << sulrats[i] 
                               << " (diff = " << diff << ", left = " << derivatives[i-1] 
                               << ", right = " << derivatives[i] << ")";
    }
}

TEST(IsorropiaSmoothnessTest, Case2DifferentiabilitySweep) {
    Isorropia::Solver solver;
    Isorropia::Input input;
    input.temp = 298.15;
    input.rh = 0.80; // standard Relative Humidity
    input.iprob = 0; // Forward solver

    double sulfate = 1e-6;
    input.w[1] = sulfate; // Sulfate
    input.w[3] = 0.5 * sulfate; // Nitrate present

    double h = 1e-5;
    std::vector<double> sulrats;
    std::vector<double> derivatives;

    // Sweep across sulrat from 0.80 to 2.20
    for (double sulrat = 0.80; sulrat <= 2.20; sulrat += 0.002) {
        input.w[2] = (sulrat - h) * sulfate; // Ammonia (left side)
        Isorropia::State state_left;
        solver.solve(input, state_left);

        input.w[2] = (sulrat + h) * sulfate; // Ammonia (right side)
        Isorropia::State state_right;
        solver.solve(input, state_right);

        double dWater_dAmmonia = (state_right.water - state_left.water) / (2.0 * h * sulfate);

        EXPECT_FALSE(std::isnan(dWater_dAmmonia)) << "NaN derivative at sulrat = " << sulrat;
        EXPECT_FALSE(std::isinf(dWater_dAmmonia)) << "Inf derivative at sulrat = " << sulrat;
        
        sulrats.push_back(sulrat);
        derivatives.push_back(dWater_dAmmonia);
    }

    // Verify derivative properties:
    // 1. Boundedness
    for (size_t i = 0; i < derivatives.size(); ++i) {
        EXPECT_LT(std::abs(derivatives[i]), 1e4) << "Derivative blew up at sulrat = " << sulrats[i] << " with value " << derivatives[i];
    }

    // 2. Smoothness/Differentiability
    for (size_t i = 1; i < derivatives.size(); ++i) {
        double diff = std::abs(derivatives[i] - derivatives[i-1]);
        EXPECT_LT(diff, 500.0) << "Discontinuity/jump in derivative found between sulrat = " 
                               << sulrats[i-1] << " and " << sulrats[i] 
                               << " (diff = " << diff << ", left = " << derivatives[i-1] 
                               << ", right = " << derivatives[i] << ")";
    }
}
