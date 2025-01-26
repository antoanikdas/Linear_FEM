#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Node.h"
#include "Element.h"

class FEMSolver {
private:
    std::vector<Node> nodes;        // All nodes
    std::vector<Element> elements;  // All elements
    Eigen::MatrixXd globalK;        // Global stiffness matrix
    Eigen::VectorXd globalF;        // Global force vector
    Eigen::VectorXd displacements;  // Displacement vector

public:
    FEMSolver(const std::vector<Node>& nodes, const std::vector<Element>& elements);

    void assembleGlobalStiffnessMatrix();
    void applyBoundaryConditions(const std::vector<int>& fixedDofs, const std::vector<double>& prescribedValues);
    void solve();

    void printResults();

    // Getter for displacements
    const Eigen::VectorXd& getDisplacements() const { return displacements; }
};
