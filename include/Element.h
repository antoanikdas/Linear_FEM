#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Node.h"

class Element {
private:
    std::vector<Node*> nodes;  // Pointer to nodes of this element
    double E, nu;              // Material properties

public:
    Element(const std::vector<Node*>& nodes, double E, double nu)
        : nodes(nodes), E(E), nu(nu) {}

    Eigen::MatrixXd calculateStiffnessMatrix() const;  // Marked as const
    const std::vector<Node*>& getNodes() const { return nodes; }  // Getter for nodes
    
    // Compute strain and stress at integration points
    std::vector<std::pair<Eigen::VectorXd, Eigen::VectorXd>> computeStrainStressAtIntegrationPoints(
        const Eigen::VectorXd& globalDisplacements) const;
};
