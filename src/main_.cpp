#include <iostream>
#include <vector>
#include "Node.h"
#include "Element.h"
#include "FEMSolver.h"

int main() {
    // Define nodes
    // std::vector<Node> nodes = {
    //     Node(0.0, 0.0, 0),  // Node 1
    //     Node(1.0, 0.0, 1),  // Node 2
    //     Node(1.0, 1.0, 2),  // Node 3
    //     Node(0.0, 1.0, 3)   // Node 4
    // };

    // Define elements
    // std::vector<Element> elements = {
    //     Element({&nodes[0], &nodes[1], &nodes[2], &nodes[3]}, 100.0, 0.25)
    // };

    // Define nodes
    std::vector<Node> nodes = {
        Node(0.0, 0.0, 0),  // Node 1
        Node(0.5, 0.0, 1),  // Node 2
        Node(1.0, 0.0, 2),  // Node 3
        Node(0.0, 0.5, 3),  // Node 4
        Node(0.5, 0.5, 4),  // Node 5
        Node(1.0, 0.5, 5),  // Node 6
        Node(0.0, 1.0, 6),  // Node 7
        Node(0.5, 1.0, 7),  // Node 8
        Node(1.0, 1.0, 8)   // Node 9
    };

    // Define elements (4 Quad4 elements)
    std::vector<Element> elements = {
        // Bottom-left element
        Element({&nodes[0], &nodes[1], &nodes[4], &nodes[3]}, 100.0, 0.25),
        // Bottom-right element
        Element({&nodes[1], &nodes[2], &nodes[5], &nodes[4]}, 100.0, 0.25),
        // Top-left element
        Element({&nodes[3], &nodes[4], &nodes[7], &nodes[6]}, 100.0, 0.25),
        // Top-right element
        Element({&nodes[4], &nodes[5], &nodes[8], &nodes[7]}, 100.0, 0.25)
    };

    // Initialize solver
    FEMSolver solver(nodes, elements);

    // Assemble global stiffness matrix
    solver.assembleGlobalStiffnessMatrix();

    // Apply boundary conditions
    std::vector<int> fixedDofs = {0, 1, 6, 12, 4, 10, 16};             // Fixed DoFs
    std::vector<double> prescribedValues = {0.0, 0.0, 0.0, 0.0, 0.01, 0.01, 0.01};  // Zero displacements
    
    
    solver.applyBoundaryConditions(fixedDofs, prescribedValues);

    // Solve the system
    solver.solve();

    // Print results
    solver.printResults();

    // Compute strain and stress for each element
    for (const auto& element : elements) {
    auto results = element.computeStrainStressAtIntegrationPoints(solver.getDisplacements());
    int gp = 1;
    for (const auto& [strain, stress] : results) {
        std::cout << "Gauss Point " << gp++ << ":\n";
        std::cout << "  Strain: " << strain.transpose() << "\n";
        std::cout << "  Stress: " << stress.transpose() << "\n";
    }
}

    return 0;
}
