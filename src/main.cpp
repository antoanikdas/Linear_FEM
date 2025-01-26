#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "Node.h"
#include "Element.h"
#include "FEMSolver.h"

// Function to read nodes, elements, and boundary conditions from a file
bool readInputFile(const std::string& filename,
                   std::vector<Node>& nodes,
                   std::vector<Element>& elements,
                   std::vector<int>& fixedDofs,
                   std::vector<int>& prescribedDofs,
                   std::vector<double>& prescribedValues) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        return false;
    }

    std::string line, section;
    while (std::getline(inputFile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        if (line == "NODES") {
            section = "NODES";
            continue;
        } else if (line == "ELEMENTS") {
            section = "ELEMENTS";
            continue;
        } else if (line == "BOUNDARY_CONDITIONS") {
            section = "BOUNDARY_CONDITIONS";
            continue;
        }

        if (section == "NODES") {
            int id;
            double x, y;
            iss >> id >> x >> y;
            nodes.emplace_back(x, y, id - 1);  // Node IDs start from 0
        } else if (section == "ELEMENTS") {
            int id, n1, n2, n3, n4;
            iss >> id >> n1 >> n2 >> n3 >> n4;
            elements.emplace_back(
                std::vector<Node*>{&nodes[n1 - 1], &nodes[n2 - 1], &nodes[n3 - 1], &nodes[n4 - 1]},
                100.0, 0.25);
        } else if (section == "BOUNDARY_CONDITIONS") {
            int type, nodeID;
            char direction;
            double value;
            iss >> type >> nodeID >> direction >> value;

            int dof = (nodeID - 1) * 2 + (direction == 'x' ? 0 : 1);
            if (type == 0) {  // Fixed DOF
                fixedDofs.push_back(dof);
            } else if (type == 1) {  // Prescribed DOF
                prescribedDofs.push_back(dof);
                prescribedValues.push_back(value);
            }
        }
    }

    inputFile.close();
    return true;
}

int main() {
    // Data containers
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<int> fixedDofs, prescribedDofs;
    std::vector<double> prescribedValues;

    // Read input data from file
    if (!readInputFile("../data/input.dat", nodes, elements, fixedDofs, prescribedDofs, prescribedValues)) {
        return 1;
    }

    // Combine fixed and prescribed boundary conditions
    fixedDofs.insert(fixedDofs.end(), prescribedDofs.begin(), prescribedDofs.end());
    std::vector<double> prescribedAllValues(fixedDofs.size(), 0.0);
    for (size_t i = 0; i < prescribedDofs.size(); ++i) {
        prescribedAllValues[fixedDofs.size() - prescribedDofs.size() + i] = prescribedValues[i];
    }

    // Initialize solver
    FEMSolver solver(nodes, elements);

    // Assemble global stiffness matrix
    solver.assembleGlobalStiffnessMatrix();

    // Apply boundary conditions
    solver.applyBoundaryConditions(fixedDofs, prescribedAllValues);

    // Solve the system
    solver.solve();

    // Write results to a file
    std::ofstream outputFile("../data/results.dat");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open results.txt for writing.\n";
        return 1;
    }

    // Write displacements to the file
    outputFile << "Nodal Displacements:\n";
    const Eigen::VectorXd& displacements = solver.getDisplacements();
    for (size_t i = 0; i < displacements.size(); i += 2) {
        outputFile << "Node " << (i / 2) + 1 << ": u_x = " << displacements[i]
                   << ", u_y = " << displacements[i + 1] << "\n";
    }

    // Write strains and stresses for each element
    outputFile << "\nElement Strain and Stress at Integration Points:\n";
    for (size_t i = 0; i < elements.size(); ++i) {
        auto results = elements[i].computeStrainStressAtIntegrationPoints(solver.getDisplacements());
        outputFile << "Element " << i + 1 << ":\n";
        int gp = 1;
        for (const auto& [strain, stress] : results) {
            outputFile << "  Gauss Point " << gp++ << ":\n";
            outputFile << "    Strain: " << strain.transpose() << "\n";
            outputFile << "    Stress: " << stress.transpose() << "\n";
        }
    }

    outputFile.close();
    std::cout << "Results written to results.txt\n";

    return 0;
}
