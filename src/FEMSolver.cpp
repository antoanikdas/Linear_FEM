#include "FEMSolver.h"
#include <iostream>

FEMSolver::FEMSolver(const std::vector<Node>& nodes, const std::vector<Element>& elements)
    : nodes(nodes), elements(elements) {
    int dof = nodes.size() * 2;  // Each node has 2 DoFs (x and y)
    globalK = Eigen::MatrixXd::Zero(dof, dof);
    globalF = Eigen::VectorXd::Zero(dof);
    displacements = Eigen::VectorXd::Zero(dof);
}

void FEMSolver::assembleGlobalStiffnessMatrix() {
    for (const auto& element : elements) {
        Eigen::MatrixXd localK = element.calculateStiffnessMatrix();

        // Get the element nodes
        const std::vector<Node*>& elementNodes = element.getNodes();
        std::vector<int> dofs = {
            elementNodes[0]->id * 2, elementNodes[0]->id * 2 + 1,
            elementNodes[1]->id * 2, elementNodes[1]->id * 2 + 1,
            elementNodes[2]->id * 2, elementNodes[2]->id * 2 + 1,
            elementNodes[3]->id * 2, elementNodes[3]->id * 2 + 1
        };

        // Assemble the local stiffness matrix into the global stiffness matrix
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                globalK(dofs[i], dofs[j]) += localK(i, j);
            }
        }
    }
}

void FEMSolver::applyBoundaryConditions(const std::vector<int>& fixedDofs, const std::vector<double>& prescribedValues) {
    for (size_t i = 0; i < fixedDofs.size(); ++i) {
        int dof = fixedDofs[i];
        double value = prescribedValues[i];

        // Adjust the force vector for nonzero boundary conditions
        for (int j = 0; j < globalK.rows(); ++j) {
            if (j != dof) {
                globalF[j] -= globalK(j, dof) * value;
            }
        }

        // Apply boundary conditions by modifying the stiffness matrix
        for (int j = 0; j < globalK.rows(); ++j) {
            globalK(dof, j) = 0.0;
            globalK(j, dof) = 0.0;
        }
        globalK(dof, dof) = 1.0;
        globalF[dof] = value;
    }
}

void FEMSolver::solve() {
    displacements = globalK.colPivHouseholderQr().solve(globalF);
}

void FEMSolver::printResults() {
    std::cout << "Displacements:" << std::endl;
    std::cout << displacements << std::endl;
}
