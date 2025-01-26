#include "Element.h"
#include <cmath>
#include <iostream>

Eigen::MatrixXd Element::calculateStiffnessMatrix() const {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(8, 8);  // 8x8 stiffness matrix

    // Plane strain material matrix
    double coeff = E / ((1 + nu) * (1 - 2 * nu));
    Eigen::MatrixXd D(3, 3);
    D << coeff * (1 - nu), coeff * nu, 0,
         coeff * nu, coeff * (1 - nu), 0,
         0, 0, coeff * (1 - 2 * nu) / 2;

    // Gauss points and weights
    std::vector<double> gauss = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    std::vector<double> weights = {1.0, 1.0};

    for (double xi : gauss) {
        for (double eta : gauss) {
            // Compute Jacobian and shape function derivatives
            Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2, 2);
            Eigen::MatrixXd dN_dxi(4, 2);
            dN_dxi << -(1 - eta) / 4.0, -(1 - xi) / 4.0,
                       (1 - eta) / 4.0, -(1 + xi) / 4.0,
                       (1 + eta) / 4.0,  (1 + xi) / 4.0,
                      -(1 + eta) / 4.0,  (1 - xi) / 4.0;

            // Compute Jacobian matrix
            for (int i = 0; i < 4; ++i) {
                J(0, 0) += dN_dxi(i, 0) * nodes[i]->x;
                J(0, 1) += dN_dxi(i, 1) * nodes[i]->x;
                J(1, 0) += dN_dxi(i, 0) * nodes[i]->y;
                J(1, 1) += dN_dxi(i, 1) * nodes[i]->y;
            }

            double detJ = J.determinant();
            Eigen::MatrixXd invJ = J.inverse();

            // Transform shape function derivatives to physical coordinates
            Eigen::MatrixXd dN_dx = dN_dxi * invJ;

            // Assemble B matrix
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);
            for (int i = 0; i < 4; ++i) {
                B(0, 2 * i)     = dN_dx(i, 0);
                B(1, 2 * i + 1) = dN_dx(i, 1);
                B(2, 2 * i)     = dN_dx(i, 1);
                B(2, 2 * i + 1) = dN_dx(i, 0);
            }

            // Add contribution to stiffness matrix
            K += B.transpose() * D * B * detJ * weights[0] * weights[1];
        }
    }

    return K;
}

std::vector<std::pair<Eigen::VectorXd, Eigen::VectorXd>> Element::computeStrainStressAtIntegrationPoints(
    const Eigen::VectorXd& globalDisplacements) const {
    // Plane strain material matrix
    double coeff = E / ((1 + nu) * (1 - 2 * nu));
    Eigen::MatrixXd D(3, 3);
    D << coeff * (1 - nu), coeff * nu, 0,
         coeff * nu, coeff * (1 - nu), 0,
         0, 0, coeff * (1 - 2 * nu) / 2;

    // Extract local displacements
    Eigen::VectorXd localDisplacements(8);
    for (int i = 0; i < 4; ++i) {
        localDisplacements(2 * i)     = globalDisplacements[2 * nodes[i]->id];
        localDisplacements(2 * i + 1) = globalDisplacements[2 * nodes[i]->id + 1];
    }

    // Gauss points
    std::vector<double> gauss = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};

    // Results at each integration point
    std::vector<std::pair<Eigen::VectorXd, Eigen::VectorXd>> results;

    for (double xi : gauss) {
        for (double eta : gauss) {
            // Compute Jacobian and shape function derivatives
            Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2, 2);
            Eigen::MatrixXd dN_dxi(4, 2);
            dN_dxi << -(1 - eta) / 4.0, -(1 - xi) / 4.0,
                       (1 - eta) / 4.0, -(1 + xi) / 4.0,
                       (1 + eta) / 4.0,  (1 + xi) / 4.0,
                      -(1 + eta) / 4.0,  (1 - xi) / 4.0;

            for (int i = 0; i < 4; ++i) {
                J(0, 0) += dN_dxi(i, 0) * nodes[i]->x;
                J(0, 1) += dN_dxi(i, 1) * nodes[i]->x;
                J(1, 0) += dN_dxi(i, 0) * nodes[i]->y;
                J(1, 1) += dN_dxi(i, 1) * nodes[i]->y;
            }

            double detJ = J.determinant();
            Eigen::MatrixXd invJ = J.inverse();
            Eigen::MatrixXd dN_dx = dN_dxi * invJ;

            // Assemble B matrix
            Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 8);
            for (int i = 0; i < 4; ++i) {
                B(0, 2 * i)     = dN_dx(i, 0);
                B(1, 2 * i + 1) = dN_dx(i, 1);
                B(2, 2 * i)     = dN_dx(i, 1);
                B(2, 2 * i + 1) = dN_dx(i, 0);
            }

            // Compute strain and stress
            Eigen::VectorXd strain = B * localDisplacements;
            Eigen::VectorXd stress = D * strain;

            // Store results
            results.emplace_back(strain, stress);
        }
    }

    return results;
}