#pragma once

class Node {
public:
    double x, y;  // Coordinates
    int id;       // Node ID
    Node(double x, double y, int id) : x(x), y(y), id(id) {}
};
