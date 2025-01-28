## **FEM Quad4 Solver**

This repository contains a finite element method (FEM) solver for 2D Quad4 elements implemented in **C++**. The program calculates displacements, strains, and stresses for a given mesh of 4-node quadrilateral elements under plane strain conditions.

---

### **Project Structure**

```
FEM_Quad4/
├── CMakeLists.txt       # Build configuration
├── src/                 # Source files
│   ├── main.cpp         
│   ├── Node.cpp         
│   ├── Element.cpp      
│   ├── FEMSolver.cpp    
├── include/             # Header files
│   ├── Node.h           
│   ├── Element.h        
│   ├── FEMSolver.h      
├── data/                 # Data files
│   ├── input.dat            
|   ├── results.dat          
└── README.md            # Project description
```

---

### **Input File Format**

The input file (`input.dat`) specifies the nodes, elements, and boundary conditions. Example:

```
# Nodes (ID, x, y)
NODES
1 0.0 0.0
2 0.5 0.0
3 1.0 0.0
4 0.0 0.5
5 0.5 0.5
6 1.0 0.5
7 0.0 1.0
8 0.5 1.0
9 1.0 1.0

# Elements (ID, Node1, Node2, Node3, Node4)
ELEMENTS
1 1 2 5 4
2 2 3 6 5
3 4 5 8 7
4 5 6 9 8

# Boundary Conditions (Type, NodeID, Direction, Value)
# Type: 0 = Fixed, 1 = Prescribed Displacement
BOUNDARY_CONDITIONS
0 1 x 0.0
0 1 y 0.0
0 4 x 0.0
0 7 x 0.0
1 3 x 0.01
1 6 x 0.01
1 9 x 0.01
```

---


### **Build Instructions**

#### Prerequisites

- **CMake** (version 3.10 or higher)
- **C++17** compiler
- **Eigen3** library (used for matrix operations)

#### Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/antoanikdas/Linear_FEM.git
   cd Linear_FEM
   ```
   Update the **Eigen3** submodule
   ```bash
   git submodule update --init
   ```
   Alternatively, 
   Clone Eigen3 library or add path in CMakeLists.txt
   ```bash
   git clone https://gitlab.com/libeigen/eigen.git
   ```
2. Build the project:
   ```bash
   mkdir build
   cd build
   cmake ..
   cmake --build .
   ```

3. Run the program:
   ```bash
   ./FEM_Quad4
   ```

---

### **Example Workflow**

1. Create an `input.dat` file with the mesh and boundary conditions (example provided above).
2. Run the solver.
3. Check `results.dat` for nodal displacements, strains, and stresses.

---

### **Contributing**

Contributions are welcome! If you encounter any bugs or have suggestions for improvements, feel free to open an issue or submit a pull request.

---

### **Acknowledgments**

This project uses the [Eigen3](https://eigen.tuxfamily.org/) library for matrix operations.
