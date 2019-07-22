#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "polyscope/surface_mesh.h"

namespace LWS {

    using namespace geometrycentral;

    double fs_greens_2d(Vector3 center, Vector3 pt);
    Vector3 fs_greens_grad_2d(Vector3 center, Vector3 pt);

    class BoundaryElement;

    class BoundaryNode {
        public:
        BoundaryNode(int id, Vector3 pos, Vector3 nor);
        BoundaryElement* element1;
        BoundaryElement* element2;
        double DualLength();
        // Nodal values for element 1
        double u;
        double v;
        Vector3 position;
        Vector3 normal;
        int nodeID;
    };

    class BoundaryElement {
        public:
        static double gauss_abscissa[4];
        static double gauss_weights[4];

        BoundaryElement(int id, BoundaryNode *n1, BoundaryNode *n2);
        BoundaryNode* node1;
        BoundaryNode* node2;
        double CoeffIntA_ia(BoundaryNode* collocation_pt, BoundaryNode* node);
        double CoeffIntB_ia(BoundaryNode* collocation_pt, BoundaryNode* node);
        
        double CoeffIntA_ia(Vector3 collocation_pt, BoundaryNode* node);
        double CoeffIntB_ia(Vector3 collocation_pt, BoundaryNode* node);

        double length();
        int elementID;

        private:
        Vector3 position(double t);
        Vector3 normal(double t);
        double basis_fn(BoundaryNode* node, double t);

        double integrand_A(Vector3 collocation_pt, BoundaryNode* node, double t);
        double integrand_A_11(Vector3 collocation_pt, BoundaryNode* node, double t);
        double integrand_B(Vector3 collocation_pt, BoundaryNode* node, double t);
        double integrand_B_11(Vector3 collocation_pt, BoundaryNode* node, double t);
    };

}