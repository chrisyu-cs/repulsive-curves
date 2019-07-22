#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "polyscope/surface_mesh.h"

#include "geometrycentral/numerical/suitesparse_utilities.h"
#include <SuiteSparseQR.hpp>
#include <cholmod.h>

#include "boundary_element.h"

#include <Eigen/Core>
#include <Eigen/Dense>

namespace LWS {

    using cholmod_dense = ::cholmod_dense;
    using namespace geometrycentral;

    class LWSBoundaryLoop {
        public:
        surface::HalfedgeMesh* mesh;
        surface::VertexPositionGeometry* geom;

        LWSBoundaryLoop(surface::HalfedgeMesh* m, surface::VertexPositionGeometry* g);
        surface::Vertex GetVertex(int index);
        BoundaryElement* GetElement(int index);
        int NumVertices();
        void ComputeGradients(double areaCoeff, double bLengthCoeff);
        void StepGradients(double h);
        bool isValid();
        void BoundaryLaplacian();
        void BoundarySaddle();
        void symbolicFactorization(cholmod_sparse* matrix);
        void numericFactorization(cholmod_sparse* matrix);
        void StepBoundaryGradient(double areaCoeff, double lenCoeff, double h);
        void ReparameterizeArc();
        void CreateNodesAndElements();

        void PrintVertices2D();
        
        size_t GlobalIndex(int index);
        int PrevVertex(int index);
        int NextVertex(int index);
        Vector3 Position(int vert);
        Vector3 BisectorNormal(int index);
        Vector3 LengthWeightedNormal(int index);
        Vector3 EdgeNormal(int index);
        Vector3 VertexTangent(int index);
        double EdgeLength(int index);
        double LoopCurvature(int index);

        void ScaleByHalf();

        void MoveAlongNormal(double h);

        void CreateBemMatrix();
        void LaplaceTestRHS();
        void AllOnesRHS();
        void SolveV();
        void MoveByV(double h);

        double ValueAtInterior(Vector3 point);

        Eigen::VectorXd boundary_u;
        Eigen::VectorXd boundary_v;
        cholmod_dense* GetFullRHS(CholmodContext *c);
        double DualLength(int x);
        double TotalLength();
        Vector3 Barycenter();

        private:
        bool valid;
        surface::BoundaryLoop boundary;
        CholmodContext context;
        std::vector<surface::Vertex> vertices;
        std::vector<size_t> loopToGlobal;
        std::vector<double> interiorAngles;

        cholmod_dense* gradients;
        cholmod_sparse* boundary_bilaplacian;
        cholmod_sparse* boundary_saddle;
        cholmod_dense* workspaceY, *workspaceE;
        cholmod_factor* symbolicFactor;
        cholmod_dense* fullRHS;

        Eigen::MatrixXd bemMatrixA;
        Eigen::MatrixXd bemMatrixB;
        Eigen::VectorXd rhs_vec;

        std::vector<BoundaryNode*> bem_nodes;
        std::vector<BoundaryElement*> bem_elements; 

        double totalLength;
        std::vector<double> arcLength;
    };

}