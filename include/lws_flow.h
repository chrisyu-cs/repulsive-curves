#pragma once

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "lws_energy.h"

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/halfedge_element_types.h"
#include "geometrycentral/numerical/suitesparse_utilities.h"

#include <SuiteSparseQR.hpp>
#include <cholmod.h>

#include <Eigen/Sparse>

#include "vertexderivatives.h"

namespace LWS {

    // for some reason VSC says every occurrence of "cholmod_dense" is
    // "ambiguous" if this isn't here
    using cholmod_dense = ::cholmod_dense;
    void addTriplet(cholmod_triplet* M, int rowInd, int colInd, double value);

    typedef Eigen::Triplet<double> EigenTriplet;
    typedef Eigen::SparseMatrix<double> EigenSparse;

    class GradientSolver {
        public:
        GradientSolver();
        ~GradientSolver();

        void setSurface(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
        void setEnergy(LWSEnergyFunction* lwsFunc);
        void computeCotanLaplacian();
        void computeSaddleMatrix();
        void symbolicFactorization(cholmod_sparse* matrix);
        void numericFactorization(cholmod_sparse* matrix);

        EigenSparse CotanLaplacianEigen();

        Vector3 vertexGradient(GradientFunction gradFunc, surface::Vertex v);
        void computeVertexGradients(GradientFunction gradFunc);
        cholmod_dense* transformGradient(cholmod_dense* rhs);
        cholmod_dense* copyGradient(cholmod_dense* rhs);
        void projectOntoNormal(cholmod_dense* gradient);
        void backprojectConstraints(cholmod_factor* factor);
        void stepProjectedGradient(EnergyFunction energyFunc, GradientFunction gradFunc, double h);

        void FixedStep(cholmod_dense* gradient, double delta);
        bool LineSearchStep(EnergyFunction energy, cholmod_dense* gradient, double initGuess);

        cholmod_sparse* laplacian;
        cholmod_sparse* bilaplacian;
        cholmod_sparse* saddleMatrix;

        CholmodContext* GetContext();

        private:
        int doublingLimit;
        void printMatrix(cholmod_sparse* matrix, std::string name);
        void printMatrix(cholmod_dense* matrix, std::string name);
        double EvaluateGradientStep(EnergyFunction energy, cholmod_dense* gradient, double delta);

        surface::HalfedgeMesh* mesh;
        surface::VertexPositionGeometry* geom;
        surface::VertexData<Vector3>* originalPositions;
        cholmod_factor* symbolicFactor;
        LWSEnergyFunction* lwsFunc;
        
        cholmod_dense* vertexGradients;
        cholmod_dense* transformedGradients;
        cholmod_dense* workspaceY, *workspaceE;

        CholmodContext context;

        double coeffA, coeffB, coeffC;
    };

    void test_derivatives(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
    void test_quantities(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
    void test_dihedral(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom);
    void step_vertex_flow(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom, EnergyFunction energyFunc, GradientFunction gradFunc, double h);

    enum FlowType {
        LWS_AREA, LWS_MEAN_CURVATURE, LWS_GAUSS_GURVATURE
    };

}
