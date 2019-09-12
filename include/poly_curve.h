#pragma once

#include "geometrycentral/utilities/vector3.h"
#include "Eigen/Sparse"

namespace LWS {
    class PolyCurve;
    using namespace geometrycentral;

    struct PointOnCurve {
        int pIndex;
        PolyCurve* curve;
        bool operator ==(const PointOnCurve &other);
        bool operator !=(const PointOnCurve &other);
        Vector3 Position();
        void SetPosition(Vector3 pos);
        Vector3 Tangent();
        double DualLength();

        int GlobalIndex();

        PointOnCurve Next();
        PointOnCurve Prev();
    };

    struct IndexedMatrix {
        Eigen::SparseMatrix<double> M;
        int fineOffset;
        int coarseOffset;
    };

    enum class MultigridMode {
        MatrixOnly,
        Barycenter,
        EdgeLengths
    };

    class MultigridOperator {
        public:
        MultigridOperator();
        int lowerSize;
        int upperSize;
        std::vector<IndexedMatrix> matrices;
        Eigen::VectorXd mapUpward(Eigen::VectorXd v, MultigridMode mode);
        Eigen::VectorXd mapDownward(Eigen::VectorXd v, MultigridMode mode);
    };

    class PolyCurve {
        public:
        std::vector<Vector3> positions;
        PolyCurve(int numVerts);
        PolyCurve(std::vector<Vector3> ps);
        int index;
        int offset;

        int NumVertices();
        int NextVertex(int v);
        int PrevVertex(int v);
        PointOnCurve GetCurvePoint(int v);

        Vector3 Position(int v);
        void SetPosition(int v, Vector3 pos);

        Vector3 BisectorNormal(int index);
        Vector3 LengthWeightedNormal(int index);
        Vector3 EdgeNormal(int index);
        Vector3 VertexTangent(int index);
        double TotalLength();
        double DualLength(int i);
        Vector3 Barycenter();
        Vector3 TotalLengthGradient(int i);

        PolyCurve* Coarsen(Eigen::SparseMatrix<double> &prolongOp, Eigen::SparseMatrix<double> &sparsifyOp);
    };

    class PolyCurveGroup {
        public:
        std::vector<PolyCurve*> curves;
        void AddCurve(PolyCurve* c);
        int NumVertices();
        void BoundingCube(Vector3 &center, double &width);
        PointOnCurve GetCurvePoint(int v); 
        Vector3 Barycenter();
        double TotalLength();
        int NextIndexInCurve(int v);
        int GlobalIndex(PointOnCurve c);

        PolyCurveGroup* Coarsen(MultigridOperator &prolongOps, MultigridOperator &sparsifyOps);
    };
}