#include "extra_potentials.h"

#include "tpe_energy_sc.h"

namespace LWS {
    // Base class doesn't do anything
    CurvePotential::CurvePotential() {}
    CurvePotential::~CurvePotential() {}
    double CurvePotential::CurrentValue(PolyCurveNetwork* curves) {
        return 0;
    }
    void CurvePotential::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient) {}

    TotalLengthPotential::TotalLengthPotential(double wt) {
        weight = wt;
    }

    double TotalLengthPotential::CurrentValue(PolyCurveNetwork* curves) {
        double len = curves->TotalLength();
        return weight * 0.5 * len * len;
    }

    void TotalLengthPotential::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient) {
        int nVerts = curves->NumVertices();
        double len = curves->TotalLength();

        Eigen::MatrixXd lenGradient;
        lenGradient.setZero(gradient.rows(), gradient.cols());

        for (int i = 0; i < nVerts; i++) {
            // Every vertex length needs to be differentiated by its neighbors
            CurveVertex* v_i = curves->GetVertex(i);
            int nNeighbors = v_i->numEdges();
            for (int j = 0; j < nNeighbors; j++) {
                CurveVertex* v_neighbor = curves->GetVertex(j);
                Vector3 lenGrad = TPESC::length_wrt_vert(v_i, v_neighbor);
                AddToRow(lenGradient, v_neighbor->GlobalIndex(), weight * lenGrad);
            }
            // Length also needs to be differentiated by itself
            Vector3 selfGrad = TPESC::length_wrt_vert(v_i, v_i);
            AddToRow(lenGradient, v_i->GlobalIndex(), weight * selfGrad);
        }

        gradient += (weight * len) * lenGradient;
    }

    // AreaPotential::AreaPotential(double wt) {
    //     weight = wt;
    // }

    


}