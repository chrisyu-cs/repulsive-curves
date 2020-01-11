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

    LengthDifferencePotential::LengthDifferencePotential(double wt) {
        weight = wt;
    }

    double LengthDifferencePotential::CurrentValue(PolyCurveNetwork* curves) {
        int nVerts = curves->NumVertices();
        double diffs2 = 0;
        for (int i = 0; i < nVerts; i++) {
            diffs2 += LenDiff(curves, i);
        }
        return 0.5 * weight * diffs2;
    }

    double LengthDifferencePotential::LenDiff(PolyCurveNetwork* curves, int i) {
        CurveVertex* v_i = curves->GetVertex(i);
        if (v_i->numEdges() != 2) {
            return 0;
        }
        CurveEdge* e1 = v_i->edge(0);
        CurveEdge* e2 = v_i->edge(1);
        double lenDiff = e2->Length() - e1->Length();
        
        return lenDiff * lenDiff;

    }

    void LengthDifferencePotential::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient) {
        int nVerts = curves->NumVertices();
        double len = curves->TotalLength();

        for (int i = 0; i < nVerts; i++) {
            // Each length difference centered at i affects vertices (i-1, i, i+1).
            double lDiff = LenDiff(curves, i);
            if (lDiff == 0) continue;

            CurveVertex* v_i = curves->GetVertex(i);
            CurveEdge* e1 = v_i->edge(0);
            CurveEdge* e2 = v_i->edge(1);
            CurveVertex* v_prev = e1->Opposite(v_i);
            CurveVertex* v_next = e2->Opposite(v_i);

            // Differentiate (length(e2) - length(e1)) by the 3 vertices
            Vector3 deriv2_next = TPESC::edge_length_wrt_vert(e2, v_next);
            Vector3 deriv2_mid = TPESC::edge_length_wrt_vert(e2, v_i);
            Vector3 deriv1_mid = -TPESC::edge_length_wrt_vert(e1, v_i);
            Vector3 deriv1_prev = -TPESC::edge_length_wrt_vert(e1, v_prev);
            
            double lenWeight = lDiff * weight;

            AddToRow(gradient, v_prev->GlobalIndex(), lenWeight * deriv1_prev);
            AddToRow(gradient, v_i->GlobalIndex(), lenWeight * (deriv1_mid + deriv2_mid));
            AddToRow(gradient, v_next->GlobalIndex(), lenWeight * deriv2_next);
        }
    }

    // AreaPotential::AreaPotential(double wt) {
    //     weight = wt;
    // }

    PinBendingPotential::PinBendingPotential(double wt) {
        weight = wt;
    }

    double PinBendingPotential::CurrentValue(PolyCurveNetwork* curves) {
        int nVerts = curves->NumVertices();
        double energy = 0;

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* v_i = curves->GetVertex(i);
            if (v_i->numEdges() != 2) continue;
            CurveEdge* e1 = v_i->edge(0);
            CurveEdge* e2 = v_i->edge(1);

            Vector3 vec1 = e1->Vector();
            Vector3 vec2 = e2->Vector();
            double pair = 1 - dot(vec1.normalize(), vec2.normalize());
            energy += pair * pair;
        }

        return 0.5 * weight * energy;
    }

    void PinBendingPotential::AddGradient(PolyCurveNetwork* curves, Eigen::MatrixXd &gradient) {
        int nVerts = curves->NumVertices();

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* v_i = curves->GetVertex(i);
            if (v_i->numEdges() != 2) continue;
            CurveEdge* e1 = v_i->edge(0);
            CurveEdge* e2 = v_i->edge(1);

            Vector3 vec1 = e1->Vector();
            Vector3 vec2 = e2->Vector();
            double pair = 1 - dot(vec1.normalize(), vec2.normalize());
            
            CurveVertex* v_prev = e1->Opposite(v_i);
            CurveVertex* v_next = e2->Opposite(v_i);

            // We want to differentiate dot(vec1_hat, vec2_hat)
            VertJacobian vec1_deriv_prev = TPESC::edge_tangent_wrt_vert(e1, v_prev);
            VertJacobian vec1_deriv_i = TPESC::edge_tangent_wrt_vert(e1, v_i);
            VertJacobian vec2_deriv_i = TPESC::edge_tangent_wrt_vert(e2, v_i);
            VertJacobian vec2_deriv_next = TPESC::edge_tangent_wrt_vert(e2, v_next);

            Vector3 deriv_prev = vec1_deriv_prev.LeftMultiply(vec2);
            Vector3 deriv_i = vec1_deriv_i.LeftMultiply(vec2) + vec2_deriv_i.LeftMultiply(vec1);
            Vector3 deriv_next = vec2_deriv_next.LeftMultiply(vec1);

            double pairWeight = pair * weight;
            AddToRow(gradient, v_prev->GlobalIndex(), pairWeight * deriv_prev);
            AddToRow(gradient, v_i->GlobalIndex(), pairWeight * deriv_i);
            AddToRow(gradient, v_next->GlobalIndex(), pairWeight * deriv_next);
        }
    }



}