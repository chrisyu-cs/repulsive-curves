#include "spatial/tpe_kdtree.h"
#include "utils.h"
#include "tpe_energy_sc.h"

#include <ctime>
#include <cstdlib> 

namespace LWS {

    inline double GetCoordFromBody(VertexBody6D body, int axis) {
        switch (axis) {
            case 0: return body.pt.position.x;
            case 1: return body.pt.position.y;
            case 2: return body.pt.position.z;
            case 3: return body.pt.tangent.x;
            case 4: return body.pt.tangent.y;
            case 5: return body.pt.tangent.z;
            default:
            std::cout << "Invalid axis passed to GetCoordFromBody: " << axis << std::endl;
            return body.pt.position.x;
        }
    }

    inline PosTan ReplaceCoord(PosTan body, int axis, double value) {
        switch (axis) {
            case 0: body.position.x = value; break;
            case 1: body.position.y = value; break;
            case 2: body.position.z = value; break;
            case 3: body.tangent.x = value; break;
            case 4: body.tangent.y = value; break;
            case 5: body.tangent.z = value; break;
            default:
            std::cout << "Invalid axis passed to ReplaceCoord: " << axis << std::endl; break;
        }
        return body;
    }

    inline int NextAxis(int axis) {
        switch (axis) {
            case 0: return 1;
            case 1: return 2;
            case 2: return 3;
            case 3: return 4;
            case 4: return 5;
            case 5: return 0;
            default: return 0;
        }
    }

    KDTreeNode3D::KDTreeNode3D(std::vector<VertexBody6D> &points, int axis, PosTan mins, PosTan maxs) {
        thresholdTheta = 0.25;
        children[0] = 0;
        children[1] = 0;

        minCoords = mins;
        maxCoords = maxs;
        splitAxis = axis;

        if (points.size() == 0) {
            isLeaf = true;
            isEmpty = true;
            totalMass = 0;
        }
        // If there's only one point, then we don't need to do anything
        // except set the fields from that one point
        else if (points.size() == 1) {
            body = points[0];
            isLeaf = true;
            isEmpty = false;
            totalMass = body.mass;
            centerOfMass = body.pt.position;
            averageTangent = body.pt.tangent;
        }
        // Otherwise we recursively split the set
        else {
            // Reserve space for splitting the points into lesser and greater
            int nPoints = points.size();
            std::vector<VertexBody6D> lesserPoints;
            lesserPoints.reserve(nPoints / 2 + 1);
            std::vector<VertexBody6D> greaterPoints;
            greaterPoints.reserve(nPoints / 2 + 1);

            // Compute the plane over which to split the points
            splitPoint = AxisSplittingPlane(points, axis, mins, maxs);

            // Split the points over the median
            for (int i = 0; i < nPoints; i++) {
                double coord = GetCoordFromBody(points[i], axis);

                if (coord <= splitPoint) {
                    lesserPoints.push_back(points[i]);
                }
                else {
                    greaterPoints.push_back(points[i]);
                }
            }

            // Compute the middle extents for the two halves
            PosTan raisedMin = ReplaceCoord(minCoords, axis, splitPoint);
            PosTan loweredMax = ReplaceCoord(maxCoords, axis, splitPoint);

            // Recursively construct children
            int nextAxis = NextAxis(axis);
            KDTreeNode3D* lesserNode = new KDTreeNode3D(lesserPoints, nextAxis, minCoords, loweredMax);
            KDTreeNode3D* greaterNode = new KDTreeNode3D(greaterPoints, nextAxis, raisedMin, maxCoords);

            children[0] = lesserNode;
            children[1] = greaterNode;

            combineValuesFromChildren();
            
            isLeaf = false;
            isEmpty = false;
        }
    }


    KDTreeNode3D::~KDTreeNode3D() {
        if (children[0]) delete children[0];
        if (children[1]) delete children[1];
    }

    double KDTreeNode3D::AxisSplittingPlane(std::vector<VertexBody6D> &points, int axis, PosTan mins, PosTan maxs) {
        int nPoints = points.size();
        int median_threshold = 99;
        int median_c = std::min(nPoints, median_threshold);
        std::vector<double> median_vals(median_c);

        // If there is a small enough number of points, we just take the true median
        if (nPoints <= median_threshold) {
            for (int i = 0; i < median_c; i++) {
                // Choose a random index and add its value
                median_vals[i] = GetCoordFromBody(points[i], axis);
            }
        }
        else {
            // Randomly subsample some number of points
            for (int i = 0; i < median_c; i++) {
                // Choose a random index and add its value
                int rand_i = rand() % nPoints;
                median_vals[i] = GetCoordFromBody(points[rand_i], axis);
            }
        }

        // Sort the subsampled points, and choose the median
        std::sort(median_vals.begin(), median_vals.end());

        double median_coord;
        int m;
        // In the even case, we need to average the two middle elements
        if (median_c % 2 == 0) {
            m = median_c / 2;
            median_coord = (median_vals[m] + median_vals[m-1]) / 2;
        }
        // In the odd case, we just pick the single median
        else {
            m = median_c / 2;
            median_coord = median_vals[m];
        }

        return median_coord;
    }

    void KDTreeNode3D::combineValuesFromChildren() {
        totalMass = children[0]->totalMass + children[1]->totalMass;
        // Compute average positions, phis, thetas
        centerOfMass = (children[0]->totalMass * children[0]->centerOfMass
            + children[1]->totalMass * children[1]->centerOfMass) / totalMass;
        averageTangent = (children[0]->totalMass * children[0]->averageTangent
            + children[1]->totalMass * children[1]->averageTangent) / totalMass;
    }

    void KDTreeNode3D::recomputeCentersOfMass(PolyCurveGroup* curves) {
        if (isEmpty) {
            totalMass = 0;
        }
        else if (isLeaf) {
            totalMass = body.mass;
            centerOfMass = body.pt.position;
            averageTangent = body.pt.tangent;
        }
        else {
            if (children[0]) children[0]->recomputeCentersOfMass(curves);
            if (children[1]) children[1]->recomputeCentersOfMass(curves);
            combineValuesFromChildren();
        }
    }

    inline double KDTreeNode3D::nodeRadius() {
        // Compute diagonal distance from corner to corner
        double diag = norm(maxCoords.position - minCoords.position);
        return diag;
    }

    bool KDTreeNode3D::shouldUseCell(Vector3 vertPos) {
        Vector3 toCenter = centerOfMass - vertPos;
        double d = norm(centerOfMass - vertPos);
        toCenter.normalize();

        // Compute max spread in the perpendicular direction

        double ratio = nodeRadius() / d;
        // TODO: take into account some tangent-related criteria?
        return ratio < thresholdTheta;
    }

    void KDTreeNode3D::accumulateVertexEnergy(double &result, PointOnCurve &i_pt,
    PolyCurveGroup* curves, double alpha, double beta) {
        if (isEmpty) {
            return;
        }
        else if (isLeaf) {
            // If this is a leaf, then it only has one vertex in it, so just use it
            PointOnCurve j = curves->GetCurvePoint(body.vertIndex1);
            // Don't sum if it's the same vertex
            if (j == i_pt) return;
            // Add contribution to energy at i from vertex j
            result += TPESC::tpe_pair(i_pt, j, alpha, beta);
        }
        else {
            if (shouldUseCell(i_pt.Position())) {
                Vector3 tangent = body.pt.tangent;
                // std::cout << "Using cell as body" << std::endl;
                tangent.normalize();
                // This cell is far enough away that we can treat it as a single body
                result += TPESC::tpe_pair_pts(i_pt.Position(), body.pt.position, tangent, i_pt.DualLength(), body.mass, alpha, beta);
            }
            else {
                // Otherwise we continue recursively traversing the tree
                for (int i = 0; i < 2; i++) {
                    if (children[i]) {
                        children[i]->accumulateVertexEnergy(result, i_pt, curves, alpha, beta);
                    }
                }
            }
        }
    }

    void KDTreeNode3D::accumulateTPEGradient(std::vector<Vector3> &gradients, PointOnCurve &i_pt,
    PolyCurveGroup* curves, double alpha, double beta) {
        if (isEmpty) {
            return;
        }
        else if (isLeaf) {
            // If this is a leaf, then it only has one vertex in it, so just use it
            PointOnCurve j = curves->GetCurvePoint(body.vertIndex1);
            // Don't sum if it's the same vertex
            if (j == i_pt) return;

            PointOnCurve i_prev = i_pt.Prev();
            PointOnCurve i_next = i_pt.Next();

            PointOnCurve j_prev = j.Prev();
            PointOnCurve j_next = j.Next();

            int iInd = curves->GlobalIndex(i_pt);
            int jInd = curves->GlobalIndex(j);

            // Differentiate both terms for previous, middle, and next
            gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(i_pt, j, alpha, beta, i_prev);
            gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(i_pt, j, alpha, beta, i_pt);
            gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(i_pt, j, alpha, beta, i_next);
            // Avoid double-counting on the reverse terms with these checks
            if (i_prev != j_prev && i_prev != j && i_prev != j_next) {
                gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(j, i_pt, alpha, beta, i_prev);
            }
            if (i_pt != j_prev && i_pt != j && i_pt != j_next) {
                gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(j, i_pt, alpha, beta, i_pt);
            }
            if (i_next != j_prev && i_next != j && i_next != j_next) {
                gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(j, i_pt, alpha, beta, i_next);
            }
        }
        else {
            if (shouldUseCell(i_pt.Position())) {
                Vector3 tangent = body.pt.tangent;
                tangent.normalize();
                // This cell is far enough away that we can treat it as a single body
                TangentMassPoint j{tangent, body.mass, body.pt.position, PointOnCurve{0, 0}, PointOnCurve{0, 0}};

                PointOnCurve i_prev = i_pt.Prev();
                PointOnCurve i_next = i_pt.Next();

                // Differentiate both terms for previous, middle, and next
                gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(i_pt, j, alpha, beta, i_prev);
                gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(j, i_pt, alpha, beta, i_prev);
                gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(i_pt, j, alpha, beta, i_pt);
                gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(j, i_pt, alpha, beta, i_pt);
                gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(i_pt, j, alpha, beta, i_next);
                gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(j, i_pt, alpha, beta, i_next);
            }
            else {
                // Otherwise we continue recursively traversing the tree
                for (int i = 0; i < 2; i++) {
                    if (children[i]) {
                        children[i]->accumulateTPEGradient(gradients, i_pt, curves, alpha, beta);
                    }
                }
            }
        }
    }


    inline void max_body(PosTan &old_max, PosTan next) {
        old_max.position = vector_max(old_max.position, next.position);
        old_max.tangent = vector_max(old_max.tangent, next.tangent);
    }
    
    inline void min_body(PosTan &old_min, PosTan next) {
        old_min.position = vector_min(old_min.position, next.position);
        old_min.tangent = vector_min(old_min.tangent, next.tangent);
    }

    inline VertexBody6D vertToBody(PolyCurveGroup* curves, PolyCurve* curve, int i) {
        PointOnCurve p = curve->GetCurvePoint(i);
        Vector3 pos = p.Position();
        Vector3 tangent = p.Tangent();

        PosTan ptan{pos, tangent};
        double mass = p.DualLength();
        int globalIndex = curves->GlobalIndex(p);

        return VertexBody6D{ptan, mass, globalIndex, -1};
    }

    KDTreeNode3D* CreateKDTreeFromCurve(PolyCurveGroup *curves) {
        std::vector<VertexBody6D> verts(curves->NumVertices());
        PosTan maxCoords = vertToBody(curves, curves->curves[0], 0).pt;
        PosTan minCoords = maxCoords;

        // Loop over all the vertices
        for (size_t i = 0; i < curves->curves.size(); i++) {
            PolyCurve* curve = curves->curves[i];
            size_t nVerts = curve->NumVertices();
            for (size_t v = 0; v < nVerts; v++) {
                VertexBody6D curBody = vertToBody(curves, curve, v);
                // Keep track of max/min extents
                max_body(maxCoords, curBody.pt);
                min_body(minCoords, curBody.pt);
                // Put vertex body into full list
                verts[curBody.vertIndex1] = curBody;
            }
        }

        std::cout << "Creating tree..." << std::endl;
        KDTreeNode3D* tree = new KDTreeNode3D(verts, 0, minCoords, maxCoords);
        std::cout << "Finished creating tree..." << std::endl;

        return tree;
    }
}

