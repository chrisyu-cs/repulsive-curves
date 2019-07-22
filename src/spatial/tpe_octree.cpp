#include "spatial/tpe_octree.h"
#include "utils.h"
#include "tpe_energy_sc.h"

namespace LWS {

    OctreeNode3D* CreateOctreeFromCurve(PolyCurveGroup* curves) {

        double bbox_width;
        Vector3 bbox_center;
        curves->BoundingCube(bbox_center, bbox_width);

        OctreeNode3D* root = new OctreeNode3D(bbox_center, bbox_width);

        for (size_t i = 0; i < curves->curves.size(); i++) {
            PolyCurve* curve = curves->curves[i];
            size_t nVerts = curve->NumVertices();
            for (size_t v = 0; v < nVerts; v++) {
                int globalIndex = curve->offset + v;
                VertexBody body{curve->Position(v), curve->VertexTangent(v), curve->DualLength(v), globalIndex};
                root->AddToSubtree(body);
            }
        }

        return root;
    }

    OctreeNode3D::OctreeNode3D(Vector3 center, double width) {
        // Set children to null in the beginning
        for (int i = 0; i < 8; i++) {
            children[i] = 0;
        }
        // Initialize empty node
        totalMass = 0;
        centerOfMass = Vector3{0, 0, 0};
        isLeaf = true;
        body.vertIndex = -1;
        // Record center and width
        nodeCenter = center;
        nodeWidth = width;

        theta = 0.25;
    }

    OctreeNode3D::~OctreeNode3D() {
        // Delete all of the children
        for (int i = 0; i < 8; i++) {
            if (children[i]) {
                delete children[i];
            }
        }
    }

    inline bool OctreeNode3D::isEmpty() {
        return body.vertIndex == -1;
    }

    void OctreeNode3D::updateCenterOfMass(Vector3 pos, Vector3 tangent, double mass) {
        // Because tangent signs don't make a difference, we could
        // make it so that all tangents point in the forward-x direction.
        // if (tangent.x < 0) tangent *= -1;

        Vector3 weightedPos = totalMass * centerOfMass + mass * pos;
        Vector3 weightedT = totalMass * averageTangent + mass * tangent;
        totalMass += mass;
        centerOfMass = weightedPos / totalMass;
        averageTangent = weightedT / totalMass;
    }

    void OctreeNode3D::recomputeCentersOfMass(PolyCurveGroup* curves) {
        centerOfMass = Vector3{0, 0, 0};
        averageTangent = Vector3{0, 0, 0};
        totalMass = 0;

        // On an internal node, recursively get centers of mass from all children
        if (!isLeaf) {
            for (int i = 0; i < 8; i++) {
                if (children[i]) {
                    children[i]->recomputeCentersOfMass(curves);
                    updateCenterOfMass(children[i]->centerOfMass, children[i]->averageTangent, children[i]->totalMass);
                }
            }
        }
        // On a leaf node, just set center to the single vertex inside
        else {
            PointOnCurve c = curves->GetCurvePoint(body.vertIndex);
            body.mass = c.DualLength();
            body.position = c.Position();
            body.tangent = c.Tangent();

            updateCenterOfMass(body.position, body.tangent, body.mass);
        }
    }

    int OctreeNode3D::findOctant(Vector3 position) {
        bool xLess = (position.x < nodeCenter.x);
        bool yLess = (position.y < nodeCenter.y);
        bool zLess = (position.z < nodeCenter.z);

        if (xLess) {
            if (yLess) {
                if (zLess) return 0; // x, y, z all less
                else return 1; // x, y less; z more
            }
            else {
                if (zLess) return 2; // x less; y more; z less
                else return 3; // x less; y more, z more
            }
        }
        else {
            if (yLess) {
                if (zLess) return 4; // x more; y, z less
                else return 5; // x more; y less; z more
            }
            else {
                if (zLess) return 6; // x, y more; z less
                else return 7; // x, y, z all more
            }
        }
    }

    Vector3 OctreeNode3D::getChildCenter(int octant) {
        double qtr = nodeWidth / 4;
        switch (octant) {
            case 0:
            return nodeCenter + Vector3{-qtr, -qtr, -qtr};
            case 1:
            return nodeCenter + Vector3{-qtr, -qtr, qtr};
            case 2:
            return nodeCenter + Vector3{-qtr, qtr, -qtr};
            case 3:
            return nodeCenter + Vector3{-qtr, qtr, qtr};
            case 4:
            return nodeCenter + Vector3{qtr, -qtr, -qtr};
            case 5:
            return nodeCenter + Vector3{qtr, -qtr, qtr};
            case 6:
            return nodeCenter + Vector3{qtr, qtr, -qtr};
            case 7:
            return nodeCenter + Vector3{qtr, qtr, qtr};
            default:
            std::cout << "Did not get a valid octant" << std::endl;
            return nodeCenter;
        }
    }

    void OctreeNode3D::createChildIfNonexistent(int childIndex) {
        // If the child hasn't been created, then create it
        if (!children[childIndex]) {
            Vector3 childCenter = getChildCenter(childIndex);
            double childWidth = nodeWidth / 2;
            children[childIndex] = new OctreeNode3D(childCenter, childWidth);
        }
    }

    void OctreeNode3D::AddToSubtree(VertexBody b) {
        Vector3 cdist = b.position - nodeCenter;
        cdist =Vector3{fabs(cdist.x), fabs(cdist.y), fabs(cdist.z)};
        if (cdist.x > nodeWidth / 2 || cdist.y > nodeWidth / 2 || cdist.z > nodeWidth / 2) {
            std::cout << "Position " << b.position << " is not in bounding box of node" << std::endl;
            std::cout << "  centered at " << nodeCenter << ", width " << nodeWidth << std::endl;
            return;
        }

        // If this node is empty, add the current body to it
        if (isEmpty()) {
            body = b;
        }
        // Otherwise, if this node is internal, recursively add to the appropriate quadrant
        else {
            if (!isLeaf) {
                // Update center of mass of this internal node
                updateCenterOfMass(b.position, b.tangent, b.mass);
                int childIndex = findOctant(b.position);
                // Create child node if necessary
                createChildIfNonexistent(childIndex);
                // Now that it's been created, add the child
                children[childIndex]->AddToSubtree(b);
            }
            // Otherwise the node is a leaf, and also has an existing body already in it
            else {
                // Whatever happens, the node won't be a leaf anymore
                isLeaf = false;
                // Since the other body was already added earlier, we only need
                // to account for the new b's center of mass.
                updateCenterOfMass(b.position, b.tangent, b.mass);
                // Get child indices
                int bIndex = findOctant(b.position);
                int cIndex = findOctant(body.position);
                // Create children if necessary
                createChildIfNonexistent(bIndex);
                createChildIfNonexistent(cIndex);
                // Recursively insert both children
                children[bIndex]->AddToSubtree(b);
                children[cIndex]->AddToSubtree(body);
            }
        }
    }

    void printTriple(PolyCurveGroup* curves, PointOnCurve c1, PointOnCurve c2, PointOnCurve wrt, bool reverse) {
        if ((curves->GlobalIndex(c1) == 2 && curves->GlobalIndex(c2) == 3) ||
            (curves->GlobalIndex(c1) == 3 && curves->GlobalIndex(c2) == 2)) {
            std::string rev = (reverse) ? "reverse" : "forward";
            std::cout << "Added (" << curves->GlobalIndex(c1) << ", "
                << curves->GlobalIndex(c2) << ") wrt " << curves->GlobalIndex(wrt) << " from " << rev << std::endl;
        }
    }

    void OctreeNode3D::accumulateVertexEnergy(double &result, PointOnCurve &i_pt,
    PolyCurveGroup* curves, double alpha, double beta) {
        if (isLeaf) {
            // If this is a leaf, then it only has one vertex in it, so just use it
            PointOnCurve j = curves->GetCurvePoint(body.vertIndex);
            // Don't sum if it's the same vertex
            if (j == i_pt) return;
            // Add contribution to energy at i from vertex j
            result += TPESC::tpe_pair(i_pt, j, alpha, beta);
        }
        else {
            double d = norm(centerOfMass - i_pt.Position());
            double ratio = nodeWidth / d;

            if (ratio < theta) {
                Vector3 tangent = body.tangent;
                // std::cout << "Using cell as body" << std::endl;
                tangent.normalize();
                // This cell is far enough away that we can treat it as a single body
                TangentMassPoint j{tangent, body.mass, body.position, PointOnCurve{0, 0}, PointOnCurve{0, 0}};
                // Add contribution to energy at i from body j
                result += TPESC::tpe_pair_pts(i_pt.Position(), j.point, i_pt.Tangent(), i_pt.DualLength(), j.mass, alpha, beta);
            }
            else {
                // Otherwise we continue recursively traversing the tree
                for (int i = 0; i < 8; i++) {
                    if (children[i]) {
                        children[i]->accumulateVertexEnergy(result, i_pt, curves, alpha, beta);
                    }
                }
            }
        }
    }

    void OctreeNode3D::accumulateTPEGradient(std::vector<Vector3> &gradients, PointOnCurve &i_pt,
    PolyCurveGroup* curves, double alpha, double beta) {
        if (isLeaf) {
            // If this is a leaf, then it only has one vertex in it, so just use it
            PointOnCurve j = curves->GetCurvePoint(body.vertIndex);
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
            double d = norm(centerOfMass - i_pt.Position());
            double ratio = nodeWidth / d;

            if (ratio < theta) {
                Vector3 tangent = body.tangent;
                // std::cout << "Using cell as body" << std::endl;
                tangent.normalize();
                // This cell is far enough away that we can treat it as a single body
                TangentMassPoint j{tangent, body.mass, body.position, PointOnCurve{0, 0}, PointOnCurve{0, 0}};

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
                for (int i = 0; i < 8; i++) {
                    if (children[i]) {
                        children[i]->accumulateTPEGradient(gradients, i_pt, curves, alpha, beta);
                    }
                }
            }
        }
    }

}