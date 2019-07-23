#include "vertexderivatives.h"
#include "utils.h"

using namespace geometrycentral;
using namespace std;

namespace LWS {

    // Computes vertex dual area as 1/3 surrounding triangle areas
    double VertexArea(surface::VertexPositionGeometry* geom, surface::Vertex vert) {
        double sum = 0;
        for (surface::Face f : vert.adjacentFaces()) {
            sum += geom->faceArea(f) / 3;
        }
        return sum;
    }

    // Computes mean curvature as sum of 1 / 4 (dihedral angle * edge length)
    double VertexMeanCurvature(surface::VertexPositionGeometry* geom, surface::Vertex vert) {
        if (vert.isBoundary()) return 0;

        double sum = 0;
        for (surface::Edge edge : vert.adjacentEdges()) {
            sum += (geom->edgeLength(edge) * DihedralAngle(geom, edge));
        }
        return sum / 4;
    }

    // Computes Gauss curvature as (2pi - angle sum)
    double VertexGaussCurvature(surface::VertexPositionGeometry* geom, surface::Vertex vert) {
        if (vert.isBoundary()) return 0;

        double sum = 0;
        for (surface::Corner c : vert.adjacentCorners()) {
            sum += geom->cornerAngle(c);
        }

        return 2 * M_PI - sum;
    }

    // Compute A - 2H + K
    double VertexTubularEnergy(surface::VertexPositionGeometry* geom, surface::Vertex vert) {
        return VertexArea(geom, vert) - 2 * VertexMeanCurvature(geom, vert) + VertexGaussCurvature(geom, vert);
    }

    // Sums E(v) over all vertices v in the given mesh.
    double MeshEnergy(surface::VertexPositionGeometry* geom, surface::HalfedgeMesh* mesh, EnergyFunction energy) {
        double sumEnergy = 0;
        for (surface::Vertex v : mesh->vertices()) {
            double e = energy(geom, v);
            sumEnergy += e;
        }
        return sumEnergy;
    }

    Vector3 GradientArea(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other) {
        // If differentiating wrt self, then need to add up all
        // derivatives of surrounding triangle areas.
        if (base == other) {
            Vector3 gradient = Vector3{0, 0, 0};

            for (const surface::Halfedge &he : base.outgoingHalfedges()) {
                // Find the direction pointing away from the opposite edge
                surface::Halfedge opposite = he.next();
                // Don't consider normals of boundary loop "faces"
                if (!he.face().isBoundaryLoop()) {
                    Vector3 oppVec = heVector(geom, he);
                    Vector3 faceNormal = geom->faceNormal(he.face());
                    Vector3 perp = cross(faceNormal, oppVec);
                    perp = perp.normalize();
                    // The gradient is in that direction, and its magnitude is the opposing length
                    gradient += norm(oppVec) * perp;
                }
            }

            return gradient / 6;
        }

        // Otherwise we have to find if the other vertex is a neighbor
        else {
            surface::Halfedge adjacentHE;
            bool found = false;

            for (const surface::Halfedge &neighborHE : base.outgoingHalfedges()) {
                // Find the half-edge between this vertex and the other one
                if (other == neighborHE.twin().vertex()) {
                    found = true;
                    adjacentHE = neighborHE;
                    break;
                }
            }

            if (!found) return Vector3{0, 0, 0};

            // Assuming we found it, add gradients from the two opposite triangles
            Vector3 gradient = Vector3{0, 0, 0};
            surface::Face f1 = adjacentHE.face();
            surface::Face f2 = adjacentHE.twin().face();

            if (!f1.isBoundaryLoop()) {
                // Find the opposite edge from the "other" vertex
                surface::Halfedge opp1 = adjacentHE.next().next();
                Vector3 oppVec1 = heVector(geom, opp1);
                Vector3 perp1 = cross(geom->faceNormal(f1), oppVec1);
                perp1 = perp1.normalize();

                gradient += norm(oppVec1) * perp1;
            }
            if (!f2.isBoundaryLoop()) {
                // Find the opposite edge from the "other" vertex
                surface::Halfedge opp2 = adjacentHE.twin().next();
                Vector3 oppVec2 = heVector(geom, opp2);
                Vector3 perp2 = cross(geom->faceNormal(f2), oppVec2);
                perp2 = perp2.normalize();

                gradient += norm(oppVec2) * perp2;
            }
            return gradient / 6;
        }
    }

    double DihedralAngle(surface::VertexPositionGeometry* geom, surface::Edge base) {
        // TODO: check sign convention of dihedral angle in base library
        return -geom->edgeDihedralAngles[base];
    }

    Vector3 GradientDihedralAngle(surface::VertexPositionGeometry* geom, surface::Edge base, surface::Vertex other) {
        if (base.isBoundary()) return Vector3{0, 0, 0};

        double length = geom->edgeLength(base);
        surface::Halfedge he = base.halfedge();

        // If the vertex is one of the endpoints of this edge, evaluate the "hard" case
        if (he.vertex() == other) {
            // Case x0 in bending energies paper
            double cot03 = geom->halfedgeCotanWeight(he.next().next());
            double cot04 = geom->halfedgeCotanWeight(he.twin().next());

            Vector3 n0 = geom->faceNormal(he.face());
            Vector3 n1 = geom->faceNormal(he.twin().face());

            return (1. / length) * (cot03 * n0 + cot04 * n1);
        }

        else if (he.twin().vertex() == other) {
            // Case x1
            double cot01 = geom->halfedgeCotanWeight(he.next());
            double cot02 = geom->halfedgeCotanWeight(he.twin().next().next());

            Vector3 n0 = geom->faceNormal(he.face());
            Vector3 n1 = geom->faceNormal(he.twin().face());

            return (1. / length) * (cot01 * n0 + cot02 * n1);
        }

        // Otherwise, if the vertex is not on the edge but is opposite
        // one of the adjacent faces, use the "easy" formula
        else {
            // Check if the vertex is in the first face
            surface::Face face1 = he.face();
            for (const surface::Vertex &faceVert : face1.adjacentVertices()) {
                if (faceVert == other) {
                    double coeff = -length / (2 * geom->faceArea(face1));
                    return coeff * geom->faceNormal(face1);
                }
            }
            // Check if the vertex is in the second face
            surface::Face face2 = he.twin().face();
            for (const surface::Vertex &faceVert : face2.adjacentVertices()) {
                if (faceVert == other) {
                    double coeff = -length / (2 * geom->faceArea(face2));
                    return coeff * geom->faceNormal(face2);
                }
            }

            // If the vertex isn't adjacent to the edge or an adjacent face,
            // then it doesn't affect the dihedral angle.
            return Vector3{0, 0, 0};
        }
    }

    Vector3 GradientDihedralNumerical(surface::VertexPositionGeometry* geom, surface::Edge base, surface::Vertex other, double h) {
        Vector3 origPos = geom->vertexPositions[other];
        double origAngle = DihedralAngle(geom, base);

        // Differentiate wrt x
        geom->vertexPositions[other] = Vector3{origPos.x + h, origPos.y, origPos.z};
        double angleX = DihedralAngle(geom, base);
        double diffX = (angleX - origAngle) / h;
        // Differentiate wrt y
        geom->vertexPositions[other] = Vector3{origPos.x, origPos.y + h, origPos.z};
        double angleY = DihedralAngle(geom, base);
        double diffY = (angleY - origAngle) / h;
        // Differentiate wrt z
        geom->vertexPositions[other] = Vector3{origPos.x, origPos.y, origPos.z + h};
        double angleZ = DihedralAngle(geom, base);
        double diffZ = (angleZ - origAngle) / h;
        
        // Restore original
        geom->vertexPositions[other] = origPos;

        return Vector3{diffX, diffY, diffZ};
    }

    Vector3 GradientMeanCurvature(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other) {
        if (base.isBoundary()) return Vector3{0, 0, 0};

        // If differentiating wrt self, add contributions from all gradients of
        // surrounding (edge lengths / dihedral angles).
        if (base == other) {
            Vector3 gradient = Vector3{0, 0, 0};

            Vector3 centerPos = geom->vertexPositions[base];

            // Compute length / dihedral angle for each neighboring edge
            for (surface::Halfedge he : base.outgoingHalfedges()) {
                Vector3 neighborPos = geom->vertexPositions[he.twin().vertex()];
                
                // Product rule: want (grad length) / dihedral + length / (grad dihedral).
                
                // First term is (grad length) * dihedral angle.
                // Gradient of length is just the edge vector away from the neighbor.
                Vector3 lengthGrad = centerPos - neighborPos;
                lengthGrad = lengthGrad.normalize();
                double dihedral = DihedralAngle(geom, he.edge());

                // Second term is length * (grad dihedral angle)
                double length = norm(centerPos - neighborPos);
                Vector3 dihedralGrad = GradientDihedralAngle(geom, he.edge(), base);

                gradient += (lengthGrad * dihedral + length * dihedralGrad);
            }
            // Return the sum
            return gradient / 4;
        }

        // Otherwise, check if vertex is a neighbor.
        else {
            surface::Halfedge adjacentHE;
            bool found = false;
            for (const surface::Halfedge &neighborHE : base.outgoingHalfedges()) {
                // Look for the halfedge connecting to the neighbor
                if (neighborHE.twin().vertex() == other) {
                    adjacentHE = neighborHE;
                    found = true;
                    break;
                }
            }

            // If we don't find it, then it isn't a neighbor and the gradient is 0.
            if (!found) return Vector3{0, 0, 0};

            // Also need the previous and next halfedges around this vertex, since
            // moving "other" will also affect those dihedral angles
            surface::Halfedge nextHE = adjacentHE.twin().next();
            // Loop around to find the previous
            surface::Halfedge previousHE;
            surface::Halfedge current = adjacentHE;
            do {
                previousHE = current;
                current = current.twin().next();
            }
            while (current != adjacentHE);

            // Differentiate dihedral angle for the previous edge
            double lengthPrev = geom->edgeLength(previousHE.edge());
            Vector3 gradPrev = lengthPrev * GradientDihedralAngle(geom, previousHE.edge(), other);

            // Differentiate dihedral angle for middle edge
            Vector3 lengthGrad = geom->vertexPositions[other] - geom->vertexPositions[base];
            lengthGrad = lengthGrad.normalize();
            double dihedral = DihedralAngle(geom, adjacentHE.edge());
            double length = geom->edgeLength(adjacentHE.edge());
            Vector3 dihedralGrad = GradientDihedralAngle(geom, adjacentHE.edge(), other);
            Vector3 gradMiddle = lengthGrad * dihedral + length * dihedralGrad;

            // Differentiate dihedral angle for next edge
            double lengthNext = geom->edgeLength(nextHE.edge());
            Vector3 gradNext = lengthNext * GradientDihedralAngle(geom, nextHE.edge(), other);

            Vector3 gradient = gradPrev + gradMiddle + gradNext;
            return gradient / 4;
        }
    }

    Vector3 GradientGaussCurvature(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other) {
        // Gradient undefined for boundary vertices
        if (base.isBoundary()) return Vector3{0, 0, 0};

        // If differentiating by self, add up all gradients of surrounding tip angles
        if (other == base) {
            Vector3 gradient = Vector3{0, 0, 0};
            Vector3 p_base = geom->vertexPositions[base];

            for (const surface::Halfedge &he : base.outgoingHalfedges()) {
                // Gradient of a tip angle = negative gradients of the other
                // two angles in the triangle wrt this vertex.
                surface::Face f = he.face();
                Vector3 N_ij = geom->faceNormal(f);

                // Assuming vertices are ordered (base, i, j), take gradient of angle i
                Vector3 p_i = geom->vertexPositions[he.next().vertex()];
                Vector3 e_i = p_i - p_base;
                double l_i = norm(e_i);
                // Move in the direction perpendicular to the edge (base, i), in the
                // direction that will widen the angle.
                Vector3 gradTheta_i = -cross(N_ij, e_i);
                gradTheta_i = gradTheta_i.normalize();
                gradTheta_i /= l_i;

                // Take gradient of angle j
                Vector3 p_j = geom->vertexPositions[he.next().next().vertex()];
                Vector3 e_j = p_j - p_base;
                double l_j = norm(e_j);
                // Same principle as derivative of angle i
                Vector3 gradTheta_j = cross(N_ij, e_j);
                gradTheta_j = gradTheta_j.normalize();
                gradTheta_j /= l_j;

                // Gradient of angle ij is the negative of the other two, since the
                // sum of angles is constant.
                gradient = gradient - gradTheta_i - gradTheta_j;
            }
            // Result is minus the sum, because Gauss curvature is 2pi - sum.
            return -gradient;
        }

        else {
            // Search for the halfedge connecting base to other
            surface::Halfedge adjacentHE;
            bool found = false;
            for (const surface::Halfedge &neighborHE : base.outgoingHalfedges()) {
                if (neighborHE.twin().vertex() == other) {
                    adjacentHE = neighborHE;
                    found = true;
                    break;
                }
            }
            // If the other vertex isn't adjacent, gradient is 0.
            if (!found) return Vector3{0, 0, 0};

            // Need to differentiate two tip angles on either side of the connecting edge.
            Vector3 e_j = heVector(geom, adjacentHE);
            double l_j = norm(e_j);

            Vector3 N_ij = geom->faceNormal(adjacentHE.twin().face());
            Vector3 N_jk = geom->faceNormal(adjacentHE.face());

            Vector3 gradTheta_ij = cross(N_ij, e_j) / (l_j * l_j);
            Vector3 gradTheta_jk = -cross(N_jk, e_j) / (l_j * l_j);
            
            // Gradient is again negated because K = 2pi - sum
            return -gradTheta_ij - gradTheta_jk;
        }
    }

    // Gradient of (A - 2H + K) is just (grad A - 2 grad H + grad K)
    Vector3 GradientTubularEnergy(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other) {
        Vector3 gradA = GradientArea(geom, base, other);
        Vector3 gradH = GradientMeanCurvature(geom, base, other);
        Vector3 gradK = GradientGaussCurvature(geom, base, other);

        return gradA - 2 * gradH + gradK;
    }

    Vector3 GradientNumerical(surface::VertexPositionGeometry* geom, EnergyFunction energy, surface::Vertex base, surface::Vertex other, double h) {
        Vector3 origPos = geom->vertexPositions[other];
        double origEnergy = energy(geom, base);

        // Differentiate wrt x
        geom->vertexPositions[other] = Vector3{origPos.x + h, origPos.y, origPos.z};
        double energyX = energy(geom, base);
        double diffX = (energyX - origEnergy) / h;
        // Differentiate wrt y
        geom->vertexPositions[other] = Vector3{origPos.x, origPos.y + h, origPos.z};
        double energyY = energy(geom, base);
        double diffY = (energyY - origEnergy) / h;
        // Differentiate wrt z
        geom->vertexPositions[other] = Vector3{origPos.x, origPos.y, origPos.z + h};
        double energyZ = energy(geom, base);
        double diffZ = (energyZ - origEnergy) / h;

        // Restore original value
        geom->vertexPositions[other] = origPos;

        Vector3 deriv = Vector3{diffX, diffY, diffZ};
        //cout << deriv << endl;
        return deriv;
    }

    Vector3 GradientAreaNumerical(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other) {
        return GradientNumerical(geom, &VertexArea, base, other, 0.0001);
    }

    Vector3 GradientMeanNumerical(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other) {
        return GradientNumerical(geom, &VertexMeanCurvature, base, other, 0.0001);
    }

    Vector3 GradientGaussNumerical(surface::VertexPositionGeometry* geom, surface::Vertex base, surface::Vertex other) {
        return GradientNumerical(geom, &VertexGaussCurvature, base, other, 0.0001);
    }
}
