#include "boundary_derivatives.h"
#include "utils.h"

using namespace geometrycentral;
using namespace std;

namespace LWS {

    void TestBoundaryQuantities(surface::HalfedgeMesh* mesh, surface::VertexPositionGeometry* geom) {
        cout << "Test boundaries" << endl;
        for (surface::Vertex v : mesh->vertices()) {
            Vector3 length = GradientBLengthNum(geom, v, v);
            Vector3 curvature = GradientBCurvatureNum(geom, v, v);
            Vector3 torsion = GradientBTorsionNum(geom, v, v);

            if (v.isBoundary()) {
                /*
                for (surface::Vertex v2 : v.adjacentVertices()) {
                    if (v2.isBoundary()) {
                        Vector3 lenNum = GradientBLengthNum(surface->geometry, v, v2);
                        Vector3 lenGrad = GradientBLength(surface->geometry, v, v2);

                        cout << "Length gradient " << lenGrad << endl;
                        cout << "Length grad num " << lenNum << endl;

                        Vector3 curvNum = GradientBCurvatureNum(surface->geometry, v, v2);
                        Vector3 curvGrad = GradientBCurvature(surface->geometry, v, v2);

                        cout << "Curvature gradient " << curvGrad << endl;
                        cout << "Curvature grad num " << curvNum << endl;
                    }
                }
                */

                HalfedgePair pair = BoundaryEdges(v);

                surface::Vertex first = v;
                Vector3 torsionNum1 = GradientBTorsionNum(geom, v, v);
                Vector3 torsion1 = GradientBTorsion(geom, v, v);
                double lenNum1 = norm(torsionNum1);
                double len1 = norm(torsion1);
                double lenRatio = lenNum1 / len1;
                torsion1.normalize();
                torsionNum1.normalize();
                double dirDot = dot(torsion1, torsionNum1);

                cout << "Dot / ratio = " << dirDot << " / " << lenRatio << endl;
            }
        }
    }

    HalfedgePair BoundaryEdges(surface::Vertex vert) {
        HalfedgePair pair;

        bool forwardFound = false;
        bool reverseFound = false;

        for (surface::Halfedge he : vert.outgoingHalfedges()) {
            // Check if this is the forward boundary halfedge
            if (!he.isInterior()) {
                pair.next = he;
                forwardFound = true;
            }
            // Check if this is the reverse boundary halfedge
            else if (!he.twin().isInterior()) {
                pair.prev = he;
                reverseFound = true;
            }
        }

        pair.valid = (forwardFound && reverseFound);
        return pair;
    }

    double VertexBoundaryLength(surface::VertexPositionGeometry* geom, surface::Vertex vert) {
        HalfedgePair pair = BoundaryEdges(vert);

        if (pair.valid) {
            return 0.5 * (geom->edgeLength(pair.next.edge()) + geom->edgeLength(pair.prev.edge()));
        }
        else return 0;
    }

    double VertexBoundaryCurvature(surface::VertexPositionGeometry* geom, surface::Vertex vert) {
        HalfedgePair pair = BoundaryEdges(vert);

        if (pair.valid) {
            Vector3 next_vec = heVector(geom, pair.next);
            Vector3 prev_vec = heVector(geom, pair.prev.twin());
            next_vec.normalize();
            prev_vec.normalize();

            // Compute the turning angle
            double angle = acos(dot(next_vec, prev_vec));
            return angle;
        }
        else return 0;
    }

    double VertexBoundaryTorsion(surface::VertexPositionGeometry* geom, surface::Vertex vert) {
        HalfedgePair pair = BoundaryEdges(vert);

        if (pair.valid) {
            Vector3 next_vec = heVector(geom, pair.next);
            Vector3 prev_vec = heVector(geom, pair.prev.twin());
            Vector3 next_next_vec = heVector(geom, pair.next.next());
            next_vec.normalize();
            prev_vec.normalize();
            next_next_vec.normalize();

            // Compute the two binormals around the vertex
            Vector3 binormal1 = cross(prev_vec, next_vec);
            Vector3 binormal2 = cross(next_vec, next_next_vec);

            if (norm2(binormal1) < 1e-5) return 0;
            if (norm2(binormal2) < 1e-5) return 0;

            binormal1.normalize();
            binormal2.normalize();

            // Compute the rotation angle between the binormals
            Vector3 forward = cross(binormal1, binormal2);
            double sinA = dot(forward, next_vec);
            return asin(sinA);
        }
        else return 0;
    }

    Vector3 GradientBLength(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other) {
        HalfedgePair pair = BoundaryEdges(vert);
        if (pair.valid) {
            // If differentiating wrt self, need to consider both side
            if (vert == other) {
                // Gradient is half the sum of unit vectors along both outgoing edges
                Vector3 forward = heVector(geom, pair.next);
                forward.normalize();
                Vector3 reverse = heVector(geom, pair.prev);
                reverse.normalize();

                return -0.5 * (forward + reverse);
            }
            // Otherwise, only consider the one edge between other and vert
            else {
                if (pair.next.twin().vertex() == other) {
                    Vector3 forward = heVector(geom, pair.next);
                    forward.normalize();
                    return forward / 2;
                }
                else if (pair.prev.twin().vertex() == other) {
                    Vector3 back = heVector(geom, pair.prev);
                    back.normalize();
                    return back / 2;
                }
                else return Vector3{0, 0, 0};
            }
        }
        else return Vector3{0, 0, 0};
    }

    Vector3 GradientBCurvature(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other) {
        HalfedgePair pair = BoundaryEdges(vert);
        if (pair.valid) {
            Vector3 next_vec = heVector(geom, pair.next);
            Vector3 prev_incoming = heVector(geom, pair.prev.twin());
            Vector3 prev_outgoing = heVector(geom, pair.prev);

            double next_length = norm(next_vec);
            double prev_length = norm(prev_incoming);

            Vector3 binormal = cross(prev_outgoing / prev_length, next_vec / next_length);

            if (norm2(binormal) < 1e-5) return Vector3{0, 0, 0};
            binormal.normalize();

            // If differentiating wrt self, consider both edges
            if (vert == other) {
                // Tangent of circle centered on previous vertex
                Vector3 prev_tangent = cross(binormal, prev_incoming);
                // Circle centered on next vertex
                Vector3 next_tangent = cross(binormal, next_vec);

                // If we move along the circumference of a circle at unit
                // speed, we change the radians at a rate of 1 / radius.

                // (On a unit circle, this is just 1, since it takes 2pi time
                // to go all the way around, and the radians change by 2pi in that time.
                // If we double the radius, it now takes 4pi time to go around, but
                // the radians remain the same, so the rate of change is halved.)
                double prev_rate = 1 / prev_length;
                double next_rate = 1 / next_length;

                prev_tangent.normalize();
                next_tangent.normalize();

                return prev_rate * prev_tangent + next_rate * next_tangent;
            }

            else {
                if (pair.next.twin().vertex() == other) {
                    Vector3 next_tangent = cross(binormal, next_vec);
                    double next_rate = 1 / next_length;
                    next_tangent.normalize();

                    return -next_rate * next_tangent;
                }
                else if (pair.prev.twin().vertex() == other) {
                    Vector3 prev_tangent = cross(binormal, prev_incoming);
                    double prev_rate = 1 / prev_length;
                    prev_tangent.normalize();

                    return -prev_rate * prev_tangent;
                }
            }
        }
        return Vector3{0, 0, 0};
    }

    Vector3 VectorToLine(Vector3 pt, Vector3 endpt1, Vector3 endpt2) {
        Vector3 lineVec = endpt1 - endpt2;
        lineVec.normalize();
        // Find the projection of the displacement onto the line
        Vector3 toLine = endpt1 - pt;
        Vector3 parallel = dot(lineVec, toLine) * lineVec;
        // Subtract that from the displacement to get the perpendicular component
        toLine = toLine - parallel;
        return toLine;
    }

    Vector3 NormalOfTriangle(Vector3 p1, Vector3 p2, Vector3 p3) {
        Vector3 e12 = p2 - p1;
        Vector3 e23 = p3 - p2;
        Vector3 normal = cross(e12, e23);
        normal.normalize();
        return normal;
    }

    Vector3 GradientBTorsion(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other) {
        HalfedgePair pair = BoundaryEdges(vert);
        if (pair.valid) {
            Vector3 middleEdge = heVector(geom, pair.next);
            Vector3 lastEdge = heVector(geom, pair.next.next());
            Vector3 firstEdge = heVector(geom, pair.prev.twin());

            surface::Vertex v1 = pair.prev.twin().vertex();
            surface::Vertex v2 = pair.prev.vertex();
            surface::Vertex v3 = pair.next.twin().vertex();
            surface::Vertex v4 = pair.next.next().twin().vertex();

            Vector3 p1 = geom->vertexPositions[v1];
            Vector3 p2 = geom->vertexPositions[v2];
            Vector3 p3 = geom->vertexPositions[v3];
            Vector3 p4 = geom->vertexPositions[v4];

            // Case for first vertex in the chain
            if (other == v1) {
                // The triangle we want to rotate is (v1, v2, v3)
                // and the one we want to rotate towards is (v2, v3, v4)
                Vector3 triNormal = NormalOfTriangle(p1, p2, p3);
                Vector3 oppNormal = NormalOfTriangle(p2, p3, p4);
                
                // If the normals flip signs, we'll need to negate the rotation at the end
                double normalSign = dot(triNormal, oppNormal);
                normalSign = normalSign / abs(normalSign);

                // We also need to correct for if the edge orientation
                // doesn't match the intended rotation axis
                Vector3 rotationAxis = cross(triNormal, oppNormal);
                rotationAxis.normalize();
                Vector3 edgeVec = p3 - p2;
                edgeVec.normalize();
                double rotationDirection = dot(rotationAxis, edgeVec);

                Vector3 perpVec = VectorToLine(p1, p2, p3);
                double radius = norm(perpVec);
                Vector3 tangentDir = cross(rotationAxis, perpVec);
                tangentDir.normalize();
                double rate = 1 / radius;

                return normalSign * rotationDirection * rate * tangentDir;
            }
            // Case for last vertex
            else if (other == v4) {
                // The triangle to rotate is (v2, v3, v4)
                // and the opposite triangle is (v1, v2, v3)
                Vector3 triNormal = NormalOfTriangle(p2, p3, p4);
                Vector3 oppNormal = NormalOfTriangle(p1, p2, p3);
                
                // Normals sign correction
                double normalSign = dot(triNormal, oppNormal);
                normalSign = normalSign / abs(normalSign);

                // Rotation axis sign correction
                Vector3 rotationAxis = cross(triNormal, oppNormal);
                rotationAxis.normalize();
                Vector3 edgeVec = p3 - p2;
                edgeVec.normalize();
                double rotationDirection = dot(rotationAxis, edgeVec);

                Vector3 perpVec = VectorToLine(p2, p3, p4);
                double radius = norm(perpVec);
                Vector3 tangentDir = cross(rotationAxis, perpVec);
                tangentDir.normalize();
                double rate = 1 / radius;

                return -1 * normalSign * rotationDirection * rate * tangentDir;
            }

            else if (other == v2) {
                Vector3 triNormal = NormalOfTriangle(p2, p3, p1);
                Vector3 oppNormal = NormalOfTriangle(p3, p1, p4);
                // Normal sign correction
                double normalSign = dot(triNormal, oppNormal);
                normalSign = normalSign / abs(normalSign);

                Vector3 rotationAxis = cross(triNormal, oppNormal);
                rotationAxis.normalize();
                Vector3 perpVec = VectorToLine(p2, p3, p1);
                double rate = 1 / norm(perpVec);

                Vector3 tangentDir = cross(rotationAxis, perpVec);
                tangentDir.normalize();
                // Rotation direction correction
                Vector3 edgeVec = p1 - p3;
                edgeVec.normalize();
                double rotationDirection = dot(rotationAxis, edgeVec);

                Vector3 deriv1 = normalSign * rotationDirection * rate * tangentDir;

                // Part for second triangle
                triNormal = NormalOfTriangle(p2, p3, p4);
                oppNormal = NormalOfTriangle(p3, p4, p1);
                // Normal sign correction
                normalSign = dot(triNormal, oppNormal);
                normalSign = normalSign / abs(normalSign);

                rotationAxis = cross(triNormal, oppNormal);
                rotationAxis.normalize();
                perpVec = VectorToLine(p2, p3, p4);
                rate = 1 / norm(perpVec);

                tangentDir = cross(rotationAxis, perpVec);
                tangentDir.normalize();
                // Rotation direction correction
                edgeVec = p4 - p3;
                edgeVec.normalize();
                rotationDirection = dot(rotationAxis, edgeVec);

                Vector3 deriv2 = -1 * normalSign * rotationDirection * rate * tangentDir;

                return deriv1 + deriv2;
            }
            else return Vector3{0, 0, 0};
        }
        else return Vector3{0, 0, 0};
    }

    Vector3 GradientBLengthNum(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other) {
        HalfedgePair pair = BoundaryEdges(vert);
        if (pair.valid) {
            return GradientNumerical(geom, &VertexBoundaryLength, vert, other, 0.0001);
        }
        else return Vector3{0, 0, 0};
    }

    Vector3 GradientBCurvatureNum(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other) {
        HalfedgePair pair = BoundaryEdges(vert);
        if (pair.valid) {
            return GradientNumerical(geom, &VertexBoundaryCurvature, vert, other, 0.0001);
        }
        else return Vector3{0, 0, 0};
    }

    Vector3 GradientBTorsionNum(surface::VertexPositionGeometry* geom, surface::Vertex vert, surface::Vertex other) {
        HalfedgePair pair = BoundaryEdges(vert);
        if (pair.valid) {
            return GradientNumerical(geom, &VertexBoundaryTorsion, vert, other, 0.0001);
        }
        else return Vector3{0, 0, 0};
    }
}