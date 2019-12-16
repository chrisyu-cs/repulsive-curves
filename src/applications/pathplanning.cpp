#include "applications/pathplanning.h"
#include <iostream>
#include <fstream>

namespace LWS {
    namespace Applications {

        bool EdgeIntersectsYPlane(CurveEdge* edge, double y, Vector3 &intersection) {
            Vector3 prev = edge->prevVert->Position();
            Vector3 next = edge->nextVert->Position();

            // Make 'prev' the lower one, and 'next' the higher one
            if (prev.y > next.y) {
                Vector3 temp = next;
                next = prev;
                prev = temp;
            }

            if (y < prev.y || y > next.y) {
                return false;
            }
            if (y == prev.y) {
                intersection = prev;
                return true;
            }
            if (y == next.y) {
                intersection = next;
                return true;
            }

            // y is somewhere between prev.y and next.y, so solve for the point of intersection
            double lengthY = next.y - prev.y;
            double offsetY = y - prev.y;
            double t = offsetY / lengthY;

            intersection = (1 - t) * prev + t * next;
            return true;
        }

        std::vector<std::vector<Vector3>> SamplePathsAlongY(PolyCurveNetwork* curves, int nSamples) {
            int nComponents = curves->NumComponents();
            int nVerts = curves->NumVertices();
            int nEdges = curves->NumEdges();

            double minY = curves->GetVertex(0)->Position().y;
            double maxY = curves->GetVertex(0)->Position().y;

            std::vector<std::unordered_set<int>> componentSets(nComponents);
            std::vector<std::vector<Vector3>> pathPositions(nComponents);

            // Group indices by which component they are in
            for (int c = 0; c < nComponents; c++) {
                int nVinC = curves->NumVerticesInComponent(c);
                for (int i = 0; i < nVinC; i++) {
                    CurveVertex* v_i = curves->GetVertexInComponent(c, i);
                    componentSets[c].insert(v_i->GlobalIndex());
                    // At the same time, get min and max Y bounds
                    minY = fmin(minY, v_i->Position().y);
                    maxY = fmax(maxY, v_i->Position().y);
                }
            }

            double step = (maxY - minY) / (nSamples - 1);

            // Consider one slice in the Y direction at a time
            for (int y = 0; y < nSamples; y++) {
                double curY = minY + step * y;
                if (y == nSamples - 1) curY = maxY;

                std::cout << "Processing slice " << y << " (height " << curY << ")" << std::endl;

                // Find all edge intersections with the current Y slice
                for (int i = 0; i < nEdges; i++) {
                    Vector3 point{0, 0, 0};
                    CurveEdge* e_i = curves->GetEdge(i);
                    // If this edge intersects with the Y plane, then that gives a point
                    // on the trajectory
                    if (EdgeIntersectsYPlane(e_i, curY, point)) {
                        // Find which component this edge is in
                        for (int c = 0; c < nComponents; c++) {
                            if (componentSets[c].count(e_i->prevVert->GlobalIndex())) {
                                // Push onto the trajectory for the path from that component
                                pathPositions[c].push_back(point);
                                break;
                            }
                        }
                    }
                }
            }

            return pathPositions;
        }

        void SampleAndWritePaths(PolyCurveNetwork* curves, int nSamples, std::string fname) {
            std::vector<std::vector<Vector3>> paths = SamplePathsAlongY(curves, nSamples);

            std::ofstream myfile(fname);
            std::vector<std::array<int, 2>> edges;
            int start = 1;

            for (size_t c = 0; c < paths.size(); c++) {
                int prev = -1;
                for (size_t i = 0; i < paths[c].size(); i++) {
                    Vector3 p = paths[c][i];
                    myfile << "v " << p.x << " " << p.y << " " << p.z << std::endl;

                    if (prev < 0) prev = start;
                    else {
                        int next = prev + 1;
                        edges.push_back({prev, next});
                        prev = next;
                    }

                }
                myfile << std::endl;
                start += paths[c].size();
            }

            for (size_t i = 0; i < edges.size(); i++) {
                myfile << "l " << edges[i][0] << " " << edges[i][1] << std::endl;
            }

            myfile.close();
        }

    }
}