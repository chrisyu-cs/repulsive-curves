#include "obstacles/mesh_obstacle.h"
#include "tpe_energy_sc.h"

namespace LWS {
    MeshObstacle::MeshObstacle(std::shared_ptr<HalfedgeMesh> m, std::shared_ptr<VertexPositionGeometry> geom) {
        mesh = m;
        geometry = geom;

    }


}