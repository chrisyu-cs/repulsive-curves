#include "args/args.hxx"
#include "json/json.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/string_cast.hpp"

#include "geometrycentral/surface/meshio.h"
#include "lws_app.h"
#include "polyscope/gl/ground_plane.h"
#include "utils.h"
#include "spacecurves/space_curve.h"

#include "spatial/tpe_kdtree.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

namespace LWS {

  Vector3 colorOfValue(double value) {
      value = clamp(value, -1.0, 1.0);
      value = (value + 1) / 2;
      return Vector3{value, 2 * value * (1 - value), 1 - value};
  }

  void LWSApp::outputFrame() {
    char buffer[5];
    std::snprintf(buffer, sizeof(buffer), "%04d", LWSOptions::frameNum);
    string fname = "frames/frame" + std::string(buffer) + ".png";
    LWSOptions::frameNum++;
    polyscope::screenshot(fname, false);
    std::cout << "Wrote screenshot to " << fname << std::endl;
  }

  void LWSApp::plotBHError(double alpha, double beta) {
    BVHNode3D* bvh = CreateBVHFromCurve(curves);
    std::vector<BHPlotData> data;
    for (int i = 0; i < curves->NumVertices(); i++) {
      bvh->testGradientSingle(data, curves->GetCurvePoint(i), curves, alpha, beta);
    }

    std::ofstream outFile("plot.csv");

    for (size_t i = 0; i < data.size(); i++) {
      std::cout << data[i].theta << ", " << data[i].error << ", " << data[i].gradientNorm << ", " << data[i].minWidth << ", " << data[i].maxWidth << std::endl;
      outFile << data[i].theta << ", " << data[i].error << ", " << data[i].gradientNorm << ", " << data[i].minWidth << ", " << data[i].maxWidth << std::endl;
    }

    std::cout << "Wrote data points to plot.csv" << std::endl;

    outFile.close();

    delete bvh;
  }

  void LWSApp::customWindow() {

    ImGuiIO& io = ImGui::GetIO();

    ImGui::Begin("Hull settings", &LWSOptions::showWindow);

    if (io.KeysDown[(int)' '] && io.KeysDownDurationPrev[(int)' '] == 0) {
      cout << "space bar" << endl;
    }
    if (ImGui::Button("Test derivatives")) {
      //TestBoundaryQuantities(mesh, geom);
    }

    if (ImGui::Button("Check MV product")) {
      tpeSolver->CompareMatrixVectorProduct();
    }

    ImGui::Checkbox("Run TPE", &LWSOptions::runTPE);
    ImGui::Checkbox("Output frames", &LWSOptions::outputFrames);

    bool buttonStepTPE = ImGui::Button("Single TPE step");
    bool buttonPlotTPE = ImGui::Button("Plot TPE gradient");

    if (LWSOptions::runTPE || buttonStepTPE) {
      if (LWSOptions::outputFrames && LWSOptions::frameNum == 0) {
        outputFrame();
      }

      bool good_step = tpeSolver->StepSobolevProjLS(true);
      UpdateCurvePositions();
      if (!good_step) {
        std::cout << "Stopped because line search could not take a step." << std::endl;
        LWSOptions::runTPE = false;
      }
      
      if (LWSOptions::outputFrames) {
        outputFrame();
      }
    }

    if (ImGui::Button("Plot BH error")) {
      plotBHError(3, 6);
    }

    if (ImGui::Button("Curve to OBJ")) {
      std::cout << "TODO" << std::endl;
    }

    ImGui::Checkbox("Use Sobalev", &LWSOptions::useSobalev);

    double delta = 0.001;

    ImGui::End();
  }

  void LWSApp::centerLoopBarycenter(PolyCurveGroup* curves) {
    Vector3 center = curves->Barycenter();

    for (PolyCurve* loop : curves->curves) {
      int nVerts = loop->NumVertices();
      for (int i = 0; i < nVerts; i++)
      loop->positions[i] = loop->positions[i] - center;
    }

    UpdateCurvePositions();
  }

  std::string LWSApp::getCurveName(int i) {
    return surfaceName + "-curve" + std::to_string(i);
  }

  void LWSApp::initSolver() {
    if (!tpeSolver) {
      // Set up solver
      tpeSolver = new TPEFlowSolverSC(curves);
    }
  }

  void LWSApp::UpdateCurvePositions() {
    // Update the positions on the space curve
    for (size_t i = 0; i < curves->curves.size(); i++)
    {
      PolyCurve* polyCurve = curves->curves[i];
      polyscope::SpaceCurve* curve = polyscope::getSpaceCurve(getCurveName(i));
      std::vector<glm::vec3> curve_vecs(curves->NumVertices());
      for (size_t i = 0; i < curve_vecs.size(); i++) {
        Vector3 v = polyCurve->Position(i);
        curve_vecs[i] = glm::vec3{v.x, v.y, v.z};
      }
      curve->updatePoints(curve_vecs);
    }
    polyscope::requestRedraw();
  }

  void LWSApp::processFileOBJ(std::string filename) {
    if (!curves) {
      curves = new PolyCurveGroup();
    }

    std::cout << "Make curves for " << filename << std::endl;

    std::tie(mesh, geom) = loadMesh(filename);
    geom->requireVertexPositions();

    VertexData<size_t> indices = mesh->getVertexIndices();

    std::vector<Vector3> positions;

    for (BoundaryLoop b : mesh->boundaryLoops()) {
      positions.clear();
      Halfedge he = b.halfedge().next();
      Halfedge start = b.halfedge().next();

      do {
        Vector3 v = geom->vertexPositions[he.vertex()];
        positions.push_back(v);
        he = he.next();
      }
      while (he != start);

      PolyCurve* pc = new PolyCurve(positions);
      curves->AddCurve(pc);
      std::cout << "Added boundary curve of length " << pc->NumVertices() << std::endl;
    }

    surfaceName = polyscope::guessNiceNameFromPath(filename);
  }

  void LWSApp::DisplayCurves() {
    std::vector<Vector3> positions;

    size_t nVerts = curves->NumVertices();

    for (size_t c = 0; c < curves->curves.size(); c++) {
      positions.clear();
      PolyCurve* curve = curves->curves[c];
      int cVerts = curve->NumVertices();
      for (int i = 0; i < cVerts; i++) {
        positions.push_back(curve->Position(i));
      }
      polyscope::registerSpaceCurve(getCurveName(c), positions, true);
    }

    centerLoopBarycenter(curves);
  }
}

LWS::LWSApp* LWS::LWSApp::instance;

bool endsWith(const std::string& str, const std::string& suffix) {
  return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void processFile(LWS::LWSApp* app, string filename) {
  // Dispatch to correct varient
  if (endsWith(filename, ".obj")) {
    app->processFileOBJ(filename);
  } else {
    cerr << "Unrecognized file type for " << filename << endl;
  }
}

void customWindow() {
  LWS::LWSApp::instance->customWindow();
}


int main(int argc, char** argv) {
  // Configure the argument parser
  args::ArgumentParser parser("A simple demo of Polyscope.\nBy "
                              "Nick Sharp (nsharp@cs.cmu.edu)",
                              "");
  args::PositionalList<string> files(parser, "files", "One or more files to visualize");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;

    std::cerr << parser;
    return 1;
  }
  // Options
  polyscope::options::autocenterStructures = false;
  // polyscope::view::windowWidth = 600;
  // polyscope::view::windowHeight = 800;
  polyscope::gl::groundPlaneEnabled = false;

  LWS::LWSApp* app = new LWS::LWSApp();
  LWS::LWSApp::instance = app;

  // Initialize polyscope
  polyscope::init();
  // Add a few gui elements
  polyscope::state::userCallback = &customWindow;

  for (string f : files) {
    processFile(LWS::LWSApp::instance, f);
  }

  std::vector<std::array<double, 3>> vertexPositions;
  std::vector<std::vector<size_t>> faceIndices;

  polyscope::SurfaceMesh* m = polyscope::registerSurfaceMesh("empty", vertexPositions, faceIndices);
  app->DisplayCurves();
  app->initSolver();

  // Show the gui
  polyscope::show();

  return 0;
}
