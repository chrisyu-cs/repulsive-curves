#include "args/args.hxx"
#include "json/json.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/string_cast.hpp"

#include "geometrycentral/surface/meshio.h"
#include "lws_app.h"
#include "polyscope/gl/ground_plane.h"
#include "utils.h"

#include "polyscope/curve_network.h"
#include "product/test_matrices.h"
#include "spatial/tpe_bvh.h"
#include "multigrid_domain.h"
#include "multigrid_hierarchy.h"

#include <random>

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

  void printMatrix(Eigen::MatrixXd &G, int precision) {
    std::cout.precision(precision);
    for (int i = 0; i < G.rows(); i++) {
      for (int j = 0; j < G.cols(); j++) {
        if (G(i,j) == 0) std::cout.precision(1);
        else std::cout.precision(precision);
        std::cout << std::fixed << G(i, j) << ((j == G.cols() - 1) ? "" : ", ");
      }
      std::cout << std::endl;
    }
  }

  void writeSparseMatrix(std::ofstream &file, Eigen::MatrixXd &G, int precision) {
    file.precision(precision);
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> values;

    for (int i = 0; i < G.rows(); i++) {
      for (int j = 0; j < G.cols(); j++) {
        if (G(i,j) != 0) {
          rows.push_back(i);
          cols.push_back(j);
          values.push_back(G(i, j));
        }
      }
    }

    for (size_t i = 0; i < rows.size(); i++) {
      file << rows[i];
      if (i < rows.size() - 1) file << ", ";
    }
    file << std::endl;

    for (size_t i = 0; i < cols.size(); i++) {
      file << cols[i];
      if (i < cols.size() - 1) file << ", ";
    }
    file << std::endl;

    for (size_t i = 0; i < values.size(); i++) {
      file << values[i];
      if (i < values.size() - 1) file << ", ";
    }
    file << std::endl;
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

    if (ImGui::Button("Print Sobolev gram")) {
      int nVerts = curves->NumVertices();
      Eigen::MatrixXd G;
      G.setZero(nVerts, nVerts);
      SobolevCurves::SobolevLengthScaled(curves, 3, 6, G, 0.1);
      printMatrix(G, 17);
    }

    if (ImGui::Button("Test coarsening")) {
      int nVerts = curves->NumVertices();
      int logNumVerts = log2(nVerts) - 4;

      // LWS::BVHNode3D* tree = CreateBVHFromCurve(curves);
      // Eigen::MatrixXd gradients;
      // gradients.setZero(nVerts, 3);
      // tpeSolver->FillGradientVectorBH(tree, gradients);
      // Eigen::VectorXd x = gradients.col(0);
      // delete tree;

      std::uniform_real_distribution<double> unif(-1, 1);
      std::default_random_engine re;
      re.seed(42);

      Eigen::VectorXd x(nVerts);
      for (int i = 0; i < nVerts; i++) {
        
        x(i) = unif(re);
      }

      // Eigen::VectorXd x;
      // x.setZero(nVerts);
      // for (int i = 0; i < nVerts; i++) {
      //   x(i) = sin(2 * M_PI * (double)i / nVerts);
      // }

      // double xMean = x.sum() / nVerts;
      // for (int i = 0; i < nVerts; i++) {
      //   x(i) -= xMean;
      // }
      // std::cout << "Mean of rhs = " << x.sum() << std::endl;
      
      // Eigen::VectorXd x;
      // x.setZero(nVerts);
      // x(nVerts / 3) = 1;
      // x(3 * nVerts / 4) = -1;

      // Eigen::VectorXd x;
      // x.setZero(nVerts);
      // x(0) = 1;

      long matrixStart = Utils::currentTimeMilliseconds();
      // Interval1DDomain* domain = new Interval1DDomain(nVerts);
      PolyCurveDomain* domain = new PolyCurveDomain(curves, 3, 6);
      Eigen::MatrixXd A = domain->GetFullMatrix();

      // printMatrix(A, 10);

      long multigridStart = Utils::currentTimeMilliseconds();

      std::cout << "Matrix assembly = " << (multigridStart - matrixStart) << " ms" << std::endl;

      // MultigridHierarchy<Interval1DDomain, DenseMatrixMult>* hierarchy =
      //  new MultigridHierarchy<Interval1DDomain, DenseMatrixMult>(domain, logNumVerts);

      MultigridHierarchy<PolyCurveDomain, DenseMatrixMult>* hierarchy =
        new MultigridHierarchy<PolyCurveDomain, DenseMatrixMult>(domain, logNumVerts);

      // for (size_t i = 0; i < hierarchy->levels.size(); i++) {
      //   Eigen::MatrixXd M = hierarchy->levels[i]->GetFullMatrix();
      //   std::string fname = "laplacian" + std::to_string(M.rows()) + ".csv";
      //   std::ofstream Mfile(fname);
      //   writeSparseMatrix(Mfile, M, 10);
      //   Mfile.close();
      //   std::cout << "Wrote to " << fname << std::endl;
      // }

      // for (size_t i = 0; i < hierarchy->prolongationOps.size(); i++) {
      //   Eigen::MatrixXd M = hierarchy->prolongationOps[i].matrices[0].M.toDense();
      //   std::string fname = "prolongation_" + std::to_string(M.rows()) + "x" + std::to_string(M.cols()) + ".csv";
      //   std::ofstream Mfile(fname);
      //   writeSparseMatrix(Mfile, M, 10);
      //   Mfile.close();
      //   std::cout << "Wrote to " << fname << std::endl;
      // }

      Eigen::VectorXd sol = hierarchy->VCycleSolve(x, MultigridMode::MatrixOnly);
      long multigridEnd = Utils::currentTimeMilliseconds();
      std::cout << "Multigrid time = " << (multigridEnd - multigridStart) << " ms" << std::endl;

      Eigen::VectorXd ref_sol = A.partialPivLu().solve(x);
      long solveEnd = Utils::currentTimeMilliseconds();
      std::cout << "Direct solve time = " << (solveEnd - multigridEnd) << " ms" << std::endl;

      Eigen::MatrixXd sols(sol.rows(), 3);
      sols.col(0) = x;
      sols.col(1) = sol;
      sols.col(2) = ref_sol.block(0, 0, nVerts, 1);

      // Eigen::BiCGSTAB<Eigen::MatrixXd, Eigen::IdentityPreconditioner> bicg(A);
      // bicg.setMaxIterations(10);
      // Eigen::VectorXd bicg_sol = bicg.solve(x);
      
      // std::cout << "Normal BiCGStab finished in " << bicg.iterations() << std::endl;
      // Eigen::VectorXd bicg_resid = (x - A * bicg_sol);
      // std::cout << "BiCGStab residual = " << bicg_resid.lpNorm<Eigen::Infinity>() << std::endl;
      // std::cout << "BiCGStab error = " << (100 * (bicg_sol - ref_sol).norm() / ref_sol.norm()) << " percent" << std::endl;

      for (int i = 0; i < sol.rows(); i++) {
        // std::cout << i << ", " << x(i) << ", " << sol(i) << ", " << ref_sol(i) << std::endl;
      }

      Eigen::VectorXd diff = sol - ref_sol;

      std::cout << "Reference norm = " << ref_sol.norm() << std::endl;
      std::cout << "Multigrid norm = " << sol.norm() << std::endl;
      std::cout << "Difference norm = " << diff.norm() << std::endl;
      std::cout << "Multigrid error from ground truth = " << 100 * (diff.norm() / ref_sol.norm()) << " percent" << std::endl;

      double finalResidual = (x - A * sol).lpNorm<Eigen::Infinity>();

      double dot = sol.normalized().dot(ref_sol.normalized());
      std::cout << "Dot product between directions = " << dot << std::endl;
      std::cout << "Multigrid relative residual = " << (finalResidual / x.lpNorm<Eigen::Infinity>()) << std::endl;

      delete hierarchy;
    }

    if (ImGui::Button("Build hierarchy")) {
      int nVerts = curves->NumVertices();
      int logNumVerts = log2(nVerts) - 4;

      PolyCurveDomain* domain = new PolyCurveDomain(curves, 3, 6);

      MultigridHierarchy<PolyCurveDomain, DenseMatrixMult>* hierarchy =
        new MultigridHierarchy<PolyCurveDomain, DenseMatrixMult>(domain, logNumVerts);
      
      for (size_t i = 1; i < hierarchy->levels.size(); i++) {
        PolyCurveGroup* g = static_cast<PolyCurveDomain*>(hierarchy->levels[i])->curves;
        DisplayCurves(g, "multigrid" + std::to_string(i));

        Eigen::VectorXd position_xs = g->GetPositionMatrix().col(0);
        Eigen::VectorXd prolonged = hierarchy->prolongationOps[i - 1].prolong(position_xs, MultigridMode::MatrixOnly);
        Eigen::VectorXd restricted = hierarchy->prolongationOps[i - 1].restrictWithPinv(prolonged, MultigridMode::MatrixOnly);

        Eigen::VectorXd diff = position_xs - restricted;

        std::cout << "Level " << i << " (RJx - x) = " << diff.norm() << std::endl;
      }
    }

    if (ImGui::Button("Output frame")) {
      outputFrame();
    }

    ImGui::Checkbox("Run TPE", &LWSOptions::runTPE);
    ImGui::Checkbox("Output frames", &LWSOptions::outputFrames);

    bool buttonStepTPE = ImGui::Button("Single TPE step");
    bool buttonPlotTPE = ImGui::Button("Plot TPE gradient");

    if (LWSOptions::runTPE || buttonStepTPE) {
      if (LWSOptions::outputFrames && LWSOptions::frameNum == 0) {
        outputFrame();
      }

      bool good_step = tpeSolver->StepSobolevLS(false);
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
      for (int i = 0; i < nVerts; i++) {
        loop->positions[i] = loop->positions[i] - center;
      }
    }

    UpdateCurvePositions();
  }

  void LWSApp::initSolver() {
    if (!tpeSolver) {
      // Set up solver
      tpeSolver = new TPEFlowSolverSC(curves);
    }
  }

  void LWSApp::UpdateCurvePositions() {
    // Update the positions on the space curve
    polyscope::CurveNetwork* curveNetwork = polyscope::getCurveNetwork(surfaceName);
    std::vector<glm::vec3> curve_vecs(curves->NumVertices());

    for (size_t i = 0; i < curves->curves.size(); i++)
    {
      PolyCurve* c = curves->curves[i];
      int nVerts = c->NumVertices();
      for (int j = 0; j < nVerts; j++) {
        Vector3 v = c->positions[j];
        curve_vecs[c->offset + j] = glm::vec3{v.x, v.y, v.z};
      }
    }

    curveNetwork->updateNodePositions(curve_vecs);
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
      Halfedge he = b.halfedge().twin();
      Halfedge start = b.halfedge().twin();

      // int i = 0;

      do {
        Vector3 v = geom->vertexPositions[he.vertex()];
        positions.push_back(v);
        // std::cout << i++ << ", " << v << std::endl;
        he = he.next();
      }
      while (he != start);

      PolyCurve* pc = new PolyCurve(positions);
      curves->AddCurve(pc);
      std::cout << "Added boundary curve of length " << pc->NumVertices() << std::endl;
    }

    surfaceName = polyscope::guessNiceNameFromPath(filename);
  }

  void LWSApp::DisplayCurves(PolyCurveGroup* curves, std::string name) {
    size_t nVerts = curves->NumVertices();
    std::vector<glm::vec3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    
    for (size_t i = 0; i < curves->curves.size(); i++) {
      PolyCurve* c = curves->curves[i];
      int nVerts = c->NumVertices();
      // Add interior edges and vertices
      for (int i = 0; i < nVerts; i++) {
        nodes.push_back(glm::vec3{c->positions[i].x, c->positions[i].y, c->positions[i].z});
        size_t v_i = c->offset + i;
        if (i > 0) {
          edges.push_back({v_i - 1, v_i});
        }
      }
      // Close loop
      edges.push_back({size_t(c->offset + nVerts - 1), size_t(c->offset)});
    }

    polyscope::registerCurveNetwork(name, nodes, edges);
    polyscope::getCurveNetwork(name)->radius *= 5;

    centerLoopBarycenter(curves);
  }

  void LWSApp::DisplayCyclicList(std::vector<Vector3> &positions, std::string name) {
    std::vector<std::array<size_t, 2>> edges;

    for (size_t i = 0; i < positions.size() - 1; i++) {
      edges.push_back({i, i+1});
    }
    edges.push_back({positions.size() - 1, 0});
    polyscope::registerCurveNetwork(name, positions, edges);
    polyscope::getCurveNetwork(name)->radius *= 5;
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
  app->DisplayCurves(app->curves, app->surfaceName);
  app->initSolver();

  // Show the gui
  polyscope::show();

  return 0;
}
