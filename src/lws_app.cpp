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
#include "multigrid/multigrid_domain.h"
#include "multigrid/multigrid_hierarchy.h"
#include "multigrid/nullspace_projector.h"

#include "poly_curve_network.h"

#include <limits>
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

    if (ImGui::Button("Print tree")) {
      int nVerts = curves->NumVertices();
      LWS::BVHNode3D* tree = CreateEdgeBVHFromCurve(curves);
      BlockClusterTree* mult = new BlockClusterTree(curves, tree, 0.5, 2, 4);

      mult->PrintData();

      tree->assignIDs();

      std::ofstream treefile;
      treefile.open("tree-indices.csv");
      tree->printIDs(treefile);
      treefile.close();

      std::ofstream admfile;
      admfile.open("admissible-clusters.csv");
      mult->PrintAdmissibleClusters(admfile);
      admfile.close();

      std::ofstream inadmfile;
      inadmfile.open("inadmissible-clusters.csv");
      mult->PrintInadmissibleClusters(inadmfile);
      inadmfile.close();
    }

    if (ImGui::Button("Test scaling")) {
      Eigen::MatrixXd initGradients;
      initGradients.setZero(curves->NumVertices(), 3);
      double initEnergy = tpeSolver->CurrentEnergyDirect();
      tpeSolver->FillGradientVectorDirect(initGradients);

      curves->positions = curves->positions * 2;

      double scaledEnergy = tpeSolver->CurrentEnergyDirect();
      Eigen::MatrixXd scaledGradients;
      scaledGradients.setZero(curves->NumVertices(), 3);
      tpeSolver->FillGradientVectorDirect(scaledGradients);
      
      Eigen::MatrixXd gradientRatio = scaledGradients;

      for (int i = 0; i < curves->NumVertices(); i++) {
        for (int j = 0; j < 3; j++) {
          gradientRatio(i, j) = scaledGradients(i, j) / initGradients(i, j);
        }
      }

      std::cout << "Gradient ratio = \n" << gradientRatio << std::endl;
      std::cout << "Energy ratio = " << scaledEnergy << " / " << initEnergy << " = " << (scaledEnergy / initEnergy) << std::endl;

    }

    if (ImGui::Button("Test multiply")) {
      std::cout << "Machine epsilon: " << std::numeric_limits<double>::epsilon() << std::endl;

      int nVerts = curves->NumVertices();
      LWS::BVHNode3D* tree = CreateEdgeBVHFromCurve(curves);
      BlockClusterTree* mult = new BlockClusterTree(curves, tree, 0.1, 2, 4);

      mult->PrintData();

      LWS::BVHNode3D* vert_tree = CreateBVHFromCurve(curves);
      Eigen::MatrixXd gradients;
      gradients.setZero(nVerts, 3);
      tpeSolver->FillGradientVectorBH(tree, gradients);
      Eigen::VectorXd x = gradients.col(0);
      std::cout << "Filled BH gradient vector" << std::endl;

      mult->TestAdmissibleMultiply(x);
      long tree_start = Utils::currentTimeMilliseconds();
      Eigen::VectorXd b_tree;
      b_tree.setZero(nVerts);
      mult->Multiply(x, b_tree);
      long tree_end = Utils::currentTimeMilliseconds();
      std::cout << "Tree time = " << (tree_end - tree_start) << " ms" << std::endl;

      Eigen::MatrixXd A;
      A.setZero(nVerts, nVerts);
      SobolevCurves::SobolevGramMatrix(curves, 2, 4, A);

      long mult_start = Utils::currentTimeMilliseconds();
      Eigen::VectorXd b_mat = A * x;
      long mult_end = Utils::currentTimeMilliseconds();
      std::cout << "Dense time = " << (mult_end - mult_start) << " ms" << std::endl;

      double norm_diff = (b_mat - b_tree).norm();
      double norm_mat = b_mat.norm();

      std::cout << "Reference norm = " << norm_mat << std::endl;
      std::cout << "Tree norm = " << b_tree.norm() << std::endl;

      std::cout << "Multiplication accuracy = " << 100 * (norm_diff / norm_mat) << " percent" << std::endl;
      delete mult;
      delete tree;
      delete vert_tree;
    }

    if (ImGui::Button("Test coarsen")) {
      int nVerts = curves->NumVertices();
      int logNumVerts = log2(nVerts) - 2;

      std::vector<MultigridOperator> ops(logNumVerts);
      std::vector<PolyCurveNetwork*> ps(logNumVerts);

      ps[0] = curves;

      for (int i = 0; i < logNumVerts - 1; i++) {
        std::cout << "Coarsening curve network of length " << ps[i]->NumVertices() << std::endl;
        ps[i+1] = ps[i]->Coarsen(ops[i]);
        DisplayCurves(ps[i + 1], "coarsened" + std::to_string(i));
      }

      for (int i = logNumVerts - 1; i > 0; i--) {
        Eigen::VectorXd coarsePos(ps[i]->NumVertices() * 3);
        coarsePos.setZero();
        MatrixIntoVectorX3(ps[i]->positions, coarsePos);

        Eigen::VectorXd finePos = ops[i - 1].prolong(coarsePos, ProlongationMode::Matrix3Only);

        Eigen::MatrixXd finePosMat(ps[i - 1]->NumVertices(), 3);
        finePosMat.setZero();
        VectorXdIntoMatrix(finePos, finePosMat);

        ps[i - 1]->positions = finePosMat;

        DisplayCurves(ps[i - 1], "prolonged" + std::to_string(i));
      }
    }

    if (ImGui::Button("Test solve")) {
      int nVerts = curves->NumVertices();
      LWS::BVHNode3D* tree = CreateBVHFromCurve(curves);
      Eigen::MatrixXd gradients;
      gradients.setZero(nVerts, 3);
      tpeSolver->FillGradientVectorBH(tree, gradients);
      delete tree;

      Eigen::VectorXd b;
      b.setZero(3 * nVerts);
      MatrixIntoVectorX3(gradients, b);

      using TestDomain = EdgeLengthNullProjectorDomain;
      TestDomain* domain = new TestDomain(curves, 2, 4, 0.5);

      Eigen::MatrixXd A = domain->GetFullMatrix();
      Eigen::VectorXd sol = domain->DirectSolve(b);

      std::cout << A << std::endl;
      std::cout << sol << std::endl;

      delete domain;
    }

    if (ImGui::Button("Test 3x saddle")) {
      // Get the TPE gradient as a test problem
      int nVerts = curves->NumVertices();
      int logNumVerts = log2(nVerts) - 4;
      std::cout << "Using " << logNumVerts << " levels" << std::endl;

      long setupStart = Utils::currentTimeMilliseconds();
      LWS::BVHNode3D* tree = CreateBVHFromCurve(curves);
      Eigen::MatrixXd gradients;
      gradients.setZero(nVerts, 3);
      tpeSolver->FillGradientVectorBH(tree, gradients);

      using TestDomain = EdgeLengthNullProjectorDomain;
      TestDomain* domain = new TestDomain(curves, 2, 4, 0.5);

      // Reshape the V x 3 matrix into a 3V x 1 vector
      Eigen::VectorXd gradientsLong(domain->NumRows());
      gradientsLong.setZero();
      MatrixIntoVectorX3(gradients, gradientsLong);

      EdgeLengthConstraint constraint(curves);
      Eigen::SparseMatrix<double> B;
      constraint.FillConstraintMatrix(B);

      Eigen::VectorXd gradientOrig = gradientsLong;

      if (domain->GetMode() == ProlongationMode::Matrix3AndProjector) {
          gradientsLong = curves->constraintProjector->ProjectToNullspace(gradientsLong);
      }

      long setupEnd = Utils::currentTimeMilliseconds();
      std::cout << "Setup time = " << (setupEnd - setupStart) << " ms" << std::endl;

      // Multigrid solve
      long multigridStart = Utils::currentTimeMilliseconds();
      MultigridHierarchy<TestDomain>* hierarchy = new MultigridHierarchy<TestDomain>(domain, logNumVerts);
      std::cout << "Created hierarchy" << std::endl;
      Eigen::VectorXd sol = hierarchy->VCycleSolve<MultigridHierarchy<TestDomain>::EigenCG>(gradientsLong);
      long multigridEnd = Utils::currentTimeMilliseconds();
      std::cout << "Multigrid time = " << (multigridEnd - multigridStart) << " ms" << std::endl;

      // std::cout << "Well-separated cluster time (total) = " << domain->tree->wellSepTime << " ms" << std::endl;
      // std::cout << "Ill-separated cluster time (total) = " << domain->tree->illSepTime << " ms" << std::endl;
      // std::cout << "Time spent on subtree traversals = " << domain->tree->traversalTime << " ms" << std::endl;

      // Direct solve
      long directStart = Utils::currentTimeMilliseconds();
      Eigen::VectorXd ref_sol = domain->DirectSolve(gradientOrig);
      long directEnd = Utils::currentTimeMilliseconds();
      std::cout << "Direct time = " << (directEnd - directStart) << " ms" << std::endl;
      Eigen::VectorXd diff = ref_sol - sol;

      Eigen::VectorXd constrValues = B * sol.block(0, 0, B.cols(), 1);
      Eigen::VectorXd ref_constrValues = B * ref_sol.block(0, 0, B.cols(), 1);

      std::cout << "Final constraint violation = " << constrValues.lpNorm<Eigen::Infinity>() << std::endl;
      std::cout << "Reference constraint violation = " << ref_constrValues.lpNorm<Eigen::Infinity>() << std::endl;

      Eigen::MatrixXd comp(sol.rows(), 2);
      comp.col(0) = sol;
      comp.col(1) = ref_sol;

      std::cout << "Reference norm = " << ref_sol.norm() << std::endl;
      std::cout << "Multigrid norm = " << sol.norm() << std::endl;
      std::cout << "Difference norm = " << diff.norm() << std::endl;
      std::cout << "Multigrid error from ground truth = " << 100 * (diff.norm() / ref_sol.norm()) << " percent" << std::endl;

      double dot = sol.normalized().dot(ref_sol.normalized());
      std::cout << "Dot product between directions = " << dot << std::endl;

      delete tree;
      delete hierarchy;
    }

    if (ImGui::Button("Output frame")) {
      outputFrame();
    }

    ImGui::Checkbox("Run TPE", &LWSOptions::runTPE);
    ImGui::Checkbox("Output frames", &LWSOptions::outputFrames);

    bool buttonStepTPE = ImGui::Button("Single TPE step");
    bool buttonPlotTPE = ImGui::Button("Plot TPE gradient");

    ImGui::Checkbox("Use multigrid", &LWSOptions::useMultigrid);

    if (LWSOptions::runTPE || buttonStepTPE) {
      if (LWSOptions::outputFrames && LWSOptions::frameNum == 0) {
        outputFrame();
      }
      bool good_step;
      if (LWSOptions::useMultigrid) {
        good_step = tpeSolver->StepSobolevLSIterative(0);
      }
      else {
        good_step = tpeSolver->StepSobolevLS(true);
      }
      UpdateCurvePositions();
      if (!good_step) {
        std::cout << "Stopped because line search could not take a step." << std::endl;
        LWSOptions::runTPE = false;
      }
      
      if (LWSOptions::outputFrames) {
        outputFrame();
      }
    }

    if (ImGui::Button("Curve to OBJ")) {
      std::cout << "TODO" << std::endl;
    }

    ImGui::Checkbox("Use Sobalev", &LWSOptions::useSobalev);

    double delta = 0.001;

    ImGui::End();
  }

  void LWSApp::centerLoopBarycenter(PolyCurveNetwork* curves) {
    Vector3 center = curves->Barycenter();
    int nVerts = curves->NumVertices();

    for (int i = 0; i < nVerts; i++) {
      CurveVertex* v = curves->GetVertex(i);
      v->SetPosition(v->Position() - center);
    }

    UpdateCurvePositions();
  }

  void LWSApp::initSolver() {
    if (!tpeSolver) {
      // Set up solver
      tpeSolver = new TPEFlowSolverSC(curves, 3, 6);
    }
  }

  void LWSApp::UpdateCurvePositions() {
    // Update the positions on the space curve
    polyscope::CurveNetwork* curveNetwork = polyscope::getCurveNetwork(surfaceName);
    size_t nVerts = curves->NumVertices();
    std::vector<glm::vec3> curve_vecs(nVerts);

    for (size_t i = 0; i < nVerts; i++)
    {
      CurveVertex* v_i = curves->GetVertex(i);
      Vector3 v = v_i->Position();
      curve_vecs[v_i->GlobalIndex()] = glm::vec3{v.x, v.y, v.z};
    }

    curveNetwork->updateNodePositions(curve_vecs);
    polyscope::requestRedraw();
  }

  void LWSApp::processFileOBJ(std::string filename) {
    if (curves) delete curves;
    std::cout << "Make curves for " << filename << std::endl;

    std::tie(mesh, geom) = loadMesh(filename);
    geom->requireVertexPositions();

    VertexData<size_t> indices = mesh->getVertexIndices();

    int nVerts = mesh->nVertices();
    std::vector<Vector3> all_positions(nVerts);
    std::vector<std::array<size_t, 2>> all_edges;

    for (Vertex v : mesh->vertices()) {
      size_t index = indices[v];
      all_positions[index] = geom->vertexPositions[v];
    }

    for (BoundaryLoop b : mesh->boundaryLoops()) {
      Halfedge he = b.halfedge().twin();
      Halfedge start = b.halfedge().twin();

      do {
        Vector3 v = geom->vertexPositions[he.vertex()];
        size_t index = indices[he.vertex()];
        size_t next_index = indices[he.next().vertex()];
        all_edges.push_back({index, next_index});
        
        he = he.next();
      }
      while (he != start);
      std::cout << "Processed boundary curve of length " << b.degree() << std::endl;
    }

    curves = new PolyCurveNetwork(all_positions, all_edges);

    surfaceName = polyscope::guessNiceNameFromPath(filename);
  }

  void LWSApp::DisplayCurves(PolyCurveNetwork* curves, std::string name) {

    std::vector<glm::vec3> nodes;
    std::vector<std::array<size_t, 2>> edges;
    
    for (int i = 0; i < curves->NumVertices(); i++) {
      // Add interior edges and vertices
      Vector3 p = curves->GetVertex(i)->Position();
      nodes.push_back(glm::vec3{p.x, p.y, p.z});
    }
    for (int i = 0; i < curves->NumEdges(); i++) {
      CurveEdge* e = curves->GetEdge(i);
      edges.push_back({(size_t)e->prevVert->GlobalIndex(), (size_t)e->nextVert->GlobalIndex()});
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

  std::cout << "Using Eigen version " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << std::endl;

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
