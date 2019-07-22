#include "spacecurves/space_curve.h"

#include "polyscope/file_helpers.h"
#include "polyscope/combining_hash_functions.h"
#include "polyscope/gl/colors.h"
#include "polyscope/gl/materials/materials.h"
#include "polyscope/gl/shaders.h"
#include "polyscope/gl/shaders/cylinder_shaders.h"
#include "polyscope/gl/shaders/sphere_shaders.h"
#include "polyscope/pick.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_count_quantity.h"

#include "stb_image.h"

#include "imgui.h"

#include <unordered_map>
#include <utility>

#include <fstream>

namespace polyscope {
    
    // Initialize statics
    const std::string SpaceCurve::structureTypeName = "Space Curve";

    SpaceCurve::~SpaceCurve() { deleteProgram(); }

    void SpaceCurve::deleteProgram() {
        if (program != nullptr) {
            delete program;
            program = nullptr;
        }
    }

    void SpaceCurve::updatePoints(const std::vector<glm::vec3>& newPositions) {
        size_t nP = nPoints();
        for (size_t i = 0; i < nP; i++) {
            // Set the base points
            points[i] = newPositions[i];
            size_t i_next = (i + 1) % nP;
            // If this isn't a closed loop, skip the segment (n-1, 0)
            if (!isClosedLoop && i_next == 0) shifted_points[i_next] = points[i_next];
            // Set the shifted points, one index later
            else shifted_points[i_next] = newPositions[i];
        }

        program->setAttribute("a_position_tail", points);
        program->setAttribute("a_position_tip", shifted_points);
    }

    void SpaceCurve::draw() {
        if (!enabled) {
            return;
        }

        if (program == nullptr) {
            prepare();
        }

        setCurveUniforms(program);

        program->draw();

        // Draw the quantities
        for (auto x : quantities) {
            x.second->draw();
        }
    }
    
    std::string SpaceCurve::typeName() {
        return "Space Curve";
    }
    
    void SpaceCurve::prepare() {
        // It not quantity is coloring the surface, draw with a default color
        program = new gl::GLProgram(&gl::PASSTHRU_CYLINDER_VERT_SHADER, &gl::CYLINDER_GEOM_SHADER,
            &gl::CYLINDER_FRAG_SHADER, gl::DrawMode::Points);

        setMaterialForProgram(*program, "wax");

        // Populate draw buffers
        fillGeometryBuffers();
    }

    void SpaceCurve::preparePick() {
        // Request pick indices
        size_t pickCount = points.size();
        size_t pickStart = pick::requestPickBufferRange(this, pickCount);

        // Create a new pick program
        safeDelete(pickProgram);
        pickProgram = new gl::GLProgram(&gl::SPHERE_COLOR_VERT_SHADER, &gl::SPHERE_COLOR_BILLBOARD_GEOM_SHADER,
                                        &gl::SPHERE_COLOR_PLAIN_BILLBOARD_FRAG_SHADER, gl::DrawMode::Points);

        // Fill an index buffer
        std::vector<glm::vec3> pickColors;
        for (size_t i = pickStart; i < pickStart + pickCount; i++) {
            glm::vec3 val = pick::indToVec(i);
            pickColors.push_back(pick::indToVec(i));
        }

        // Store data in buffers
        pickProgram->setAttribute("a_position", points);
        pickProgram->setAttribute("a_color", pickColors);
    }

    void SpaceCurve::fillGeometryBuffers() {
        program->setAttribute("a_position_tail", points);
        program->setAttribute("a_position_tip", shifted_points);
    }

    // Helper to set uniforms
    void SpaceCurve::setCurveUniforms(gl::GLProgram* p) {

        glm::mat4 viewMat = getModelView();
        p->setUniform("u_modelView", glm::value_ptr(viewMat));

        glm::mat4 projMat = view::getCameraPerspectiveMatrix();
        p->setUniform("u_projMatrix", glm::value_ptr(projMat));

        p->setUniform("u_radius", pointRadius * state::lengthScale);
        p->setUniform("u_color", pointColor);
    }

    // Helper to set uniforms
    void SpaceCurve::setPickUniforms(gl::GLProgram* p) {

        glm::mat4 viewMat = getModelView();
        p->setUniform("u_modelView", glm::value_ptr(viewMat));

        glm::mat4 projMat = view::getCameraPerspectiveMatrix();
        p->setUniform("u_projMatrix", glm::value_ptr(projMat));
        p->setUniform("u_pointRadius", pointRadius);
    }

    void SpaceCurve::buildCustomOptionsUI() {}

    void SpaceCurve::buildPickUI(size_t localPickID) {
        ImGui::TextUnformatted(("#" + std::to_string(localPickID) + "  ").c_str());
        ImGui::SameLine();
        ImGui::TextUnformatted(to_string(points[localPickID]).c_str());

        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Indent(20.);

        // TODO: quantities?

        ImGui::Indent(-20.);
    }

    void SpaceCurve::buildCustomUI() {
        ImGui::PushID(name.c_str()); // ensure there are no conflicts with identically-named labels

        if (ImGui::TreeNode(name.c_str())) {

            // Print stats
            ImGui::Text("# points: %lld", static_cast<long long int>(points.size()));


            ImGui::Checkbox("Enabled", &enabled);
            ImGui::SameLine();
            ImGui::ColorEdit3("Point color", (float*)&pointColor, ImGuiColorEditFlags_NoInputs);
            ImGui::SameLine();

            // Options popup
            if (ImGui::Button("Options")) {
            ImGui::OpenPopup("OptionsPopup");
            }
            if (ImGui::BeginPopup("OptionsPopup")) {

            // Quantities
            if (ImGui::MenuItem("Clear Quantities")) {
                std::cout << "TODO: removeAllQuantities()" << std::endl;
            }
            if (ImGui::MenuItem("Write points to file")) writePointsToFile();


            ImGui::EndPopup();
            }

            ImGui::SliderFloat("Point Radius", &pointRadius, 0.0, .1, "%.5f", 3.);

            // Draw the quantities
            for (auto x : quantities) {
                x.second->drawUI();
            }

            ImGui::TreePop();
        }
        ImGui::PopID();
    }

    void SpaceCurve::drawPick() {
        if (!enabled) {
            return;
        }
        // Set uniforms
        //setPickUniforms(pickProgram);
        //pickProgram->draw();
    }

    double SpaceCurve::lengthScale() {
        // Measure length scale as twice the radius from the center of the bounding box
        auto bound = boundingBox();
        glm::vec3 center = 0.5f * (std::get<0>(bound) + std::get<1>(bound));

        double lengthScale = 0.0;
        for (glm::vec3& p : points) {
            lengthScale = std::max(lengthScale, (double)glm::length2(p - center));
        }

        return 2 * std::sqrt(lengthScale);
    }

    std::tuple<glm::vec3, glm::vec3> SpaceCurve::boundingBox() {

        glm::vec3 min = glm::vec3{1, 1, 1} * std::numeric_limits<float>::infinity();
        glm::vec3 max = -glm::vec3{1, 1, 1} * std::numeric_limits<float>::infinity();

        for (glm::vec3& rawP : points) {
            glm::vec3 p = glm::vec3(objectTransform * glm::vec4(rawP, 1.0));
            min = componentwiseMin(min, p);
            max = componentwiseMax(max, p);
        }

        return std::make_tuple(min, max);
    }

    void SpaceCurve::writePointsToFile(std::string filename) {
        if (filename == "") {
            filename = promptForFilename();
            if (filename == "") {
            return;
            }
        }

        std::cout << "Writing point cloud " << name << " to file " << filename << std::endl;

        std::ofstream outFile(filename);
        outFile << "#Polyscope point cloud " << name << std::endl;
        outFile << "#displayradius " << (pointRadius * state::lengthScale) << std::endl;

        for (size_t i = 0; i < points.size(); i++) {
            outFile << points[i] << std::endl;
        }

        outFile.close();
    }

    void SpaceCurve::writePointsToOBJ(std::string filename) {
        if (filename == "") {
            filename = promptForFilename();
            if (filename == "") {
            return;
            }
        }

        std::cout << "Writing OBJ of  " << name << " to file " << filename << std::endl;

        std::ofstream outFile(filename);

        for (size_t i = 0; i < points.size(); i++) {
            outFile << "v " << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
        }

        outFile << "f ";
        for (size_t i = 0; i < points.size(); i++) {
            outFile << (i + 1) << " ";
        }
        outFile << std::endl;


        outFile.close();
    }

    void SpaceCurve::addCurveQuantity(CurveQuantity* quantity) {
        std::shared_ptr<CurveQuantity> ptr;
        ptr.reset(quantity);
        addCurveQuantity(ptr);
    }

    void SpaceCurve::toggleIsClosed(bool isClosed) {
        isClosedLoop = isClosed;
        updatePoints(points);
    }

    void SpaceCurve::addCurveQuantity(std::shared_ptr<CurveQuantity> quantity) {
        // Delete old if in use
        bool wasEnabled = false;
        if (quantities.find(quantity->name) != quantities.end()) {
            wasEnabled = quantities[quantity->name]->isEnabled();
            removeQuantity(quantity->name);
        }

        // Store
        quantities[quantity->name] = quantity;

        // Re-enable the quantity if we're replacing an enabled quantity
        if (wasEnabled) {
            quantity->enable();
        }
    }

    std::shared_ptr<CurveQuantity> SpaceCurve::getCurveQuantity(std::string name, bool errorIfAbsent) {
        // Check if exists
        if (quantities.find(name) == quantities.end()) {
            if (errorIfAbsent) {
            polyscope::error("No quantity named " + name + " registered");
            }
            return nullptr;
        }
        return quantities[name];
    }

    void SpaceCurve::removeQuantity(std::string name) {
        if (quantities.find(name) == quantities.end()) {
            return;
        }

        std::shared_ptr<CurveQuantity> q = quantities[name];
        quantities.erase(name);
        /*
        if (activeSurfaceQuantity == q.get()) {
            clearActiveSurfaceQuantity();
        }
        */
    }

    void SpaceCurve::removeAllQuantities() {
        while (quantities.size() > 0) {
            removeQuantity(quantities.begin()->first);
        }
    }


    CurveQuantity::CurveQuantity(std::string name_, SpaceCurve* mesh_) : name(name_), parent(mesh_) {}
    CurveQuantity::~CurveQuantity() {}

    void CurveQuantity::buildVertexInfoGUI(size_t vInd) {}
    void CurveQuantity::buildFaceInfoGUI(size_t fInd) {}
    void CurveQuantity::buildEdgeInfoGUI(size_t eInd) {}
    void CurveQuantity::buildHalfedgeInfoGUI(size_t heInd) {}

    bool CurveQuantity::isEnabled() { return enabled; }

    void CurveQuantity::enable() { enabled = true; }
    void CurveQuantity::disable() { enabled = false; }

    void CurveQuantity::setEnabled(bool newEnabled) {
        if (enabled == false && newEnabled == true) {
            enable();
        } else if (enabled == true && newEnabled == false) {
            disable();
        }
    }
}
