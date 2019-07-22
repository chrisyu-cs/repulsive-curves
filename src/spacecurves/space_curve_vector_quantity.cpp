#include "spacecurves/space_curve_vector_quantity.h"

#include "polyscope/gl/materials/materials.h"
#include "polyscope/gl/shaders.h"
#include "polyscope/gl/shaders/vector_shaders.h"
#include "polyscope/utilities.h"
#include "polyscope/file_helpers.h"

#include <fstream>

namespace polyscope {

    CurveVectorQuantity::CurveVectorQuantity(std::string name, SpaceCurve* mesh_, CurveElement definedOn_, VectorType vectorType_)
    : CurveQuantity(name, mesh_), vectorType(vectorType_), definedOn(definedOn_) {

    // Don't forget to call finishConstructing() in children classes!
    }

    CurveVectorQuantity::~CurveVectorQuantity() {
        safeDelete(program);
    }

    void CurveVectorQuantity::finishConstructing() {
        // Create a mapper (default mapper is identity)
        if (vectorType == VectorType::AMBIENT) {
            mapper.setMinMax(vectors);
        } else {
            mapper = AffineRemapper<glm::vec3>(vectors, DataType::MAGNITUDE);
        }

        // Default viz settings
        if (vectorType != VectorType::AMBIENT) {
            lengthMult = .02;
        } else {
            lengthMult = 1.0;
        }
        radiusMult = .0005;
        vectorColor = getNextUniqueColor();
    }

    void CurveVectorQuantity::draw() {
        if (!enabled) return;

        if (program == nullptr) prepare();

        // Set uniforms
        glm::mat4 viewMat = parent->getModelView();
        program->setUniform("u_modelView", glm::value_ptr(viewMat));

        glm::mat4 projMat = view::getCameraPerspectiveMatrix();
        program->setUniform("u_projMatrix", glm::value_ptr(projMat));

        program->setUniform("u_radius", radiusMult * state::lengthScale);
        program->setUniform("u_color", vectorColor);

        if (vectorType == VectorType::AMBIENT) {
            program->setUniform("u_lengthMult", 1.0);
        } else {
            program->setUniform("u_lengthMult", lengthMult * state::lengthScale);
        }

        program->draw();
    }

    void CurveVectorQuantity::prepare() {
        program = new gl::GLProgram(&gl::PASSTHRU_VECTOR_VERT_SHADER, &gl::VECTOR_GEOM_SHADER, &gl::SHINY_VECTOR_FRAG_SHADER,
                                    gl::DrawMode::Points);

        // Fill buffers
        std::vector<glm::vec3> mappedVectors;
        for (glm::vec3& v : vectors) {
            mappedVectors.push_back(mapper.map(v));
        }

        program->setAttribute("a_vector", mappedVectors);
        program->setAttribute("a_position", vectorRoots);

        setMaterialForProgram(*program, "wax");
    }

    void CurveVectorQuantity::drawUI() {


    if (ImGui::TreeNode((name + " (" + getCurveElementTypeName(definedOn) + " vector)").c_str())) {
        ImGui::Checkbox("Enabled", &enabled);
        ImGui::SameLine();
        ImGui::ColorEdit3("Color", (float*)&vectorColor, ImGuiColorEditFlags_NoInputs);
        ImGui::SameLine();


        // === Options popup
        if (ImGui::Button("Options")) {
        ImGui::OpenPopup("OptionsPopup");
        }
        if (ImGui::BeginPopup("OptionsPopup")) {
        if (ImGui::MenuItem("Write to file")) writeToFile();
        ImGui::EndPopup();
        }


        // Only get to set length for non-ambient vectors
        if (vectorType != VectorType::AMBIENT) {
        ImGui::SliderFloat("Length", &lengthMult, 0.0, .1, "%.5f", 3.);
        }

        ImGui::SliderFloat("Radius", &radiusMult, 0.0, .1, "%.5f", 3.);

        { // Draw max and min magnitude
        ImGui::TextUnformatted(mapper.printBounds().c_str());
        }

        drawSubUI();

        ImGui::TreePop();
    }
    }

    void CurveVectorQuantity::drawSubUI() {}

    void CurveVectorQuantity::writeToFile(std::string filename) {

    if (filename == "") {
        filename = promptForFilename();
        if (filename == "") {
        return;
        }
    }

    std::cout << "Writing surface curve quantity " << name << " to file " << filename << std::endl;

    std::ofstream outFile(filename);
    outFile << "#Vectors written by polyscope from Surface Curve Quantity " << name << std::endl;
    outFile << "#displayradius " << (radiusMult * state::lengthScale) << std::endl;
    outFile << "#displaylength " << (lengthMult * state::lengthScale) << std::endl;

    for (size_t i = 0; i < vectors.size(); i++) {
        if (glm::length(vectors[i]) > 0) {
        outFile << vectorRoots[i] << " " << vectors[i] << std::endl;
        }
    }

    outFile.close();
    }

    // ========================================================
    // ==========           Vertex Vector            ==========
    // ========================================================

    CurveVertexVectorQuantity::CurveVertexVectorQuantity(std::string name, std::vector<glm::vec3> vectors_,
                                                            SpaceCurve* mesh_, VectorType vectorType_)
        : CurveVectorQuantity(name, mesh_, CurveElement::VERTEX, vectorType_), vectorField(vectors_) {

        for (size_t v = 0; v < parent->nPoints(); v++) {
            vectorRoots.push_back(parent->points[v]);
            vectors.push_back(vectorField[v]);
        }

        finishConstructing();
    }

    void CurveVertexVectorQuantity::buildVertexInfoGUI(size_t iV) {
        ImGui::TextUnformatted(name.c_str());
        ImGui::NextColumn();

        std::stringstream buffer;
        buffer << vectorField[iV];
        ImGui::TextUnformatted(buffer.str().c_str());

        ImGui::NextColumn();
        ImGui::NextColumn();
        ImGui::Text("magnitude: %g", glm::length(vectorField[iV]));
        ImGui::NextColumn();
    }

}
