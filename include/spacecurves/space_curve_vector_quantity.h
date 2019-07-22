#pragma once

#include "space_curve.h"

namespace polyscope {
    class CurveVectorQuantity : public CurveQuantity {
        public:
        CurveVectorQuantity(std::string name, SpaceCurve* mesh_, CurveElement definedOn_,
                                VectorType vectorType_ = VectorType::STANDARD);

        virtual ~CurveVectorQuantity() override;

        virtual void draw() override;
        virtual void drawUI() override;

        // Allow children to append to the UI
        virtual void drawSubUI();

        // Do work shared between all constructors
        void finishConstructing();

        // === Members
        const VectorType vectorType;
        std::vector<glm::vec3> vectorRoots;
        std::vector<glm::vec3> vectors;
        float lengthMult; // longest vector will be this fraction of lengthScale (if not ambient)
        float radiusMult; // radius is this fraction of lengthScale
        glm::vec3 vectorColor;
        CurveElement definedOn;

        // The map that takes values to [0,1] for drawing
        AffineRemapper<glm::vec3> mapper;

        void writeToFile(std::string filename = "");

        // GL things
        void prepare();
        gl::GLProgram* program = nullptr;
    };


    class CurveVertexVectorQuantity : public CurveVectorQuantity {
        public:
        CurveVertexVectorQuantity(std::string name, std::vector<glm::vec3> vectors_, SpaceCurve* mesh_,
                                    VectorType vectorType_ = VectorType::STANDARD);

        std::vector<glm::vec3> vectorField;
        virtual void buildVertexInfoGUI(size_t vInd) override;
    };
    
    template <class T>
    void SpaceCurve::addVertexVectorQuantity(std::string name, const T& vectors, VectorType vectorType) {

        validateSize(vectors, nPoints(), "vertex vector quantity " + name);

        std::shared_ptr<CurveVectorQuantity> q = std::make_shared<CurveVertexVectorQuantity>(
            name, standardizeVectorArray<glm::vec3, T, 3>(vectors), this, vectorType);
        addCurveQuantity(q);
    }
}
