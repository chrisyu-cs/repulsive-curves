#pragma once

#include <vector>

#include "polyscope/affine_remapper.h"
#include "polyscope/color_management.h"
#include "polyscope/gl/gl_utils.h"
#include "polyscope/polyscope.h"
#include "polyscope/standardize_data_array.h"
#include "polyscope/structure.h"

namespace polyscope {

    enum class CurveElement { VERTEX = 0, EDGE };

    inline std::string getCurveElementTypeName(CurveElement type) {
    switch (type) {
        case CurveElement::VERTEX:
            return "vertex";
        case CurveElement::EDGE:
            return "edge";
        }
        throw std::runtime_error("broken");
    }

    class SpaceCurve;

    // Data defined on a space curve
    class CurveQuantity {
        public:
        // Base constructor which sets the name
        CurveQuantity(std::string name, SpaceCurve* mesh);
        virtual ~CurveQuantity() = 0;

        // Draw the quantity on the surface Note: for many quantities (like scalars)
        // this does nothing, because drawing happens in the mesh draw(). However
        // others (ie vectors) need to be drawn.
        virtual void draw() = 0;

        // Draw the ImGUI ui elements
        virtual void drawUI() = 0;

        // Build GUI info about this element
        virtual void buildVertexInfoGUI(size_t vInd);
        virtual void buildFaceInfoGUI(size_t fInd);
        virtual void buildEdgeInfoGUI(size_t eInd);
        virtual void buildHalfedgeInfoGUI(size_t heInd);

        virtual void enable();
        virtual void disable();
        bool isEnabled();
        void setEnabled(bool newEnabled);

        // === Member variables ===
        const std::string name;
        SpaceCurve* const parent;

        protected:
        bool enabled = false; // should be set by enable() and disable()
    };

    class SpaceCurve : public QuantityStructure<SpaceCurve> {
        public:
        // === Member functions ===

        // Construct a new point cloud structure
            
        template <class T>
        SpaceCurve(std::string name, const T& points_, bool isClosed)
            : QuantityStructure<SpaceCurve>(name),
            points(standardizeVectorArray<glm::vec3, 3, T>(points_)),
            shifted_points(points.size())
        {
            isClosedLoop = isClosed;
            initialBaseColor = getNextUniqueColor();
            pointColor = initialBaseColor;

            // Copy the points into the shifted list
            size_t nP = points.size();
            for (size_t i = 0; i < nP; i++) {
                size_t i_next = (i + 1) % nP;
                // If this isn't a closed loop, skip the segment (n-1, 0)
                if (!isClosedLoop && i_next == 0) shifted_points[i_next] = points[i_next];
                else shifted_points[i_next] = points[i];
            }

            prepare();
            preparePick();
        }

        ~SpaceCurve();

        void updatePoints(const std::vector<glm::vec3>& newPositions);

        // Render the the structure on screen
        virtual void draw() override;

        // Do setup work related to drawing, including allocating openGL data
        void prepare();
        void preparePick();

        // Build the imgui display
        virtual void buildCustomUI() override;
        virtual void buildPickUI(size_t localPickID) override;
        virtual void buildCustomOptionsUI() override;

        virtual std::string typeName() override;
        
        void setCurveUniforms(gl::GLProgram* p);
        void setPickUniforms(gl::GLProgram* p);

        // Render for picking
        virtual void drawPick() override;

        // A characteristic length for the structure
        virtual double lengthScale() override;

        // Axis-aligned bounding box for the structure
        virtual std::tuple<glm::vec3, glm::vec3> boundingBox() override;

        // The points of this curve, in order
        std::vector<glm::vec3> points;
        std::vector<glm::vec3> shifted_points;
        size_t nPoints() const { return points.size(); }

        // Misc data
        bool enabled = true;
        static const std::string structureTypeName;

        // Small utilities
        void deleteProgram();
        void writePointsToFile(std::string filename = "");
        void writePointsToOBJ(std::string filename = "");

        // === Quantity-related

        // general form
        void addCurveQuantity(CurveQuantity* quantity); // will be deleted internally when appropriate
        void addCurveQuantity(std::shared_ptr<CurveQuantity> quantity);
        std::shared_ptr<CurveQuantity> getCurveQuantity(std::string name, bool errorIfAbsent = true);

        void removeQuantity(std::string name);
        void removeAllQuantities();

        // Vectors (expect vector array, inner type must be indexable with correct dimension)
        template <class T>
        void addVertexVectorQuantity(std::string name, const T& vectors, VectorType vectorType = VectorType::STANDARD);

        void toggleIsClosed(bool isClosed);

        private:
        // Quantities
        std::map<std::string, std::shared_ptr<CurveQuantity>> quantities;

        // Visualization parameters
        glm::vec3 initialBaseColor;
        glm::vec3 pointColor;
        float pointRadius = 0.005;

        void fillGeometryBuffers();
        bool isClosedLoop;
        
        // Drawing related things
        gl::GLProgram* program = nullptr;
        gl::GLProgram* pickProgram = nullptr;
    };

    template <class T>
    void registerSpaceCurve(std::string name, const T& points, bool isClosed, bool replaceIfPresent = true) {
        SpaceCurve* s = new SpaceCurve(name, points, isClosed);
        bool success = registerStructure(s);
        if (!success) delete s;
    }

    // Shorthand to get a point cloud from polyscope
    inline SpaceCurve* getSpaceCurve(std::string name = "") {
        return dynamic_cast<SpaceCurve*>(getStructure(SpaceCurve::structureTypeName, name));
    }
}

#include "space_curve_vector_quantity.h"
