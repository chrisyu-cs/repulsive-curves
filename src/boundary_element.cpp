#include "boundary_element.h"

namespace LWS {

    using namespace geometrycentral;
    
    double fs_greens_2d(Vector3 center, Vector3 pt) {
        double r = norm(center - pt);
        return (-1. / (2 * M_PI)) * log(r);
    }

    Vector3 fs_greens_grad_2d(Vector3 center, Vector3 pt) {
        // Gradient of the length r (wrt pt) is pointing outward from the center
        Vector3 unit = pt - center;
        double r = norm(unit);
        unit.normalize();
        // Gradient of log(r) = 1/r * (grad r)
        double p2r = -1 / (2 * M_PI * r);
        return p2r * unit;
    }

    BoundaryNode::BoundaryNode(int id, Vector3 pos, Vector3 nor) {
        nodeID = id;
        position = pos;
        normal = nor;
    }

    double BoundaryNode::DualLength() {
        return (element1->length() + element2->length()) / 2;
    }


    BoundaryElement::BoundaryElement(int id, BoundaryNode *n1, BoundaryNode *n2) {
        elementID = id;
        node1 = n1;
        node2 = n2;
    }

    double BoundaryElement::gauss_abscissa[4] = {-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526};
    double BoundaryElement::gauss_weights[4] = { 0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538};

    Vector3 BoundaryElement::position(double t) {
        return (1 - t) * node1->position + t * node2->position;
    }

    Vector3 BoundaryElement::normal(double t) {
        // TODO: spherical interpolation?
        Vector3 n = (1 - t) * node1->normal + t * node2->normal;
        n.normalize();
        return n;
    }

    double BoundaryElement::basis_fn(BoundaryNode* node, double t) {
        if (node == node1) {
            return 1 - t;
        }
        else if (node == node2) {
            return t;
        }
        else return 0;
    }

    double BoundaryElement::integrand_A(Vector3 collocation_pt, BoundaryNode* node, double t) {
        // Assumes that t is in [0, 1]
        Vector3 eval_pos = position(t);
        Vector3 eval_normal = normal(t);
        Vector3 greens_grad = fs_greens_grad_2d(collocation_pt, eval_pos);
        double greens_nderiv = dot(greens_grad, eval_normal);
        double basis = basis_fn(node, t);
        return basis * greens_nderiv;
    }

    double BoundaryElement::integrand_A_11(Vector3 collocation_pt, BoundaryNode* node, double t) {
        // Rescale t from [-1, 1] to [0, 1].
        double scaled = (t + 1) / 2;
        return integrand_A(collocation_pt, node, scaled);
    }

    double BoundaryElement::CoeffIntA_ia(Vector3 collocation_pt, BoundaryNode* node) {
        double sum = 0;
        for (int i = 0; i < 4; i++) {
            double t = BoundaryElement::gauss_abscissa[i];
            double w = BoundaryElement::gauss_weights[i];
            sum += w * integrand_A_11(collocation_pt, node, t);
        }
        return (length() / 2) * sum;
    }

    double BoundaryElement::length() {
        return norm(node1->position - node2->position);
    }

    double BoundaryElement::CoeffIntA_ia(BoundaryNode* collocation_pt, BoundaryNode* node) {
        return CoeffIntA_ia(collocation_pt->position, node);
    }

    double BoundaryElement::integrand_B(Vector3 collocation_pt, BoundaryNode* node, double t) {
        // Assumes that t is in [0, 1]
        Vector3 eval_pos = position(t);
        double greens = fs_greens_2d(collocation_pt, eval_pos);
        double basis = basis_fn(node, t);
        return basis * greens;
    }

    double BoundaryElement::integrand_B_11(Vector3 collocation_pt, BoundaryNode* node, double t) {
        // Rescale t from [-1, 1] to [0, 1].
        double scaled = (t + 1) / 2;
        return integrand_B(collocation_pt, node, scaled);
    }

    double BoundaryElement::CoeffIntB_ia(Vector3 collocation_pt, BoundaryNode* node) {
        double sum = 0;
        for (int i = 0; i < 4; i++) {
            double t = BoundaryElement::gauss_abscissa[i];
            double w = BoundaryElement::gauss_weights[i];
            sum += w * integrand_B_11(collocation_pt, node, t);
        }
        return (length() / 2) * sum;
    }

    double BoundaryElement::CoeffIntB_ia(BoundaryNode* collocation_pt, BoundaryNode* node) {
        return CoeffIntB_ia(collocation_pt->position, node);
    }
}