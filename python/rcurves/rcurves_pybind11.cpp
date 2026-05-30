#include <pybind11/pybind11.h>

#include "lws_app.h"
#include "geometrycentral/utilities/vector3.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

namespace LWS {
class PyLWSApp : public LWSApp {
public:
    // inherit the constructors
    using LWSApp::LWSApp;

    void AddMeshObstacle(const std::string& objName,
                         const Vector3& center,
                         double p,
                         double weight) override {
        PYBIND11_OVERLOAD(
            void,
            LWSApp,
            AddMeshObstacle,
            objName,
            center,
            p,
            weight
            );
    }

    void AddPlaneObstacle(const Vector3& center,
                          const Vector3& normal,
                          double p,
                          double weight) override {
        PYBIND11_OVERLOAD(
            void,
            LWSApp,
            AddPlaneObstacle,
            center,
            normal,
            p,
            weight
            );
    }

    void AddSphereObstacle(const Vector3& center, double radius) override {
        PYBIND11_OVERLOAD(
            void,
            LWSApp,
            AddSphereObstacle,
            center,
            radius
            );
    }
};
}


PYBIND11_MODULE(rcurves_python, m) {
    m.doc() = "Python bindings for rcurves";

    py::class_<geometrycentral::Vector3>(m, "RCVector3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &geometrycentral::Vector3::x)
        .def_readwrite("y", &geometrycentral::Vector3::y)
        .def_readwrite("z", &geometrycentral::Vector3::z)
        .def("__repr__", [](const geometrycentral::Vector3& v) {
            return "RCVector3(" + std::to_string(v.x) + ", " +
                   std::to_string(v.y) + ", " +
                   std::to_string(v.z) + ")";
        });

    py::class_<LWS::LWSApp, LWS::PyLWSApp>(m, "LWSApp")
        .def(py::init<>())
        .def("AddMeshObstacle",
             &LWS::LWSApp::AddMeshObstacle,
             py::arg("objName"),
             py::arg("center"),
             py::arg("p"),
             py::arg("weight"))
        .def("AddPlaneObstacle",
             &LWS::LWSApp::AddPlaneObstacle,
             py::arg("center"),
             py::arg("normal"),
             py::arg("p"),
             py::arg("weight"))
        .def("AddSphereObstacle",
             &LWS::LWSApp::AddSphereObstacle,
             py::arg("center"),
             py::arg("radius"));

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}