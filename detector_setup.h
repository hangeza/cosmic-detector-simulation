#pragma once

#include <vector>

#include "geometry_types.h"

/** @brief DetectorSetup - class for managing geometric objects of type ExtrudedObject
 * This class stores a list of ExtrudedObject detector objects for convenience.
 * An arbitrary number of detector objects with individual alignments can be supplied to the
 * constructor. Additionaly, a global rotation of the whole setup can be set.
 * @note The first detector of the vector which is supplied to the constructor is set as
 * reference detector. The assignment of the ref detector can be changed with the
 * DetectorSetup::ref_detector() method
*/
class DetectorSetup {
public:
    struct TriggerClass {
        unsigned number { 1u };
        bool invert { false };
        std::function<bool(bool, bool)> logic_function { std::logical_or<bool> {} };
    };
    DetectorSetup() = default;
    DetectorSetup(const DetectorSetup& other);
    DetectorSetup(DetectorSetup&& other);
    DetectorSetup(const std::vector<ExtrudedObject>& detectorlist, const ExtrudedObject& ref_volume = ExtrudedObject::invalid_volume());

    auto detectors() -> std::vector<ExtrudedObject>& { return m_detectors; }
    auto detectors() const -> const std::vector<ExtrudedObject>& { return m_detectors; }
    void add_detector(const ExtrudedObject& det);
    void set_ref_volume(const ExtrudedObject& ref_volume) { m_ref_volume = ref_volume; }
    auto ref_volume() const -> const ExtrudedObject& { return m_ref_volume; }
    auto bounding_box() const -> std::pair<Point, Point>;
    void rotate(const Vector& rot_axis, double rot_angle);
    void reset_rotation();

    auto intersection(const Line& path) const -> std::vector<LineSegment>;

private:
    auto get_largest_bounding_box() const -> std::pair<Point, Point>;
    std::vector<ExtrudedObject> m_detectors {};
    ExtrudedObject m_ref_volume {};
    std::string m_name {};
};
