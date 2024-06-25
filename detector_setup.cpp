#include "detector_setup.h"
#include "algebra_types.h"
#include "algebra_utils.h"
#include "geometry_types.h"

#include <iostream>

DetectorSetup::DetectorSetup(const std::vector<ExtrudedObject>& detectorlist, const ExtrudedObject& ref_volume)
    : m_detectors(detectorlist)
{
    if (ref_volume == ExtrudedObject::invalid_volume()) {
        // invalid ref volume...setting to 2*bounding box
        auto bounds { this->get_largest_bounding_box() };
        //auto bounds { this->bounding_box() };
        bounds.first *= 1.1;
        bounds.second *= 1.1;
        m_ref_volume = std::move(ExtrudedObject(bounds));
    } else {
        m_ref_volume = ref_volume;
    }
}

DetectorSetup::DetectorSetup(const DetectorSetup& other)
    : m_detectors(other.m_detectors)
    , m_ref_volume(other.m_ref_volume)
    , m_name(other.m_name)
{
}

void DetectorSetup::add_detector(const ExtrudedObject& det)
{
    m_detectors.push_back(ExtrudedObject { det });
    //m_detectors.emplace_back( ExtrudedObject { det } );
}

auto DetectorSetup::bounding_box() const -> std::pair<Point, Point>
{
    Vector min_coordinates {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()
    };
    Vector max_coordinates {
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest()
    };

    for (const auto detector : m_detectors) {
        auto bb = detector.bounding_box();
        auto minbound = bb.first;
        auto maxbound = bb.second;

        if (minbound[0] < min_coordinates[0]) {
            min_coordinates[0] = minbound[0];
        }
        if (maxbound[0] > max_coordinates[0]) {
            max_coordinates[0] = maxbound[0];
        }
        if (minbound[1] < min_coordinates[1]) {
            min_coordinates[1] = minbound[1];
        }
        if (maxbound[1] > max_coordinates[1]) {
            max_coordinates[1] = maxbound[1];
        }
        if (minbound[2] < min_coordinates[2]) {
            min_coordinates[2] = minbound[2];
        }
        if (maxbound[2] > max_coordinates[2]) {
            max_coordinates[2] = maxbound[2];
        }
    }

    return std::make_pair<Point, Point>(std::move(min_coordinates), std::move(max_coordinates));
}

auto DetectorSetup::get_largest_bounding_box() const -> std::pair<Point, Point>
{
    Vector min_coordinates {
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()
    };
    Vector max_coordinates {
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest()
    };

    auto current_rot_axis { R3::Base::X };
    while (current_rot_axis != R3::NullVec) {
        for (auto rot_angle { 0. }; rot_angle < 2 * pi(); rot_angle += pi() / 16) {
            auto rotated_setup { *this };
            rotated_setup.rotate(current_rot_axis, rot_angle);
            auto bb = rotated_setup.bounding_box();
            auto minbound = bb.first;
            auto maxbound = bb.second;

            if (minbound[0] < min_coordinates[0]) {
                min_coordinates[0] = minbound[0];
            }
            if (maxbound[0] > max_coordinates[0]) {
                max_coordinates[0] = maxbound[0];
            }
            if (minbound[1] < min_coordinates[1]) {
                min_coordinates[1] = minbound[1];
            }
            if (maxbound[1] > max_coordinates[1]) {
                max_coordinates[1] = maxbound[1];
            }
            if (minbound[2] < min_coordinates[2]) {
                min_coordinates[2] = minbound[2];
            }
            if (maxbound[2] > max_coordinates[2]) {
                max_coordinates[2] = maxbound[2];
            }
        }
        if (current_rot_axis == R3::Base::X) {
            current_rot_axis = R3::Base::Y;
        } else if (current_rot_axis == R3::Base::Y) {
            current_rot_axis = R3::Base::Z;
        } else
            current_rot_axis = R3::NullVec;
    }

    return std::make_pair<Point, Point>(std::move(min_coordinates), std::move(max_coordinates));
}

void DetectorSetup::rotate(const Vector& rot_axis, double rot_angle)
{
    if (inEpsilon(rot_angle))
        return;
    for (auto& detector : m_detectors) {
        //Point pos { detector.position() };
        //pos = ::rotate(pos, rot_axis, rot_angle);
        //detector.set_position(pos);
        detector.add_rotation(rot_axis, rot_angle);
    }
}
