#include "geometry_types.h"
#include "algebra_utils.h"
#include "matrix.h"

#include <iostream>
#include <limits>
#include <algorithm>


bool operator==(const ExtrudedObject& lhs, const ExtrudedObject& rhs)
{
    bool equal { false };
    if (lhs.m_vertices.size() != rhs.m_vertices.size()) return false;
    if (lhs.m_thickness != rhs.m_thickness) return false;
    if (!isFuzzySame(lhs.m_position, rhs.m_position)) return false;
    if (!std::equal(lhs.m_vertices.cbegin(), lhs.m_vertices.cend(), rhs.m_vertices.cbegin(), [](const Point& a, const Point& b) { return isFuzzySame(a,b);})) return false;
    if (lhs.m_rotation_matrix != rhs.m_rotation_matrix) return false;
    return true;
}

bool operator!=(const ExtrudedObject& lhs, const ExtrudedObject& rhs) 
{
    return !(lhs == rhs); 
}

auto Line::operator()(double t) const -> Point
{
    return p + t * q;
}

auto Line::distance(const Point& point) const -> double
{
    return norm((point - p) - q * q * (point - p));
}

auto Line::generate(Point p0, double theta, double phi) -> Line
{
    return {
        p0,
        { std::sin(theta) * std::cos(phi),
            std::sin(theta) * std::sin(phi),
            std::cos(theta) }
    };
}

void Line::rotate(const Vector& rot_axis, double rot_angle)
{
    p = { ::rotate(p, rot_axis, rot_angle) };
    q = { ::rotate(q, rot_axis, rot_angle) };
}

auto LineSegment::length() const -> double
{
    return norm(line(t_start) - line(t_end));
}

auto Plane::distance(const Point& point) const -> double
{
    return ((point - p) * normal).sum();
}

auto Plane::intersection(const Line& line) const -> Point
{
    if (isFuzzySame(p, line.p)) {
        // reference points are identical
        return p;
    }
    double v1 { ((p - line.p) * normal).sum() };
    double v2 { (line.q * normal).sum() };

    if (inEpsilon(v1)) {
        //std::cerr << "Error in Plane::intersection(const Line&): line is contained entirely in plane v1="<<v1<<" v2="<<v2<<" this->p="<<p<<" line.p="<<line.p<<" this->normal()="<<normal()<<" (p-line.p)="<<p-line.p<<" (p-line.p)*normal="<<(p-line.p)*normal()<<"\n";
        return p;
    }

    if (inEpsilon(v2)) {
        //std::cerr << "Error in Plane::intersection(const Line&): line is not intersecting the plane\n";
        throw no_intersection("Plane::intersection(const Line&): line is not intersecting the plane");
        return {};
    }
    double t { v1 / v2 };
    return line(t);
}

void Plane::rotate(const Vector& rot_axis, double rot_angle)
{
    ::rotate(normal, rot_axis, rot_angle);
}

void Plane::rotate(const matrix2d<double>& rot_matrix)
{
    //    std::cout<<"orig plane: p="<<p<<" norm="<<normal<<"\n";
    Plane new_plane { rot_matrix * p, rot_matrix * normal };
    //    std::cout<<"new plane:  p="<<new_plane.p<<" norm="<<new_plane.normal<<"\n";
    *this = new_plane;
}

ExtrudedObject::ExtrudedObject()
    : m_vertices( {} )
    , m_position( {} )
    , m_thickness( 0. )
    , m_planes( {} )
    , m_rotation_matrix(R3::Identity)
{
}

ExtrudedObject::ExtrudedObject(const ExtrudedObject& other)
    : m_vertices(other.m_vertices)
    , m_position(other.m_position)
    , m_thickness(other.m_thickness)
    , m_planes(other.m_planes)
    , m_rotation_matrix(other.m_rotation_matrix)
{
}


ExtrudedObject& ExtrudedObject::operator=(const ExtrudedObject& other) 
{
    m_vertices = other.m_vertices;
    m_position = other.m_position;
    m_thickness = other.m_thickness;
    m_planes = other.m_planes;
    m_rotation_matrix = other.m_rotation_matrix;
    return *this;
}

ExtrudedObject::ExtrudedObject(ExtrudedObject&& other)
    : m_vertices(std::move(other.m_vertices))
    , m_position(std::move(other.m_position))
    , m_thickness(std::move(other.m_thickness))
    , m_planes(std::move(other.m_planes))
    , m_rotation_matrix(std::move(other.m_rotation_matrix))
{
//    std::cout<<"ExtrudedObject::ExtrudedObject(ExtrudedObject&&) 1"<<std::endl;
    //m_planes = getPlanes();
//    std::cout<<"ExtrudedObject::ExtrudedObject(ExtrudedObject&&) 2"<<std::endl;
}

ExtrudedObject::ExtrudedObject(const std::vector<Point>& vertices, const Point& position, double thickness)
    : m_vertices(vertices)
    , m_position(position)
    , m_thickness(thickness)
{
    m_planes = getPlanes();
}

ExtrudedObject::ExtrudedObject(const Point& position, double radius, double thickness, std::size_t nr_vertices)
    : m_position(position)
    , m_thickness(thickness)
{
    for (std::size_t n { 0 }; n < nr_vertices; ++n) {
        const double angle { twopi() * n / nr_vertices };
        Point p {
            m_position[0] + radius * std::cos(angle),
            m_position[1] + radius * std::sin(angle)
        };
        m_vertices.emplace_back(std::move(p));
    }
    m_planes = getPlanes();
}

ExtrudedObject::ExtrudedObject(const std::pair<Point, Point>& bounding_box)
{
    const auto minbound { bounding_box.first };
    const auto maxbound { bounding_box.second };
    
    m_vertices.push_back( { minbound[0], minbound[1] } );
    m_vertices.push_back( { maxbound[0], minbound[1] } );
    m_vertices.push_back( { maxbound[0], maxbound[1] } );
    m_vertices.push_back( { minbound[0], maxbound[1] } );
    
    m_position = { (minbound + maxbound)/2 };
    m_position[2] = minbound[2];
    m_thickness = maxbound[2] - minbound[2];
    
    m_planes = getPlanes();
}

const auto ExtrudedObject::position() const -> Point
{
    return m_position;
}

void ExtrudedObject::set_position(const Point& new_pos)
{
    m_position = new_pos;
    //m_planes = getPlanes();
}

auto ExtrudedObject::thickness() const -> double
{
    return m_thickness;
}

auto ExtrudedObject::contains(const Point& point) const -> bool
{
    for (const auto& plane : m_planes) {
        double dist { plane.distance(point) };
        if (dist > DEFAULT_EPSILON)
            return false;
    }
    return true;
}

auto ExtrudedObject::intersection(const Line& path) const -> LineSegment
{
    std::vector<Point> hitpoints {};
    for (const auto& plane : m_planes) {
        Point hitpoint {};
        try {
            hitpoint = { plane.intersection(path) };
        } catch (std::exception& e) {
            continue;
        }
        hitpoints.push_back(std::move(hitpoint));
    }

    if (hitpoints.empty()) {
        // no intersections with volume found
        return LineSegment {};
    }
    auto it = hitpoints.begin();
    while (it != hitpoints.end()) {
        if (contains(*it)) {
            ++it;
        } else
            it = hitpoints.erase(it);
    }
    if (hitpoints.size() > 1) {
        it = hitpoints.begin();
        while (it != hitpoints.end()) {
            auto it2 = std::next(it);
            while (it2 != hitpoints.end()) {
                if (isFuzzySame(*it, *it2)) {
                    it2 = hitpoints.erase(it2);
                } else
                    ++it2;
            }
            ++it;
        }
    }
    if (hitpoints.size() == 2) {
        return LineSegment {
            { hitpoints[0], hitpoints[1] - hitpoints[0] },
            0., 1.
        };
    } else if (hitpoints.size() > 2) {
        std::cerr << "ExtrudedObject::intersection(const Line&): strange nr of intersection points: " << hitpoints.size() << "\n";
    }
    return LineSegment {};
}

auto ExtrudedObject::getPlanes() const -> std::vector<Plane>
{
    std::vector<Plane> planes {};
    if (m_vertices.size() < 3) {
        throw std::runtime_error("Error in ExtrudedObject::getPlanes(): insufficient number of vertices (" + std::to_string(m_vertices.size()) + ")!");
        std::cerr << "Error in ExtrudedObject::getPlanes(): insufficient number of vertices (" << m_vertices.size() << ")!\n";
        return planes;
    }
    for (auto vertex { m_vertices.begin() };
         vertex != std::prev(m_vertices.end());
         ++vertex) {
        Point p0 { Point { (*vertex)[0], (*vertex)[1], 0. } };
        Point p1 { Point { (*std::next(vertex))[0], (*std::next(vertex))[1], 0. } };
        Point p2 { Point { p0[0], p0[1], m_thickness } };
        Plane plane { p0, cross_product(p1 - p0, p2 - p0) };
        plane.rotate(m_rotation_matrix);
        plane.p += m_position;
        planes.push_back(std::move(plane));
    }
    Point p0 { Point { (*std::prev(m_vertices.end()))[0], (*std::prev(m_vertices.end()))[1], 0. } };
    Point p1 { Point { (*m_vertices.begin())[0], (*m_vertices.begin())[1], 0. } };
    Point p2 { Point { p0[0], p0[1], m_thickness } };
    Plane plane { p0, cross_product(p1 - p0, p2 - p0) };
    plane.rotate(m_rotation_matrix);
    plane.p += m_position;
    planes.push_back(std::move(plane));

    Plane bot_plane { R3::NullVec, { 0., 0., -1. } };
    bot_plane.rotate(m_rotation_matrix);
    bot_plane.p += m_position;
    planes.push_back(std::move(bot_plane));

    Plane top_plane { { 0., 0., m_thickness }, { 0., 0., 1. } };
    top_plane.rotate(m_rotation_matrix);
    top_plane.p += m_position;
    planes.push_back(std::move(top_plane));

    return planes;
}

auto ExtrudedObject::get_vertices() const -> std::vector<Point>
{
    std::vector<Point> points {};
    if (m_vertices.size() < 3) {
        throw std::runtime_error("Error in ExtrudedObject::getVertices(): insufficient number of vertices (" + std::to_string(m_vertices.size()) + ")!");
        std::cerr << "Error in ExtrudedObject::getVertices(): insufficient number of vertices (" << m_vertices.size() << ")!\n";
        return points;
    }
    for (const auto vertex : m_vertices) {
        Point p1 { Point { vertex[0], vertex[1], 0. } };
        Point p2 { Point { vertex[0], vertex[1], m_thickness } };
        p1 = m_rotation_matrix * p1;
        p2 = m_rotation_matrix * p2;
        p1 += m_position;
        p2 += m_position;
        points.push_back(std::move(p1));
        points.push_back(std::move(p2));
    }
    return points;
}

auto ExtrudedObject::bounding_box() const -> std::pair<Point, Point>
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

    auto vertices { get_vertices() };
    
    for (const auto vertex : vertices) {
        if (vertex[0] < min_coordinates[0]) {
            min_coordinates[0] = vertex[0];
        } else if (vertex[0] > max_coordinates[0]) {
            max_coordinates[0] = vertex[0];
        }
        if (vertex[1] < min_coordinates[1]) {
            min_coordinates[1] = vertex[1];
        } else if (vertex[1] > max_coordinates[1]) {
            max_coordinates[1] = vertex[1];
        }
        if (vertex[2] < min_coordinates[2]) {
            min_coordinates[2] = vertex[2];
        } else if (vertex[2] > max_coordinates[2]) {
            max_coordinates[2] = vertex[2];
        }
    }
    
    if (min_coordinates[0] > max_coordinates[0])
        std::swap(min_coordinates[0], max_coordinates[0]);
    if (min_coordinates[1] > max_coordinates[1])
        std::swap(min_coordinates[1], max_coordinates[1]);
    if (min_coordinates[2] > max_coordinates[2])
        std::swap(min_coordinates[2], max_coordinates[2]);
    return std::make_pair<Point, Point>(std::move(min_coordinates), std::move(max_coordinates));
}

void ExtrudedObject::add_rotation(const Vector& rot_axis, double rot_angle)
{
//    std::cout<<"matrix before rot:\n";
//    std::cout<<m_rotation_matrix;
    Point pos { position() };
    pos = ::rotate(pos, rot_axis, rot_angle);
    set_position(pos);
    matrix2d<double> K { 3,
        { 0., -rot_axis[2], rot_axis[1],
            rot_axis[2], 0., -rot_axis[0],
            -rot_axis[1], rot_axis[0], 0. } };
    matrix2d R { R3::Identity + std::sin(rot_angle) * K + (1. - std::cos(rot_angle)) * (K * K) };
    m_rotation_matrix = R * m_rotation_matrix;
//    std::cout<<"matrix after rot:\n";
//    std::cout<<m_rotation_matrix;
    m_planes = getPlanes();
}

auto ExtrudedObject::get_rotation_matrix() -> const matrix2d<double>&
{
    return m_rotation_matrix;
}

void ExtrudedObject::reset_rotation_matrix()
{
    m_rotation_matrix = R3::Identity;
}
