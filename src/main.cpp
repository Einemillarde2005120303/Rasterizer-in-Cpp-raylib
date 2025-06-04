#include <raylib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <tuple>
#include <cmath>
#include <algorithm>

#define M_PI 3.14159265358979323846

class Point_2d {
    public:
    float x, y;

    Point_2d(float x = 0, float y = 0) {
        this->x = x;
        this->y = y;
    }
};

class Point_3d {
    public:
    float x, y, z;

    Point_3d(float x = 0, float y = 0, float z = 0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Point_3d operator+(const Point_3d& other) const {
        return Point_3d(x + other.x, y + other.y, z + other.z);
    }

    Point_3d operator-(const Point_3d& other) const {
        return Point_3d(x - other.x, y - other.y, z - other.z);
    }

    Point_3d operator*(float scalar) const {
        return Point_3d(x * scalar, y * scalar, z * scalar);
    }

    Point_3d cross(const Point_3d& other) const {
        return Point_3d(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }

    Point_3d operator/(float scalar) const {
        if (scalar == 0) {
            throw std::runtime_error("Division by zero error.");
        }
        return Point_3d(x / scalar, y / scalar, z / scalar);
    }

    float length() const {
        return sqrt(x * x + y * y + z * z);
    }

    Point_3d normalize() const {
        float len = this->length();
        if (len == 0) {
            throw std::runtime_error("Cannot normalize a zero-length vector.");
        }
        return Point_3d(x / len, y / len, z / len);
    }

    float dot(const Point_3d& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    bool operator==(const Point_3d& other) const {
        return abs(x - other.x) < 0.00001f && abs(y - other.y) < 0.00001f && abs(z - other.z) < 0.00001f;
    }

    bool operator!=(const Point_3d& other) const {
        return !(*this == other);
    }
};

class Utils {
    public: 
        static std::vector<std::string> split_str(const std::string& str, char delimiter);
        static std::string read_file(std::string path);
        static std::tuple<std::vector<std::array<float, 3>>, std::vector<std::array<float, 3>>, std::vector<std::array<float, 2>>, std::vector<std::array<std::array<int, 3>, 3>>> read_obj(std::string path);
        static Point_3d ray_equation(Point_3d p, Point_3d d, float t);
        static float ray_intersects_triangle(Point_3d a, Point_3d b, Point_3d c, Point_3d p, Point_3d d);
};

std::string Utils::read_file(std::string path) {
    std::ifstream file(path);
    if (!file) {
        std::cerr << "Error opening file: \"" << path << "\"\n";
        return "";
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

std::vector<std::string> Utils::split_str(const std::string& str, char delimiter) {
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        result.push_back(token);
    }

    return result;
}

// Foolproof this read_obj function (idk how tho)#
// MrurBo, if you're reading this, HELP.
std::tuple<std::vector<std::array<float, 3>>, std::vector<std::array<float, 3>>, std::vector<std::array<float, 2>>, std::vector<std::array<std::array<int, 3>, 3>>> Utils::read_obj(std::string path) {
    std::string file_str = Utils::read_file(path);
    if (file_str.empty()) return std::tuple<std::vector<std::array<float, 3>>, std::vector<std::array<float, 3>>, std::vector<std::array<float, 2>>, std::vector<std::array<std::array<int, 3>, 3>>>{};

    std::vector<std::string> file_arr = Utils::split_str(file_str, '\n');
    std::vector<std::array<float, 3>> v, vn;
    std::vector<std::array<float, 2>> vt;
    std::vector<std::array<std::array<int, 3>, 3>> f;

    for (const std::string& line : file_arr) {
        std::vector<std::string> comps = Utils::split_str(line, ' ');

        if (comps[0] == "v" && comps.size() >= 4) {
            v.push_back({std::stof(comps[1]), std::stof(comps[2]), std::stof(comps[3])});
        } else if (comps[0] == "vn" && comps.size() >= 4) {
            vn.push_back({std::stof(comps[1]), std::stof(comps[2]), std::stof(comps[3])});
        } else if (comps[0] == "vt" && comps.size() >= 3) {
            vt.push_back({std::stof(comps[1]), std::stof(comps[2])});
        } else if (comps[0] == "f" && comps.size() >= 4) {
            std::array<std::array<int, 3>, 3> face;
            for (int i = 0; i < 3; ++i) {
                std::vector<std::string> vertexData = Utils::split_str(comps[i + 1], '/');
                face[i] = {std::stoi(vertexData[0]), std::stoi(vertexData[1]), std::stoi(vertexData[2])};
            }
            f.push_back(face);
        }
    }

    return std::make_tuple(v, vn, vt, f);
}

Point_3d Utils::ray_equation(Point_3d p, Point_3d d, float t) {
    return p + (d * t);
}

float Utils::ray_intersects_triangle(Point_3d a, Point_3d b, Point_3d c, Point_3d p, Point_3d d) {
    d = d.normalize();

    Point_3d e1 = a - b;
    Point_3d e2 = b - c;
    Point_3d e3 = c - a;

    Point_3d n = e1.cross(e2);

    if (std::abs(n.dot(d)) < 1e-5) return -1;
    float dist = 0 - (n.dot(p) - n.dot(a)) / n.dot(d);

    Point_3d i = Utils::ray_equation(p, d, dist);

    Point_3d v0 = i - a;
    Point_3d v1 = i - b;
    Point_3d v2 = i - c;

    Point_3d m0 = e1.cross(v0).normalize();
    Point_3d m1 = e2.cross(v1).normalize();
    Point_3d m2 = e3.cross(v2).normalize();

    if (m0 == m1 && m1 == m2) {
        return dist;
    } else {
        return -1;
    }
}

class Cam {
    public:
    inline static float x = 0;
    inline static float y = 0;
    inline static float z = 0;
    inline static float roth = 0;
    inline static float rotv = 0;

    inline static const float movespeed = 2;
    inline static const float turnspeed = 2;
    inline static const float fov = 60;
    inline static const int screen_width = 960;
    inline static const int screen_height = 540;
    static std::vector<Point_3d> default_rays;
};

std::vector<Point_3d> Cam::default_rays = []() {
    std::vector<Point_3d> array = {};
    const float fov_increment = tan(Cam::fov / 2) * 2;
    for (size_t i = 0; i < Cam::screen_width; ++i) {
        for (size_t j = 0; j < Cam::screen_height; ++j) {
            array.push_back(Point_3d((i - round(Cam::screen_width / 2)) * fov_increment, (j - round(Cam::screen_height / 2)) * fov_increment, 1));
        }
    }
    return array;
}();

std::vector<Point_3d> calc_rays() {
    const Point_3d ihat = Point_3d(sin(Cam::roth + (M_PI / 2)), cos(Cam::roth + (M_PI / 2)), 0);
    const Point_3d jhat = Point_3d(sin(Cam::roth) * cos(Cam::rotv + (M_PI / 2)), cos(Cam::roth) * cos(Cam::rotv + (M_PI / 2)), sin(Cam::rotv + (M_PI / 2)));
    const Point_3d khat = Point_3d(sin(Cam::roth) * cos(Cam::rotv), cos(Cam::roth) * cos(Cam::rotv), sin(Cam::rotv));

    std::vector<Point_3d> rays(Cam::screen_width * Cam::screen_height);

    std::transform(Cam::default_rays.begin(), Cam::default_rays.end(), rays.begin(), [ihat, jhat, khat](Point_3d v) {
        return Point_3d(
            v.dot(Point_3d(ihat.x, jhat.x, khat.x)),
            v.dot(Point_3d(ihat.y, jhat.y, khat.y)),
            v.dot(Point_3d(ihat.z, jhat.z, khat.z))
        );
    });

    return rays;
}

int main() {
    InitWindow(Cam::screen_width, Cam::screen_height, "Rasteriser In C++ (raylib)");
    SetTargetFPS(60);


    while (!WindowShouldClose()) {
        // Update


        // Draw
        BeginDrawing();
        ClearBackground(BLACK);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}