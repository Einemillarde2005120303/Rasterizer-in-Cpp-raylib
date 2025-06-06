#include <raylib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <random>

#define PI 3.14159265358979323846f

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
            return Point_3d(
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x
            );
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
                return Point_3d(0, 0, 0);
            }
            return *this / len;
        }

        Point_3d transform(const Point_3d& ihat, const Point_3d& jhat, const Point_3d& khat) const {
            return Point_3d(
                x * ihat.x + y * jhat.x + z * khat.x,
                x * ihat.y + y * jhat.y + z * khat.y,
                x * ihat.z + y * jhat.z + z * khat.z
            );
        }

        float dot(const Point_3d& other) const {
            return x * other.x + y * other.y + z * other.z;
        }

        bool operator==(const Point_3d& other) const {
            return abs(x - other.x) <= 1e-6 && abs(y - other.y) <= 1e-6 && abs(z - other.z) <= 1e-6;
        }

        bool operator!=(const Point_3d& other) const {
            return !(*this == other);
        }
};

class Object {
    public:
        inline static std::vector<Object> objects = {};

        std::vector<Point_3d> v, vn;
        std::vector<Point_2d> vt;
        std::vector<std::array<std::array<int, 3>, 3>> f;
        Point_3d ihat, jhat, khat, origin;

        Object(std::vector<Point_3d> v = {}, std::vector<Point_3d> vn = {}, std::vector<Point_2d> vt = {}, std::vector<std::array<std::array<int, 3>, 3>> f = {}) {
            this->v = v;
            this->vn = vn;
            this->vt = vt;
            this->f = f;
            this->origin = Point_3d();
            this->ihat = Point_3d(1, 0, 0);
            this->jhat = Point_3d(0, 1, 0);
            this->khat = Point_3d(0, 0, 1);
            objects.push_back(*this);
        }

        Object move(Point_3d p) {
            Object obj = *this;
            obj.origin = p;
            return obj;
        }

        Object rotate(float angle, char axis) const {
            Object obj = *this;
            const float sina = sin(angle);
            const float cosa = cos(angle);
            if (axis == 'x') {
                const Point_3d ihat_ = Point_3d(1, 0, 0);
                const Point_3d jhat_ = Point_3d(0, cosa, sina);
                const Point_3d khat_ = Point_3d(0, -sina, cosa);

                obj.ihat = ihat.transform(ihat_, jhat_, khat_);
                obj.jhat = jhat.transform(ihat_, jhat_, khat_);
                obj.khat = khat.transform(ihat_, jhat_, khat_);
            } else if (axis == 'y') {
                const Point_3d ihat_ = Point_3d(cosa, 0, -sina);
                const Point_3d jhat_ = Point_3d(0, 1, 0);
                const Point_3d khat_ = Point_3d(sina, 0, cosa);

                obj.ihat = ihat.transform(ihat_, jhat_, khat_);
                obj.jhat = jhat.transform(ihat_, jhat_, khat_);
                obj.khat = khat.transform(ihat_, jhat_, khat_);
            } else if (axis == 'z') {
                const Point_3d ihat_ = Point_3d(cosa, sina, 0);
                const Point_3d jhat_ = Point_3d(-sina, cosa, 0);
                const Point_3d khat_ = Point_3d(0, 0, 1);

                obj.ihat = ihat.transform(ihat_, jhat_, khat_);
                obj.jhat = jhat.transform(ihat_, jhat_, khat_);
                obj.khat = khat.transform(ihat_, jhat_, khat_);
            } else {
                throw std::invalid_argument("Invalid axis for rotation. Use 'x', 'y', or 'z'.");
            }
            return obj;
        }

        Object scale(float scalar, char axis) const {
            Object obj = *this;
            if (axis == 'x') {
                obj.ihat = ihat * scalar;
            } else if (axis == 'y') {
                obj.jhat = jhat * scalar;
            } else if (axis == 'z') {
                obj.khat = khat * scalar;
            } else {
                throw std::invalid_argument("Invalid axis for scaling. Use 'x', 'y', or 'z'.");
            }
            return obj;
        }

        Object shear(float shear_factor, char axis, char by_axis) const {
            Object obj = *this;
            if (axis == by_axis) {
                throw std::invalid_argument("Cannot shear an axis by itself.");
            }
            if (axis == 'x') {
                if (by_axis == 'y') {
                    obj.ihat = ihat + jhat * shear_factor;
                } else if (by_axis == 'z') {
                    obj.ihat = ihat + khat * shear_factor;
                } else {
                    throw std::invalid_argument("Invalid by_axis for shearing x. Use 'y' or 'z'.");
                }
            } else if (axis == 'y') {
                if (by_axis == 'x') {
                    obj.jhat = jhat + ihat * shear_factor;
                } else if (by_axis == 'z') {
                    obj.jhat = jhat + khat * shear_factor;
                } else {
                    throw std::invalid_argument("Invalid by_axis for shearing y. Use 'x' or 'z'.");
                }
            } else if (axis == 'z') {
                if (by_axis == 'x') {
                    obj.khat = khat + ihat * shear_factor;
                } else if (by_axis == 'y') {
                    obj.khat = khat + jhat * shear_factor;
                } else {
                    throw std::invalid_argument("Invalid by_axis for shearing z. Use 'x' or 'y'.");
                }
            } else {
                throw std::invalid_argument("Invalid axis for shearing. Use 'x', 'y', or 'z'.");
            }
            return obj;
        }

        std::vector<Point_3d> get_v() const {
            std::vector<Point_3d> vertices;
            for (const Point_3d& vertex : v) {
                vertices.push_back(vertex.transform(ihat, jhat, khat) + origin);
            }
            return vertices;
        }

        std::vector<Point_3d> get_vn() const {
            std::vector<Point_3d> normals;
            for (const Point_3d& normal : vn) {
                normals.push_back(normal.transform(ihat, jhat, khat));
            }
            return normals;
        }
};

class Utils {
    public: 
        static std::vector<std::string> split_str(const std::string& str, char delimiter);
        static std::string read_file(std::string path);
        static std::tuple<std::vector<Point_3d>, std::vector<Point_3d>, std::vector<Point_2d>, std::vector<std::array<std::array<int, 3>, 3>>> read_obj(std::string path);
        static Point_3d ray_equation(Point_3d p, Point_3d d, float t);
        static float ray_intersects_triangle(Point_3d a, Point_3d b, Point_3d c, Point_3d p, Point_3d d);
};

std::string Utils::read_file(std::string path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Could not open file: \"" + path + "\"");
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

// Foolproof this read_obj function (idk how tho)
// MrurBo, if you're reading this, HELP.
std::tuple<std::vector<Point_3d>, std::vector<Point_3d>, std::vector<Point_2d>, std::vector<std::array<std::array<int, 3>, 3>>> Utils::read_obj(std::string path) {
    std::string file_str = Utils::read_file(path);
    if (file_str.empty()) return std::tuple<std::vector<Point_3d>, std::vector<Point_3d>, std::vector<Point_2d>, std::vector<std::array<std::array<int, 3>, 3>>>{};

    std::vector<std::string> file_arr = Utils::split_str(file_str, '\n');
    std::vector<Point_3d> v, vn;
    std::vector<Point_2d> vt;
    std::vector<std::array<std::array<int, 3>, 3>> f;

    for (const std::string& line : file_arr) {
        std::vector<std::string> comps = Utils::split_str(line, ' ');

        if (comps[0] == "v" && comps.size() >= 4) {
            v.push_back(Point_3d(std::stof(comps[1]), std::stof(comps[2]), std::stof(comps[3])));
        } else if (comps[0] == "vn" && comps.size() >= 4) {
            vn.push_back(Point_3d(std::stof(comps[1]), std::stof(comps[2]), std::stof(comps[3])));
        } else if (comps[0] == "vt" && comps.size() >= 3) {
            vt.push_back(Point_2d(std::stof(comps[1]), std::stof(comps[2])));
        } else if (comps[0] == "f" && comps.size() >= 4) {
            std::array<std::array<int, 3>, 3> face;
            for (int i = 0; i < 3; i++) {
                std::vector<std::string> vertexData = Utils::split_str(comps[i + 1], '/');
                face[i] = {std::stoi(vertexData[0]) - 1, std::stoi(vertexData[1]) - 1, std::stoi(vertexData[2]) - 1};
            }
            f.push_back(face);
        }
    }

    return std::make_tuple(v, vn, vt, f);
}

Point_3d Utils::ray_equation(Point_3d p, Point_3d d, float t) {
    return p + (d.normalize() * t);
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

    return m0 == m1 && m1 == m2 ? dist : -1;
}

class Cam {
    public:
        inline static float x = 0;
        inline static float y = 0;
        inline static float z = 0;
        inline static float roth = 0;
        inline static float rotv = 0;

        inline static const float movespeed = 6 / 60.0f;
        inline static const float turnspeed = PI / 45;
        inline static const float fov = 90 * (PI / 180);
        inline static const int screen_width = 192;
        inline static const int screen_height = 108;
};

std::vector<Point_3d> calc_rays() {
    const float aspect = static_cast<float>(Cam::screen_width) / Cam::screen_height;

    const Point_3d forward = Point_3d(
        sin(Cam::roth) * cos(Cam::rotv),
        sin(Cam::rotv),
        cos(Cam::roth) * cos(Cam::rotv)
    );
    const Point_3d right = Point_3d(
        sin(Cam::roth + PI / 2),
        0,
        cos(Cam::roth + PI / 2)
    );
    const Point_3d up = right.cross(forward);

    float half_fov = tan(Cam::fov / 2);

    std::vector<Point_3d> rays;
    for (int y = 0; y < Cam::screen_height; y++)
    for (int x = 0; x < Cam::screen_width; x++) {
        float px = (x / (Cam::screen_width / 2.0f) - 1.0f) * half_fov * aspect;
        float py = (1.0f - y / (Cam::screen_height / 2.0f)) * half_fov;

        rays.push_back(forward + right * px + up * py);
    }

    return rays;
}

std::vector<Color> rand_colors() {
    std::vector<Color> colors;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    for (const Object& obj : Object::objects) {
        for (const auto& face : obj.f) {
            unsigned char r = static_cast<unsigned char>(dist(gen));
            unsigned char g = static_cast<unsigned char>(dist(gen));
            unsigned char b = static_cast<unsigned char>(dist(gen));

            colors.push_back({
                static_cast<unsigned char>(r),
                static_cast<unsigned char>(g),
                static_cast<unsigned char>(b),
                255
            });
        }
    }
    return colors;
}

void draw_world(std::vector<Color> colors) {
    std::vector<Point_3d> dirs = calc_rays();
    Point_3d p = Point_3d(Cam::x, Cam::y, Cam::z);

    for (int i = 0; i < dirs.size(); i++) {
        Point_3d d = dirs[i];
        Point_2d screen_pos = Point_2d(i % Cam::screen_width, i / Cam::screen_width);

        int closest_face_idx = -1;
        float min_dist = INFINITY;

        for (const Object& obj : Object::objects) {
            std::vector<Point_3d> v = obj.get_v();

            for (int j = 0; j < obj.f.size(); j++) {
                const auto& face = obj.f[j];
                Point_3d a = v[face[0][0]];
                Point_3d b = v[face[1][0]];
                Point_3d c = v[face[2][0]];
                float dist = Utils::ray_intersects_triangle(a, b, c, p, d);
                if (dist > 0 && dist < min_dist) {
                    min_dist = dist;
                    closest_face_idx = j;
                }

            }
        }

        if (closest_face_idx != -1) DrawPixel(screen_pos.x, screen_pos.y, colors[closest_face_idx]);
    }
}

int main() {
    InitWindow(Cam::screen_width, Cam::screen_height, "Rasteriser In   ++(raylib)");
    SetTargetFPS(60);

    auto file_data = Utils::read_obj("assets/Cube.obj");
    auto [v, vn, vt, f] = file_data;
    Object(v, vn, vt, f);
    Object::objects[0] = Object::objects[0].move(Point_3d(0, 0, 4));

    std::vector<Color> colors = rand_colors();

    while (!WindowShouldClose()) {
        // Update
        float deltaTime = GetFrameTime() * 60;

        Cam::roth += (IsKeyDown(KEY_RIGHT) - IsKeyDown(KEY_LEFT)) * Cam::turnspeed * deltaTime;
        Cam::rotv -= (IsKeyDown(KEY_UP) - IsKeyDown(KEY_DOWN)) * Cam::turnspeed * deltaTime;
        Cam::rotv = std::clamp(Cam::rotv, -PI / 2, PI / 2);

        float sinh = sin(Cam::roth);
        float cosh = cos(Cam::roth);
        float sinhp90 = sin(Cam::roth + PI / 2);
        float coshp90 = cos(Cam::roth + PI / 2);

        Cam::x += (IsKeyDown(KEY_W) - IsKeyDown(KEY_S)) * Cam::movespeed * sinh * deltaTime;
        Cam::x += (IsKeyDown(KEY_D) - IsKeyDown(KEY_A)) * Cam::movespeed * sinhp90 * deltaTime;
        Cam::z += (IsKeyDown(KEY_W) - IsKeyDown(KEY_S)) * Cam::movespeed * cosh * deltaTime;
        Cam::z += (IsKeyDown(KEY_D) - IsKeyDown(KEY_A)) * Cam::movespeed * coshp90 * deltaTime;

        Cam::y += (IsKeyDown(KEY_Q) - IsKeyDown(KEY_E)) * Cam::movespeed * deltaTime;


        // Draw
        BeginDrawing();
        ClearBackground(BLACK);

        draw_world(colors);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}