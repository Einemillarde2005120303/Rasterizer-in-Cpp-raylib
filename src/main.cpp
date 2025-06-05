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
        inline static std::vector<Object> objects = std::vector<Object>{};

        std::vector<Point_3d> v, vn;
        std::vector<Point_2d> vt;
        std::vector<std::array<std::array<int, 3>, 3>> f;
        Point_3d ihat, jhat, khat, origin;

        Object(std::vector<Point_3d> v = std::vector<Point_3d>{}, std::vector<Point_3d> vn = std::vector<Point_3d>{}, std::vector<Point_2d> vt = std::vector<Point_2d>{}, std::vector<std::array<std::array<int, 3>, 3>> f = std::vector<std::array<std::array<int, 3>, 3>>{}) {
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

        Object rotate(float angle, char axis1, char axis2) const {
            Object obj = *this;
            if (axis1 == 'x' && axis2 == 'y' || axis1 == 'y' && axis2 == 'x') {
                obj.ihat = Point_3d(
                    ihat.x * cos(angle) - ihat.y * sin(angle),
                    ihat.x * sin(angle) + ihat.y * cos(angle),
                    ihat.z
                );
                obj.jhat = Point_3d(
                    jhat.x * cos(angle) - jhat.y * sin(angle),
                    jhat.x * sin(angle) + jhat.y * cos(angle),
                    jhat.z
                );
            } else if (axis1 == 'y' && axis2 == 'z' || axis1 == 'z' && axis2 == 'y') {
                obj.jhat = Point_3d(
                    jhat.x,
                    jhat.y * cos(angle) - jhat.z * sin(angle),
                    jhat.y * sin(angle) + jhat.z * cos(angle)
                );
                obj.khat = Point_3d(
                    khat.x,
                    khat.y * cos(angle) - khat.z * sin(angle),
                    khat.y * sin(angle) + khat.z * cos(angle)
                );
            } else if (axis1 == 'x' && axis2 == 'z' || axis1 == 'z' && axis2 == 'x') {
                obj.ihat = Point_3d(
                    ihat.x,
                    ihat.y * cos(angle) - ihat.z * sin(angle),
                    ihat.y * sin(angle) + ihat.z * cos(angle)
                );
                obj.khat = Point_3d(
                    khat.x,
                    khat.y * cos(angle) - khat.z * sin(angle),
                    khat.y * sin(angle) + khat.z * cos(angle)
                );
            } else {
                throw std::invalid_argument("Invalid axes for rotation. Use 'x', 'y', or 'z'.");
            }
        }

        Object scale(float scalar, char axis) const {
            Object obj = *this;
            if (axis == 'x') {
                obj.ihat = Point_3d(ihat.x * scalar, ihat.y * scalar, ihat.z * scalar);
            } else if (axis == 'y') {
                obj.jhat = Point_3d(jhat.x * scalar, jhat.y * scalar, jhat.z * scalar);
            } else if (axis == 'z') {
                obj.khat = Point_3d(khat.x * scalar, khat.y * scalar, khat.z * scalar);
            } else {
                throw std::invalid_argument("Invalid axis for scaling. Use 'x', 'y', or 'z'.");
            }
            return obj;
        }

        Object shear(float shear_factor, char axis1, char axis2) const {
            Object obj = *this;
            if (axis1 == 'x' && axis2 == 'y' || axis1 == 'y' && axis2 == 'x') {
                obj.ihat = Point_3d(
                    ihat.x + shear_factor * jhat.x,
                    ihat.y + shear_factor * jhat.y,
                    ihat.z + shear_factor * jhat.z
                );
                obj.jhat = Point_3d(
                    jhat.x + shear_factor * ihat.x,
                    jhat.y + shear_factor * ihat.y,
                    jhat.z + shear_factor * ihat.z
                );
            } else if (axis1 == 'y' && axis2 == 'z' || axis1 == 'z' && axis2 == 'y') {
                obj.jhat = Point_3d(
                    jhat.x + shear_factor * khat.x,
                    jhat.y + shear_factor * khat.y,
                    jhat.z + shear_factor * khat.z
                );
                obj.khat = Point_3d(
                    khat.x + shear_factor * jhat.x,
                    khat.y + shear_factor * jhat.y,
                    khat.z + shear_factor * jhat.z
                );
            } else if (axis1 == 'x' && axis2 == 'z' || axis1 == 'z' && axis2 == 'x') {
                obj.ihat = Point_3d(
                    ihat.x + shear_factor * khat.x,
                    ihat.y + shear_factor * khat.y,
                    ihat.z + shear_factor * khat.z
                );
                obj.khat = Point_3d(
                    khat.x + shear_factor * ihat.x,
                    khat.y + shear_factor * ihat.y,
                    khat.z + shear_factor * ihat.z
                );
            } else {
                throw std::invalid_argument("Invalid axes for shearing. Use 'x', 'y', or 'z'.");
            }
            return obj;
        }

        std::vector<Point_3d> get_v() const {
            std::vector<Point_3d> vertices;
            for (const Point_3d& vertex : v) {
                vertices.push_back(Point_3d(
                    origin.x + vertex.dot(Point_3d(ihat.x, jhat.x, khat.x)),
                    origin.y + vertex.dot(Point_3d(ihat.y, jhat.y, khat.y)),
                    origin.z + vertex.dot(Point_3d(ihat.z, jhat.z, khat.z))
                ));
            }
            return vertices;
        }

        std::vector<Point_3d> get_vn() const {
            std::vector<Point_3d> normals;
            for (const Point_3d& normal : vn) {
                normals.push_back(Point_3d(
                    normal.dot(Point_3d(ihat.x, jhat.x, khat.x)),
                    normal.dot(Point_3d(ihat.y, jhat.y, khat.y)),
                    normal.dot(Point_3d(ihat.z, jhat.z, khat.z))
                ));
            }
            return normals;
        }

        std::vector<Point_2d> get_vt() const {
            std::vector<Point_2d> texture_coords;
            for (const Point_2d& coord : vt) {
                texture_coords.push_back(Point_2d(
                    coord.x * ihat.x + coord.y * jhat.x,
                    coord.x * ihat.y + coord.y * jhat.y
                ));
            }
            return texture_coords;
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
            for (int i = 0; i < 3; ++i) {
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

        inline static const float movespeed = 8 / 60.0f;
        inline static const float turnspeed = PI / 60;
        inline static const float fov = 60 * (PI / 180);
        inline static const int screen_width = 160;
        inline static const int screen_height = 90;
        static std::vector<Point_3d> default_rays;
};

std::vector<Point_3d> Cam::default_rays = []() {
    std::vector<Point_3d> rays;
    float aspect_ratio = static_cast<float>(Cam::screen_width) / Cam::screen_height;
    float half_fov = tan(Cam::fov / 2);

    for (int j = 0; j < Cam::screen_height; ++j) {
        for (int i = 0; i < Cam::screen_width; ++i) {
            float x = (2.0f * (i + 0.5f) / Cam::screen_width - 1) * aspect_ratio * half_fov;
            float y = (1 - 2.0f * (j + 0.5f) / Cam::screen_height) * half_fov;
            rays.push_back(Point_3d(x, y, 1).normalize());
        }
    }
    return rays;
}();

std::vector<Point_3d> calc_rays() {
    const Point_3d forward = Point_3d(
        sin(Cam::roth) * cos(Cam::rotv),
        sin(Cam::rotv),
        cos(Cam::roth) * cos(Cam::rotv)
    );
    const Point_3d right = Point_3d(
        sin(Cam::roth + (PI / 2)),
        0,
        cos(Cam::roth + (PI / 2))
    );
    const Point_3d up = Point_3d(
        sin(Cam::roth) * cos(Cam::rotv + (PI / 2)),
        sin(Cam::rotv + (PI / 2)),
        cos(Cam::roth) * cos(Cam::rotv + (PI / 2))
    );

    std::vector<Point_3d> rays;
    for (const Point_3d& ray : Cam::default_rays) {
        Point_3d transformed_ray = right * ray.x + up * ray.y + forward * ray.z;
        rays.push_back(transformed_ray);
    }

    return rays;
}

// Temporary (thanks ChatGPT)
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

            unsigned char max_component = std::max({r, g, b});

            if (max_component == 0) {
                colors.push_back({255, 255, 255, 255});
                continue;
            }

            float scale = 255.0f / max_component;
            colors.push_back({
                static_cast<unsigned char>(r * scale),
                static_cast<unsigned char>(g * scale),
                static_cast<unsigned char>(b * scale),
                255
            });
        }
    }
    return colors;
}

void draw_world(std::vector<Color> colors) {
    std::vector<Point_3d> dirs = calc_rays();
    Point_3d p = Point_3d(Cam::x, Cam::y, Cam::z);

    for (int i = 0; i < dirs.size(); ++i) {
        Point_3d d = dirs[i];
        Point_2d screen_pos = Point_2d(i % Cam::screen_width, i / Cam::screen_width);

        int closest_face_idx = -1;
        float min_dist = INFINITY;

        for (const Object& obj : Object::objects) {
            std::vector<Point_3d> v = obj.get_v();

            for (int j = 0; j < obj.f.size(); ++j) {
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
    InitWindow(Cam::screen_width, Cam::screen_height, "Rasteriser In C++ (raylib)");
    SetTargetFPS(60);

    auto file_data = Utils::read_obj("assets/Cube.obj");
    auto [v, vn, vt, f] = file_data;
    Object(v, vn, vt, f);
    Object::objects[0] = Object::objects[0].move(Point_3d(0, 0, 5));

    std::vector<Color> colors = rand_colors();

    while (!WindowShouldClose()) {
        // Update
        float deltaTime = GetFrameTime() * 60;

        Cam::roth += (IsKeyDown(KEY_RIGHT) - IsKeyDown(KEY_LEFT)) * Cam::turnspeed * deltaTime;
        Cam::rotv += (IsKeyDown(KEY_UP) - IsKeyDown(KEY_DOWN)) * Cam::turnspeed * deltaTime;
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

        std::cout << GetFPS() << " FPS" << "\n";

        // Draw
        BeginDrawing();
        ClearBackground(BLACK);

        draw_world(colors);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}