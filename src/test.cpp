#include "./utils.cpp"

int main() {
    const auto [v, vn, vt, f] = Utils::read_obj("../assets/Cube.obj");

    std::cout << "Vertices:\n";
    for (const auto& vertex : v) {
        std::cout << "(" << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << ")\n";
    }

    std::cout << "\nVertex Normals:\n";
    for (const auto& normal : vn) {
        std::cout << "(" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")\n";
    }

    std::cout << "\nTexture Coordinates:\n";
    for (const auto& tex : vt) {
        std::cout << "(" << tex[0] << ", " << tex[1] << ")\n";
    }

    std::cout << "\nFaces:\n";
    for (const auto& face : f) {
        std::cout << "(" << face[0] << ", " << face[1] << ", " << face[2] << ")\n";
    }

    return 0;
}