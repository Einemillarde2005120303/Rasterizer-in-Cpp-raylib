#include "./utils.cpp"

int main() {
    const auto [v, vn, vt, f] = Utils::read_obj("../assets/Cube.obj");

    for (const auto& vertex : v) {
        std::cout << "v " << vertex[0] << " " << vertex[1] << " " << vertex[2] << "\n";
    }
    for (const auto& normal : vn) {
        std::cout << "vn " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
    }
    for (const auto& tex : vt) {
        std::cout << "vt " << tex[0] << " " << tex[1] << "\n";
    }
    for (const auto& face : f) {
        std::cout << "f " << face[0][0] << "/" << face[0][1] << "/" << face[0][2] << " " << face[1][0] << "/" << face[1][1] << "/" << face[1][2] << " " << face[2][0] << "/" << face[2][1] << "/" << face[2][2] << "\n";
    }

    return 0;
}