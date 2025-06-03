#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <tuple>

class Utils {
    public: 
        static std::vector<std::string> split_str(const std::string& str, char delimiter);
        static std::string read_file(std::string path);
        static std::tuple<std::vector<std::array<float, 3>>, std::vector<std::array<float, 3>>, std::vector<std::array<float, 2>>, std::vector<std::array<std::array<int, 3>, 3>>> read_obj(std::string path);
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

std::tuple<std::vector<std::array<float, 3>>, std::vector<std::array<float, 3>>, std::vector<std::array<float, 2>>, std::vector<std::array<std::array<int, 3>, 3>>> Utils::read_obj(std::string path) {
    std::string file_str = Utils::read_file(path);
    if (file_str.empty()) return std::tuple<std::vector<std::array<float, 3>>, std::vector<std::array<float, 3>>, std::vector<std::array<float, 2>>, std::vector<std::array<std::array<int, 3>, 3>>>{};
;

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