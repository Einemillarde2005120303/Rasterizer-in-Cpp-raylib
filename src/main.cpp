#include <raylib.h>
#include <cmath>
#include <array>
#include <algorithm>

class Point {
    public:
        float x;
        float y;

        Point(float x, float y) {
            this->x = x;
            this->y = y;
        }
};

float tri_area(Point a, Point b, Point c) {
    return std::abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) * 0.5;
}

void draw_triangle(Point a, Point b, Point c) {
    const std::array<float, 3> xs = {a.x, b.x, c.x};
    const std::array<float, 3> ys = {a.y, b.y, c.y};

    const float min_x = *std::min_element(xs.begin(), xs.end());
    const float min_y = *std::min_element(ys.begin(), ys.end());
    const float max_x = *std::max_element(xs.begin(), xs.end());
    const float max_y = *std::max_element(ys.begin(), ys.end());

    const float A = tri_area(a, b, c);
    float A1, A2, A3, t;
    unsigned char red, green, blue;

    for (size_t y = min_y; y <= max_y; ++y) {
        for (size_t x = min_x; x <= max_x; ++x) {
            A1 = tri_area(Point(x, y), b, c);
            A2 = tri_area(a, Point(x, y), c);
            A3 = tri_area(a, b, Point(x, y));

            if (std::abs((A1 + A2 + A3) - A) < 0.1f) {
                t = 255 / A;

                red = (unsigned char)std::clamp(t * A1, 0.0f, 255.0f);
                green = (unsigned char)std::clamp(t * A2, 0.0f, 255.0f);
                blue = (unsigned char)std::clamp(t * A3, 0.0f, 255.0f);
                Color color = { red, green, blue, 255 };

                DrawPixel(x, y, color);
            }
        }
    }
}

int main() {
    InitWindow(960, 540, "RGB Triangle");
    SetTargetFPS(60);


    while (!WindowShouldClose()) {
        // Update


        // Draw
        BeginDrawing();
        ClearBackground(BLACK);
        draw_triangle(Point(480, 10), Point(180, 520), Point(780, 520));
        EndDrawing();
    }

    CloseWindow();
    return 0;
}