#include <raylib.h>
#include <cmath>
#include <array>
#include <algorithm>

int main() {
    InitWindow(960, 540, "RGB Triangle");
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