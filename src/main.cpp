#include <raylib.h>

int main() {
    InitWindow(960, 540, "Rasteriser In C++ (raylib)");
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