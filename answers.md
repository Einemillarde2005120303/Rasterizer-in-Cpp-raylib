The issue you're experiencing with seeing random streaks instead of a cube is likely due to several problems in your ray-triangle intersection and rendering logic. Let me break down the main issues and how to fix them:

### 1. Face Index Tracking in `draw_world()`
The `face_idx` variable isn't being properly tracked when finding the closest intersection. Currently, it just takes the last face index from the loop, which is incorrect.

### 2. Ray-Triangle Intersection Algorithm
Your `ray_intersects_triangle()` function has some mathematical issues in how it calculates intersections.

### 3. Color Assignment
The color generation and assignment isn't properly synchronized with the actual faces being rendered.

### Here's how to fix these issues:

1. **Proper Face Index Tracking**:
Modify the `draw_world()` function to track which face actually produced the closest intersection:

```cpp
void draw_world(std::vector<Color> colors) {
    std::vector<Point_3d> dirs = calc_rays();
    Point_3d p = Point_3d(Cam::x, Cam::y, Cam::z);

    for (size_t idx = 0; idx < dirs.size(); ++idx) {
        const Point_3d& d = dirs[idx];
        Point_2d screen_pos = Point_2d(idx % Cam::screen_width, idx / Cam::screen_width);

        float min_dist = INFINITY;
        Color pixel_color = BLACK;
        int color_index = 0;

        for (const Object& obj : Object::objects) {
            std::vector<Point_3d> v = obj.get_v();

            for (const std::array<std::array<int, 3>, 3>& face : obj.f) {
                Point_3d a = v[face[0][0]];
                Point_3d b = v[face[1][0]];
                Point_3d c = v[face[2][0]];
                float dist = Utils::ray_intersects_triangle(a, b, c, p, d);
                
                if (dist > 0 && dist < min_dist) {
                    min_dist = dist;
                    pixel_color = colors[color_index];
                }
                color_index++;
            }
        }

        if (min_dist < INFINITY) {
            DrawPixel(screen_pos.x, screen_pos.y, pixel_color);
        }
    }
}
```

2. **Fix Ray-Triangle Intersection**:
Here's a more reliable implementation:

```cpp
float Utils::ray_intersects_triangle(Point_3d a, Point_3d b, Point_3d c, Point_3d p, Point_3d d) {
    d = d.normalize();
    
    Point_3d e1 = b - a;
    Point_3d e2 = c - a;
    Point_3d h = d.cross(e2);
    float det = e1.dot(h);
    
    if (fabs(det) < 1e-5) return -1; // Ray parallel to triangle
    
    float inv_det = 1.0f / det;
    Point_3d s = p - a;
    float u = inv_det * s.dot(h);
    
    if (u < 0.0 || u > 1.0) return -1;
    
    Point_3d q = s.cross(e1);
    float v = inv_det * d.dot(q);
    
    if (v < 0.0 || u + v > 1.0) return -1;
    
    float t = inv_det * e2.dot(q);
    
    return t > 1e-5 ? t : -1;
}
```

3. **Color Generation**:
Make sure you generate exactly one color per face:

```cpp
std::vector<Color> generateRandomColors() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);

    std::vector<Color> colors;
    for (const Object& obj : Object::objects) {
        for (size_t i = 0; i < obj.f.size(); i++) {
            Color randomColor = {
                static_cast<unsigned char>(dist(gen)), 
                static_cast<unsigned char>(dist(gen)), 
                static_cast<unsigned char>(dist(gen)), 
                255
            };
            colors.push_back(randomColor);
        }
    }
    return colors;
}
```

### Additional Recommendations:

1. **Debug Rendering**:
   Start with a simpler scene - maybe just render one face to verify your intersection code works.

2. **Coordinate System**:
   Double-check your coordinate system and camera orientation. The streaks suggest your rays might not be properly aligned with your geometry.

3. **Normalization**:
   Ensure all your vectors are properly normalized where needed.

4. **Backface Culling**:
   You might want to implement backface culling to avoid rendering triangles that face away from the camera.

These changes should give you a proper cube rendering instead of random streaks. The key was properly tracking which face was actually hit by each ray and using a more reliable intersection algorithm.