#include <iostream>
#include <vector>
#include <thread>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

Color TraceRay(const Ray& ray, const std::vector<LightSource>& light_sources, const Hittable& scene, int trace_depth);
Color Shade(const std::vector<LightSource>& light_sources,
    const Hittable& hittable_collection,
    const HitRecord& hit_record,
    int trace_depth);


int main() {
    // TODO: Set your workdir (absolute path) here. Don't forget the last slash
    const std::string work_dir("/Users/jibeombae/3B (HKU)/comp3271-graphics/A3/PA3_Release/");

    // Construct scene
    Scene scene(work_dir, "scene/teapot_area_light.toml");
    const Camera& camera = scene.camera_;
    int width = camera.width_;
    int height = camera.height_;

    std::vector<unsigned char> image(width * height * 4, 0);
    float* image_buffer = new float[width * height * 3];
    for (int i = 0; i < width * height * 3; ++i) {
		image_buffer[i] = 0.f;
	}

    int spp = 128;
    int NUM_THREAD = 1;

    
    float progress = 0.f;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Color color(0.f, 0.f, 0.f);
            for (int i = 0; i < spp; i++)
            {
                float bias_x = get_random_float();
                float bias_y = get_random_float();
                Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
                color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
            }
            color /= float(spp);
            int idx = 3 * ((height - y - 1) * width + x);
            image_buffer[idx] += color.r;
            image_buffer[idx + 1] += color.g;
            image_buffer[idx + 2] += color.b;


            float curr_progress = float(x * height + y) / float(height * width);
            if (curr_progress > progress + 0.05f) {
                progress += 0.05f;
                std::cout << "Progress (thread " << 1 << "/" << NUM_THREAD << "): " << progress << std::endl;
            }
        }
    }

    // copy from image_buffer to image
    for (int i = 0; i < width * height; ++i) {
        for (int j = 0; j < 3; ++j) {
			image[4 * i + j] = (uint8_t)(glm::min(image_buffer[3 * i + j], 1.f - 1e-5f) * 256.f);
		}
		image[4 * i + 3] = 255;
	}
    

    std::vector<unsigned char> png;
    unsigned error = lodepng::encode(png, image, width, height);
    lodepng::save_file(png, work_dir + "outputs/teapot_area_light.png");
}


Vec toWorld(const Vec& localRay, const Vec& N) {
    // TODO: add your code here.
    // This function transforms a vector from local coordinate to world coordinate.
}

Vec SampleHemisphere(const HitRecord& hit_record)
{
    // TODO: add your code here.
    // This function randomly samples a direction on the hemisphere.
    // It will calls toWorld() to transform the local coordinate to world coordinate.
}



Color Shade(const std::vector<LightSource>& light_sources,
            const Hittable& hittable_collection,
            const HitRecord& hit_record,
            int trace_depth) {
    // TODO: Add your code here.
    Color color(0.f, 0.f, 0.f);

    return color;
}

Color TraceRay(const Ray& ray,
               const std::vector<LightSource>& light_sources,
               const Hittable& hittable_collection,
               int trace_depth) {
    // TODO: Add your code here.
    HitRecord record;
    Color color(0.0f, 0.0f, 0.0f);

    if(hittable_collection.Hit(ray, &record)) {
        // if there's a hit, calculate the color
        color = Shade(light_sources, hittable_collection, record, trace_depth);
    }

    return color;
}