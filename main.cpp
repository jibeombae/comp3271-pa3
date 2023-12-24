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


// Function to render a specific range of pixels in the image
void runThread(int start_x, int end_x, int width, int height, const Camera& camera,
                 const Scene& scene, int spp, float* image_buffer) {
    for (int x = start_x; x < end_x; x++) {
        for (int y = 0; y < height; y++) {
            Color color(0.f, 0.f, 0.f);
            for (int i = 0; i < spp; i++) {
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
        }
    }
}

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

    int spp = 256;
    int NUM_THREAD = 128;

    int pixels_per_thread = width / NUM_THREAD;

    // Vector to store thread objects
    std::vector<std::thread> threads;

    // Launch threads
    for (int i = 0; i < NUM_THREAD; ++i) {
        int start_x = i * pixels_per_thread;
        int end_x = (i + 1) * pixels_per_thread;
        if (i == NUM_THREAD - 1) {
            // Last thread should handle any remaining pixels
            end_x = width;
        }
        // Split the image pixels and make multiple threads work on it
        threads.emplace_back(runThread, start_x, end_x, width, height, std::ref(camera), std::ref(scene), spp, image_buffer);
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
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
    Vec C, B;

    // Check if N is far from the y-axis
    if (std::abs(N.x) > std::abs(N.y)) {
        C = glm::normalize(Vec(N.z, 0, -N.x));
    } else {
        // Handle the case when N is close to the y-axis
        C = glm::normalize(Vec(0, -N.z, N.y));
    }

    // Calculate B using cross product
    B = glm::normalize(glm::cross(C, N));

    // Create the transformation matrix T
    float T[3][3] = {
        {C.x, N.x, B.x},
        {C.y, N.y, B.y},
        {C.z, N.z, B.z}
    };

    // Multiply the localRay vector by the transformation matrix T
    float x = T[0][0] * localRay.x + T[0][1] * localRay.y + T[0][2] * localRay.z;
    float y = T[1][0] * localRay.x + T[1][1] * localRay.y + T[1][2] * localRay.z;
    float z = T[2][0] * localRay.x + T[2][1] * localRay.y + T[2][2] * localRay.z;

    return Vec(x, y, z);
}

Vec SampleHemisphere(const HitRecord& hit_record)
{
    // TODO: add your code here.
    // This function randomly samples a direction on the hemisphere.
    // It will calls toWorld() to transform the local coordinate to world coordinate.
    
    float theta = static_cast<float>(2.0 * M_PI * get_random_float());
    float phi = static_cast<float>(0.5 * M_PI * get_random_float());

    // Direction vector in local space
    Vec localDirection(
        sin(phi) * cos(theta),
        cos(phi),
        sin(phi) * sin(theta)
    );

    // Transform the direction to world space
    return toWorld(localDirection, hit_record.normal);
}



Color Shade(const std::vector<LightSource>& light_sources,
            const Hittable& hittable_collection,
            const HitRecord& hit_record,
            int trace_depth) {
    // TODO: Add your code here.
    Color color(0.f, 0.f, 0.f);

    color = hit_record.material.k_a * hit_record.material.ambient;

    if (trace_depth < kMaxTraceDepth){
        // Pick a random direction from here and keep going.
        Ray reflected_ray;
        reflected_ray.o = hit_record.position;

        // Sample random dir on hemisphere
        reflected_ray.d = SampleHemisphere(hit_record);
        const float pdf = 1 / M_PI;

        // Compute the BRDF for this ray
        float nwi = glm::dot(reflected_ray.d, hit_record.normal);
        if (nwi < 0.00001) nwi = 0.0;
        Color BRDF = hit_record.material.k_d * hit_record.material.diffuse / M_PI;

        // Recursively trace reflected light sources.
        Color incoming_rad = TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth + 1);

        // Apply the Rendering Equation here.
        color += hit_record.material.emission + (incoming_rad * BRDF * nwi / pdf);

    }

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