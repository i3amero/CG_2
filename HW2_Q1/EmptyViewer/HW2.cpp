#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>
#include <limits>
#include <cmath>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;

// Set screen resolution
const int WIDTH = 512, HEIGHT = 512;//Resolution that assigned in HW1(512 * 512)
const float INF = std::numeric_limits<float>::infinity();

// Light position (Q1 requirement: point light at -4,4,-3)
const vec3 lightPos = vec3(-4.0f, 4.0f, -3.0f);
// Light color (white) and intensity (assume 1.0)
const vec3 lightColor = vec3(1.0f, 1.0f, 1.0f);

// Material structure
struct Material {
    vec3 ka;        // Ambient reflectance
    vec3 kd;        // Diffuse reflectance
    vec3 ks;        // Specular reflectance
    float shininess;// Specular exponent (power)
};

// Ray Class
struct Ray {
    vec3 origin, direction;
    Ray(const vec3& o, const vec3& d) : origin(o), direction(normalize(d)) {}
};

// Structure to store intersection information
struct HitInfo {
    bool hit;
    float t;             // Intersection distance
    vec3 normal;         // Normal at intersection
    Material material;   // Material of the intersected object
    vec3 hitPoint;       // Intersection point
};

//Surface Class (Abstract Class)
struct Surface {
    Material material;
    virtual bool intersect(const Ray& ray, float& t) const = 0;
    virtual vec3 getNormal(const vec3& hitPoint) const = 0;
};

// Plane class
struct Plane : public Surface {
    vec3 normal;
    float d;

    Plane(const vec3& n, float d, const Material& m) {
        normal = normalize(n);
        this->d = d;
        material = m;
    }

    bool intersect(const Ray& ray, float& t) const override {
        float denom = dot(normal, ray.direction);
        if (fabs(denom) < 1e-6f) return false; // Almost parallel
        t = -(dot(normal, ray.origin) + d) / denom;
        return (t >= 0.0f);
    }

    vec3 getNormal(const vec3& /*hitPoint*/) const override {
        // The plane has the same normal everywhere
        return normal;
    }
};

// Sphere class
struct Sphere : public Surface {
    vec3 center;
    float radius;

    Sphere(const vec3& c, float r, const Material& m) {
        center = c;
        radius = r;
        material = m;
    }

    bool intersect(const Ray& ray, float& t) const override {
        vec3 oc = ray.origin - center;
        float a = dot(ray.direction, ray.direction);
        float b = 2.0f * dot(oc, ray.direction);
        float c = dot(oc, oc) - radius * radius;
        float discriminant = b * b - 4 * a * c;

        if (discriminant < 0.0f) return false;

        float sqrtD = sqrtf(discriminant);
        float t0 = (-b - sqrtD) / (2.0f * a);
        float t1 = (-b + sqrtD) / (2.0f * a);

        // Use the smaller positive t; if not positive, use the larger one
        if (t0 > 0.0f) t = t0;
        else if (t1 > 0.0f) t = t1;
        else return false;

        return true;
    }

    vec3 getNormal(const vec3& hitPoint) const override {
        return normalize(hitPoint - center);
    }
};

// Camera class
struct Camera {
    vec3 eye, u, v, w;// Location of Camera and Direction Vector
    float l, r, b, t, d;// Viewport Setting

    Camera() : eye(0, 0, 0), u(1, 0, 0), v(0, 1, 0), w(0, 0, 1),
        l(-0.1f), r(0.1f), b(-0.1f), t(0.1f), d(0.1f) {}

    Ray generateRay(float i, float j) {
        float px = l + (r - l) * (i + 0.5f) / WIDTH;
        float py = b + (t - b) * (j + 0.5f) / HEIGHT;
        vec3 direction = -d * w + px * u + py * v;
        return Ray(eye, direction);
    }
};

// Scene class
struct Scene {
    std::vector<Surface*> objects;

    void addObject(Surface* obj) {
        objects.push_back(obj);
    }

    // Return the intersection info with the closest object
    bool intersect(const Ray& ray, HitInfo& hitInfo) const {
        float closestT = INF;
        bool isHit = false;
        Surface* hitObject = nullptr;

        for (auto& obj : objects) {
            float t;
            if (obj->intersect(ray, t) && t < closestT) {
                closestT = t;
                hitObject = obj;
                isHit = true;
            }
        }

        if (isHit && hitObject) {
            hitInfo.hit = true;
            hitInfo.t = closestT;
            hitInfo.hitPoint = ray.origin + ray.direction * closestT;
            hitInfo.normal = hitObject->getNormal(hitInfo.hitPoint);
            hitInfo.material = hitObject->material;
        }
        return isHit;
    }
};

// Determine if there's a shadow: if a shadow ray from the intersection point
// to the light hits another object first, it's in shadow
bool inShadow(const Scene& scene, const vec3& point, const vec3& lightDir) {
    // A small offset to avoid floating-point errors
    vec3 shadowOrigin = point + lightDir * 0.001f;
    Ray shadowRay(shadowOrigin, lightDir);

    HitInfo temp;
    if (scene.intersect(shadowRay, temp)) {
        // If t is not too small and is less than the distance to the light, it's in shadow
        float distToHit = length(temp.hitPoint - point);
        float distToLight = length(lightPos - point);
        if (distToHit < distToLight) {
            return true;
        }
    }
    return false;
}

// Function for Phong shading calculation
vec3 phongShading(const Scene& scene, const Ray& ray) {
    HitInfo hit;
    vec3 color(0.0f);

    // If there's an intersection with the scene, compute Phong illumination
    if (scene.intersect(ray, hit)) {
        // 1) Ambient
        vec3 ambient = hit.material.ka * lightColor;

        // 2) Check shadows
        vec3 L = normalize(lightPos - hit.hitPoint);
        bool shadowed = inShadow(scene, hit.hitPoint, L);

        // 3) Diffuse, Specular
        vec3 diffuse(0.0f), specular(0.0f);
        if (!shadowed) {
            float nl = dot(hit.normal, L);
            if (nl > 0.0f) {
                diffuse = hit.material.kd * nl * lightColor;
            }

            // Reflection vector R = reflect(-L, N)
            // V is the camera direction (i.e., -ray.direction)
            vec3 V = normalize(-ray.direction);
            vec3 R = reflect(-L, hit.normal);
            float rv = dot(R, V);
            if (rv > 0.0f && length(hit.material.ks) > 0.0f) {
                specular = hit.material.ks * pow(rv, hit.material.shininess);
            }
        }

        color = ambient + diffuse + specular;
        // Clamp color to [0,1] for display
        color = clamp(color, 0.0f, 1.0f);
    }

    return color;
}

// Function that Renders Scene
void render(const Scene& scene, Camera& camera) {
    // Use triple size for RGB storage
    std::vector<unsigned char> image(WIDTH * HEIGHT * 3, 0);

    for (int j = 0; j < HEIGHT; ++j) {
        for (int i = 0; i < WIDTH; ++i) {
            Ray ray = camera.generateRay((float)i, (float)j);
            vec3 color = phongShading(scene, ray);

            // Convert final color to 0~255
            int idx = (j * WIDTH + i) * 3;
            image[idx + 0] = (unsigned char)(color.r * 255.99f);
            image[idx + 1] = (unsigned char)(color.g * 255.99f);
            image[idx + 2] = (unsigned char)(color.b * 255.99f);
        }
    }

    // Draw pixels with GL_RGB
    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, image.data());
}

int main() {
    if (!glfwInit()) return -1;
    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "CG HW1 - Phong Ray Tracer", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glewInit();

    // Camera setup
    Camera camera;

    // Create the scene
    Scene scene;

    // -------------------------------------------------------
    // Set material for each object (example values from question)
    // -------------------------------------------------------
    // Plane P: ka=(0.2,0.2,0.2), kd=(1,1,1), ks=(0,0,0), specular power=0
    Material planeMat;
    planeMat.ka = vec3(0.2f, 0.2f, 0.2f);
    planeMat.kd = vec3(1.0f, 1.0f, 1.0f);
    planeMat.ks = vec3(0.0f, 0.0f, 0.0f);
    planeMat.shininess = 1.0f; // Not used if ks=0, but set to 1

    // Sphere S1: ka=(0.2,0,0), kd=(1,0,0), ks=(0,0,0), power=0
    Material s1Mat;
    s1Mat.ka = vec3(0.2f, 0.0f, 0.0f);
    s1Mat.kd = vec3(1.0f, 0.0f, 0.0f);
    s1Mat.ks = vec3(0.0f, 0.0f, 0.0f);
    s1Mat.shininess = 1.0f;

    // Sphere S2: ka=(0.2,0.2,0), kd=(0.5,0.5,0.5), ks=(0.5,0.5,0.5), power=32
    Material s2Mat;
    s2Mat.ka = vec3(0.2f, 0.2f, 0.0f);
    s2Mat.kd = vec3(0.5f, 0.5f, 0.5f);
    s2Mat.ks = vec3(0.5f, 0.5f, 0.5f);
    s2Mat.shininess = 32.0f;

    // Sphere S3: ka=(0,0,0.2), kd=(0,0,1), ks=(0,0,0), power=0
    Material s3Mat;
    s3Mat.ka = vec3(0.0f, 0.0f, 0.2f);
    s3Mat.kd = vec3(0.0f, 0.0f, 1.0f);
    s3Mat.ks = vec3(0.0f, 0.0f, 0.0f);
    s3Mat.shininess = 1.0f;

    // -------------------------------------------------------
    // Add objects (same coordinates as HW1 example)
    // -------------------------------------------------------
    // Plane: normal=(0,1,0), d=2 => y=-2 plane
    scene.addObject(new Plane(vec3(0, 1, 0), 2.0f, planeMat));

    // Sphere S1: center=(-4,0,-7), radius=1
    scene.addObject(new Sphere(vec3(-4, 0, -7), 1.0f, s1Mat));

    // Sphere S2: center=(0,0,-7), radius=2
    scene.addObject(new Sphere(vec3(0, 0, -7), 2.0f, s2Mat));

    // Sphere S3: center=(4,0,-7), radius=1
    scene.addObject(new Sphere(vec3(4, 0, -7), 1.0f, s3Mat));

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Rendering
        render(scene, camera);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
