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
#include <cstdlib>
#include <ctime>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;

const int WIDTH = 512;
const int HEIGHT = 512;
const float INF = std::numeric_limits<float>::infinity();
const vec3 lightPos = vec3(-4.0f, 4.0f, -3.0f);
const vec3 lightColor = vec3(1.0f);
const float GAMMA = 2.2f;
const int SAMPLES = 4; // Reduce for testing speed

struct Material {
    vec3 ka, kd, ks;
    float shininess;
};

struct Ray {
    vec3 origin, direction;
    Ray(const vec3& o, const vec3& d) : origin(o), direction(normalize(d)) {}
};

struct HitInfo {
    bool hit;
    float t;
    vec3 normal, hitPoint;
    Material material;
};

struct Surface {
    Material material;
    virtual bool intersect(const Ray& ray, float& t) const = 0;
    virtual vec3 getNormal(const vec3& hitPoint) const = 0;
};

struct Plane : public Surface {
    vec3 normal;
    float d;
    Plane(const vec3& n, float d, const Material& m) {
        normal = normalize(n); this->d = d; material = m;
    }
    bool intersect(const Ray& ray, float& t) const override {
        float denom = dot(normal, ray.direction);
        if (fabs(denom) < 1e-6f) return false;
        t = -(dot(normal, ray.origin) + d) / denom;
        return t >= 0.0f;
    }
    vec3 getNormal(const vec3&) const override { return normal; }
};

struct Sphere : public Surface {
    vec3 center;
    float radius;
    Sphere(const vec3& c, float r, const Material& m) {
        center = c; radius = r; material = m;
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
        t = (t0 > 0.0f) ? t0 : ((t1 > 0.0f) ? t1 : -1);
        return t > 0;
    }
    vec3 getNormal(const vec3& hitPoint) const override {
        return normalize(hitPoint - center);
    }
};

struct Camera {
    vec3 eye, u, v, w;
    float l, r, b, t, d;
    Camera() {
        eye = vec3(0, 0, 0);
        w = normalize(vec3(0, 0, 1));     // Z- 쪽 보기
        u = vec3(1, 0, 0);
        v = vec3(0, 1, 0);
        l = -1.0f;
        r = 1.0f;
        b = -1.0f;
        t = 1.0f;
        d = 1.0f;
    }

    Ray generateRay(float i, float j, float offsetX = 0.5f, float offsetY = 0.5f) {
        float px = l + (r - l) * (i + offsetX) / WIDTH;
        float py = b + (t - b) * (j + offsetY) / HEIGHT;
        vec3 dir = -d * w + px * u + py * v;
        return Ray(eye, dir);
    }
};

struct Scene {
    std::vector<Surface*> objects;
    void addObject(Surface* obj) { objects.push_back(obj); }
    bool intersect(const Ray& ray, HitInfo& hitInfo) const {
        float closestT = INF;
        bool isHit = false;
        for (auto obj : objects) {
            float t;
            if (obj->intersect(ray, t) && t < closestT) {
                closestT = t;
                hitInfo.hit = true;
                hitInfo.t = t;
                hitInfo.hitPoint = ray.origin + t * ray.direction;
                hitInfo.normal = obj->getNormal(hitInfo.hitPoint);
                hitInfo.material = obj->material;
                isHit = true;
            }
        }
        return isHit;
    }
};

bool inShadow(const Scene& scene, const vec3& point, const vec3& lightDir) {
    Ray shadowRay(point + 0.001f * lightDir, lightDir);
    HitInfo temp;
    if (scene.intersect(shadowRay, temp)) {
        float distToHit = length(temp.hitPoint - point);
        float distToLight = length(lightPos - point);
        return distToHit < distToLight;
    }
    return false;
}

vec3 phongShading(const Scene& scene, const Ray& ray) {
    HitInfo hit;
    vec3 color(0);
    if (scene.intersect(ray, hit)) {
        vec3 ambient = hit.material.ka * lightColor;
        vec3 L = normalize(lightPos - hit.hitPoint);
        bool shadowed = inShadow(scene, hit.hitPoint, L);
        vec3 diffuse(0), specular(0);
        if (!shadowed) {
            float nl = dot(hit.normal, L);
            if (nl > 0) diffuse = hit.material.kd * nl * lightColor;
            vec3 R = reflect(-L, hit.normal);
            vec3 V = normalize(-ray.direction);
            float rv = dot(R, V);
            if (rv > 0 && length(hit.material.ks) > 0)
                specular = hit.material.ks * pow(rv, hit.material.shininess);
        }
        color = clamp(ambient + diffuse + specular, 0.0f, 1.0f);
    }
    return color;
}

void render(const Scene& scene, Camera& camera) {
    std::vector<unsigned char> image(WIDTH * HEIGHT * 3, 0);
    std::srand(static_cast<unsigned>(std::time(0)));
    for (int j = 0; j < HEIGHT; ++j) {
        for (int i = 0; i < WIDTH; ++i) {
            vec3 color(0);
            for (int s = 0; s < SAMPLES; ++s) {
                float offsetX = static_cast<float>(rand()) / RAND_MAX;
                float offsetY = static_cast<float>(rand()) / RAND_MAX;
                Ray ray = camera.generateRay((float)i, (float)j, offsetX, offsetY);
                color += phongShading(scene, ray);
            }
            color /= float(SAMPLES);
            color = clamp(color, 0.0f, 1.0f);
            color.r = pow(color.r, 1.0f / GAMMA);
            color.g = pow(color.g, 1.0f / GAMMA);
            color.b = pow(color.b, 1.0f / GAMMA);
            int idx = (j * WIDTH + i) * 3;
            image[idx + 0] = static_cast<unsigned char>(color.r * 255.99f);
            image[idx + 1] = static_cast<unsigned char>(color.g * 255.99f);
            image[idx + 2] = static_cast<unsigned char>(color.b * 255.99f);
        }
    }
    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, image.data());
}

int main() {
    if (!glfwInit()) return -1;
    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "Ray Tracer", NULL, NULL);
    if (!window) return -1;
    glfwMakeContextCurrent(window);
    glewInit();
    glViewport(0, 0, WIDTH, HEIGHT);
    glClearColor(0, 0, 0, 1);

    Camera camera;
    Scene scene;

    Material planeMat = { vec3(0.2f), vec3(1.0f), vec3(0.0f), 1.0f };
    Material s1Mat = { vec3(0.2f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f), vec3(0.0f), 1.0f };
    Material s2Mat = { vec3(0.0f, 0.2f, 0.0f), vec3(0.5f), vec3(0.5f), 32.0f };
    Material s3Mat = { vec3(0.0f, 0.0f, 0.2f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f), 1.0f };

    scene.addObject(new Plane(vec3(0, 1, 0), 2.0f, planeMat));
    scene.addObject(new Sphere(vec3(-4, 0, -7), 1.0f, s1Mat));
    scene.addObject(new Sphere(vec3(0, 0, -7), 2.0f, s2Mat));
    scene.addObject(new Sphere(vec3(4, 0, -7), 1.0f, s3Mat));

    // Render only once
    glClear(GL_COLOR_BUFFER_BIT);
    render(scene, camera);
    glfwSwapBuffers(window);

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}