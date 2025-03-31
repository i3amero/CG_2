202011378 차현준 Computer Graphics HW2 Repository<br/>

--HW1 Description-<br/>
This task does not require a separate library installation.<br/>
Run the sln files in the HW2_Q1, HW2_Q2, and HW2_Q3 folders as VS22, and then run the source file HW2.cpp.
<br/>

Q1 – Basic Ray Tracer with Object Intersection<br/>
What was implemented:<br/>
Ray-sphere and ray-plane intersection<br/>
Camera ray generation<br/>
Scene composed of three spheres and one plane<br/>
Monochrome image: white pixel if a ray hits any object, black otherwise<br/><br/>
Explanation:<br/>
In Q1, I implemented a basic ray tracer that checks whether a ray intersects with any object in the scene. The scene includes a plane and three spheres. For each pixel, a ray is generated from the camera, and if the ray hits an object, the pixel is set to white; otherwise, it is set to black. This exercise focuses solely on geometric intersection without any lighting or shading.
<br/><br/>
Q2 – Ray Tracer with Phong Illumination<br/>
What was implemented:<br/>
Phong shading model (ambient + diffuse + specular)
<br/>
A point light source<br/>
Gamma correction<br/>
Color rendering based on material properties<br/>

Explanation:<br/>
In Q2, I extended the ray tracer by incorporating the Phong illumination model. Each object is assigned a material with ambient, diffuse, and specular properties. A point light source is placed in the scene, and lighting is computed using the Phong reflection formula. Gamma correction is applied to convert linear color to display color. This step introduces realistic lighting and surface shading.
<br/><br/>
Q3 – Ray Tracer with Anti-Aliasing<br/>
What was implemented:<br/>
Stochastic sampling (64 rays per pixel)<br/>
Box filter (averaging color across samples)<br/>
Anti-aliasing to smooth jagged edges<br/>
Retained Phong shading and gamma correction<br/><br/>
Explanation:<br/>
In Q3, I added anti-aliasing to the ray tracer to reduce jagged edges and improve visual quality. For each pixel, 64 rays are randomly sampled within the pixel area. The color results are averaged using a box filter. This produces smoother transitions and more realistic rendering. The Phong shading model and gamma correction remain active, but now with improved image quality thanks to supersampling.
<br/><br/>

![image](https://github.com/user-attachments/assets/ecaeefaf-3a28-4c3e-b17f-7c572754a526)
<br/>
HW2 Q1 Image
<br/>
![image](https://github.com/user-attachments/assets/e783d878-568c-47f1-904b-d3a92939a05c)
<br/>
HW2 Q2 Image
<br/>
![image](https://github.com/user-attachments/assets/7ddd2572-9eab-4c03-8989-6e2bc93a1d65)
<br/>
HW2 Q3 Image
