#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "write_png.h"

// simple distance function raymarcher
// see also: http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

// useful functions
void normalize(double &x, double &y, double &z) {
    double invsqrt = 1.0/std::sqrt(x*x+y*y+z*z);
    x *= invsqrt; y *= invsqrt; z *= invsqrt;
}
double clamp(double x, double a, double b) {
    return x<a?a:x>b?b:x;
}
double smoothstep(double x) {
    x = clamp(x, 0, 1);
    return (3 - 2*x)*x*x;
}
// smooth min and max
double smin(double a, double b, double width) {
    double diff = (a-b)/width;
    return a*smoothstep(0.5-diff)+b*smoothstep(0.5+diff);
}
double smax(double a, double b, double width) {
    double diff = (a-b)/width;
    return a*smoothstep(0.5+diff)+b*smoothstep(0.5-diff);
}
// distance functions for various shapes
double sphere(double x, double y, double z, double cx, double cy, double cz, double r) {
    x-=cx; y-=cy; z-=cz;
    return std::sqrt(x*x+y*y+z*z) - r;
}
double octahedron(double x, double y, double z, double cx, double cy, double cz, double r) {
    x-=cx; y-=cy; z-=cz;
    return (std::abs(x)+std::abs(y)+std::abs(z))*0.57735026919 - r;
}
double box(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1, double r) {
    x-=clamp(x, x0, x1); y-=clamp(y, y0, y1); z-=clamp(z, z0, z1);
    return std::sqrt(x*x+y*y+z*z) - r;
}
double bevel_box(double x, double y, double z, double x0, double x1, double y0, double y1, double z0, double z1, double r) {
    x-=clamp(x, x0, x1); y-=clamp(y, y0, y1); z-=clamp(z, z0, z1);
    return (std::abs(x)+std::abs(y)+std::abs(z))*0.57735026919 - r;
}

// distance function for non reflecting surfaces
double diffuse(double x, double y, double z) {
    double d = std::sqrt(x*x+z*z);
    return std::min(
            std::min(
                sphere(x,y,z, 3,0,3, 0.7),
                sphere(x,y,z, -1.5f,0,0, 1.0)
            ),
            std::min(
                y+2+0.5*(std::cos(d)*exp(-0.01*d*d)),
                bevel_box(x,y,z, -4,-3,1,2,2,3, 0.4)
            )
        );
}
// distance function for reflecting surfaces
double mirror(double x, double y, double z) {
    return std::min(
            octahedron(x,y,z, 0,3,3, 1.5),
            box(x,y,z, 4,5,2,3,3,10, 0.5)
        );
}

struct Material {
    double (*distance)(double, double, double); // distance function
    double color[3]; // surface color
    double specular; // specular strength
    double diffuse;  // diffuse strength
    double mirror;   // reflection strength
};
// material definitions
Material materials[2] = {
    {diffuse, {1,0.5,0.5}, 0.9, 0.1, 0.0},
    {mirror, {0.2,0.2,0.2}, 0.1, 0.2, 1.0}
};
int materialcount = 2;
// point light sources {x,y,z, r,g,b}
double lights[][6] = {
    {300,100,-300, 100000,100000,100000},
    {10,10,0, 100,100,100}
};
int lightcount = 2;
// distance function of all materials
double all(double x, double y, double z) {
    double mindist = 1.e10;
    for(int i = 0;i<materialcount;++i) {
        double dist = materials[i].distance(x,y,z);
        if(dist<mindist) mindist = dist;
    }
    return mindist;
}
// approximates the gradient of the distance function with finite differences
void grad(double x, double y, double z, double &dx, double &dy, double &dz) {
    double eps = 1.e-5f*(std::abs(x)+std::abs(y)+std::abs(z));
    double inv2eps = 0.5f/eps;
    dx = (all(x+eps,y,z)-all(x-eps,y,z))*inv2eps;
    dy = (all(x,y+eps,z)-all(x,y-eps,z))*inv2eps;
    dz = (all(x,y,z+eps)-all(x,y,z-eps))*inv2eps;
}
// d -= n*dot(d,n), vector reflection
void reflect(double &dx, double &dy, double &dz, double nx, double ny, double nz) {
    double dot = dx*nx+dy*ny+dz*nz;
    dx -= 2*dot*nx; dy -= 2*dot*ny; dz -= 2*dot*nz;
}
// this gets called if a ray "hits" the sky
void sky(double dx, double dy, double dz, double &r, double &g, double &b) {
    (void)dx; (void)dy; (void)dz;
    r += 0x87/255.0; g += 0xCE/255.0; b += 0xEB/255.0;
    //~ r += 0; g += 0; b += 0;
}
// partial trace function for shadow calculation
double trace(double &x, double &y, double &z, double dx, double dy, double dz, double max_dist = 1.e6) {
    const int max_depth = 16*1024;
    const double precision = 1.e-10f;

    double d_tot = 0;
    double d = all(x,y,z);
    int i = 0;
    for(;i<max_depth && d>precision && d<max_dist;++i) {
        d_tot += d;
        x += (0.8*d)*dx; y += (0.8*d)*dy; z += (0.8*d)*dz;
        d = all(x,y,z);
    }
    return std::min(d_tot, max_dist);
}
// full recursive trace function 
void trace(double &x, double &y, double &z, double dx, double dy, double dz, double &r, double &g, double &b, int recursion = 10) {
    if(recursion == 0) return;

    const int max_depth = 16*1024;
    const double max_dist = 1.e6;
    const double precision = 1.e-10f;
    
    // iterate along ray until we hit a surface or exceed limits
    double d = all(x,y,z);
    int i = 0;
    for(;i<max_depth && d>precision && d<max_dist;++i) {
        x += (0.8*d)*dx; y += (0.8*d)*dy; z += (0.8*d)*dz;
        d = all(x,y,z);
    }

    if(i == max_depth || d>=max_dist) { // didn't hit a surface -> sky
        sky(dx, dy, dz, r, g, b);
    } else {
        double nx, ny, nz;
        grad(x, y, z, nx, ny, nz);
        normalize(nx, ny, nz);
        
        // figure out which material we hit
        int minindex = 0;
        double mindist = materials[minindex].distance(x,y,z);;
        for(int j = 1;j<materialcount;++j) {
            double dist = materials[j].distance(x,y,z);
            if(dist<mindist) {
                mindist = dist;
                minindex = j;
            }
        }

        // add lights
        for(int i = 0;i<lightcount;++i) {
            double lx = lights[i][0]-x;
            double ly = lights[i][1]-y;
            double lz = lights[i][2]-z;
            double dist2 = lx*lx+ly*ly+lz*lz;
            double dist = std::sqrt(dist2);
            double attenuation = 1.0/dist2;
            double bright = 0;
            
            // lambert diffuse brdf
            normalize(lx, ly, lz);
            double dot = lx*nx+ly*ny+lz*nz;
            bright += materials[minindex].diffuse*std::max(0.0, dot);

            // phong specular brdf
            double dx2 = dx, dy2 = dy, dz2 = dz;
            reflect(dx2, dy2, dz2, nx, ny, nz);
            double dot2 = dx2*lx+dy2*ly+dz2*lz;
            bright += materials[i].specular*std::pow(std::max(0.0, dot2), 64);

            // check if there is an object obstructing the light
            double x1 = x+10.0*precision*lx, y1 = y+10.0*precision*ly, z1 = z+10.0*precision*lz;
            double free_dist = trace(x1, y1, z1, lx, ly, lz, dist);
            bright *= dist-free_dist>1-2*precision?0:1;
            bright *= attenuation;

            r += materials[minindex].color[0]*bright*lights[i][3];
            g += materials[minindex].color[1]*bright*lights[i][4];
            b += materials[minindex].color[2]*bright*lights[i][5];
        }
        // if we hit a mirror surface continue tracing in the reflected direction
        if(materials[minindex].mirror>0) {
            double dx2 = dx, dy2 = dy, dz2 = dz;
            reflect(dx2, dy2, dz2, nx, ny, nz);
            x+=10.0*precision*dx2; y+=10.0*precision*dy2; z+=10.0*precision*dz2;
            double r1 = 0, g1= 0, b1 = 0;
            trace(x,y,z, dx2,dy2,dz2, r1,g1,b1, recursion-1);
            r1 *= materials[minindex].mirror*materials[minindex].color[0];
            g1 *= materials[minindex].mirror*materials[minindex].color[1];
            b1 *= materials[minindex].mirror*materials[minindex].color[2];
            r += r1; g += g1; b += b1;
        }
    }
}
// supersampling offsets
double samples[][2] = {
    {0,0},
    {0.33,0.33},{-0.33,0.33},{-0.33,-0.33},{0.33,-0.33},
    {0.33,0},{-0.33,0},{0,0.33},{0,-0.33},
};

double gamma(double x) {
    return std::pow(x, 1.0/2.2);
}

int main() {
    int width = 640, height = 480; // output resolution
    int samplecount = 9; // 1 = no antialiasing, 9 is higher quality than 5 (only 1,5,9 make sense)

    double aspect = double(width)/height;
    std::vector<unsigned char> data(4*width*height);
    
    // dynamic scheduling since the workload is highly irregular
    #pragma omp parallel for schedule(dynamic, 1)
    for(int x = 0;x<width;++x) {
        for(int y = 0;y<height;++y) {
            double r=0.0, g=0.0, b=0.0;
            for(int s = 0;s<samplecount;++s) {
                double fx = 0, fy = 0, fz = -5;
                double dx = (2.0*(x+samples[s][0])/width-1.0)*aspect;
                double dy = 1.0-2.0*(y+samples[s][1])/height;
                double dz = 1.5; // changing this changes the field of view
                normalize(dx, dy, dz);
                trace(fx, fy, fz, dx, dy, dz, r, g, b);
            }
            r/=samplecount; g/=samplecount; b/=samplecount;
            
            // actually relevant part for image writing:
            // in this case we write doubles which might not be in [0,1]
            // but the color data we write is one byte per channel
            // so we have to clamp and rescale
            int index = 4*(y*width+x);
            data[index+0] = clamp(gamma(r), 0.0, 1.0)*0xFFu; // red
            data[index+1] = clamp(gamma(g), 0.0, 1.0)*0xFFu; // green
            data[index+2] = clamp(gamma(b), 0.0, 1.0)*0xFFu; // blue
            data[index+3] = 0xFFu; // alpha (opacity)
        }
    }

    write_png("fun.png", width, height, &data[0]);

    return 0;
}
