#include <iostream>
#include <vector>

#include "write_png.h"

int main() {
    int width = 640, height = 480; // output resolution

    std::vector<unsigned char> data(4*width*height);

    for(int x = 0;x<width;++x) {
        for(int y = 0;y<height;++y) {
            int index = 4*(y*width+x);
            data[index+0] = x/float(width)*0xFFu;       // red
            data[index+1] = y/float(height)*0xFFu;      // green
            data[index+2] = (1-x/float(width))*0xFFu;   // blue
            data[index+3] = 0xFFu;                      // alpha (opacity)
        }
    }

    write_png("boring.png", width, height, &data[0]);

    return 0;
}
