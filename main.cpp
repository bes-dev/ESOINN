#include <iostream>
#include <fstream>
#include <random>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "ESOINN.h"
#include <boost/numeric/ublas/vector.hpp>

typedef boost::numeric::ublas ub;

int main() {
    soinn::ESOINN model(4, 4, 30);
    std::ifstream ifs;
    ifs.open("input.txt");
    ifs.open("iris.txt");
    ub::vector<double> input(4);
    int num;
    int tmp;
    for(int t = 0; t < 150; t++) {
        for(int i = 0; i < 4; i++) {
            ifs>>tmp;
            input(i) = double(tmp);
        }
        model.process(input);
        std::cout<<t<<std::endl;
    }
    model.classify();
    std::cout<<"number of classes: "<<model.getNumberOfClasses()<<std::endl;
    std::cout<<"number of vertices: "<<model.getNumberOfVertices()<<std::endl;
}
