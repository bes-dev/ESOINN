#include <iostream>
#include <vector>
#include <string>
#include "ESOINN.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <random>
#include <ctime>
#include <fstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace boost::numeric;
using namespace soinn;

struct Options
{
    std::string filename_load;
    std::string filename_save;
    std::string img;
    bool save_file;
    bool load_file;
    bool visible;
    int sl_iterations;
};

int options(int argc, char ** argv, Options& opts)
{
    po::options_description desc("General options");
    desc.add_options()("help,H", "Produce help message.");
    desc.add_options()("load,L", po::value<std::string>(&opts.filename_load)->default_value("soinn.xml"), "file name with a model for load");
    desc.add_options()("image,I", po::value<std::string>(&opts.img)->default_value("img.png"), "the name of the image file to the classes");
    desc.add_options()("save,S", po::value<std::string>(&opts.filename_save)->default_value("soinn.xml"), "file name with a model for save");
    desc.add_options()("fload,O", po::value<bool>(&opts.load_file)->default_value(0), "load model?(flag)");
    desc.add_options()("fsave,C", po::value<bool>(&opts.save_file)->default_value(1), "save model?(flag)");
    desc.add_options()("iteration,T", po::value<int>(&opts.sl_iterations)->default_value(3000), "the number of iterations for training the second layer");
    desc.add_options()("visible,V", po::value<bool>(&opts.visible)->default_value(1), "visible mode(1 - on; 0 - off)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    return 0;
}

void draw(cv::Mat &image, Graph graph)
{
    const cv::Scalar colors[6] = {cv::Scalar(255, 0, 0), cv::Scalar(0, 255, 0), cv::Scalar(0, 0, 255), cv::Scalar(255, 255, 0), cv::Scalar(0, 255, 255), cv::Scalar(255, 0, 255)};
    EdgeIterator edge_begin, edge_end;
    boost::tie(edge_begin, edge_end) = boost::edges(graph);
    for(EdgeIterator i = edge_begin; i != edge_end; i++)
    {
        Vertex vertex_s = boost::source(*i, graph);
        Vertex vertex_t = boost::target(*i, graph);
        cv::Point p1(graph[vertex_s].weight[0], graph[vertex_s].weight[1]);
        cv::Point p2(graph[vertex_t].weight[0], graph[vertex_t].weight[1]);
        cv::line(image, p1, p2, colors[graph[vertex_s].classId%6]);
    }
    VertexIterator vertex_begin, vertex_end;
    boost::tie(vertex_begin, vertex_end) = boost::vertices(graph);
    for(VertexIterator i = vertex_begin; i != vertex_end; i++)
    {
		int id = std::abs(graph[*i].classId)%6;
        cv::circle(image, cv::Point(graph[*i].weight[0], graph[*i].weight[1]), 3, colors[std::abs(graph[*i].classId)%6]);
    }
}


int main(int argc, char **argv)
{
	std::srand(0);
    Options opts;
    if( options(argc, argv, opts))
    {
        return 1;
    }
    soinn::ESOINN model;
    ublas::vector<double> a(2);
    ublas::vector<double> x1(2), x2(2);
    x1(0) = 320;
    x1(1) = 220;
    x2(0) = 320;
    x2(1) = 240;

    model.init(x1, x2);
    cv::Mat img = cv::imread(opts.img, 0);
    if(!opts.load_file)
    {
        std::vector<cv::Point> points;
        for(int i = 0; i < img.size().width; ++i)
        {
            for(int j = 0; j < img.size().height; ++j)
            {
                if(img.at<char>(cv::Point(i, j)) != 0)
                {
                    points.push_back(cv::Point(i, j));
                }
            }
        }

        int size = points.size();
        for(int i = 0; i < 10000; ++i)
        {
            std::cout<<i<<"\n";
            int id = std::rand() % size;
            a(0) = points[id].x;
            a(1) = points[id].y;
            model.addSignal(a);
        }
        model.classify();
		cv::Mat img1(img.size(), CV_32FC3);

        draw(img1, model.getGraph());
        cv::imshow("ESOINN", img1);
        cv::waitKey();

    }

    return 0;
}
