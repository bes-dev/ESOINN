/*
 * author: Sergei Belousov aka BeS
 * email: belbes122@yandex.ru
 */
#ifndef ESOINN_H
#define ESOINN_H

#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>

namespace soinn {

struct VertexProperties {
    boost::numeric::ublas::vector<double> weight;
    int classId;
    double density;
    int numberOfSignals;
    double S;
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version) {
        ar & BOOST_SERIALIZATION_NVP(weight);
        ar & BOOST_SERIALIZATION_NVP(classId);
        ar & BOOST_SERIALIZATION_NVP(density);
        ar & BOOST_SERIALIZATION_NVP(numberOfSignals);
        ar & BOOST_SERIALIZATION_NVP(S);
    }
};

struct EdgeProperties {
    int age;
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version) {
        ar & BOOST_SERIALIZATION_NVP(age);
    }
};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, boost::property<boost::vertex_index_t, size_t, VertexProperties>, EdgeProperties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;
typedef boost::graph_traits<Graph>::vertices_size_type VerticesSizeType;
typedef std::map<Vertex, int> ComponentMap;

class ESOINNException
{
public:
    ESOINNException(std::string message):message(message){}
    ~ESOINNException(){}
    inline std::string getMessage() {
        return message;
    }
private:
    std::string message;
};

class ESOINN {
public:
    ESOINN(int dim = 2, int ageMax = 30, int iterationThreshold = 50, double c1 = 0.001, double c2 = 1.0);
    ~ESOINN();
    void setParams(int dim = 2, int ageMax = 30, int iterationThreshold = 50, double c1 = 0.001, double c2 = 1.0);
    void process(const boost::numeric::ublas::vector<double>& inputSignal);
    void classify();
    Graph getGraph();
    int getNumberOfClasses();
    int getNumberOfVertices();
    boost::numeric::ublas::vector<double> getCenterOfCluster(int classId);
    VertexProperties getBestMatch(boost::numeric::ublas::vector<double>& inputSignal);
    void save(std::string filename);
    void load(std::string filename);
    void clear();
private:
    int dim;
    Graph graph;
    int ageMax;
    int iterationCount;
    int iterationThreshold;
    int numberOfClasses;
    double c1, c2;

    void addSignal(const boost::numeric::ublas::vector<double>& inputSignal);
    std::pair<Vertex,Vertex> findWinners(const boost::numeric::ublas::vector<double>& inputSignal);
    bool isWithinThreshold(const boost::numeric::ublas::vector<double>& inputSignal, Vertex& firstWinner, Vertex& secondWinner);
    double getSimilarityThreshold(const Vertex& vertex);
    void incrementEdgesAge(Vertex& vertex);
    bool needAddEdge(Vertex& firstWinner, Vertex& secondWinner);
    bool needMergeClasses(Vertex &a, Vertex &b);
    void mergeClasses(int A, int B);
    double meanDensity(int classId);
    double maxDensity(int classId);
    double densityThershold(double mean, double max);
    double meanDistance(Vertex& vertex);
    void updateDensity(Vertex& vertex);
    void updateWeights(Vertex& firstWinner, const boost::numeric::ublas::vector<double> &inputSignal);
    void deleteOldEdges();
    void updateClassLabels();
    void markClasses();
    void partitionClasses();
    void markAdjacentVertices(Vertex &vertex, int cID);
    void deleteNoiseVertex();
    double distance(const boost::numeric::ublas::vector<double> &x, const boost::numeric::ublas::vector<double> &y);
private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version) {
        ar & BOOST_SERIALIZATION_NVP(dim);
        ar & BOOST_SERIALIZATION_NVP(ageMax);
        ar & BOOST_SERIALIZATION_NVP(iterationCount);
        ar & BOOST_SERIALIZATION_NVP(iterationThreshold);
        ar & BOOST_SERIALIZATION_NVP(numberOfClasses);
        ar & BOOST_SERIALIZATION_NVP(c1);
        ar & BOOST_SERIALIZATION_NVP(c2);
        ar & BOOST_SERIALIZATION_NVP(graph);
    }
};
}

#endif // ESOINN_H
