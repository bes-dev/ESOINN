/*
 * author: Sergei Belousov aka BeS
 * email: belbes122@yandex.ru
 */
#include "ESOINN.h"

#include <fstream>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/iteration_macros.hpp>

#define E1(t) 1./t
#define E2(t) 1./(100*t)

using namespace soinn;
using namespace boost::numeric;

ESOINN::ESOINN(int dim, int ageMax, int iterationThreshold, double c1, double c2): dim(dim), ageMax(ageMax), iterationThreshold(iterationThreshold), c1(c1), c2(c2){
}

ESOINN::~ESOINN() {
}

Graph ESOINN::getGraph() {
    return graph;
}

void ESOINN::process(const boost::numeric::ublas::vector<double> &inputSignal) {
    if(inputSignal.size() != dim) {
        throw ESOINNException(std::string("Incorrect dimension of input signal in ESOINN::addSignal()."));
    } else {
        addSignal(inputSignal);
    }
}

void ESOINN::addSignal(const boost::numeric::ublas::vector<double> &inputSignal) {
    if(boost::num_vertices(graph) < 2) {
        Vertex vertex = boost::add_vertex(graph);
        graph[vertex].weight = ublas::vector<double>(inputSignal);
        graph[vertex].classId = -1;
        graph[vertex].density = 0.;
        graph[vertex].numberOfSignals = 0;
        graph[vertex].S = 0;
        return;
    }
    Vertex firstWinner, secondWinner;
    boost::tie(firstWinner, secondWinner) = findWinners(inputSignal);
    if(!isWithinThreshold(inputSignal, firstWinner, secondWinner)){
        Vertex vertex = boost::add_vertex(graph);
        graph[vertex].weight = ublas::vector<double>(inputSignal);
        graph[vertex].classId = -1;
        graph[vertex].density = 0.;
        graph[vertex].numberOfSignals = 0;
        graph[vertex].S = 0;
        return;
    }
    incrementEdgesAge(firstWinner);
    if(needAddEdge(firstWinner, secondWinner)) {
        Edge e = boost::add_edge(firstWinner, secondWinner, graph).first;
        graph[e].age = 0;
    } else {
        boost::remove_edge(firstWinner, secondWinner, graph);
    }
    updateDensity(firstWinner);
    updateWeights(firstWinner, inputSignal);
    deleteOldEdges();
    if(iterationCount % iterationThreshold == 0) {
        updateClassLabels();
    }
    iterationCount++;
}

double ESOINN::distance(const boost::numeric::ublas::vector<double> &x, const boost::numeric::ublas::vector<double> &y){
    return ublas::norm_2( x - y );
}

std::pair<Vertex,Vertex> ESOINN::findWinners(const boost::numeric::ublas::vector<double> &inputSignal) {
    Vertex firstWinner = NULL;
    Vertex secondWinner = NULL;
    double firstWinnerDistance = std::numeric_limits<double>::max();
    double secondWinnerDistance = std::numeric_limits<double>::max();
    VertexIterator current, end;
    boost::tie(current, end) = boost::vertices(graph);
    for(; current != end; current++) {
        double dist = distance(inputSignal, graph[*current].weight);
        if(dist < firstWinnerDistance) {
            secondWinner = firstWinner;
            secondWinnerDistance = firstWinnerDistance;
            firstWinner = *current;
            firstWinnerDistance = dist;
        } else if(dist < secondWinnerDistance) {
            secondWinner = *current;
            secondWinnerDistance = dist;
        }
    }
    return std::pair<Vertex,Vertex>(firstWinner, secondWinner);
}

bool ESOINN::isWithinThreshold(const boost::numeric::ublas::vector<double>& inputSignal, Vertex& firstWinner, Vertex& secondWinner) {
    if(distance(inputSignal, graph[firstWinner].weight) > getSimilarityThreshold(firstWinner)) {
        return false;
    }
    if(distance(inputSignal, graph[secondWinner].weight) > getSimilarityThreshold(secondWinner)) {
        return false;
    }
    return true;
}

double ESOINN::getSimilarityThreshold(const Vertex& vertex) {
    double dist = 0.0;
    if(!boost::out_degree(vertex, graph)) {
        dist = std::numeric_limits<double>::max();
        VertexIterator current, end;
        boost::tie(current, end) = boost::vertices(graph);
        for(; current != end; current++) {
            if(*current != vertex) {
                double distCurrent = distance(graph[vertex].weight, graph[*current].weight);
                if(distCurrent < dist) {
                    dist = distCurrent;
                }
            }
        }
    } else {
        dist = std::numeric_limits<double>::min();
        AdjacencyIterator current, end;
        boost::tie(current, end) = boost::adjacent_vertices(vertex, graph);
        for(; current != end; current++) {
            double distCurrent = distance(graph[vertex].weight, graph[*current].weight);
            if(distCurrent > dist) {
                dist = distCurrent;
            }
        }
    }
    return dist;
}

void ESOINN::incrementEdgesAge(Vertex& vertex) {
    OutEdgeIterator current, end;
    boost::tie(current, end) = boost::out_edges(vertex, graph);
    for(; current != end; current++) {
        graph[*current].age++;
    }
}

bool ESOINN::needAddEdge(Vertex& firstWinner, Vertex &secondWinner) {
    if(graph[firstWinner].classId == -1 || graph[secondWinner].classId == -1) {
        return true;
    } else if(graph[firstWinner].classId == graph[secondWinner].classId) {
        return true;
    } else if(graph[firstWinner].classId != graph[secondWinner].classId && needMergeClasses(firstWinner, secondWinner)) {
        return true;
    }
    return false;
}

bool ESOINN::needMergeClasses(Vertex &a, Vertex &b) {
    int A = graph[a].classId;
    double meanA = meanDensity(A);
    double maxA = maxDensity(A);
    double thresholdA = densityThershold(meanA, maxA);
    int B = graph[b].classId;
    double meanB = meanDensity(B);
    double maxB = maxDensity(B);
    double thresholdB = densityThershold(meanB, maxB);
    double minAB = std::min(graph[a].density, graph[b].density);
    if(minAB > thresholdA * maxA && minAB > thresholdB * maxB) {
        return true;
    }
    return false;
}

void ESOINN::mergeClasses(int A, int B) {
    int classId = std::min(A, B);
    VertexIterator current, end;
    boost::tie(current, end) = boost::vertices(graph);
    for(; current != end; current++) {
        if(graph[*current].classId == A || graph[*current].classId == B) {
            graph[*current].classId = classId;
        }
    }
}

double ESOINN::meanDensity(int classId) {
    if(classId == -1) return 0.0;
    int n = 0;
    double density = 0.0;
    VertexIterator current, end;
    boost::tie(current, end) = boost::vertices(graph);
    for(; current != end; current++) {
        if(graph[*current].classId == classId) {
            n++;
            density += graph[*current].density;
        }
    }
    density *= 1./double(n);
    return density;
}

double ESOINN::maxDensity(int classId) {
    double density = std::numeric_limits<double>::min();
    VertexIterator current, end;
    boost::tie(current, end) = boost::vertices(graph);
    for(; current != end; current++) {
        if(graph[*current].density > density && graph[*current].classId == classId) {
            density = graph[*current].density;
        }
    }
    return density;
}

double ESOINN::densityThershold(double mean, double max) {
    double threshold;
    if(2.0 * mean >= max) {
        threshold = 0.0;
    } else if(3.0 * mean >= max && max > 2.0 * mean) {
        threshold = 0.5;
    } else {
        threshold = 1.0;
    }
    return threshold;
}

void ESOINN::updateDensity(Vertex& vertex) {
    double mDistance = meanDistance(vertex);
    graph[vertex].numberOfSignals++;
    graph[vertex].S += 1./((1 + mDistance)*(1 + mDistance));
    graph[vertex].density = graph[vertex].S/double(graph[vertex].numberOfSignals);
}

void ESOINN::updateWeights(Vertex& firstWinner, const boost::numeric::ublas::vector<double> &inputSignal) {
    graph[firstWinner].weight += E1(graph[firstWinner].numberOfSignals) * (inputSignal - graph[firstWinner].weight);
    AdjacencyIterator current, end;
    boost::tie(current, end) = boost::adjacent_vertices(firstWinner, graph);
    for(; current != end; current++) {
        graph[*current].weight += E2(graph[firstWinner].numberOfSignals) * (inputSignal - graph[*current].weight);
    }
}

double ESOINN::meanDistance(Vertex& vertex) {
    double mDistance = 0.0;
    int m = 0;
    VertexIterator current, end;
    boost::tie(current, end) = boost::vertices(graph);
    for(; current != end; current++) {
        mDistance += distance(graph[vertex].weight, graph[*current].weight);
        m++;
    }
    mDistance *= 1./double(m);
    return mDistance;
}

void ESOINN::deleteOldEdges() {
    EdgeIterator current, end;
    boost::tie(current, end) = boost::edges(graph);
    EdgeIterator next = current;
    for(;current != end; current = next) {
        next ++;
        if(graph[*current].age > ageMax) {
            Vertex vertexS = boost::source(*current, graph);
            Vertex vertexT = boost::target(*current, graph);
            boost::remove_edge(*current, graph);
        }
    }
}

void ESOINN::updateClassLabels() {
    markClasses();
    partitionClasses();
    deleteNoiseVertex();
}

void ESOINN::markClasses() {
    std::list<VertexIterator> vertexList;
    VertexIterator begin, end;
    boost::tie(begin, end) = boost::vertices(graph);
    for(VertexIterator current = begin; current != end; current++) {
        graph[*current].classId = -1;
        vertexList.push_back(current);
    }
    vertexList.sort([&](VertexIterator &a, VertexIterator &b) -> bool {
        if(graph[*a].density > graph[*b].density) return true;
        return false;
    });
    int classCount = 0;
    for(std::list<VertexIterator>::iterator current = vertexList.begin(); current != vertexList.end(); current++) {
        if(graph[**current].classId == -1) {
            graph[**current].classId = classCount;
            markAdjacentVertices(**current, classCount++);
        }
    }
}

void ESOINN::partitionClasses() {
    EdgeIterator current, end;
    boost::tie(current, end) = boost::edges(graph);
    EdgeIterator next = current;
    for(;current != end; current = next) {
        next ++;
        Vertex vertexS = boost::source(*current, graph);
        Vertex vertexT = boost::target(*current, graph);
        if(graph[vertexS].classId != graph[vertexT].classId) {
            if(needMergeClasses(vertexS, vertexT)) {
                mergeClasses(graph[vertexS].classId, graph[vertexT].classId);
            } else {
                boost::remove_edge(*current, graph);
            }
        }
    }
}

void ESOINN::markAdjacentVertices(Vertex &vertex, int cID) {
    AdjacencyIterator current, end;
    boost::tie(current, end) = boost::adjacent_vertices(vertex, graph);
    for(; current != end; current++){
        if(graph[*current].classId == -1 && graph[*current].density < graph[vertex].density) {
            graph[*current].classId = cID;
            Vertex v = *current;
            markAdjacentVertices(v, cID);
        }
    }
}

void ESOINN::deleteNoiseVertex() {
    VertexIterator begin, end;
    boost::tie(begin, end) = boost::vertices(graph);
    VertexIterator next = begin;
    for(VertexIterator current = begin; current != end; current = next) {
        next++;
        double mean = meanDensity(graph[*current].classId);
        if((boost::out_degree(*current, graph) == 2 && graph[*current].density < c1* mean) ||
                (boost::out_degree(*current, graph) == 1 && graph[*current].density < c2* mean) ||
                (boost::out_degree(*current, graph) == 0)) {
            boost::clear_vertex(*current, graph);
            boost::remove_vertex(*current, graph);
        }
    }
}

void ESOINN::classify() {
    deleteNoiseVertex();
    size_t index = 0;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        boost::put(boost::vertex_index, graph, v, index++);
    }
    ComponentMap component;
    boost::associative_property_map<ComponentMap> componentMap(component);
    numberOfClasses = connected_components(graph, componentMap);
    BGL_FORALL_VERTICES(v, graph, Graph) {
        graph[v].classId = boost::get(componentMap, v);
    }
}

void ESOINN::save(std::string filename) {
    std::ofstream ofs(filename.c_str());
    boost::archive::xml_oarchive oa(ofs);
    oa << BOOST_SERIALIZATION_NVP(this);
}

void ESOINN::load(std::string filename) {
    clear();
    std::ifstream ifs(filename.c_str());
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(*const_cast<ESOINN*>(this));
}

void ESOINN::clear() {
    graph.clear();
    numberOfClasses = 0;
}

void ESOINN::setParams(int dim, int ageMax, int iterationThreshold, double c1, double c2) {
    this->dim = dim;
    this->ageMax = ageMax;
    this->iterationThreshold = iterationThreshold;
    this->c1 = c1;
    this->c2 = c2;
}

int ESOINN::getNumberOfClasses() {
    return numberOfClasses;
}

int ESOINN::getNumberOfVertices() {
    return boost::num_vertices(graph);
}

boost::numeric::ublas::vector<double> ESOINN::getCenterOfCluster(int classId) {
    double density = -1;
    Vertex center;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(graph[v].classId == classId && graph[v].density > density) {
            center = v;
            density = graph[center].density;
        }
    }
    return graph[center].weight;
}

VertexProperties ESOINN::getBestMatch(boost::numeric::ublas::vector<double>& inputSignal) {
    Vertex firstWinner = NULL;
    double firstWinnerDistance = std::numeric_limits<double>::max();
    VertexIterator current, end;
    boost::tie(current, end) = boost::vertices(graph);
    for(; current != end; current++) {
        double dist = distance(inputSignal, graph[*current].weight);
        if(dist < firstWinnerDistance) {
            firstWinner = *current;
            firstWinnerDistance = dist;
        }
    }
    return graph[firstWinner];
}
