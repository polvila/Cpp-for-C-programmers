//
//  homework2.cpp
//  Homework2
//
//  Created by Pol Vilà Saló on 22/3/16.
//  Copyright © 2016 Pol Vilà Saló. All rights reserved.
//
//  This program implements a Monte Carlo simulation that calculates the average shortest path in a graph.
//  The graphs are generated randomly. The shortest path is found by the Dijkstra's algorithm.
//
//  The data structure used to store the graph is an edge list (with vertex, weigh pairs) because a graph with
//  densities between 20% and 40%, I think that it's a sparse graph, for this reason, it's worst to use
//  a connectivity matrix.
//
//  The vertices of the graph use the vector position as id.
//
//  Note: This should be compiled with C++11 enable!

#include <iostream>
#include <vector>
#include <random>

using namespace std;

class PairVertexValue
{
public:
    // This class stores the vertex, value pairs
    //
    // @param vertex id of the vertex
    // @param value weight of the edge
    PairVertexValue(unsigned int vertex, double value) : vertex_(vertex), value_(value) {
        if(value < 0)
            throw invalid_argument("This value must be positive.");
    }

    unsigned int getVertex_() const {
        return vertex_;
    }

    double getValue_() const {
        return value_;
    }

    void setVertex_(unsigned int vertex_) {
        PairVertexValue::vertex_ = vertex_;
    }

    void setValue_(double value_) {
        PairVertexValue::value_ = value_;
    }

    bool operator<(PairVertexValue element) {
        return value_ < element.value_;
    }
    bool operator>(PairVertexValue element) {
        return value_ > element.value_;
    }

    bool operator== (const PairVertexValue& element) const {
        return vertex_ == element.getVertex_() && value_ == element.getValue_();
    }

    PairVertexValue clone() const{
        return PairVertexValue(getVertex_(),getValue_());
    }

private:
    unsigned int vertex_;
    double value_;
};

class Graph
{
public:
    // Uses a vector (each position is a vertex) of vectors (each position is a vertex, weight pair) as adjacency
    // structure. This class implements an undirected graph.
    //
    // @param num_vertices number of vertices
    Graph(unsigned int num_vertices) : v_(num_vertices){
        for (int i = 0; i < num_vertices; i++)
            adj_.push_back(vector<PairVertexValue>());       //Initialize the graph_ with empty edge vectors
    }

    // @return the number of vertices in the graph
    unsigned int V() const { return v_; }

    // @return the number of edges in the graph
    unsigned int E() const { return e_; }

    // Tests whether there is an edge from node x to node y
    bool adjacent(unsigned int x, unsigned int y) const {
        for(vector<PairVertexValue>::size_type i = 0; i != adj_[x].size(); i++)
            if(adj_[x][i].getVertex_() == y)
                return true;
        return false;
    }

    // Lists all vertex y (and its weight) such that there is an edge from x to y
    vector<PairVertexValue> neighbors(unsigned int x) const {
        vector<PairVertexValue> result;
        for(vector<PairVertexValue>::size_type i = 0; i != adj_[x].size(); i++)
            result.push_back(adj_[x][i].clone());
        return result;
    }

    // Adds to G the edge from x to y, if it is not there
    void add_edge(unsigned int x, unsigned int y, double weight) {
        if(!adjacent(x, y)){
            adj_[x].push_back(PairVertexValue(y, weight));
            adj_[y].push_back(PairVertexValue(x, weight));
            e_++;
        }
    }

    // Removes the edge from x to y, if it is there
    void delete_edge(unsigned int x, unsigned int y) {
        for(vector<PairVertexValue>::iterator it = adj_[x].begin(); it != adj_[x].end(); it++)
            if(it->getVertex_() == y)
                adj_[x].erase(it);
    }

    // @return the value associated to the edge (x,y)
    double get_edge_value(unsigned int x, unsigned int y) const {
        for(vector<PairVertexValue>::size_type i = 0; i != adj_[x].size(); i++)
            if(adj_[x][i].getVertex_() == y)
                return adj_[x][i].getValue_();
        return numeric_limits<double>::max();
    }

    // Sets the value associated to the edge (x,y) to v
    void set_edge_value(unsigned int x, unsigned int y, double v) {
        for(vector<PairVertexValue>::size_type i = 0; i != adj_[x].size(); i++)
            if(adj_[x][i].getVertex_() == y)
                adj_[x][i].setValue_(v);
    }

private:
    const unsigned int v_;
    unsigned int e_;
    vector< vector<PairVertexValue> > adj_;
};

// PriorityQueue is to always have access to the vertex with the next shortest link in the shortest
// path calculation at the top of the queue.
// This class implements a min heap to get the next min distance vertex at the top of the list.
// Reference: https://en.wikipedia.org/wiki/Binary_heap
class PriorityQueue
{
public:

    // Removes the top element of the queue
    void minPriority() {
        int pointer = 0;
        minHeap[0] = minHeap.back();
        minHeap.pop_back();
        while(true){
            int left_node = pointer*2+1;
            int right_node = pointer*2+2;
            int next_node;
            if(left_node >= minHeap.size()) break;
            next_node = left_node;
            if(right_node < minHeap.size() && minHeap[right_node] < minHeap[left_node])
                next_node = right_node;
            swap(minHeap[pointer], minHeap[next_node]);
            pointer = next_node;
        }
    }

    // Insert queue_element into queue
    void insert(PairVertexValue element) {
        minHeap.push_back(element);
        int pointer = size()-1;
        int father = (pointer-1)/2;
        while(pointer > 0 && minHeap[father] > minHeap[pointer]){
            swap(minHeap[father], minHeap[pointer]);
            pointer = father;
            father = (pointer-1)/2;
        }
    }

    // @return the top element of the queue
    PairVertexValue top() {
        return minHeap.front();
    }

    // @return the number of queue_elements
    int size() {
        return static_cast<int>(minHeap.size());
    }

private:
    vector<PairVertexValue> minHeap;
};

// This class implements the mechanics of Dijkstra’s algorithm
class ShortestPath
{
public:
    // This class permits apply the Dijkstra algorithm to the given graph
    //
    // @param graph The graph that you want to apply Dijkstra algorithm.
    ShortestPath(Graph graph) :
            graph_(graph),
            min_distances_(graph_.V(), numeric_limits<double>::max()) { }

    // This method implements the dijkstra algorithm
    // Reference: http://en.wikipedia.org/wiki/Dijkstra's_algorithm
    //
    // @return the path cost associated with the shortest path between u and w
    double path_size(unsigned int u, unsigned int w) {
        PriorityQueue priorityQueue;

        priorityQueue.insert(PairVertexValue(u, 0));

        while(priorityQueue.size()>0) {
            //Gets the top node
            PairVertexValue node = priorityQueue.top();
            //Delete the top node in the priorityQueue
            priorityQueue.minPriority();

            if(node.getValue_() < min_distances_[node.getVertex_()])
                min_distances_[node.getVertex_()] = node.getValue_();

            vector<PairVertexValue> neighbours = graph_.neighbors(node.getVertex_());       //Gets the u neighbours
            for(int i = 0; i != neighbours.size(); i++) {
                double dist = min_distances_[node.getVertex_()] +
                        graph_.get_edge_value(node.getVertex_(), neighbours[i].getVertex_());
                if(dist < min_distances_[neighbours[i].getVertex_()]) {
                    priorityQueue.insert(PairVertexValue(neighbours[i].getVertex_(), dist));
                }
            }
        }

        return min_distances_[w];
    }

    vector<double> get_min_distances_(){
        return min_distances_;
    }

private:
    Graph graph_;
    vector<double> min_distances_;
};

class MonteCarloSimulation
{
public:
    // Calculates the average distance of the shortest path of the all trials
    //
    // @param trials number of trials to run the simulation
    // @param density of the randomly generated graph
    // @param vertex number of vertices of the generetad graph
    // @param min_weight
    // @param max_weight    range of the random weight
    MonteCarloSimulation(int trials = 20, double density = 0.2, unsigned int vertex = 50, int min_weight = 1, int max_weight = 10) :
            density_(density), vertex_(vertex), e_(time(0)), weight_distribution_(min_weight, max_weight),
            edge_distribution_(0.0, 1.0) {

        double accum = 0;
        for (int i = 0; i < trials; i++)
            accum += simulate();
        average_ = accum / trials;
    }

    // @return the average distance of the shortest path
    double get_average() const {
        return average_;
    }

protected:
    // Creates a new random graph and @return the shortest path of that graph
    double simulate() {
        Graph graph = create_graph();
        ShortestPath shortestPath(graph);
        int counter = 0;
        double accum = 0;
        shortestPath.path_size(0, vertex_-1);  //Calculates the distances from vertex 0 to all other connected vertices
        vector<double> distances = shortestPath.get_min_distances_();
        for (int j = 1; j < vertex_; j++) {
            if (distances[j] < numeric_limits<double>::max()) {             //If there is a possible shortest path
                counter++;
                accum += distances[j];
            }
        }

        return accum / counter;
    }

    // Constructs the randomized graph and @returns it
    Graph create_graph() {
        Graph graph(vertex_);

        for (int i = 0; i < vertex_ - 1; i++) {
            for (int j = i + 1; j < vertex_; j++) {
                if (edge_distribution_(e_) < density_ && i != j)    // If the random gives a value less than the density
                    // and to avoid vertex loops
                    graph.add_edge(i, j, weight_distribution_(e_));
            }

        }
        return graph;
    }

private:
    unsigned int vertex_;
    double density_;
    double average_;
    default_random_engine e_;
    uniform_int_distribution<int> weight_distribution_;
    uniform_real_distribution<double> edge_distribution_;
};

int main() {
    cout << "Doing this homework I have learned the use of keyword const at the end of the method sentence to\n"
                    " specify that the method can't modify the object. The structure of a class and its comments is\n"
                    " a big achievement. Also I have learned to instantiate objects like this\n"
                    " PairVertexValue(y, weight). And the use of default_random_engine instead of rand() and both\n"
                    " uniform_int_distribution (for ints) and uniform_real_distribution (for real values).\n"
                    " Moreover, I have put my first exception with C++ when the distance received by parameter\n"
                    " at PairVertexValue constructor is negative and I throw the invalid_argument exception.\n"
                    " Also I have written unsigned int instead of int to emphasize that only positives ints are\n"
                    " allowed. To differentiate between class attribute and the variable passed by parameter in\n"
                    " the class constructor, I have seen that It added a underscore at the final of the class\n"
                    " attribute name. What I have achieved is a big control about the std::vector\n"
                    " (push_back(), erase()...) as well as the use of iterators. Additionally, I have modified\n"
                    " some operators for the PairVertexValue class and I have added a clone() method to avoid\n"
                    " modify the original graph.\n"
                    " And finally, I want to explain that this homework has been made with an IDE for C++ named\n"
                    " Clion (for MacOSX), It has helped much in my work. However, I hear new recommendations\n"
                    " of others IDEs for MacOSX. " << endl;

    cout << "\n********Program********" << endl;

    cout << "\nEnter number of trials: ";
    int trials;
    cin >> trials;

    MonteCarloSimulation simulation1(trials, 0.2);     //Density of 0.2
    MonteCarloSimulation simulation2(trials, 0.4);     //Density of 0.4

    cout << "The average shortest path (Density 20%) is: " << simulation1.get_average() << endl;
    cout << "The average shortest path (Density 40%) is: " << simulation2.get_average() << endl;

    return 0;
}