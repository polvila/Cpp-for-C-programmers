//
//  main.cpp
//  Homework3
//
//  Created by Pol Vilà Saló on 22/3/16.
//  Copyright © 2016 Pol Vilà Saló. All rights reserved.
//
//  This programs implements a Prim's algorith to get the minimum spanning tree (MST) of a graph, that can be get by file.
//  It can be used with the data file example of the submission webpage.

#include <iostream>
#include <vector>
#include <fstream>
#include <cfloat>

using namespace std;

const string PATH = "data1.txt";

class PairVertexValue {
public:
    // This class stores the vertex, value pairs
    //
    // @param vertex id of the vertex
    // @param value weight of the edge
    PairVertexValue(unsigned int vertex, double value) : vertex_(vertex), value_(value) {
        if (value < 0)
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

    bool operator==(const PairVertexValue &element) const {
        return vertex_ == element.getVertex_() && value_ == element.getValue_();
    }

    PairVertexValue clone() const {
        return PairVertexValue(getVertex_(), getValue_());
    }

private:
    unsigned int vertex_;
    double value_;
};


class Graph {
public:
    // Uses a vector (each position is a vertex) of vectors (each position is a vertex, weight pair) as adjacency
    // structure. This class implements an undirected graph.
    //
    // @param num_vertices number of vertices
    Graph(unsigned int num_vertices) {
        construct(num_vertices);
    }

    // Other constructor that permits create a graph from a file
    //
    // @param file_path of the input graph
    //
    // ----------- FILE FORMAT -------------
    //  mum_vertices
    //  vertex1 vertex2 weight1
    //  vertex3 vertex4 weight2
    //  ...
    Graph(string file_path) {

        ifstream data_file(file_path);
        istream_iterator<int> start(data_file), end;
        vector<unsigned int> data(start, end);        //vector contains each int in each position

        construct(*data.begin());   //Initializes a graph with the num_vertices contained in the first line of the file

        for (auto elem = data.begin() + 1; elem != data.end(); elem++) {
            unsigned int i = *elem;
            elem++;
            unsigned int j = *elem;
            elem++;
            double w = static_cast<double >(*elem);

            add_edge(i, j, w);
        }
    }

    //Method that initializes the arrays of adjacency
    void construct(unsigned int num_vertices) {
        v_ = num_vertices;
        for (int i = 0; i < num_vertices; i++)
            adj_.push_back(vector<PairVertexValue>());       //Initialize the graph_ with empty edge vectors
    }

    // @return the number of vertices in the graph
    unsigned int V() const { return v_; }

    // @return the number of edges in the graph
    unsigned int E() const { return e_; }

    // Tests whether there is an edge from node x to node y
    bool adjacent(unsigned int x, unsigned int y) const {
        for (auto i = 0; i != adj_[x].size(); i++)
            if (adj_[x][i].getVertex_() == y)
                return true;
        return false;
    }

    // Lists all vertex y (and its weight) such that there is an edge from x to y
    vector<PairVertexValue> neighbors(unsigned int x) const {
        vector<PairVertexValue> result;
        for (auto i = 0; i != adj_[x].size(); i++)
            result.push_back(adj_[x][i].clone());
        return result;
    }

    // Adds to G the edge from x to y, if it is not there
    void add_edge(unsigned int x, unsigned int y, double weight) {
        if (!adjacent(x, y)) {
            adj_[x].push_back(PairVertexValue(y, weight));
            adj_[y].push_back(PairVertexValue(x, weight));
            e_++;
        }
    }

    // Removes the edge from x to y, if it is there
    void delete_edge(unsigned int x, unsigned int y) {
        for (auto it = adj_[x].begin(); it != adj_[x].end(); it++)
            if (it->getVertex_() == y)
                adj_[x].erase(it);
    }

    // @return the value associated to the edge (x,y)
    double get_edge_value(unsigned int x, unsigned int y) const {
        for (auto i = 0; i != adj_[x].size(); i++)
            if (adj_[x][i].getVertex_() == y)
                return adj_[x][i].getValue_();
        return numeric_limits<double>::max();
    }

    // Sets the value associated to the edge (x,y) to v
    void set_edge_value(unsigned int x, unsigned int y, double v) {
        for (auto i = 0; i != adj_[x].size(); i++)
            if (adj_[x][i].getVertex_() == y)
                adj_[x][i].setValue_(v);
    }


    // Implementation of Prim Minimum Spanning Tree Algorithm.
    // @return the minimum spanning tree of the graph
    Graph getMST() {

        vector<double> dist(v_, DBL_MAX);   //array of distances from the source to each vertex, init: double max
        vector<unsigned int> edges(v_);    //array indicating for a given vertex, which vertex in the tree it is closest to
                                            // (parent)
        vector<bool> mstSet(v_, false);     //list of finished/unfinished edges

        dist[0] = 0;    //For including the first vertex in MST
        // The MST will have v_ vertices
        for (unsigned int i = 0; i < v_ - 1; i++) {

            double min_dist = DBL_MAX;
            int min_index;

            // To find the vertex with the shortest edge, from unfinished edges set
            for (int v = 0; v < v_; v++) {
                if (mstSet[v] == false && dist[v] < min_dist) {
                    min_dist = dist[v];
                    min_index = v;
                }
            }

            mstSet[min_index] = true;   //Add the vertex to the MST

            // Update the distances and the vertex parents with the min distances if the values are less than before
            // and the new vertex are from unfinished edges set
            for (auto neighbor : neighbors(i)) {
                if (neighbor.getValue_() < dist[neighbor.getVertex_()] && mstSet[neighbor.getVertex_()] == false) {
                    dist[neighbor.getVertex_()] = neighbor.getValue_();
                    edges[neighbor.getVertex_()] = i;

                }
            }
        }

        Graph g = Graph(v_);
        // Transform the MST information to Graph to print it
        for (unsigned int i = 1; i < v_; i++) {
            g.add_edge(i, edges[i], dist[i]);
        }

        return g;
    }

    // @return the total of add all the edges weights
    double getTotalCost() {
        double result = 0.0;
        for (auto i = 0; i != v_; i++) {
            for (auto j = 0; j != adj_[i].size(); j++) {
                result += adj_[i][j].getValue_();
            }
        }
        return result / 2;
    }

    // To output a graph. Friend method to use the private attributes
    friend ostream &operator<<(ostream &output, const Graph &G) {
        for (int i = 0; i < G.v_; i++) {
            for (auto j = 0; j != G.adj_[i].size(); j++) {
                output << "(" << i << ", " << G.adj_[i][j].getVertex_() << ") -> "
                << G.adj_[i][j].getValue_() << endl;
            }
        }
        return output;
    }

private:
    unsigned int v_;
    unsigned int e_;
    vector<vector<PairVertexValue> > adj_;
};

int main() {

    Graph g = Graph(PATH);

    Graph MST = g.getMST();

    cout << MST;

    cout << "Total cost of MST: " << MST.getTotalCost();

    return 0;
}