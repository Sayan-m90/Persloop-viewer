//
//  main.cpp
//  dijkstra
//
//  Created by Mahmut Bulut on 11/11/13.
//  Copyright (c) 2013 Mahmut Bulut. All rights reserved.
//

#include <unordered_map>
#include <vector>
#include <limits>
#include <algorithm>
#include <iostream>

using namespace std;

class Graph
{
    
    
public:
    unordered_map<int,  unordered_map<int, float>> vertices;
    
    void add_vertex2(int name,  unordered_map<int, float>& edges)
    {
        // Insert the connected nodes in unordered map
        unordered_map<int, float> allreadygot;
        allreadygot.insert(vertices[name].begin(),vertices[name].end());
        allreadygot.insert(unordered_map<int, float>::value_type(1,edges[1]));
        vertices.insert(unordered_map<int,  unordered_map<int, float>>::value_type(name, allreadygot));
    }

    void add_vertex(int name, int dest, float val)// unordered_map<int, float>& edges)
    {
        // Insert the connected nodes in unordered map
        //  unordered_map<int, float> allreadygot = edges;
        // unordered_map<int, float> allreadygot;
        // allreadygot.insert(vertices[name].begin(),vertices[name].end());
        // allreadygot[dest] = val;//.insert(unordered_map<int, float>::value_type(dest,val));//.insert(unordered_map<int, float>::value_type(dest,val));
        // vertices.insert(unordered_map<int,  unordered_map<int, float>>::value_type(name, allreadygot));
        vertices[name][dest] = val;
    }
    
    vector<int> shortest_path(int start, int finish)
    {
        // Second arguments -> distances
        // Find the smallest distance in the already in closed list and push it in -> previous
        unordered_map<int, float> distances;
        unordered_map<int, int> previous;
        vector<int> nodes; // Open list
        vector<int> path; // Closed list
        
        auto comparator = [&] (int left, int right) { return distances[left] > distances[right]; };
        
        for (auto& vertex : vertices)
        {
            if (vertex.first == start)
            {
                distances[vertex.first] = 0;
            }
            else
            {
                distances[vertex.first] = numeric_limits<float>::max();
            }
            
            nodes.push_back(vertex.first);
            push_heap(begin(nodes), end(nodes), comparator);
        }
        
        while (!nodes.empty())
        {
            pop_heap(begin(nodes), end(nodes), comparator);
            int smallest = nodes.back();
            nodes.pop_back();
            
            std::cout << "Open list: ";
            for( std::vector<int>::iterator i = nodes.begin(); i != nodes.end(); ++i)
                std::cout << *i << ' ';
            std::cout << std::endl;
            
            if (smallest == finish)
            {
                while (previous.find(smallest) != end(previous))
                {
                    path.push_back(smallest);
                    smallest = previous[smallest];
                    std::cout << "Closed list: ";
                    for( std::vector<int>::iterator i = path.begin(); i != path.end(); ++i)
                        std::cout << *i << ' ';
                    std::cout << std::endl;
                }
                
                break;
            }
            
            if (distances[smallest] == numeric_limits<float>::max())
            {
                break;
            }
            
            for (auto& neighbor : vertices[smallest])
            {
                float alt = distances[smallest] + neighbor.second;
                if (alt < distances[neighbor.first])
                {
                    distances[neighbor.first] = alt;
                    previous[neighbor.first] = smallest;
                    make_heap(begin(nodes), end(nodes), comparator);
                }
            }
        }
        
        return path;
    }
};


int notmain()
{
    // cout << "@author: Mahmut Bulut" << endl << endl;
    int seq = 0;
    int init_node = 1;
    int dest_node = 7;
    
    Graph g;

    g.add_vertex(1, 3, 4);
    g.add_vertex(1, 2, 1.0);
    g.add_vertex(3, 7, 2);
    g.add_vertex(3, 4, 4);
    g.add_vertex(2, 5, 2);
    
    // // // g.add_vertex(4, {});
    g.add_vertex(5, 4, 3);
    g.add_vertex(6, 3, 1);
    g.add_vertex(6, 7, 4);
    // // g.add_vertex(1, {{2, 0.8}, {3, 4}});
    g.add_vertex(7, 5, 5);
    cout << "As initial node: " << init_node << endl;
    cout << "As goal node: " << dest_node << endl;
    
    for (int vertex : g.shortest_path(init_node, dest_node))
    {
        cout << "Solution path from goal sequence : " << seq << " Node : " << vertex << endl;
        seq++;
    }
    
    cout << "Solution path from goal sequence : " << seq << " Node : " << init_node << endl;
    
return 0;
}
