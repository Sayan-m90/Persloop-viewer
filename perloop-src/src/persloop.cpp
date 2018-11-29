///////////////////////////////////////////////////////////////////////////////
//
// THIS SOFTWARE IS PROVIDED "AS-IS". THERE IS NO WARRANTY OF ANY KIND.
// NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE
// FOR ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY.
//
// Copyright (c) 2010 Jyamiti Research Group.
// CS&E Department of the Ohio State University, Columbus, OH.
// All rights reserved.
//
// Author: Sayan Mandal
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <boost/geometry.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <map>
#include <list>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_set>
#include <string>
#include <sys/stat.h>

#ifndef INFINITY
#define INFINITY std::numeric_limits< double >::infinity()
#endif // INFINITY

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/timer.hpp>
#include <boost/progress.hpp>

// #include <CGAL/Cartesian.h>

#include<ParseCommand.h>
#include <Legal.h>
#include "SimplicialComplex.h"
#include "trygraph.h"

using namespace boost;
using namespace std;
namespace bg = boost::geometry;
typedef std::vector<std::vector<int>> higherOrder;
typedef std::vector<std::vector<std::vector<float>>> higherPoint;

extern std::vector<int> complexSizes;
extern std::vector<int> accumulativeSizes;
extern int accumulativeSize;
extern float fThreshold;
extern vector<float> vecFiltrationScale;
extern SimplicialTree<bool> domain_complex;
//extern SimplicialTree<bool> range_complex;
//
//timer
extern std::clock_t start, timer1, timer2;
extern double dFuncTimeSum;
extern double dInsertTime;
extern double dCollapseTime;
extern int max_dimension;
extern int collapseCount;
extern int smallCount;
std::map<int,higherOrder> deathTracker(higherOrder thisLoop,  int birthTimethisLoop, std::map<int, higherOrder> birthOfLoops);

void printHigherOrder(higherOrder ho){
  // typedef std::vector<std::vector<int>> higherOrder;
  for(int i=0;i<ho.size();i++){
    for(int j=0;j<ho[i].size();j++)
      cout<<ho[i][j]<<" ";
    cout<<"<-->";
  }
  cout<<"\n";
}

void printHigherPoint(higherPoint ho){
  // typedef std::vector<std::vector<int>> higherOrder;
  for(int i=0;i<ho.size();i++){
    for(int j=0;j<ho[i].size();j++)
      cout<<"("<<ho[i][j][0]<<","<<ho[i][j][1]<<","<<ho[i][j][2]<<")";
    cout<<"<-->";
  }
  cout<<"\n";
}


higherOrder modifylastloop(higherOrder lastLoop){
  higherOrder temp;
  for(int i=0;i<lastLoop.size();i++){
    // for(int j=0;j<lastLoop[i].size();j++){
      if(lastLoop[i].size()!=2)
        {cout<<"error"; exit(0);}
      // cout<<lastLoop[i][0]<<" "<<lastLoop[i][1];
      // cout<<"-+-";
    }

  temp.push_back(lastLoop[0]);
  lastLoop.erase(lastLoop.begin());
  int count = 1;
  int fullsize = lastLoop.size();
  while(count<=fullsize){
    int a = temp[temp.size()-1][1]; //last element in last edge of vector, last vertex;
    int k;
    for(k=0;k<lastLoop.size();k++)
      if(lastLoop[k][1]==a ||lastLoop[k][0]==a)
        break;
      else if(k==lastLoop.size()-1)
      {
        cout<<"Loop mismatch: "<<a<<" ";
        for (int it = 0; it < temp.size(); ++it)
          cout<<temp[it][0]<<" "<<temp[it][1]<<"|";
        cout<<"\n";
        for (int it = 0; it < lastLoop.size(); ++it)
          cout<<lastLoop[it][0]<<" "<<lastLoop[it][1]<<"|";
        exit(0);
      }

    if(lastLoop[k][1]==a)   //need to reverse the edge
      { 
        std::vector<int> tttt;
        tttt.push_back(lastLoop[k][1]);
        tttt.push_back(lastLoop[k][0]);
        temp.push_back(tttt);
        
      }
    else if(lastLoop[k][0]==a){
      // cout<<"higherOrder size 5:"<<lastLoop.size();
      // getchar(); 
      temp.push_back(lastLoop[k]);
    }
    else
    {
      cout<<"Loop mismatch"; exit(0);
    }

      // cout<<"higherOrder size 3:"<<lastLoop.size();
  // getchar();

    lastLoop.erase(lastLoop.begin() + k);
    // cout<<"2";
    count++;
  // getchar();

  }
  return temp;
  
}

int simpersPart(std::vector<int>  &born,std::vector<int>  &dead, std::string simpers_file){

  ifstream ff(simpers_file.c_str());
    // cout<<"simpers part";
    if( ff.good()==false)
    {
        cout<<"simpers file "<<simpers_file<<" does not exist.";
        exit(0);
    }

    int dim, fborn;
    std::string fdead;
    int count=0;
    while(!ff.eof()){
    count ++;
      char sLine[256]="";
      ff.getline(sLine, 256);
      if(sLine==""||strlen(sLine)==0)
        return 0;
      // cout<<sLine<<",,,, "<<"\n";
      stringstream ss;
      ss.str(sLine);
      
      ss >> dim;  ss>>fborn; ss>> fdead;
      
      if(dim==0)
        continue;
      if(dim==2)
        return 0;
      born.push_back(fborn);
      // cout<<fborn<<"fb "<<fdead<<"|";
      if(fdead == "inf")
        dead.push_back(-1);
      else
          dead.push_back(stoi(fdead));

    }

    return 0;
}

void loopPrintingSingle(std::map<int,higherPoint> vloopH, int birthtime, std::string filtration_file){
  string  loops_folder;
  loops_folder = filtration_file.substr(0,filtration_file.size()-4)+"loops/";
  boost::filesystem::path dir(loops_folder.c_str());
  // boost::filesystem::create_directory(dir);
  // if(boost::filesystem::create_directory(dir)==false)
    // {cout<<"No directory";  exit(0);}

    std::vector<int> edge;
    // int k =  vloopH.begin()->first;
    higherOrder vInd;
    int count = 0;
    // higherPoint vloop = vloop.begin()->second;
  // {
    std::string file = loops_folder+std::to_string(birthtime)+".off";
    cout<<"OFF folder:"<<file<<"\n";
    ofstream ofloop(file.c_str());
    int vertexcount = 0, edgecount = 0, edgecountfirst = 0, fi=0;

  for(std::map<int,higherPoint>::iterator iter = vloopH.begin(); iter != vloopH.end(); ++iter)
  {
    higherPoint here = iter->second;

    // printHigherPoint(here);
    // cout<<"\n";
    // vertexcount += (here.size()*4);    //size of an individual loop
    edgecount += here.size();
    // if(fi==0)
      // edgecountfirst = here.size();   //first loop is of different color
    // fi++;
  }
    
  ofloop<<"OFF"<<std::endl<<(edgecount*4)<<" "<<edgecount<<" 0\n";    

  for(std::map<int,higherPoint>::iterator iter = vloopH.begin(); iter != vloopH.end(); ++iter)
  {
    higherPoint vloop_local = iter->second;
    for(int l2=0;l2<vloop_local.size();l2++){
    ofloop<<vloop_local[l2][0][0]<<" "<<vloop_local[l2][0][1]<<" "<<vloop_local[l2][0][2]<<"\n";
    ofloop<<vloop_local[l2][1][0]<<" "<<vloop_local[l2][1][1]<<" "<<vloop_local[l2][1][2]<<"\n";
    ofloop<<vloop_local[l2][0][0]<<" "<<vloop_local[l2][0][1]<<" "<<vloop_local[l2][0][2]<<"\n";
    ofloop<<vloop_local[l2][1][0]<<" "<<vloop_local[l2][1][1]<<" "<<vloop_local[l2][1][2]<<"\n";

    edge.push_back(count++);
    edge.push_back(count++);
    edge.push_back(count++);
    edge.push_back(count++);
    vInd.push_back(edge);
    edge.clear();
    if(iter==vloopH.begin())  //first cycle that was born actually
      fi++;

    }
  } 

  for(int ed=0;ed<vInd.size();ed++)
    {
      if(ed<fi)//to make first cycle a different color
        ofloop<<"4 "<<vInd[ed][0]<<" "<<vInd[ed][1]<<" "<<vInd[ed][2]<<" "<<vInd[ed][3]<<" 1.0 1.0 0.0"<<std::endl;
      else
        ofloop<<"4 "<<vInd[ed][0]<<" "<<vInd[ed][1]<<" "<<vInd[ed][2]<<" "<<vInd[ed][3]<<" 1.0 1.0 0.0"<<std::endl;
    }
    ofloop.close();

}

void loopPrinting( std::map<int,higherPoint> vloop_all, std::string filtration_file){

//higherpoint: vector, vector , vector float
  string loops_folder;
  loops_folder = filtration_file.substr(0,filtration_file.size()-4)+"loops/";
  boost::filesystem::path dir(loops_folder.c_str());
  // boost::filesystem::create_directory(dir);
  cout<<"OFF inside, loop size: "<<vloop_all.size();

  higherOrder vInd;
  
  for(std::map<int,higherPoint>::iterator iter = vloop_all.begin(); iter != vloop_all.end(); ++iter)
  {

    std::vector<int> edge;
    int k =  iter->first;
    higherPoint vloop = iter->second;
  // {
    std::string file = loops_folder+std::to_string(k)+".off";
    cout<<"OFF folder:"<<file<<"\n";
    ofstream ofloop(file.c_str());
    ofloop<<"OFF"<<std::endl<<vloop.size()*4<<" "<<vloop.size()<<" 0\n";
    int count = 0;
    for(int l2=0;l2<vloop.size();l2++){
    ofloop<<vloop[l2][0][0]<<" "<<vloop[l2][0][1]<<" "<<vloop[l2][0][2]<<"\n";
    ofloop<<vloop[l2][1][0]<<" "<<vloop[l2][1][1]<<" "<<vloop[l2][1][2]<<"\n";
    ofloop<<vloop[l2][0][0]<<" "<<vloop[l2][0][1]<<" "<<vloop[l2][0][2]<<"\n";
    ofloop<<vloop[l2][1][0]<<" "<<vloop[l2][1][1]<<" "<<vloop[l2][1][2]<<"\n";

    edge.push_back(count++);
    edge.push_back(count++);
    edge.push_back(count++);
    edge.push_back(count++);
    vInd.push_back(edge);
    edge.clear();
    // ofloop<<"#######################################################\n";

    }

    for(int ed=0;ed<vInd.size();ed++)
    {
      ofloop<<"4 "<<vInd[ed][0]<<" "<<vInd[ed][1]<<" "<<vInd[ed][2]<<" "<<vInd[ed][3]<<" 1.0 1.0 0.0"<<std::endl;
    }

    ofloop.close();
    vInd.clear();

  }
}//end of loopPrinting

int main(int argc, char *argv[] )
{
  int dimensions, noPoints, currentEdgeSize = 0, barcode_count = 0;


  string filtration_file, simpers_file;
  
  ParseCommand(argc, argv,  simpers_file, filtration_file);
  
  string loops_folder = filtration_file.substr(0,filtration_file.size()-4)+"loops/";
  boost::filesystem::path dir(loops_folder.c_str());
  boost::filesystem::create_directory(dir);
  ifstream ff(filtration_file.c_str());
  std::map<int, higherOrder> birthOfLoops;  //int: birth time, higherOrder: edges in the loop
  std::map<int, higherPoint> birthOfLoopCoordinates; //int: birth time, higherOrder: coord of edges in loop

  
  std::vector<bg::model::point<float, 3, bg::cs::cartesian>> vPoint;
  std::vector<int> v1EdgeFilt;  //edges in the filtration, start and end vertices
  std::vector<int> v2EdgeFilt;
  int generic_counter = 0;
  int li1,li2;  // vertices of the last edge which was inserted, which led to birth of a cycle
  std::vector<int>  vborn;
  std::vector<int>  vdead;
  vector<int> interm;
  float scalecount = 0;
  bool bDeathTimeOrder = true;
  bool bTimeStamp = false;
  float fMaxScale = 0.0;
  int noInserstions=0;
  int noTrinserstions = 0;
  int totinsCount = 0;
// std::cout<<"IP";
// getchar();
  simpersPart(vborn,vdead,simpers_file);

  if(ff.good()==false)
    {
        cout<<"Point+ Filtration file "<<filtration_file<<" does not exist.";
        exit(0);
    }

  
  ff >> dimensions ;      ff >> noPoints;
  cout<<"Dimming: "<<dimensions<<" #Pt: "<<noPoints<<" \n";

  if(dimensions!=3 && dimensions!=0 ){
    cout<<"Dimension should be 3";
    exit(0);
  }
  // int name[noPoints];
  for ( int itp=0; itp < noPoints; itp++ )
  {     
    float xc,yc,zc;
    ff >> xc; ff >> yc; ff >> zc;
    // cout<<"xc: "<<xc<<" yc: "<<yc<<" zc:"<<zc<<" itp:"<<itp<<" noPoints: "<<noPoints<<endl;
    // getchar();
    bg::model::point<float, 3, bg::cs::cartesian> pointbuff(xc,yc,zc);
    vPoint.push_back(pointbuff);
    
    scalecount+=1;
    filtration_step += 1;
    timer2 = std::clock();
    vecFiltrationScale.push_back(scalecount);
    interm.push_back(itp);

    domain_complex.ElementaryInsersion(interm);
    interm.clear();
    dFuncTimeSum += (std::clock() - timer2);

    complexSizes.push_back(domain_complex.ComplexSize());
    accumulativeSizes.push_back(domain_complex.accumulativeSimplexSize);

  }

    while (!ff.eof()) 
    {
      int indf;
      string infgetter;
      char sLine[256]="";
      


      ff.getline(sLine, 256);
      // cout<<sLine<<" sLine"<<endl;
      // getchar();
      if(strlen(sLine)==0||sLine[0]=='c')
        {cout<<"continuing";  continue;}

      stringstream ss;
      string buff;
      ss.str(sLine);

      if(sLine[0]=='#'){

        ss >> buff; //gets the #
        ss >> infgetter;
        if(infgetter!="inf")
          indf = atoi(infgetter.c_str());
        else
          {std::cout<<"Time should not be inf"; exit(0);}
        // cout<<sLine<<" indf: "<<indf<<" ";
        if(totinsCount+noPoints!=indf)
        {cout<<"Mismatch in insertion:"<<totinsCount<<" "<<totinsCount+noPoints<<"|"<<indf; exit(0); }
        // getchar();
        auto itf = std::find(vdead.begin(), vdead.end(), indf);

        // ******************** DEAD PART *********************
        if(itf != vdead.end() ){

          cout<<"Pers Loop Dead:sL:"<<sLine<<"\n"; 
          
          int index = std::distance(vdead.begin(), itf);
          int lid = vborn[index]; //loop_index_which_died
          cout<<"Loop born at: "<<vborn[index]<<", died at: "<<vdead[index]<<"\n";
          higherOrder lwd = birthOfLoops[lid]; // actual loop which died
          // cout<<"Actual loop death: ";
          // printHigherOrder(lwd);
          // cout<<"\n";
          map<int,higherOrder> ref = deathTracker(lwd, lid, birthOfLoops);;
          ref.insert(std::pair<int, higherOrder>(lid, lwd));//
          map<int,higherPoint> vCombDead;
          
          bool multiple = false;
          for(map<int,higherOrder>::iterator itref = ref.begin();itref!=ref.end();itref++){
            higherPoint loopbornnowCoord;
            
            higherOrder ho = itref->second;     //one cycle
              // cout<<"Birth: "<<itref->first<<" || ";

            // printHigherOrder(ho);
            // cout<<"\n";
            for(int ho_ind=0; ho_ind < ho.size(); ho_ind++)   //through each edge of the cycle
            {
              vector<float> buffpointloop;
              vector<vector<float>> buffpointloop2;

              if(ho[ho_ind].size()!=2)
              {
                cout<<"Ho size not 2. Error"; exit(0);
              }
              int li1 = ho[ho_ind][0]; //two vertices of participant edge in a cycle
              int li2 = ho[ho_ind][1];
              buffpointloop.push_back(vPoint[li1].get<0>()); 
              buffpointloop.push_back(vPoint[li1].get<1>()); 
              buffpointloop.push_back(vPoint[li1].get<2>());
              buffpointloop2.push_back(buffpointloop);

              buffpointloop.clear();
              buffpointloop.push_back(vPoint[li2].get<0>()); 
              buffpointloop.push_back(vPoint[li2].get<1>()); 
              buffpointloop.push_back(vPoint[li2].get<2>());
              buffpointloop2.push_back(buffpointloop);
              
              loopbornnowCoord.push_back(buffpointloop2);
              buffpointloop2.clear();
            }
            vCombDead[itref->first] = loopbornnowCoord;  
            if(itref!=ref.begin())
              multiple = true;
            loopbornnowCoord.clear();

          }

          lwd = modifylastloop(lwd);
          cout<<"Loop born at: "<<vborn[index]<<", died at: "<<vdead[index]<<"\n";
          // printHigherOrder(lwd);
          birthOfLoops.erase(lid);
          continue;
        
        }
        // ******************* BORN PART ************************
        else if (std::find(vborn.begin(), vborn.end(), indf) != vborn.end()){
            // cout<<sLine<<" Start of born part"<< v1EdgeFilt.size()<<"\n";
          typedef adjacency_list < listS, vecS, directedS,
            no_property, property < edge_weight_t, int > > graph_t;
          typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
          typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
          typedef pair<int, int> Edge;
  

            int sizehere = v1EdgeFilt.size();
               Graph g;
               generic_counter = 0;
            v1EdgeFilt.erase(v1EdgeFilt.begin()+v1EdgeFilt.size()-1);
            v2EdgeFilt.erase(v2EdgeFilt.begin()+v2EdgeFilt.size()-1);

            
            
            higherOrder loopbornnow;
            higherPoint loopbornnowCoord;
            vector<int> buffloop;
            vector<float> buffpointloop;
            vector<vector<float>> buffpointloop2;

            // *********************** Create shortest path METHOD 1  *************************//
            if(v1EdgeFilt.size()>300000){//{
            int count=0;
            vector<int>::iterator ib = v1EdgeFilt.begin();
            vector<int>::iterator id = v2EdgeFilt.begin();
            
            
            for(; ib!=v1EdgeFilt.end(); ib++,id++){

              g.add_vertex(*ib, *id, boost::geometry::distance(vPoint[*ib],vPoint[*id]));
              g.add_vertex(*id, *ib, boost::geometry::distance(vPoint[*id],vPoint[*ib]));

            }

            vector<int> vertt = g.shortest_path(li1, li2);

            for (int ii=0;ii<vertt.size()-1;ii++)
            {              
              buffloop.push_back(vertt[ii]);
              buffloop.push_back(vertt[ii+1]);
              loopbornnow.push_back(buffloop);
              
              buffpointloop.push_back(vPoint[vertt[ii]].get<0>()); 
              buffpointloop.push_back(vPoint[vertt[ii]].get<1>()); 
              buffpointloop.push_back(vPoint[vertt[ii]].get<2>());
              buffpointloop2.push_back(buffpointloop);

              buffpointloop.clear();
              buffpointloop.push_back(vPoint[vertt[ii+1]].get<0>()); 
              buffpointloop.push_back(vPoint[vertt[ii+1]].get<1>()); 
              buffpointloop.push_back(vPoint[vertt[ii+1]].get<2>());
              buffpointloop2.push_back(buffpointloop);
              
              loopbornnowCoord.push_back(buffpointloop2);
              buffloop.clear();
              buffpointloop.clear();
              buffpointloop2.clear();

            }
              buffloop.push_back(vertt[vertt.size()-1]);
              buffloop.push_back(li1);
              loopbornnow.push_back(buffloop);
              
              buffpointloop.push_back(vPoint[vertt.size()-1].get<0>()); 
              buffpointloop.push_back(vPoint[vertt.size()-1].get<1>()); 
              buffpointloop.push_back(vPoint[vertt.size()-1].get<2>());
              buffpointloop2.push_back(buffpointloop);

              buffpointloop.clear();
              buffpointloop.push_back(vPoint[li1].get<0>()); 
              buffpointloop.push_back(vPoint[li1].get<1>()); 
              buffpointloop.push_back(vPoint[li1].get<2>());
              buffpointloop2.push_back(buffpointloop);
              
              loopbornnowCoord.push_back(buffpointloop2);
              // cout<<"vert: "<<vertt[vertt.size()-1]<<","<<vertt[0]<<endl;
              // getchar();
              buffloop.clear();
              buffpointloop.clear();
              buffpointloop2.clear();  
            }
            // *********************** END Create shortest path METHOD 1  *************************//
            // *********************** Create shortest path METHOD 2 BOOST  *************************//
            else{
            
            Edge edge_array[v1EdgeFilt.size()*2];
            float weights[v1EdgeFilt.size()*2];//
            // cout<<"chance";
            // getchar();
            vector<int>::iterator ib = v1EdgeFilt.begin();
            vector<int>::iterator id = v2EdgeFilt.begin();            
            generic_counter = 0;

            for(; ib!=v1EdgeFilt.end(); ib++,id++){

              edge_array[generic_counter].first = *ib;
              edge_array[generic_counter].second = *id;
              weights[generic_counter] = boost::geometry::distance(vPoint[*ib],vPoint[*id]);
              generic_counter++;
              edge_array[generic_counter].first = *id;
              edge_array[generic_counter].second = *ib;
              weights[generic_counter] = boost::geometry::distance(vPoint[*id],vPoint[*ib]);
              generic_counter++;              
            }
            const int num_nodes = noPoints+1;
            
            int num_arcs = sizeof(edge_array) / sizeof(Edge);
            //graph created from the list of edges
            graph_t grph(edge_array, edge_array + num_arcs, weights, num_nodes);
            //create the property_map from edges to weights
            property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, grph);
            //create vectors to store the predecessors (p) and the distances from the root (d)
            vector<vertex_descriptor> pred(num_vertices(grph));
            vector<float> dist_root(num_vertices(grph));
            //create a descriptor for the source node
            vertex_descriptor srcv = vertex(li1, grph);
            //evaluate dijkstra on graph g with source s, predecessor_map p and distance_map d
            //note that predecessor_map(..).distance_map(..) is a bgl_named_params<P, T, R>, so a named parameter
            dijkstra_shortest_paths(grph, srcv, predecessor_map(&pred[0]).distance_map(&dist_root[0]));
            graph_traits < graph_t >::vertex_iterator vi, vs, vback, vj, vend;
            tie(vback, vend) = vertices(grph);
            tie(vi, vend) = vertices(grph);
            tie(vj, vend) = vertices(grph);
            tie(vs, vend) = vertices(grph);

          // if(indf==10132)
            // {cout<<"\n got here 6 "<<generic_counter<<" "<<noPoints<<"|"<<noInserstions<<"|"<<noTrinserstions<<"|"<<totinsCount; getchar();}

            vi = vi + li2;
            vj = vj + li1;
            // higherOrder loopbornnow;
            // higherPoint loopbornnowCoord;
            vector<int> buffloop;
            vector<float> buffpointloop;
            vector<vector<float>> buffpointloop2;
            int count=0;
            while(*vi != *vj){

              buffloop.push_back(*vi);
              buffloop.push_back(pred[*vi]);
              if(*vi==pred[*vi])
              {
                cout<<"Wrong graph construction"<<endl;
                exit(0);
              }
              loopbornnow.push_back(buffloop);
              
              buffpointloop.push_back(vPoint[*vi].get<0>()); 
              buffpointloop.push_back(vPoint[*vi].get<1>()); 
              buffpointloop.push_back(vPoint[*vi].get<2>());
              buffpointloop2.push_back(buffpointloop);

              buffpointloop.clear();
              buffpointloop.push_back(vPoint[pred[*vi]].get<0>()); 
              buffpointloop.push_back(vPoint[pred[*vi]].get<1>()); 
              buffpointloop.push_back(vPoint[pred[*vi]].get<2>());
              buffpointloop2.push_back(buffpointloop);
              
              loopbornnowCoord.push_back(buffpointloop2);
              buffloop.clear();
              buffpointloop.clear();
              buffpointloop2.clear();


              // cout<<*vi<<" "<<li1<<" "<<li2<<" "<<*vj<<" pv: "<<pred[*vi]<<" pv: "<<pred[*vj]<<"--";
              // getchar();
              vi = vs + pred[*vi];
            }

            }
            // *********************** End Create shortest path METHOD 2 BOOST  *************************//
            buffloop.push_back(li2);
            buffloop.push_back(li1);
            loopbornnow.push_back(buffloop);

            buffpointloop.push_back(vPoint[li2].get<0>()); 
            buffpointloop.push_back(vPoint[li2].get<1>()); 
            buffpointloop.push_back(vPoint[li2].get<2>());
            buffpointloop2.push_back(buffpointloop);

            buffpointloop.clear();
            buffpointloop.push_back(vPoint[li1].get<0>()); 
            buffpointloop.push_back(vPoint[li1].get<1>()); 
            buffpointloop.push_back(vPoint[li1].get<2>());
            buffpointloop2.push_back(buffpointloop);
            
            loopbornnowCoord.push_back(buffpointloop2);
            buffloop.clear();
            buffpointloop.clear();
            buffpointloop2.clear();
            //***********************************************************
            birthOfLoops.insert(std::pair<int, higherOrder>(indf, loopbornnow));
            birthOfLoopCoordinates.insert(std::pair<int, higherPoint>(indf, loopbornnowCoord));
            v1EdgeFilt.push_back(li1);
            v2EdgeFilt.push_back(li2);
            currentEdgeSize++;

            if(sizehere!=v1EdgeFilt.size())
              {cout<<"mismatch "<<sizehere<<" "<<v1EdgeFilt.size(); exit(0);}
            // *********************** END Create shortest path BOOST *************************//
            continue;
        }

        else continue;
      }// end of #
      
      
      ss.str(sLine);
      char ch;
      int index;
      vector<int> simplex;
      ss >> ch;
      if(ch!='i'){
        cout<<"Not an insertion: "<<sLine;
        exit(0);
      }

      totinsCount++;
      while (ss >> index)
        simplex.push_back(index);

      sort(simplex.begin(),simplex.end());
      if(simplex[0]>simplex[1])
        {cout<<"Error: simplices unsorted"; exit(0);}
      scalecount+=1;
      filtration_step += 1;
      timer2 = std::clock();

      vecFiltrationScale.push_back(scalecount);
      domain_complex.ElementaryInsersion(simplex);
      // simplex.clear();
        
      complexSizes.push_back(domain_complex.ComplexSize());
      accumulativeSizes.push_back(domain_complex.accumulativeSimplexSize);
      dFuncTimeSum += (std::clock() - timer2);
    
    
      if(simplex.size()>2)
       continue;

      // noInserstions++;
      v1EdgeFilt.push_back(simplex[0]);
      v2EdgeFilt.push_back(simplex[1]);
      li1 = simplex[0];
      li2 = simplex[1];
    }
    
  
  // cout<<"birthOfLoopCoordinates size:"<<birthOfLoopCoordinates.size();
  loopPrinting( birthOfLoopCoordinates, filtration_file);


  return EXIT_SUCCESS;
}
