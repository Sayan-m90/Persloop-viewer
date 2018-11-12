/*
(c) 2015 Fengtao Fan, Dayu Shi
*/
#include "SimplicialComplex.h"
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

std::vector<std::unordered_map<int, pair<int, int>>> persistences;
int filtration_step;
int time_in_each_filtration_step;
SimplicialTree<bool> domain_complex;
SimplicialTree<bool> range_complex;
vector<unordered_set<int> > homo_info;
std::vector<int> complexSizes;
std::vector<int> accumulativeSizes;
float fThreshold;
vector<float> vecFiltrationScale;
int collapseCount = 0;
int smallCount = 0;
// TreeLoopTracker_ptr tlt_ptr;	//original loops of short loop
// TreeLoopTracker_ptr tlt_ptr_filt;	// loops of each filtration: stores loop at i-th step
int max_dimension = 3;

//timer
std::clock_t start, timer1, timer2;
double dFuncTimeSum;
double dInsertTime;
double dCollapseTime;

// outer vector: loop, inner vector: edges forming the loop
typedef std::vector<std::vector<int>> higherOrder;
// int findRoot(int a);

void WritePersistence(const char* pFileName, vector<unordered_set<int> > &homo_info) {
	std::ofstream ofile;
	ofile.open(pFileName, std::ifstream::out);
	//
	std::stringstream sstr(std::stringstream::in | std::stringstream::out);
	if (ofile.is_open())
	{
		sstr << "Ranks of the persistent image of input simplicial map in all dimensions" << endl;
		for (int i = 0; i < homo_info.size(); ++i) {
			if (!homo_info[i].empty()) {
				sstr << "Dim " << i << ": " << homo_info[i].size() << " <";
				for (unordered_set<int>::iterator sIter = homo_info[i].begin(); sIter != homo_info[i].end(); ++sIter) {
					sstr << (sIter != homo_info[i].begin() ? ", " : "") << *sIter;
				}
				sstr << ">" << endl;
			}
			else {
				sstr << "Dim " << i << ": " << homo_info[i].size() << endl;
			}
		}
		//ofile << sstr.rdbuf();
		ofile.write(sstr.str().c_str(), sstr.str().size());
		//
		ofile.close();
		sstr.clear();
	}
	else
	{
		std::cout << "Can NOT open file " << pFileName << std::endl;
		exit(0);
	}
}


void ComputingPersistenceForSimplicialMap(const char* file_name_of_domain_complex, bool is_domain_complex_with_annotation,
											const char* file_name_of_range_complex,
											const char* file_name_of_simplicial_map,
											bool is_save_range_complex_with_annotation = false,
											const char* new_range_complex_file_name = NULL) 
{
	time_in_each_filtration_step = 0;

	if (filtration_step == 0)
	{
		if (is_domain_complex_with_annotation) {
			domain_complex.ReadComplexWithAnnotation(file_name_of_domain_complex);
			domain_complex.SnapshotHomologicalFeatures(homo_info);
			///////////////first step in filtration, check born time of each homology class
			for (int i = 0; i < homo_info.size(); ++i)
			{
				std::unordered_map<int, pair<int, int>> vecEmp;
				persistences.push_back(vecEmp);
				for (::unordered_set<int>::iterator it = homo_info[i].begin(); it != homo_info[i].end(); ++it)
				{
					persistences[i][*it] = std::make_pair(0, -1);
				}
			}
		}
		else {
			domain_complex.ReadComplex(file_name_of_domain_complex);
		}
		filtration_step += 1;
		complexSizes.push_back(domain_complex.ComplexSize());
	}
	else
	{
		domain_complex.clearMemory();
		domain_complex = range_complex;
	}

	//////////perform simplicial collapse
	vector<pair<int, int> > vertex_map;
	domain_complex.ReadSimplicialMap(file_name_of_simplicial_map, vertex_map);
	unordered_map<int, int> updated_vertex_map;
	domain_complex.PerformSimplicialCollapse(vertex_map, updated_vertex_map);
	/////////perform simplicial insertions
	// add the subcomplex of the range_complex which is the preimage of the simplicial map first
	range_complex.InitializeByRenamingIncomingComplex(domain_complex, updated_vertex_map);
	// add the rest of the range_complex
	range_complex.AddRemainingSimpliciesFromFile(file_name_of_range_complex);
	if (is_save_range_complex_with_annotation) 
	{
		range_complex.WriteComplexWithAnnotation(new_range_complex_file_name);
	}
	complexSizes.push_back(range_complex.ComplexSize());
	accumulativeSizes.push_back(domain_complex.accumulativeSimplexSize);
	return;
}

void ComputingPersistenceForSimplicialMapElementary(vector<vector<int> > simplexIns, const char* file_name_of_domain_complex, bool is_domain_complex_with_annotation,
													vector<string>& vecElemOpers,
													bool is_save_range_complex_with_annotation = false,
													const char* new_range_complex_file_name = NULL) 
{
	
	// cout<<domain_complex.ComplexSize();
	if (filtration_step == 0)
	{
		if (file_name_of_domain_complex != '\0')
		{
			if (is_domain_complex_with_annotation) {
				domain_complex.ReadComplexWithAnnotation(file_name_of_domain_complex);
				domain_complex.SnapshotHomologicalFeatures(homo_info);
				///////////////first step in filtration, check born time of each homology class
				for (int i = 0; i < homo_info.size(); ++i)
				{
					std::unordered_map<int, pair<int, int>> vecEmp;
					persistences.push_back(vecEmp);
					for (::unordered_set<int>::iterator it = homo_info[i].begin(); it != homo_info[i].end(); ++it)
					{
						persistences[i][*it] = std::make_pair(0, -1);
					}
				}
			}
			else {
				domain_complex.ReadComplex(file_name_of_domain_complex);
			}
		}
		complexSizes.push_back(domain_complex.ComplexSize());
		accumulativeSizes.push_back(domain_complex.accumulativeSimplexSize);
		filtration_step += 1;
	}

	bool flag = false;
	for(int count=0;count<simplexIns.size();count++){

		std::vector<int> interm = simplexIns.at(count);
		timer1 = std::clock();

		domain_complex.ElementaryInsersion(interm);
		

		interm.clear();
		dInsertTime += std::clock() - timer1;
	}

	if (is_save_range_complex_with_annotation) {
		domain_complex.WriteComplexWithAnnotation(new_range_complex_file_name);
	}
	complexSizes.push_back(domain_complex.ComplexSize());
	accumulativeSizes.push_back(domain_complex.accumulativeSimplexSize);
	filtration_step += 1;
	return;
}


ListNodePtr AddBoundary(std::vector<std::vector<int>> simplex_vertices, std::vector<int> vPointMap, SimplicialTree<bool> domain_complex){
// For each cycle
// outer vector: stores edges
	vector<SimplicialTreeNode_ptr> cyclic_edges;	
		cout<<"got here as well 1 ("<<simplex_vertices[0][0]<<","<<simplex_vertices[0][1]<<")(";
		cout<<vPointMap[simplex_vertices[0][0]]<<","<<vPointMap[simplex_vertices[0][1]]<<")\n";	
	int sizes = 1;

	for(int i=0;i<simplex_vertices.size();i++){
		for(int j=0;j<simplex_vertices[i].size();j++){
			simplex_vertices[i][0] = domain_complex.findRoot(simplex_vertices[i][0]);
			simplex_vertices[i][1] = domain_complex.findRoot(simplex_vertices[i][1]);
			// int idx = simplex_vertices[i][0];
			// while(vPointMap[idx]!=idx)
			// 	idx = vPointMap[idx];
			// int ibx = simplex_vertices[i][1];
			// while(vPointMap[ibx]!=ibx)
			// 	ibx = vPointMap[ibx];
			// if(idx==ibx)
			// 	simplex_vertices.erase(simplex_vertices.begin()+i);

					// cout<<simplex_vertices[i][j]<<" ";
			// while(vPointMap[simplex_vertices[i][0]] != simplex_vertices[i][0])
				// simplex_vertices[i][0] = vPointMap[simplex_vertices[i][0]];
			// while(vPointMap[simplex_vertices[i][1]] != simplex_vertices[i][1])
				// simplex_vertices[i][1] = vPointMap[simplex_vertices[i][1]];
			
			// if(simplex_vertices[i][0] == simplex_vertices[i][1])
				// {
					// simplex_vertices.erase(simplex_vertices.begin() + i);
					// continue;
				// }
			if(simplex_vertices[i][0] > simplex_vertices[i][1]){
				int temp = simplex_vertices[i][0];
				simplex_vertices[i][0] = simplex_vertices[i][1];
				simplex_vertices[i][1] = temp;
			}
		}
		// simplex_vertices
		// cout<<"|";
	}
	// if(simplex_vertices[0][0]==8830 && simplex_vertices[0][1]==8838 )
		// simplex_vertices[0][1] = 8832;
	// cout<<"got here as well";
	// getchar();
	// cout<<"Loop Sum:\n";
	int i=0;
	// bool emptyFlag = false;
	SimplicialTreeNode_ptr simplex ;

	do{

	// emptyFlag = false;
	simplex = domain_complex.find(simplex_vertices[i]);	//get one edge

		cout<<"i val ("<<simplex_vertices[i][0]<<","<<simplex_vertices[i][1]<<") "<<i<<" "<<simplex_vertices.size()<<"\n";
		i++;
	}while(simplex==NULL && i<simplex_vertices.size());
	
		// cout<<vPointMap[simplex_vertices[0][0]]<<","<<vPointMap[simplex_vertices[0][1]]<<")";
// getchar();
	// if(simplex!=NULL)//->index_in_filtration)
		cout<<"inside print";
	// getchar();
	if(i==simplex_vertices.size())
		{cout<<"inside print2";
		// getchar();
		// ListNodePtr p = domain_complex.find_annotation(simplex);
		// cout<<"inside print3";
		// getchar();
			ListNodePtr sum;
			// sum->next=sum;
			return sum;
		}
			// cout<<"print"<<simplex->index_in_filtration;
		// getchar();
	ListNodePtr p = domain_complex.find_annotation(simplex);

ListNodePtr sum = domain_complex.annotations[sizes]->DeepCopyAnnotationColumn(domain_complex.find_annotation(simplex));
	int dead_bit = -1;
// }
cout<<"out of it all\n";
// getchar();
// 	// getchar();
	
// 		cout<<"got here as well 3";
// 	getchar();
	
// 			cout<<"got here as well 4";
// 	getchar();
	// simplex.clear();
	for(;i<simplex_vertices.size();i++){
		// emptyFlag= false;
	// SimplicialTreeNode_ptr simplex ;
	cout<<"got here iter: "<<simplex_vertices[i][0]<<","<<simplex_vertices[i][1]<<"|"<<i<<endl;
	// cout<<domain_complex.findRoot(simplex_vertices[i][0])<<","<<domain_complex.findRoot(simplex_vertices[i][1])<<endl;
	// getchar();
	simplex = domain_complex.find(simplex_vertices[i]);	//get one edge
	if(simplex==NULL)
		continue;
	// cout<<"got here iter 2.";
	// getchar();
	ListNodePtr tempLNP = domain_complex.annotations[sizes]->DeepCopyAnnotationColumn(domain_complex.find_annotation(simplex));
	cout<<"got here iter 3.";
	// getchar();
	dead_bit = domain_complex.annotations[sizes]->sum_two_annotation_with_changed_dst(sum, tempLNP);
	}
	

	// ListNodePtr trav(sum->next);
	// 	cout << "[";
	// 	while (trav != sum) {
	// 		cout << "<" << trav->row << " " << 1 << ">" << (trav->next == sum ? "" : " ");
	// 		trav = trav->next;
	// 	}
	// 	cout << "]";

		// cout<<"Dead bit:"<<dead_bit<<" sum row: "<<sum->next->row<<"\n";
		// if(sum->next==sum)
			// cout<<"reduced to zero. cycle is dead\n";
		// else cout<<"cycle persist\n";
		return sum;
}//end of addboundary

int LowReturner(ListNodePtr sum){

	if(sum->next==sum)
		return sum->row;
	ListNodePtr trav(sum->next);

	while (trav != sum) {
		if(trav->next==sum)
			return trav->row;
			trav = trav->next;
		}
		return -1;
}


//input: vector of vector: cycle indices: 
//output: the low of the annotation of the cycle
// int lowfromHigherOrder(higherOrder &simplex_vertices){
	
// 	SimplicialTreeNode_ptr simplex ;
// 	int sizes = 1;
// 	simplex = domain_complex.find(simplex_vertices.at(0));	//get one edge
// 	ListNodePtr p = domain_complex.find_annotation(simplex);
// 	ListNodePtr sum = domain_complex.annotations[sizes]->DeepCopyAnnotationColumn(domain_complex.find_annotation(simplex));
// 	int dead_bit = -1;
// 	// simplex.clear();
// 	for(int i=1;i<simplex_vertices.size();i++){
// 	// SimplicialTreeNode_ptr simplex ;
// 	simplex = domain_complex.find(simplex_vertices.at(i));	//get one edge
// 	ListNodePtr tempLNP = domain_complex.annotations[sizes]->DeepCopyAnnotationColumn(domain_complex.find_annotation(simplex));
// 	dead_bit = domain_complex.annotations[sizes]->sum_two_annotation_with_changed_dst(sum, tempLNP);
// 	}
// 	int num=0;
// 	ListNodePtr trav(sum->next);
// 	while (trav != sum) {
// 		num += pow(2,(trav->row));
// 		// cout<<"row: "<<trav->row<<" ";
// 	if(trav->next==sum)
// 		{
// 			// cout<<"\n num: "<<num; 
// 			return num;
// 		}
// 		trav = trav->next;
// 	}
// 	return -1;
	
// }

// void CheckBoundaryBirthOfLoops(std::map<int, higherOrder> birthOfLoops){
// 	std::map<int, higherOrder>::iterator itbl;
// 	for(itbl = birthOfLoops.begin(); itbl != birthOfLoops.end(); itbl++){
// 		ListNodePtr inter = AddBoundary(itbl->second);
// 		if(inter->next==inter){
// 			std::cout<<"Cycle in original complex shouldn't have died here."<<itbl->first<<"\n";
// 			higherOrder buff = itbl->second;
// 			for(int ind=0;ind<buff.size();ind++){
// 				cout<<buff[ind][0]<<" "<<buff[ind][1]<<"<->";
// 			}
// 			exit(0);
// 		}
// 	}

// }

//takes the loop, if its annotation is zero, it has dies. if not, then a combination has dies. Find which combination is dead.
std::map<int,higherOrder> deathTracker(higherOrder thisLoop,  int birthTimethisLoop, std::vector<int> vPointMap, std::map<int, higherOrder> birthOfLoops, SimplicialTree<bool> domain_complex){

	std::vector<ListNodePtr> listofLoops;
	std::map<int, higherOrder>::iterator itbl;
	// std::map<int, higherOrder>::iterator ithr;
	std::map<int, higherOrder> vloop;
	int sizes = 1;
	int dead_bit = 0;
	std::map<int, higherOrder> reducedLoops; //will contain vertices mapped under collapse
	// reducedLoops.insert(birthOfLoops.begin(), birthOfLoops.end());;
	// cout<<"got here in search of death";
	// getchar();

	ListNodePtr newloop = AddBoundary(thisLoop, vPointMap, domain_complex);
	cout<<"got here in search of death as well";
	// getchar();
		// cout<<"got here";
	// getchar();
	if(newloop==NULL || newloop->next==newloop){
			std::cout<<"Cycle dead. No combination required"<<"\n";
			// getchar();
			vloop.insert(std::pair<int, higherOrder>(birthTimethisLoop, thisLoop));
			return vloop;
		}
		//*********** Check to see if loop 
	// for(itbl = birthOfLoops.begin(); itbl != birthOfLoops.end(); itbl++){
	// 	higherOrder buffLoop = itbl->second;
	// 	higherOrder innerLoop;
	// 	std::vector<int> buffedge;
	// 	for(int ind=0;ind< buffLoop.size();ind++){
	// 		int indloop0 = buffLoop[ind][0];
	// 		int indloop1 = buffLoop[ind][1];
	// 		if (vPointMap[indloop0] == indloop0)//This vertex has not been collapsed
	// 			buffedge.push_back(indloop0);
	// 		else{								//Push whatever this vertex is collapsed to
	// 			while(vPointMap[indloop0]!=indloop0)
	// 				indloop0 = vPointMap[indloop0];
	// 			buffedge.push_back(indloop0);

	// 		}

	// 		// if (vPointMap[indloop1] == indloop1)//This vertex has not been collapsed
	// 			// buffedge.push_back(indloop1);
	// 		// else{								//Push whatever this vertex is collapsed to
	// 			// while(vPointMap[indloop1]!=indloop1)
	// 				// indloop1 = vPointMap[indloop1];
	// 			// buffedge.push_back(indloop1);
	// 		// }
	// 		innerLoop.push_back(buffedge);
	// 	}
	// 	reducedLoops[itbl->first] = innerLoop;	//new reduced map
	// }	

	vloop.insert(std::pair<int, higherOrder>(birthTimethisLoop, thisLoop));
	for(itbl = birthOfLoops.begin(); itbl != birthOfLoops.end(); itbl++){
		if(itbl->first > birthTimethisLoop)
			continue;
		ListNodePtr inter = AddBoundary(itbl->second, vPointMap, domain_complex);
		if(inter==NULL)
		{
			std::cout<<"Cycle dead. This should be alive "<<itbl->first<<"\n";
			getchar();
			higherOrder buff = itbl->second;
			for(int ind=0;ind<buff.size();ind++){
				cout<<buff[ind][0]<<" "<<buff[ind][1]<<"|";
				int one = domain_complex.findRoot(buff[ind][0]);
				int two = domain_complex.findRoot(buff[ind][1]);
				cout<<one<<","<<two<<"|(";
				int bg=buff[ind][0];
				while(bg!=vPointMap[bg])
					bg = vPointMap[bg];
				cout<<bg<<",";

				bg=buff[ind][1];
				while(bg!=vPointMap[bg])
					bg = vPointMap[bg];
				cout<<bg<<")<-->";
			}
			exit(0);
		}

		if(inter->next==inter){
			// simp_weight.erase(itsw); 
			std::cout<<"Cycle Single. "<<itbl->first<<","<<inter->row<<"\n";
			higherOrder buff = itbl->second;
			for(int ind=0;ind<buff.size();ind++){
				cout<<buff[ind][0]<<" "<<buff[ind][1]<<"|";
				int one = domain_complex.findRoot(buff[ind][0]);
				int two = domain_complex.findRoot(buff[ind][1]);
				cout<<one<<","<<two<<"|(";
				int bg=buff[ind][0];
				while(bg!=vPointMap[bg])
					bg = vPointMap[bg];
				cout<<bg<<",";

				bg=buff[ind][1];
				while(bg!=vPointMap[bg])
					bg = vPointMap[bg];
				cout<<bg<<")<-->";
			}
			// exit(0);
			getchar();
			// continue;
		}
		
		if(LowReturner(newloop)==LowReturner(inter)){
			vloop.insert(std::pair<int, higherOrder>(itbl->first, itbl->second));
			dead_bit = domain_complex.annotations[sizes ]->sum_two_annotation_with_changed_dst(newloop, inter);
		}
		if(dead_bit == -1){	// this is the combination which got killed				
			cout<<"comb which got killed";
			return vloop;
		}

	}
	cout<<"No combination of loops are dead. Error!";
	// exit(0);
	return vloop;

}
// Takes in the current loop and set containing all loops and sees if this current one is independant, if so then this was born
// bool bornTracker(higherOrder thisLoop, std::map<int, higherOrder> birthOfLoops){
// 	// Convert loops in complex to simpers data structure
// 	std::vector<ListNodePtr> listofLoops;
// 	std::map<int, higherOrder>::iterator itbl;
// 	int countitr = 0;
// 	int sizes = 1;

// 	for(itbl = birthOfLoops.begin(); itbl != birthOfLoops.end(); itbl++){
// 		ListNodePtr inter = AddBoundary(itbl->second);
// 		if(inter->next==inter){
// 			// simp_weight.erase(itsw); 
// 			std::cout<<"Cycle in original complex shouldn't have died here."<<itbl->first<<"\n";
// 			higherOrder buff = itbl->second;
// 			for(int ind=0;ind<buff.size();ind++){
// 				cout<<buff[ind][0]<<" "<<buff[ind][1]<<"<->";
// 			}
// 			exit(0);			
// 		}

// 		listofLoops.push_back(inter);
// 		countitr++;
// 	}
// 	ListNodePtr newloop = AddBoundary(thisLoop);
// 	// cout<<"This is the problem";
// 						// getchar();

// 	// Go through each loop check independance with this one to see if it is linearly independant
// 	for(int itl=0;itl<listofLoops.size();itl++){
// 		int proxydeadbit = domain_complex.annotations[sizes]->sum_two_annotation_with_changed_dst(listofLoops[itl], newloop);
// 		if(proxydeadbit==-1)
// 			return false;
// 			// {						cout<<"no it aint";
// 						// getchar();return false;}
// 	}
// 							// cout<<"no it aint";
// 						// getchar();
// 	return true;

// }


// int independantCycleCalculate(std::multimap<float, higherOrder> &simp_weight){

// 	std::vector<ListNodePtr> listofLoops;
// 	std::vector<ListNodePtr> aliveLoops;
// 	std::vector<int> aliveLoopsIndex;
// 	std::vector<int> deadLoopsIndex;
// 	bool nullflag = true;
// 	int countitr = 0;
// 	std::multimap<float, higherOrder>::iterator itsw;

// 	for(itsw = simp_weight.begin(); itsw != simp_weight.end(); itsw++){
// 		ListNodePtr inter = AddBoundary(itsw->second);
// 		if(inter->next==inter){
// 			std::cout<<"Cycle dead poor soul.\n";
// 			return countitr;
// 		}

// 		listofLoops.push_back(inter);
// 		countitr++;
// 	}


// 	int dead_bit;
// 	int sizes = 2;
// 	itsw = simp_weight.begin();
// 	bool dead = false;
// 	int itl;
// 	for(itl=0;itl<listofLoops.size();itl++){	// go through all loops in the basis
// 		// std::cout<<"Loop "<<(itl)<<" :";
// 		dead = false;
// 		int alivesize = aliveLoops.size();
// 		for(int ita=alivesize-1; ita>=0;--ita){ //all previous, loops already alive, in matrix, check current
// 			if(LowReturner(aliveLoops[ita])==LowReturner(listofLoops[itl]))// if low is same, add with previous
// 				dead_bit = domain_complex.annotations[sizes ]->sum_two_annotation_with_changed_dst(listofLoops[itl], aliveLoops[ita]);
// 			if(dead_bit == -1){	// not part of the basis class, know which one has died
// 				deadLoopsIndex.push_back(itl);
// 				std::cout<<"died poor soul. \n************CHECK this part*************\n";	
// 				dead = true;
// 				// simp_weight.erase(itsw); 
// 				return (itl);
// 			}
// 		}
// 		//std::cout<<"is alive\n";
// 		if(dead == false)
// 			{	aliveLoopsIndex.push_back(itl);}
// 		else{
// 				std::cout<<"Died poor soul. Error: should not have come here. \n";	
// 				exit(0);
// 			}
// 			itsw++;
// 	}
// 	int returner = 0;
// 	// cout<<"No cycle present";
// 	if(dead == false)
// 		return -1;
// 	else
// 		return itl;
// }


// int independantCycleCalculate2(std::map<float,int> weights, std::vector<std::vector<std::vector<int>>> higher_simplex){
// 	// std::sort(weights.begin(),weights.end());
// 	std::map<float, int> :: iterator itr;
// 	// cout<<"cycle count:"<<higher_simplex.size();//<<"weights";
// 	// for(itr = weights.begin(); itr != weights.end(); ++itr)
// 	// 	cout<<"\nFirst: "<<itr->first<<" ,"<<itr->second<<"\n";


// 	std::vector<ListNodePtr> listofLoops;
// 	std::vector<ListNodePtr> aliveLoops;
// 	std::vector<int> aliveLoopsIndex;
// 	std::vector<int> deadLoopsIndex;
// 	bool nullflag = true;
// 	int countitr = 0;
// 	// cout<<"INDEP: higher_simp size:"<<higher_simplex.size()<<"weights:"<<weights.size()<<"\n";
// 	for(itr = weights.begin(); itr != weights.end(); ++itr){
// 		ListNodePtr inter = AddBoundary(higher_simplex[itr->second],vPointMap);
// 		if(inter->next!=inter)
// 			{std::cout<<"Cycle dead poor soul.\n"; return itr->second;}
// 		listofLoops.push_back(inter);
// 		countitr++;
// 	}
// 	// cout<<"loop size:"<<listofLoops.size()<<" ";
// 	// if(nullflag == false ){
// 	// 	std::cout<<"Cycle dead poor soul.\n";
// 	// 	return itr;
// 	// }
// 	std::vector<ListNodePtr>::iterator itl;
// 	int dead_bit;
// 	int sizes = 2;
// 	for(int itl=0;itl<listofLoops.size();itl++){	// go through all loops in the basis
// 		std::cout<<"Loop "<<(itl)<<" :";
// 		bool dead = false;
// 		int alivesize = aliveLoops.size();
// 		for(int ita=alivesize-1; ita>=0;--ita){ //all previous, loops already alive, in matrix, check current
// 			if(LowReturner(aliveLoops[ita])==LowReturner(listofLoops[itl]))// if low is same, add with previous
// 				dead_bit = domain_complex.annotations[sizes ]->sum_two_annotation_with_changed_dst(listofLoops[itl], aliveLoops[ita]);
// 			if(dead_bit == -1){	// not part of the basis class, know which one has died
// 				deadLoopsIndex.push_back(itl);
// 				std::cout<<"died poor soul\n";	
// 				dead = true;
// 				return (itl);
// 			}
// 		}
// 		if(dead == false)
// 			{std::cout<<"is alive\n";	aliveLoopsIndex.push_back(itl);}
// 		else
// 			{std::cout<<"died poor soul\n";	return (itl+1);}
// 	}
// 	int returner = 0;
// 	// cout<<returner<<"\n";
// 	return returner;

// }


bool CheckBoundary (std::vector<std::vector<int>> simplex_vertices){

	vector<SimplicialTreeNode_ptr> boundaries;
	SimplicialTreeNode_ptr simplex ;
	std::vector<int> sizes;
	for(int i=0;i<simplex_vertices.size();i++){
		// SimplicialTreeNode_ptr 
		simplex = domain_complex.find(simplex_vertices.at(i));
		
		boundaries.push_back(simplex);
		sizes.push_back((simplex_vertices.at(i)).size());
		// std::cout<<"Size: "<<(simplex_vertices.at(i)).size();
		// simplex.clear();
	}

	ListNodePtr sum = domain_complex.annotations[sizes.at(0) ]->DeepCopyAnnotationColumn(domain_complex.find_annotation(boundaries.front()));
	int dead_bit = -1;
	for (int i = 1; i < boundaries.size(); ++i){
		if (boundaries[i] == NULL){
			cerr << "one of the sub-simplecies of";
			for (int k = 0; k < simplex_vertices.size(); ++k)
				cerr << " " << simplex_vertices[k][0]<<" "<<simplex_vertices[k][1];
			cerr << " is missing.\n";
			exit(0);
		}
		ListNodePtr tempLNP = domain_complex.find_annotation(boundaries[i]);
		dead_bit = domain_complex.annotations[sizes.at(i) ]->sum_two_annotation_with_changed_dst(sum, tempLNP);
	}
	if(dead_bit== -1)
		return false;
	else return true;
	// return true;
}