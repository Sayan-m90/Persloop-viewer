/*
(c) 2015 Fengtao Fan, Dayu Shi
*/
#include "SimplicialComplex.h"
#include <fstream>
#include <sstream>
#include <string>

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

int max_dimension = 3;

//timer
std::clock_t start, timer1, timer2;
double dFuncTimeSum;
double dInsertTime;
double dCollapseTime;

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

void ComputingPersistenceForSimplicialMapElementary(const char* file_name_of_domain_complex, bool is_domain_complex_with_annotation,
													vector<string>& vecElemOpers,
													bool is_save_range_complex_with_annotation = false,
													const char* new_range_complex_file_name = NULL) 
{
	if (filtration_step == 0)
	{
		//accumulativeSize = 0;

		if (file_name_of_domain_complex[0] != '\0')
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
	//read and perform elementary operations
	for (int i = 0; i < vecElemOpers.size(); ++i)
	{
		stringstream iss;
		iss.str(vecElemOpers[i]);
		string sOperation;
		iss >> sOperation;
		vector<int> simplex;
		vector<int> removeVertices;
		int preserveVertex;
		if (sOperation.compare("i") == 0)
		{
			//assume no simplex is added twice
			int index;
			while (iss >> index)
			{
				simplex.push_back(index);
			}
			//timer
			timer1 = std::clock();
			//elementary insertion
			domain_complex.ElementaryInsersion(simplex);
			simplex.clear();
			dInsertTime += std::clock() - timer1;
		}
		else if (sOperation.compare("c") == 0)
		{
			string s;
			istringstream issConvert;
			int removeLabel;
			while (iss >> s)
			{
				if (s.compare("t") == 0)
					break;
				issConvert.str(s);
				issConvert >> removeLabel;
				removeVertices.push_back(removeLabel);
				s.clear();
				issConvert.clear();
				issConvert.str("");
			}
			iss >> preserveVertex;
			//timer
			timer1 = std::clock();
			//elementary collapse
			for (int i = 0; i < removeVertices.size(); ++i)
			{
				collapseCount++;
				if (removeVertices[i] > preserveVertex) smallCount++;
				domain_complex.ElementaryCollapse(removeVertices[i], preserveVertex);

			}
			dCollapseTime += std::clock() - timer1;
			removeVertices.clear();
		}
	}
	if (is_save_range_complex_with_annotation) {
		domain_complex.WriteComplexWithAnnotation(new_range_complex_file_name);
	}
	complexSizes.push_back(domain_complex.ComplexSize());
	accumulativeSizes.push_back(domain_complex.accumulativeSimplexSize);
	filtration_step += 1;
	return;
}