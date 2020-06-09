/*
(c) 2015 Fengtao Fan, Dayu Shi
*/

#include "SimplicialComplex.h" 

#include <ctime>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_set>
#include <fstream>

#include <boost/program_options.hpp>

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

void colorVertex(std::unordered_set<int> setVer, string filename, string outputFile)
{
	std::ifstream ifs(filename);
	std::ofstream ofs(outputFile);
	std::ifstream ifsReindex("reindex.txt");
	int numV, numF;
	string ss;
	ifs >> ss;
	ifs >> numV;
	ifs >> numF;	//numF
	ifs >> ss;	//numE
	//setVer has the new indices
	//reading reindex info
	std::unordered_map<int, int> mapRe;
	for (int i = 0; i < numV; ++i)
	{
		int iNew, iOld;
		ifsReindex >> iNew >> iOld;
		mapRe[iOld] = iNew;
	}
	ifsReindex.close();
	ifsReindex.clear();
	//writing OFF file with vertex colors
	ofs << "COFF" << endl;
	ofs << numV << " " << numF << " 0" << std::endl;
	for (int i = 0; i < numV; ++i)
	{
		ifs >> ss;
		ofs << ss << " ";
		ifs >> ss;
		ofs << ss << " ";
		ifs >> ss;
		ofs << ss << " ";
		//convert old index to new index
		if (setVer.find(mapRe[i]) != setVer.end())
		//if (i == 250 || i == 4665)
		{
			ofs << "1.0 0.0 0.0" << std::endl;
		}
		else
		{
			ofs << "0.1 0.1 0.1" << std::endl;
		}
	}
	char line[256];
	ifs.getline(line, 256);
	for (int i = 0; i < numF; ++i)
	{
		ifs.getline(line, 256);
		ofs << line << endl;
	}
	ifs.close();
	ifs.clear();
	ofs.close();
	ofs.clear();
	return;
}

bool barcodeCompare(const pair<int, int>& a, const pair<int, int>& b)
{
	float a1, a2, b1, b2;
	a1 = vecFiltrationScale[a.first];
	a2 = a.second == -1 ? vecFiltrationScale.back() : vecFiltrationScale[a.second];
	b1 = vecFiltrationScale[b.first];
	b2 = b.second == -1 ? vecFiltrationScale.back() : vecFiltrationScale[b.second];
	return a2 - a1 < b2 - b1;

	//death time order
	/*if (a.second != -1 && b.second != -1)
		return (a.second < b.second);
	if (a.second == -1 && b.second != -1)
		return false;
	if (a.second != -1 && b.second == -1)
		return true;
	return a.first < b.first;*/

	////random order
	//srand(time(NULL));
	//int m = rand();
	//int n = rand();
	//return m < n;
}

bool barcodeCompareWithTS(const std::pair<int, pair<int, int>>& a, const std::pair<int, pair<int, int>>& b)
{
	return barcodeCompare(a.second, b.second);
}

void ComputingPersistenceForSimplicialMap(const char* file_name_of_domain_complex,
	bool is_domain_complex_with_annotation,
	const char* file_name_of_range_complex,
	const char* file_name_of_simplicial_map,
	bool is_save_range_complex_with_annotation = false,
	const char* new_range_complex_file_name = NULL);

void ComputingPersistenceForSimplicialMapElementary(const char* file_name_of_domain_complex,
	bool is_domain_complex_with_annotation,
	vector<string>& vecElemOpers,
	bool is_save_range_complex_with_annotation = false,
	const char* new_range_complex_file_name = NULL);

/***********************************************/
char * strLicense = "THIS SOFTWARE IS PROVIDED \"AS-IS\". THERE IS NO WARRANTY OF ANY KIND. "
"NEITHER THE AUTHORS NOR THE OHIO STATE UNIVERSITY WILL BE LIABLE FOR "
"ANY DAMAGES OF ANY KIND, EVEN IF ADVISED OF SUCH POSSIBILITY. \n"
"\n"
"This software was developed (and is copyrighted by) the Jyamiti group at "
"The Ohio State University. Please do not redistribute this software. "
"This program is for academic research use only. This software uses the "
"Boost library (www.boost.org) "
"which is covered under their own licenses.\n"
"\n"
"The Boost library's license "
"(which applies to the Boost library ONLY and NOT to this program itself) is "
"as follows:\n"
"\n"
"LICENSE\n"
"---------------------------------------------------------------------------\n"
"Boost Software License - Version 1.0 - August 17th, 2003\n"
"\n"
"Permission is hereby granted, free of charge, to any person or organization "
"obtaining a copy of the software and accompanying documentation covered by "
"this license (the \"Software\") to use, reproduce, display, distribute, "
"execute, and transmit the Software, and to prepare derivative works of the "
"Software, and to permit third-parties to whom the Software is furnished to "
"do so, all subject to the following: \n"
"\n"
"The copyright notices in the Software and this entire statement, including "
"the above license grant, this restriction and the following disclaimer, "
"must be included in all copies of the Software, in whole or in part, and "
"all derivative works of the Software, unless such copies or derivative "
"works are solely in the form of machine-executable object code generated by "
"a source language processor. \n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, "
"FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT "
"SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE "
"FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, "
"ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER "
"DEALINGS IN THE SOFTWARE. \n"
"---------------------------------------------------------------------------\n";
/**********************************************************************/
bool ParseCommand(int argc, char** argv,
	std::string &input_domain_complex_file_name,
	std::string &input_range_complex_file_folder,
	std::string &input_simplicial_map_file_folder,
	std::string &output_range_complex_with_annotation_file_name,
	std::string &output_persistence_file_name,
	bool &is_input_domain_complex_with_annotation,
	bool &is_output_range_complex_with_annotation,
	bool &is_elementary,
	bool &is_generator,
	float &fThres,
	int &iFiltrationSize,
	int &dim)
{
	try
	{
		/* Define the program options description
		*/
		namespace po = boost::program_options;
		po::options_description desc("Simpers Usage");
		desc.add_options()
			(",h", "Help information;")
			(",l", "License information;")
			(",e", po::value<bool>(&is_elementary)->default_value(true), "The flag indicating input simplicial maps are in elementary mode or general mode (default value: true(elementary) )")
			(",n", po::value<int>(&iFiltrationSize), "Number of input simplicial maps (for general mode only)")
			(",d", po::value<std::string>(&input_domain_complex_file_name)->default_value(""), "The file name for the initial simplicial complex (optional for elementary mode: do not need to specify this parameter if input simplicial complex is empty)")
			("dflag", po::value<bool>(&is_input_domain_complex_with_annotation)->default_value(false), "Optional: The flag indicating to read intial simplicial complex with available annotations or not (default value: false)")
			(",m", po::value<std::string>(&input_simplicial_map_file_folder)->required(), "The file name (elementary mode) or folder name (general mode) for input simplicial maps (see more details in readme file)")
			(",r", po::value<std::string>(&input_range_complex_file_folder), "The folder for input simplicial complexes in the range of each simplicial map (for general mode only)")
			(",s", po::value<std::string>(&output_persistence_file_name)->default_value("pers"), "The file name for the output persistence barcodes of input simplicial maps (default value: \"pers\")")
			(",t", po::value<float>(&fThres)->default_value(0), "The threshold for denoising the output persistence barcodes (default value: 0)")
			("gflag", po::value<bool>(&is_generator)->default_value(false), "Optional: The flag indicating to compute Dim-1 persistence generators or not (default value: false, only supported in elementary mode)")
			("rflag", po::value<bool>(&is_output_range_complex_with_annotation)->default_value(false), "Optional: The flag indicating to save resulting simplicial complex or not (default value: false)")
			("rfile", po::value<std::string>(&output_range_complex_with_annotation_file_name), "Optional: The file name for saving the resulting simplicial complex with annotations (default value: false)")
			("dim", po::value<int>(&dim)->default_value(2), "The maximum dimension user needs in the output persistence barcodes (default value: 2)");
		// Parser map
		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, desc), vm);

			//
			if (vm.count("-h"))
			{
				std::cout << desc << std::endl;
				exit(EXIT_SUCCESS);
			}
			//
			if (vm.count("-l"))
			{
				std::cout << strLicense << std::endl;
				exit(EXIT_SUCCESS);
			}
			//
			po::notify(vm);
		}
		catch (boost::program_options::required_option& e)
		{
			std::cerr << "ERROR: " << e.what() << std::endl;
			return false;
		}
		catch (boost::program_options::error& e)
		{
			std::cerr << "ERROR: " << e.what() << std::endl;
			return false;
		}
	}
	catch (std::exception& e)
	{
		std::cerr << "Unhandled Exception reached the top of main: "
			<< e.what() << ", application will now exit" << std::endl;
		return false;

	}
	return true;
}
/*********************/
int main(int argc, char **argv)
{
	start = std::clock();
	std::string input_domain_complex_file_name;
	std::string input_range_complex_file_folder;
	std::string input_simplicial_map_file;		//file name for elementary mode and folder name for genaral mode
	std::string output_range_complex_with_annotation_file_name;
	std::string output_persistence_file_name;
	bool is_input_domain_complex_with_annotation;
	bool is_output_range_complex_with_annotation;
	bool is_elementary;
	bool is_generator;
	bool bDeathTimeOrder = true;
	bool bTimeStamp = false;
	//indicate if persistence barcode is ordered according to death time
	int numFiltration;
	string sDelimiter("#");
	vecFiltrationScale.push_back(0);
	float fMaxScale = 0.0;
	if (ParseCommand(argc, argv, input_domain_complex_file_name,
		input_range_complex_file_folder,
		input_simplicial_map_file,
		output_range_complex_with_annotation_file_name,
		output_persistence_file_name,
		is_input_domain_complex_with_annotation,
		is_output_range_complex_with_annotation,
		is_elementary,
		is_generator,
		fThreshold,
		numFiltration,
		max_dimension))
	{
		if (!is_elementary)
		{
			domain_complex.bGenerator = false;
			if (input_range_complex_file_folder.empty())
			{
				std::cout << "In general simplicial map mode, parameter -r is required." << endl;
				return 0;
			}
			for (int i = 1; i <= numFiltration; i++)
			{

				string sRC = input_range_complex_file_folder;
				string sMap = input_simplicial_map_file;
				stringstream ss;
				ss << i;
				sRC += "/";
				sRC += ss.str();
				sRC += ".txt";
				sMap += "/";
				sMap += ss.str();
				sMap += ".txt";
				//timer
				timer2 = std::clock();
				if (i == numFiltration)
				{
					filtration_step = i;
					ComputingPersistenceForSimplicialMap(input_domain_complex_file_name.c_str(), is_input_domain_complex_with_annotation,
						sRC.c_str(), sMap.c_str(), is_output_range_complex_with_annotation,
						output_range_complex_with_annotation_file_name.c_str());
				}
				else
				{
					if (i == 1)
						filtration_step = 0;
					else
						filtration_step = i;
					ComputingPersistenceForSimplicialMap(input_domain_complex_file_name.c_str(), is_input_domain_complex_with_annotation,
						sRC.c_str(), sMap.c_str(), false,
						output_range_complex_with_annotation_file_name.c_str());
				}
				dFuncTimeSum += (std::clock() - timer2);
			}
		}
		else
		{
			domain_complex.bGenerator = is_generator;
			if (is_generator)
				bTimeStamp = true;
			ifstream ifs_m;
			ifs_m.open(input_simplicial_map_file);
			if (!ifs_m.is_open())
			{
				std::cout << "can't open file!" << endl;
				return 0;
			}
			filtration_step = 0;
			while (!(ifs_m.eof()))
			{
				vector<string> vecOpers;
				string sOper;
				char sLine[256];
				bool bValid = false;
				while (ifs_m.getline(sLine, 256))
				{
					if (sLine[0] == '#')
					{
						bValid = true;
						stringstream ss;
						ss.str(sLine);
						string sSharp;
						float fScale;
						ss >> sSharp;
						ss >> fScale;
						vecFiltrationScale.push_back(fScale);
						break;
					}
					sOper = sLine;
					vecOpers.push_back(sOper);
				}
				if (bValid)
				{
					//timer
					timer2 = std::clock();
					ComputingPersistenceForSimplicialMapElementary(input_domain_complex_file_name.c_str(),
						is_input_domain_complex_with_annotation,
						vecOpers, false,
						output_range_complex_with_annotation_file_name.c_str());
					dFuncTimeSum += (std::clock() - timer2);
				}
			}
			ifs_m.close();
			ifs_m.clear();
			if (is_output_range_complex_with_annotation)
				domain_complex.WriteComplexWithAnnotation(output_range_complex_with_annotation_file_name.c_str());
		}

		
		//output persistence barcode and generator
		ofstream ofs_dgm(output_persistence_file_name);
		std::vector<int> h1TimeStamp;
		for (int i = 0; i < persistences.size(); ++i)
		{
			std::vector<std::pair<int, int>> vecBarcodes;
			std::vector<std::pair<int, std::pair<int, int>>> vecBarcodesWithTS;
			for (std::unordered_map<int, pair<int, int>>::iterator it = persistences[i].begin(); it != persistences[i].end(); ++it)
			{
				if (bDeathTimeOrder)
				{
					if (bTimeStamp)
						vecBarcodesWithTS.push_back(std::make_pair(it->first, std::make_pair((it->second).first, (it->second).second)));
					else
						vecBarcodes.push_back(std::make_pair((it->second).first, (it->second).second));
				}
				else
				{
					ofs_dgm << i << " ";
					if (bTimeStamp)
						ofs_dgm << it->first << " ";
					ofs_dgm << vecFiltrationScale[(it->second).first] << " ";
					if ((it->second).second == -1)
						ofs_dgm << "inf" << endl;
					else
						ofs_dgm << vecFiltrationScale[(it->second).second] << endl;
				}
			}
			if (bDeathTimeOrder)
			{
				if (bTimeStamp)
				{
					std::stable_sort(vecBarcodesWithTS.begin(), vecBarcodesWithTS.end(), barcodeCompareWithTS);
					for (int k = 0; k < vecBarcodesWithTS.size(); ++k)
					{
						ofs_dgm << i << " ";
						ofs_dgm << vecBarcodesWithTS[k].first << " ";
						ofs_dgm << vecFiltrationScale[vecBarcodesWithTS[k].second.first] << " ";
						if (vecBarcodesWithTS[k].second.second == -1)
							ofs_dgm << "inf" << endl;
						else
							ofs_dgm << vecFiltrationScale[vecBarcodesWithTS[k].second.second] << endl;
						if (i == 1) {
							h1TimeStamp.push_back(vecBarcodesWithTS[k].first);
						}
					}
				}
				else
				{
					std::stable_sort(vecBarcodes.begin(), vecBarcodes.end(), barcodeCompare);
					for (int k = 0; k < vecBarcodes.size(); ++k)
					{
						ofs_dgm << i << " ";
						ofs_dgm << vecFiltrationScale[vecBarcodes[k].first] << " ";
						if (vecBarcodes[k].second == -1)
							ofs_dgm << "inf" << endl;
						else
							ofs_dgm << vecFiltrationScale[vecBarcodes[k].second] << endl;
					}
				}
			}
		}
		ofs_dgm.close();
		ofs_dgm.clear();

		//find the maximum scale
		/*for (int k = 0; k < vecFiltrationScale.size(); ++k)
		{
			if (vecFiltrationScale[k] > fMaxScale)
				fMaxScale = vecFiltrationScale[k];
		}*/
		fMaxScale = vecFiltrationScale.back();

		//output timing result
		dFuncTimeSum = dFuncTimeSum / (double)CLOCKS_PER_SEC;
		ofstream ofsTime("Timing.txt");
		ofsTime << "Total time: " << (std::clock() - start) / (double)CLOCKS_PER_SEC << "s" << endl;
		std::cout << "Total time: " << (std::clock() - start) / (double)CLOCKS_PER_SEC << "s" << endl;
		std::cout << "Insert time: " << (dInsertTime) / (double)CLOCKS_PER_SEC << "s" << endl;
		std::cout << "Collapse time: " << (dCollapseTime) / (double)CLOCKS_PER_SEC << "s" << endl;
		ofsTime << "Simpers computation time: " << (dInsertTime + dCollapseTime) / (double)CLOCKS_PER_SEC << "s" << endl;
		std::cout << "Simpers computation time: " << (dInsertTime + dCollapseTime) / (double)CLOCKS_PER_SEC << "s" << endl;

		ofstream ofs_size("size.txt");
		int maxCS = 0;
		ofs_size << "Scales : Complex Sizes" << endl;
		for (int i = 0; i < vecFiltrationScale.size(); ++i)
		{
			ofs_size << vecFiltrationScale[i] << " : " << complexSizes[i] << endl;
			if (complexSizes[i] > maxCS)
				maxCS = complexSizes[i];
		}
		ofs_size << "Cumulative Complex size: " << accumulativeSizes.back() << endl;
		ofs_size << endl;
		ofs_size.close();
		ofs_size.clear();

		ofsTime << "Max Complex Size: " << maxCS << endl;
		std::cout << "Max Complex Size: " << maxCS << endl;
		ofsTime << "Cumulative Complex Size: " << accumulativeSizes.back() << endl;
		std::cout << "Cumulative Complex Size: " << accumulativeSizes.back() << endl;
		ofsTime << "Max Scale: " << fMaxScale << endl;
		std::cout << "Max Scale: " << fMaxScale << endl;
		ofsTime.close();
		ofsTime.clear();

		if (domain_complex.bGenerator &&  bTimeStamp)
		{
			ofstream ofs_gen("generator.txt");
			std::vector<int> path;
			int tsGen;
			while (true) {
				cout << "Do you want to see the generators for top K longest persistence bars ('Y') or manually pick one ('N') or quit ('Q')? : ";
				string choice;
				cin >> choice;
				if (choice == "Y" || choice == "y") {
					int k;
					cout << "Please enter K value (only output at most K generators): ";
					cin >> k;
					ofs_gen << "----------------------------------------------" << endl;
					ofs_gen << "Top " << k << " longest generators: " << endl;
					int n = h1TimeStamp.size();
					for (int i = n - 1; i >= n - k && i >= 0; i--) {
						int tsGen = h1TimeStamp[i];
						path = domain_complex.gen1[tsGen];
						ofs_gen << "Generator for ID " << tsGen << ":" << endl;
						cout << "Generator for ID " << tsGen << ":" << endl;
						for (int i = 0; i < path.size(); ++i) {
							ofs_gen << path[i] << " ";
							cout << path[i] << " ";
						}
						ofs_gen << endl << endl;
						cout << endl << endl;
					}
				}
				else if (choice == "N" || choice == "n" ) {
					ofs_gen << "----------------------------------------------" << endl;
					ofs_gen << "Manually pick generators: " << endl;
					do {
						cout << "Which generator to show? (enter generator ID or -1 to exit): ";
						cin >> tsGen;
						if (tsGen == -1) break;
						if (domain_complex.gen1.count(tsGen) == 0) {
							cout << "Invalid ID" << endl << endl;
							continue;
						}
						path = domain_complex.gen1[tsGen];
						ofs_gen << "Generator for ID " << tsGen << ":" << endl;
						cout << "Generator for ID " << tsGen << ":" << endl;
						for (int i = 0; i < path.size(); ++i) {
							ofs_gen << path[i] << " ";
							cout << path[i] << " ";
						}
						ofs_gen << endl << endl;
						cout << endl << endl;
					} while (tsGen != -1);
				}
				else if (choice == "Q" || choice == "q") break;
				else {
					cout << "Invalid input, please retry." << endl;
				}
			}
			ofs_gen.close();
			ofs_gen.clear();
		}

	}
	return 0;
}
