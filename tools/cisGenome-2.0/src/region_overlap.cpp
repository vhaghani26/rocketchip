/*******************************
* Ma, Wenxiu
* updated 2008.3.27
* to calculate the pairwise overlap among multiple regions
********************************/


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <map>

using namespace std;

#define INPUT_LINE_MAX 65535

#define WINDOWS_COMPILER

string itoa(int num){	
	// convert int to string;
	std::string buf;
	//buf.clear();
	if(num == 0 )
		return "0";
	while(num!=0)
	{
		buf.append(std::string(1,'0'+num%10));
		num=num/10;
	}
	// convert the buf in the reorder
	char c1;
	for(int temp= 0; temp<(int)buf.size()/2; temp ++)
	{	
		c1 = buf[temp];
		buf[temp] = buf[buf.size()-1- temp];
		buf[buf.size()-1- temp] = c1;
	}
	return buf;
};

double atod(string num){
	//convert string to double;
	bool neg = false;
	if( num[0] == '-'){
		neg = true;
		num = num.substr(1, string::npos);
		}
	string::size_type dotIndex = num.find('.');
	string intPart, fracPart;
	if(dotIndex == string::npos){
		//integar number
		intPart = num;
		fracPart.clear();
		}
	else{
		//real number
		intPart = num.substr(0, dotIndex);
		fracPart = num.substr(dotIndex+1, string::npos);
		}
	double result = atoi(intPart.c_str());
	double fraction = atoi(fracPart.c_str());
	for ( int i = 0; i < (int) fracPart.size(); i++){
		fraction = fraction/10;
		}
	result = result + fraction;
	if(neg)
		return -result;
	else
		return result;
	};

class peak{
public:
	string id;
	int int_id;
	string chr;
	int start;
	int end;
	string strand;
	void print2file(ofstream& ofs){
		ofs<<id<<"\t"
			<<chr<<"\t"
			<<start<<"\t"
			<<end<<"\t"
			<<strand<<endl;
		};
};

inline bool operator<(const peak &p1, const peak &p2){
	return (p1.chr.compare(p2.chr)<0) || ((p1.chr.compare(p2.chr)==0)&& (p1.start < p2.start));
	};

void write_vector_peak_2file(vector<peak>& peaks, string filename){
	ofstream myout(filename.c_str());
	for(int i=0; i<(int)peaks.size();i++){
		myout<<peaks[i].id<<"\t"
			<<peaks[i].chr<<"\t"
			<<itoa(peaks[i].start)<<"\t"
			<<itoa(peaks[i].end)<<"\t"
			<<peaks[i].strand<<endl;
		}
	myout.close();	
}

int main(int argc, char *argv[]){
	string usage = "Usage:	./region_overlap -n N -i peak1.cod [peak2.cod ... peakN.cod] -o OUTPUT_title [-p min_overlap_rate]";
	if (argc < 7){
		cout <<usage<<endl;
		return -1;
		}
	
	vector<vector<peak> > peaks;
	int number_input_peak_file = 0;
	vector<string> input_file_names;
	string output_file_title;
	double min_overlap_rate  = 0;
	bool min_overlap = false;
	vector<map<int, int> > maps_index2sorted;
	
	
	int i,j;


	for(i=1;i<argc;i++){
		if(string(argv[i]) == "-n")
			number_input_peak_file= atoi(argv[++i]);
		else if(string(argv[i]) == "-i"){
			for(j=0;j<number_input_peak_file;j++)
				input_file_names.push_back(argv[++i]);
			}
		else if(string(argv[i]) == "-p"){
			min_overlap = true;
			min_overlap_rate= atod(argv[++i]);
			}
		else if(string(argv[i]) == "-o")
			output_file_title= argv[++i];
		else{
			cout<<"Error!Wrong command line parameters!"<<endl
				<<usage<<endl;
			return -2;
			}
		}

	if(number_input_peak_file<=1||(int)input_file_names.size()!=number_input_peak_file||min_overlap_rate<0){
		cout<<"Error!Command line parameters do not match!"<<endl
			<<usage<<endl;
		return -3;
		}

	peaks.resize(number_input_peak_file);
	maps_index2sorted.resize(number_input_peak_file);
		
	//get input;
	for(i=0;i<(int)input_file_names.size();i++){
		ifstream myfin;
		myfin.open(input_file_names[i].c_str());
		j=0;
		if(!myfin){
			cout<<"Cant open file: "<<input_file_names[i]<<endl;
			return -3;
			}
		char line[INPUT_LINE_MAX];
		while(myfin.getline(line, INPUT_LINE_MAX)){
			peak newelement;
			istringstream iss(line);
			string datavalue;
			iss>> datavalue; //id
			if(datavalue.size() == 0){
				cout<< "Error! Wrong Input File Format! " <<endl;
				return -3;
				}
			if(datavalue.substr(0,1)=="#")
				continue;
			newelement.id = datavalue;

			iss>> datavalue; //chr
			if(datavalue.size() == 0){
				cout<< "Error! Wrong Input File Format! " <<endl;
				return -3;
				}
			newelement.chr= datavalue;

			iss>> datavalue; //start
			if(datavalue.size() == 0){
				cout<< "Error! Wrong Input File Format! " <<endl;
				return -3;
				}
			newelement.start = atoi(datavalue.c_str());

			iss>> datavalue; //end
			if(datavalue.size() == 0){
				cout<< "Error! Wrong Input File Format! " <<endl;
				return -3;
				}
			newelement.end = atoi(datavalue.c_str());

			iss>> datavalue; //strand
			if(datavalue!="+" && datavalue!="-" &&
				datavalue!="F" && datavalue!="f" &&
				datavalue!="R" && datavalue!="r")
				datavalue="+";
			newelement.strand = datavalue;

			newelement.int_id = j++;
			peaks[i].push_back(newelement);

			}
		myfin.close();
		cout<<input_file_names[i]<< " loaded..."<<endl;
		cout<<"\tnumber of regions = "<<peaks[i].size()<<endl;;
		}

	for(i=0; i<(int)peaks.size();i++){
		sort(peaks[i].begin(),peaks[i].end());
		cout<<input_file_names[i]<< " sorted..."<<endl;
		for(j=0;j<(int)peaks[i].size();j++){
			maps_index2sorted[i].insert(make_pair(peaks[i][j].int_id, j));
			}
		}

	cout<<"Counting overlaps... "<<endl;	

	vector<vector<int> > peaks_overlap_counts;
	peaks_overlap_counts.resize(peaks.size());
	for(i=0;i<(int)peaks.size();i++){
		peaks_overlap_counts[i].clear();
		for(j=0;j<(int)peaks.size();j++){	
		peaks_overlap_counts[i].push_back(0);
			}
		}
	
	vector< vector< vector< vector<int> > > > peak1_overlap;
	vector< vector< vector< vector<int> > > > peak1_50_overlap;
	peak1_overlap.resize(peaks.size());
	peak1_50_overlap.resize(peaks.size());
	
	ofstream myfout;
	string output_file;
	int i_iterator, j_iterator;
	for(i=0;i<(int)peaks.size();i++){
		peak1_overlap[i].resize(peaks.size());
		peak1_50_overlap[i].resize(peaks.size());
		for(j=0;j<(int)peaks.size();j++){
			if(i==j)
				continue;
			peak1_overlap[i][j].clear();
			peak1_50_overlap[i][j].clear();
			peak1_overlap[i][j].resize(peaks[i].size());
			peak1_50_overlap[i][j].resize(peaks[i].size());
			i_iterator =0;
			j_iterator=0;
		 	while(i_iterator<(int)peaks[i].size() &&  j_iterator<(int)peaks[j].size()){
				if(peaks[i][i_iterator].chr.compare(peaks[j][j_iterator].chr)<0 ) {
					i_iterator++;
					continue;
					}
				if(peaks[j][j_iterator].chr.compare(peaks[i][i_iterator].chr)<0 ){
					j_iterator++;
					continue;
					}
				if(peaks[i][i_iterator].end <= peaks[j][j_iterator].start){
					i_iterator++;
					continue;
					}
				if(peaks[i][i_iterator].start >= peaks[j][j_iterator].end){
					j_iterator++;
					continue;
					}
				int size1=peaks[i][i_iterator].end - peaks[i][i_iterator].start+1;
				int size2=peaks[j][j_iterator].end - peaks[j][j_iterator].start+1;
				int min_size = size1;
				if (size2<size1)
					min_size = size2;
//				int max_size = size1+size2-min_size;
				int overlap_size;
				if(peaks[i][i_iterator].end >= peaks[j][j_iterator].end){
					if(peaks[i][i_iterator].start <=peaks[j][j_iterator].start){
						overlap_size = size2;
						peak1_overlap[i][j][i_iterator].push_back(j_iterator);
						peak1_50_overlap[i][j][i_iterator].push_back(j_iterator);
						}
					else{
						overlap_size = peaks[j][j_iterator].end-peaks[i][i_iterator].start;
						if(min_overlap && overlap_size>= (int)min_size*min_overlap_rate){
							peak1_50_overlap[i][j][i_iterator].push_back(j_iterator);
							}
						peak1_overlap[i][j][i_iterator].push_back(j_iterator);
						}
					j_iterator++;
					}
				else{
					if(peaks[i][i_iterator].start >=peaks[j][j_iterator].start){
						overlap_size = size1;
						peak1_overlap[i][j][i_iterator].push_back(j_iterator);
						peak1_50_overlap[i][j][i_iterator].push_back(j_iterator);
						}
					else{
						overlap_size = peaks[i][i_iterator].end-peaks[j][j_iterator].start;
						if(min_overlap && overlap_size>= (int)min_size*min_overlap_rate){
							peak1_50_overlap[i][j][i_iterator].push_back(j_iterator);
							}
						peak1_overlap[i][j][i_iterator].push_back(j_iterator);
						}
					i_iterator++;
					}
				}
			}
		myfout.close();
		}

	ofstream myfout2;
	for(i=0;i<(int)peaks.size();i++){
		output_file= output_file_title + string(".overlap")+string(itoa(i+1));
		myfout2.open(output_file.c_str());
		for(i_iterator = 0; i_iterator<(int)peaks[i].size();i_iterator++){
			map<int,int>::iterator map_it;
			map_it = maps_index2sorted[i].find(i_iterator);
			if(map_it ==maps_index2sorted[i].end()){
				cout<<"Interval Error!"<<endl;
				return -4;
				}
			int true_iterator = map_it->second;
//			cout<<true_iterator<<endl;
			for(j=0;j<(int)peaks.size();j++){
				if(i==j)
					continue;
				if(!min_overlap){
					if((int)peak1_overlap[i][j][true_iterator].size()> 0){
						peaks_overlap_counts[i][j] += 1;
						myfout2<<"1\t";
						}
					else 
						myfout2<<"0\t";
					}
				else{
					if((int)peak1_50_overlap[i][j][true_iterator].size()> 0){
						peaks_overlap_counts[i][j] += 1;
						myfout2<<"1\t";
						}
					else 
						myfout2<<"0\t";
					}
				}
			myfout2<<endl;
			}
		myfout2.close();
		}
	
	for(i=0;i<(int)peaks.size();i++){
		peaks_overlap_counts[i][i] = 0;
		for(i_iterator = 0; i_iterator<(int)peaks[i].size();i_iterator++){
			bool true_self = true;
			for(j=0;j<(int)peaks.size();j++){
				if(i == j)
					continue;
				if(!min_overlap){
					if((int)peak1_overlap[i][j][i_iterator].size()>0){
						true_self = false;					
						break;
						}
					}
				else{
					if((int)peak1_50_overlap[i][j][i_iterator].size()>0){
						true_self = false;					
						break;
						}
					}
				}
			if(true_self){
				peaks_overlap_counts[i][i] +=1;
				}
			}
		}
	
	ofstream myfout3;
	output_file= output_file_title + string(".counts");
	myfout3.open(output_file.c_str());
	for(i=0;i<(int)peaks.size();i++){
		for(j=0;j<(int)peaks.size();j++){
			myfout3<<peaks_overlap_counts[i][j]<<"\t";
			}
		myfout3<<endl;
		}
	myfout3.close();
	
	return 0;
	
}
