/*
 * TRACE v1.0 - Topological Ring and Additive-Coordinated Cage Explorer
 * Copyright (C) 2025 Jun Wei Hsu, Shiang-Tai Lin
 * COMET Laboratory, National Taiwan University

 * This file is part of TRACE.

 * TRACE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.

 * TRACE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with TRACE. If not, see https://www.gnu.org/licenses/.
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <limits>
#include <omp.h>
#include <time.h>
#include <unordered_map>
#include <unordered_set>
#include <climits>
#include <parallel/algorithm>
using namespace std;
//========================================class and structure=====================================================
double get_elapsed_time(timespec start, timespec end) {
        return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}
struct point {
	double x,y,z;
};
struct hbond {
	string dtype,atype;
	int d,a,dh;
	double th;//theta_cut
};
struct ring {
        vector<int>  vertex;
        vector<pair<int,int> >  edge;
	vector<pair<int,int> > sort_edge;
	vector<double> IA;//interior angle
};
struct polyhedron{
	vector<int>  vertex;
	vector<pair<int,int> >  edge; 
	vector<vector<int> > face;
	string type;
	int total;
	int perfect;//0 incomplete 1 complete non EFSC 2 perfect complete cage
	vector<int> guest_idx;
	int found;
	vector<pair<int,double> > add_list;
	int cluster;
	point cp;
	int color_tier;
        polyhedron() : found(0),total(0),cluster(-1) {}
};
struct pair_hash {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const {
        size_t seed = 0;
        seed ^= hash<T1>{}(p.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash<T2>{}(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};
struct VectorHash {
	size_t operator()(const vector<int>& v) const {
		hash<int> hasher;
		size_t seed = 0;
		for (int i : v) {
			seed ^= hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};
class cage
{
public:
	// class constructor & destructor
	cage();
	~cage();
	//input.txt file & output.txt file
	ifstream ifs,ifs_water,ifs_guest,ifs_add,ifs_hbond;
	ofstream ofs_summary,ofs_detail,ofs_visual_gro,ofs_occupancy,ofs_cluster,ofs_visual_index,ofs_crystallinity,ofs_ring_count,ofs_ring_detail;
	//variables
	point **H2O_molecule,*guest_molecule,**add_molecule,box,box_2;
	int H2O,H2O_atom;//H2O molecule
	int guest,guest_atom;//guest molecule
	int add,add_atom;//additive molecule
	vector<string> add_name;
	vector<int> add_atoms;//the atoms that participate in the hydrogen bond
	unordered_map<pair<int, int>, int, pair_hash> adj_map;
	vector<vector<vector<vector<int> > > > grid;
	vector<vector<int> > H2O_adjlist,cage_adjlist; 
	vector<int> isguest;
	vector<int> H2O_color;//use for vmd visualize
	int nt;
	int initial_frame,end_frame;
	int fx,fy,fz;
	bool grid_flag;
	vector<hbond> hbond_H2O_add,hbond_add_add;
	double dt,rcut,theta,cos_theta,acut;//water-water rcut water-additive acut
	double IA_tolerance,DA_tolerance;
	int file_format1,file_format2,file_format3,max_ring;//determine gro file format
	bool cal_guest,cal_additive,cal_box;
	string cl,rm;
	vector<ring > allring;
	//cage and filled cage types count
	int cage512,cage51262,cage51263,cage51264,cage51265,cage51266,cage51268,cage4151062,cage4151063,cage4151064,cage4151065,cage425864,cage425863,cage425862,cage425861,cage435663,cage435664,other;
	//function
	void build_grid();
	//input and output function
	void input_guest(int atom_number,point *molecule,int molecule_num);//input geust molecule
	bool input_H2O(int atom_number,point **molecule,int molecule_num);//input all molecule information of H2O (all atom)
	void input_add(int atom_number,point **molecule,int molecule_num);
	bool input_hbond();
	void cal_and_output(int i);
	void print1D(vector<int> vec);
        void print2D(vector<vector<int> > vec);
	void printP(polyhedron p,ofstream& ofs_summary);
	void printP2(polyhedron p);
	void write_ring_summary(int i);
	void classify_and_process_cages(int i,int& fill,vector<polyhedron> &cage,vector<polyhedron> &SEC_cage,vector<polyhedron> &non_SEC_cage,vector<polyhedron> &IC_cage,vector<int> &cage_form);
	void process_and_write_clusters(int i,const vector<polyhedron> &SEC_cage,const vector<polyhedron> &non_SEC_cage,const vector<polyhedron> &IC_cage);
	void write_cage_details(const int& i,const vector<int> &cage_form,const vector<polyhedron> &SEC_cage,const vector<polyhedron> &non_SEC_cage,const vector<polyhedron> &IC_cage);
	void write_cage_summary(int& i,int& fill,int& SEC_cage_size, int& non_SEC_cage_size, int& IC_cage_size);	
	void output_visual(int i,const vector<polyhedron> &SEC_cage,const vector<polyhedron> &non_SEC_cage,const vector<polyhedron> &IC_cage);
	//function of calulation
	double PBC_dist(point a,point b);//cal the distance between a molecule and b in pbc
	double PBC_angle(point a,point b,point c);
	point PBC_vector(point a, point b);
        point PBC_average(const vector<point>& p);
	point PBC_average_P(const polyhedron& P);
	pair<bool,int> define_molecule(string* type, int *atom_number,string filename,int &molecule_num,bool additive = false);//define the numbers of atom per molecule and molecule type
	void ring_dfs(int istart,int& ispath,unordered_set<vector<int>, VectorHash>& path,ring& local_ring,vector<int>& visited,point n1,int level,int max_level,vector<ring> &local_rings);
	pair<bool,double> check_planarity(int& level,ring& ring,point n1);
        bool check_interior_angle(int& level,ring& ring);
        bool check_isring(int& level,ring& ring);
	void find_ring();
	void find_ring_group(const vector<ring >& ring_groups,unordered_map<pair<int, int>, vector<ring* >, pair_hash >& grouped);//when rings share the same edge count as same group
	void ring_edge(ring &r);
	void find_cup(const vector<ring >& ring_groups,const unordered_map<pair<int, int>, vector<ring*>, pair_hash>& grouped,vector<polyhedron*>& cup_group);//fully coordinated rings count as a cup
	void find_cup_group(const vector<polyhedron*>& cup_groups,unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>& grouped);//the bottom of cup is key and the cup is key value
	pair<bool, vector<polyhedron*> > containecup(const vector<int>& face1,const vector<int>& face2,unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>& grouped);
	bool contain_face(const vector<int>& vec, const vector<vector<int> >& vecs);
	bool share_edge(const vector<pair<int,int> >& vec1, const vector<pair<int,int> >& vec2);
	void find_cage(vector<polyhedron*>& cup_group,vector<polyhedron>& cage,unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>& grouped,int& type);
	//find cage cluster
        void built_cage_matrix(vector<polyhedron>& ALLcage);
	void cage_dfs(int i, vector<int>& cluster, vector<int>& visited);
	void find_cluster(vector<vector<int> >& groups,const int& n);
	//other function
	string intToString(int number);
	void built_hbond_map();
	void define_cage(polyhedron &P,int iframe,int perfect,int &fill);
	string mapToString(const map<int, int>& face_count);
        int SEC_cage_identify(polyhedron &P);
        point cross_product(point v1, point v2);
        double dihedral(point n1, point n2);
	void init();
        //recursive funtion
	void lateral_ring_combine(
        	vector<vector<ring*> >& lateral_ring_grouped_vec,
        	vector<ring*>& current_combination,
        	vector<vector<ring*> >& all_combinations,
        	size_t idx,bool& iscup);
	void lateral_cup_combine(const unordered_map<vector<int>, vector<polyhedron*> , VectorHash>& lateral_cup_grouped,
             	vector<polyhedron*>& current_combination,
             	vector<vector<polyhedron*> >& all_combinations,
             	unordered_map<vector<int>, vector<polyhedron*> , VectorHash>::const_iterator it,
	     	unordered_map<vector<int>, vector<polyhedron*> , VectorHash>::const_iterator end,int& type);
	void find_cage_ICO(int layer, int max_layer,
             	unordered_set<vector<int>, VectorHash>& uni_face,unordered_set<int>& uni_vertice,
             	vector<vector<int> >& remaining_face,
             	vector<vector<int> >& pre_layer_face,
             	unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>& grouped,int& iscage,
		vector<polyhedron>& cage,int& type);
};
void cage::init(){
	//Initialize memory
	max_ring--;
	H2O_molecule=new point*[H2O];
	for (int i = 0; i < H2O; ++i) {
		H2O_molecule[i] = new point[3];
    	}
	if(cal_guest){
		guest_molecule=new point[guest];
        	isguest.resize(guest,0);
	}
	if(cal_additive){
		add_molecule=new point*[add];
		for (int i = 0; i < add; ++i) {
                	add_molecule[i] = new point[add_atom];
        	}
	}
	H2O_adjlist.resize(H2O + add);
	H2O_color.resize(H2O,10);
}
cage::cage()
{

};
cage::~cage()
{

};
//==============================class and structure==========================================================
//==================================main code================================================================
int main(int argc, char* argv[])
{
	cout << "\033[35m"
		<< "Please ensure that the H2O.gro file is formatted as follows:\n"
		<< "  - The first atom must be O (oxygen).\n"
		<< "  - If present, H1 and H2 atoms must immediately follow as the second and third atoms, respectively.\n"
		<< "  - O       : 1-point water\n"
		<< "  - OHH     : 3-point water\n"
		<< "  - OHHM    : 4-point water\n"
		<< "  - OHHLL   : 5-point water\n\n"
		<< "The guest molecule uses the first atom as its center of mass.\n"
		<< "If you want to include multiple types of guest molecules, convert them to single-point representations\n"
		<< "using gmx make_ndx,and then merge them into a single .gro file.\n"
		<< "\033[93m"
		<< "This strict ordering is required for correct identification and processing of water molecules.\n"
		<< "\033[0m";
	cout << "\033[34m"
	     	<< "\n# Format of detail_cage.txt\n"
	     	<< "1. F is the number of polygonal faces (rings), E is the number of edges, and V is the number of vertices in the cage.\n"
	     	<< "2. t is the theoretical number of vertices, computed assuming all faces are isolated (i.e., no shared edges or vertices).\n"
	     	<< "3. Cage labels:\n"
	     	<< "     - #cage n   : SEC        (standard edge-saturated cages)\n"
	    	<< "     - @cage n   : non-SEC    (nonstandard edge-saturated cages)\n"
	     	<< "     - !cage n   : IC         (incomplete cage)\n"
	     	<< "     - #a-cage n : a-SEC      (additive-coordinated standard edge-saturated cages)\n"
                << "     - @a-cage n : a-non-SEC  (additive-coordinated nonstandard edge-saturated cages)\n"
                << "     - !a-cage n : a-IC       (additive-coordinated incomplete cage)\n"
		<< "     where `n` indicates the cluster ID the cage belongs to.\n"
	     	<< "4. Molecule indices:\n"
	    	<< "     - Water     : 1 to Nw\n"
	     	<< "     - Additive  : Nw+1 to Nw+Na\n"
	     	<< "     - Note: molecule indices do not wrap at 99999; they continue (e.g., 100000, 100001, ...).\n"
	     	<< "5. section `CP [ x y z ]` lists the center point of cage.\n"
		<< "6. Section `V [ ... ]` lists indices (1 to Nw+Na) of water and additive molecules in the cage.\n"
		<< "7. Section `g [ ... ]` lists guest molecule indices (1 to Ng) occupying the cage.\n"
		<< "8. Section `a [ ... ]` lists up to 3 additives closest to the cage center, formatted as: add_index:distance (nm).\n"
		<< "(Note: Nw = number of water molecules, Ng = number of guest molecules, Na = number of additive molecules)\n"
		<< "(Note: Indices are based on the input sequence, not on .gro file residue IDs)\n"
		<< "\033[0m" << endl;
		cout << "\033[33m";
		cout << "(Required parameter)\n"
		<< "-w   Ex: H2O.gro              H2O file\n"
		<< "\033[0m" << endl;
		cout << "\033[32m(Optional parameters)\n"
		<< "-g   Ex: CO2.gro                    Guest molecule file\n"
		<< "-a   Ex: add.gro                    Additive molecule file\n"
		<< "-h   Ex: hbond.txt                  Hydrogen bond definition file\n"
		<< "-b   initial frame                  Starting frame index (default: 0)\n"
		<< "-e   end frame                      Ending frame index\n"
		<< "-si  shift frame                    Frame shift interval (default: 0)\n"
		<< "-r   cutoff radius (nm)             H-bond O–O cutoff distance (default: 0.36)\n"
		<< "-th  cutoff angle (°)               H-bond O–O–H angle threshold (default: 35; <0 disables angle check)\n"
		<< "-DA  dihedral tolerance (°)         Tolerance for dihedral angle (default: 90)\n"
		<< "-IA  Interior angle tolerance (°)   Tolerance for interior angle  angle (default 20, applied to rings ≥7, should be <30)\n"
		<< "-mr  max ring size (6–12)           Maximum n-membered ring to consider (default: 10)\n"
		<< "-cl  yes | no                       Include ICs in cage cluster detection (default: yes)\n"
		<< "-nb  nx ny nz                       Number of sub-boxes along x, y, z (default: box_length / 1.2 nm, rounded)\n"
		<< "-nt  threads                        Number of threads to use (default: 1)\n"
		<< "(note relative or absolute path accepted)"
		<< "\033[0m" << endl;
		
	//cout << "\033[32m[Note] We recommend using -nb to partition the box. Sub-box edge lengths of ~1.2 nm in x, y, and z are effective and do not affect accuracy, especially for medium to large systems (>10,000 molecules).\033[0m" << endl;

	
	cage cage;
	string type;
	string H2O_filename,Guest_filename,add_filename,hbond_filename;
	string nbox;
	int shift_frame=0;
	//default
	cage.rcut=0.36;//O-O distance
	cage.acut=0.36;
	cage.theta=35;
	hbond_filename="hbond.txt";
	string fname_summary       = "cage.txt";
	string fname_detail        = "detail_cage.txt";
	string fname_visual_gro    = "visual.gro";
	string fname_visual_index  = "visual_index.txt";
	string fname_occupancy     = "occupancy.txt";
	string fname_cluster       = "cluster.txt";
	string fname_crystallinity = "crystallinity.txt";
	string fname_ring_count    = "ring.txt";
	string fname_ring_detail   = "ring_detail.txt";
	cage.DA_tolerance=90;
        cage.IA_tolerance=20;
	cage.cal_guest=false;
	cage.cal_additive=false;
	cage.cl="yes";
	cage.rm="no";
	cage.initial_frame=0,
	cage.end_frame=INT_MAX;
	cage.max_ring=10;
	cage.nt=1;
	cage.fx=1;
	cage.fy=1;
	cage.fz=1;
	cage.cal_box=true;
	cage.grid_flag=false;
	//argument parameter
	for (int i = 1; i < argc; i++) {
		//======These parameters are required inputs========
        	if (string(argv[i]) == "-w" && i + 1 < argc) {
            		H2O_filename = argv[i + 1];
            		i++; //next parameter
        	}
		//======These parameters are optional========
		else if (string(argv[i]) == "-g" && i + 1 < argc) {
                        Guest_filename = argv[i + 1];
                        i++; //next parameter
			cage.cal_guest=true;
                }
		else if (string(argv[i]) == "-a" && i + 1 < argc) {
                        add_filename = argv[i + 1];
                        i++; //next parameter
			cage.cal_additive=true;
                }
                else if (string(argv[i]) == "-h" && i + 1 < argc) {
                        hbond_filename = argv[i + 1];
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-b" && i + 1 < argc) {
                        cage.initial_frame = atoi(argv[i + 1]);
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-e" && i + 1 < argc) {
                        cage.end_frame = atoi(argv[i + 1]);
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-cl" && i + 1 < argc) {
                        cage.cl = argv[i + 1];
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-si" && i + 1 < argc) {
                        shift_frame = atoi(argv[i + 1]);
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-mr" && i + 1 < argc) {
                        cage.max_ring = atoi(argv[i + 1]);
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-nt" && i + 1 < argc) {
                        cage.nt = atoi(argv[i + 1]);
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-nb" && i + 3 < argc) {
			nbox = string(argv[i + 1]) + " " + argv[i + 2] + " " + argv[i + 3];
			//cage.cal_box = true;
			i += 3;
			cage.grid_flag=true;
		}
		else if (string(argv[i]) == "-r" && i + 1 < argc) {
                        cage.rcut = atof(argv[i + 1]);
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-th" && i + 1 < argc) {
                        cage.theta = atof(argv[i + 1]);
                        i++; //next parameter
                }
		else if (string(argv[i]) == "-DA" && i + 1 < argc) {
                        cage.DA_tolerance = atof(argv[i + 1]);
                        i++; //next parameter
                }
                else if (string(argv[i]) == "-IA" && i + 1 < argc) {
                        cage.IA_tolerance = atof(argv[i + 1]);
                        i++; //next parameter
                }
        }
	omp_set_num_threads(cage.nt);
        //=============================input=======================
	//H2O molecule
	ifstream infile;
	infile.open(H2O_filename.c_str());
	if (infile){ 
		pair<bool,int> format;
		cout<<H2O_filename<<" opened successfully!"<<endl;
		format=cage.define_molecule(&type,&cage.H2O_atom,H2O_filename,cage.H2O);
		cage.file_format1=format.second;
		if(!format.first)return 0;
		if (cage.H2O_atom == 2) {
			cerr << "\033[1;31m" << "Error: H2O GRO file contains only 2 atoms, which is invalid." << "\033[0m" << endl;
			return 0;
		}
		if (cage.H2O_atom == 1)
			cage.theta=-1; 
	}
	else{
		cerr<<"\033[1;31m"<<"No H2O input file found!"<< "\033[0m"<<endl;return 0;
	}
	infile.close();
	//additive
	if(cage.cal_additive){
		infile.open(add_filename.c_str());
		if (infile){
			pair<bool,int> format;
			cout<< add_filename <<" opened successfully!"<<endl;
			format=cage.define_molecule(&type,&cage.add_atom,add_filename,cage.add,cage.cal_additive);
			cage.file_format3=format.second;
			if(!format.first)return 0;
		}
		else{
			cerr<<"\033[1;31m"<<"No additive input file found!"<< "\033[0m"<<endl;return 0;		
			return 0;
		}
	}
        else{ 
		cout<<"Additive molecules will not be included in the calculations."<<endl<<endl;
		cage.add=0;
        }
        infile.close();
	//guest molecule
	if(cage.cal_guest){
		infile.open(Guest_filename.c_str());
		if (infile){
			pair<bool,int> format;
			cout<< Guest_filename <<" opened susscesful!"<<endl;
			format=cage.define_molecule(&type, &cage.guest_atom,Guest_filename,cage.guest);
			cage.file_format2=format.second;
			if(!format.first)return 0;
		}
		else{
                        cerr<<"\033[1;31m"<<"No guest input file found!"<< "\033[0m"<<endl;return 0;
                        return 0;
                }
	}
	else{
		cout<<"Guest molecules will not be included in the calculations."<<endl<<endl;
	}
	infile.close();
	if(cage.cl=="no" || cage.cl=="yes")
                cout<<"Include ICs in cage cluster detection: "<<cage.cl<<endl;
        else{
                cerr <<"\033[1;31m"<< "Error: -cl option should be yes or no!"<< "\033[0m" << endl;
                return 0;
        }
	if(cage.cal_additive){
		infile.open(hbond_filename.c_str());
		if (infile){	
			cage.ifs_hbond.open(hbond_filename.c_str());
			cout<<endl<<hbond_filename<<" opened successfully!"<<endl;
			if(cage.input_hbond())
				cout<<"Build Hbond successfully!"<<endl;
			else
				return 0;
		}
		else{
			cerr<<"\033[1;31m"<<"No hbond input file found!"<< "\033[0m"<<endl;
			return 0;	
		}
		cage.ifs_hbond.close();
	}
	infile.close();
	
	if(cage.max_ring<6 || cage.max_ring>12){
		cerr<<"\033[1;31m"<<"-mr The maximum size of n-membered rings to be considered in the calculations (valid input range: 6–12)."<< "\033[0m"<<endl;	
		return 0;
	}
	if(cage.IA_tolerance>30){
		cerr<<"\033[1;31m"<<"Error: IA tolerance must not exceed 30°."<< "\033[0m"<<endl;
		cage.IA_tolerance=30;
	}
	cout<<"\nHbond criteria:\n  rcut:              "<<cage.rcut<<" nm\n  angle threshold:   "<<cage.theta<<"°"<<endl;
        cout<<"  DA tolerance:      "<<cage.DA_tolerance<<"°"<<endl;
	cout<<"  IA tolerance:      " << cage.IA_tolerance << "° (applied to rings with size >= 7)" << endl;
        cout<<"  Max ring size:     "<<cage.max_ring<<endl;
        if(cage.theta<0)
                cout << "Skipping hydrogen bond angle consideration for this calculation." << endl;
	if(cage.cal_box){
		istringstream iss(nbox);
		iss>>cage.fx>>cage.fy>>cage.fz;
		if(cage.fx<1 || cage.fy<1 || cage.fz<1){
			cerr<<"-nbox inputs should be at least than 1!"<<endl;
			cage.grid_flag=true;
			//return 0;
		}
		if(cage.grid_flag)
			cout << "\033nBox: "<<cage.fx<<"X"<<cage.fy<<"X"<<cage.fz<<"\033[0m" << endl;
		else
			cout<<"Automatically calculating grid bounds"<<endl;
	}
	cage.init();
	cout << "Using " << omp_get_max_threads() << " thread(s)." << endl;
	//==================calculate and output=======================
		//open all file
	cage.ifs_water.open(H2O_filename.c_str());
	cage.ofs_summary.open(fname_summary.c_str());
	cage.ofs_detail.open(fname_detail.c_str());
	cage.ofs_visual_gro.open(fname_visual_gro.c_str());
	cage.ofs_visual_index.open(fname_visual_index.c_str());
	cage.ofs_cluster.open(fname_cluster.c_str());
	cage.ofs_crystallinity.open(fname_crystallinity.c_str());
	cage.ofs_ring_count.open(fname_ring_count.c_str());
	cage.ofs_ring_detail.open(fname_ring_detail.c_str());
        if(cage.cal_additive){
                cage.ifs_add.open(add_filename.c_str());
        }
	if(cage.cal_guest){
		cage.ifs_guest.open(Guest_filename.c_str());
		cage.ofs_occupancy.open(fname_occupancy.c_str());
		cage.ofs_occupancy<<"#ALL represents the total number of SEC, non-SEC, and IC types; ALL_F is the number of those filled by guests.";
		cage.ofs_occupancy<<"\n#vac is the vacancy rate, defined as 1 - (ALL_F / ALL)."<<endl;
                cage.ofs_occupancy<<"#Frame  vac     ALL     ALL_F   "<<endl;
        }
		//open all file
        cage.ofs_summary<<"#Cage type format ==> example: 4(1)5(10)6(2)==> 1,10,2 (Shows only some common SECs)"<<endl;
	cage.ofs_summary<<"#The others column includes SEC types that are not individually listed above."<<endl;
        cage.ofs_summary<<"#Standard edge-saturated cages  (SECs) nonstandard edge-saturated cages (non-SECs) Incomplete Cage (IC) others(SECs) "<<endl;
        cage.ofs_summary<<"#Frame  cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    cage    "<<endl;
        cage.ofs_summary<<"#       SEC     non-SEC IC      others  0,12,0  0,12,2  0,12,3  0,12,4  0,12,5  0,12,6  0,12,8  1,10,2  1,10,3  1,10,4  1,10,5  2,8,1   2,8,2   2,8,3   2,8,4   3,6,3   3,6,4  "<<endl;
        cage.ofs_cluster<<"#Frame cluster 1 (size) ... cluster n (size)"<<endl;
	cage.ofs_crystallinity<<"#frame 1:H2O ... Nw:H2O ... Na:add  Crystallinity (The number of cages each molecule participates)"<<endl;
        cage.ofs_ring_count<<"#frame 1:H2O ... Nw:H2O ... Na:add (The number of 4rings,5rings,6rings each molecule participates)"<<endl;
        cout<<"Start the calculation ..."<<endl;
	bool end=true;
	int i=0;//frame index
	//shift frame
	i+=shift_frame;cage.initial_frame+=shift_frame;cage.end_frame+=shift_frame;
	timespec start, eend;
	clock_gettime(CLOCK_MONOTONIC, &start);
	while(1){
		if(i>cage.end_frame)break;
		end=cage.input_H2O(cage.H2O_atom,cage.H2O_molecule,cage.H2O);
		if(end)break;
		if(cage.cal_additive){
        		cage.input_add(cage.add_atom,cage.add_molecule,cage.add);
		}
        	if(cage.cal_guest)
                	cage.input_guest(cage.guest_atom,cage.guest_molecule,cage.guest);
		if(i>=cage.initial_frame){
			cage.build_grid();
			cage.cal_and_output(i);
		}
		else{
                        //cout<<"Skip frame "<<i<<"..."<<endl;
                        i++;
                        continue;
                }	
		i++;
	}
	clock_gettime(CLOCK_MONOTONIC, &eend);
        cout << "Time consumed: " << get_elapsed_time(start, eend) << " sec" << endl;
	cage.ifs_water.close();
        cage.ifs_guest.close();
        cage.ifs_add.close();
        cage.ofs_summary.close();
        cage.ofs_detail.close();
        cage.ofs_visual_gro.close();
	cage.ofs_cluster.close();
	cage.ofs_visual_index.close();
        cage.ofs_crystallinity.close();
        cage.ofs_ring_count.close();
	cage.ofs_ring_detail.close();
	if(cage.cal_guest)
                cage.ofs_occupancy.close();
}
//=============================function of calculation==================================================
void cage::build_grid(){
	//timespec start, end;
        //clock_gettime(CLOCK_MONOTONIC, &start);
	if (grid_flag){
		if(box.x/fx<rcut || box.y/fy<rcut || box.z/fz<rcut){
			cerr << "\033[31mGrid is too small (<rcut)! Automatically calculating grid bounds.\033[0m" << endl;
			grid_flag=false;
                }	
	}
	if (!grid_flag){
		fx=round(box.x/1.2);
		fy=round(box.y/1.2);
		fz=round(box.z/1.2);
		if(fx<1)
			fx=1;
		if(fy<1)
			fy=1;
		if(fz<1)
			fz=1;
	}
	//cout<<box.x<<" "<<box.y<<" "<<box.z<<endl;
	//cout<<fx<<" "<<fy<<" "<<fz<<endl;
	grid.resize(fx);
	for(size_t i = 0; i < grid.size(); ++i) {
		grid[i].resize(fy);
		for(size_t j = 0; j < grid[i].size(); ++j) {
			grid[i][j].resize(fz);
			for(size_t k = 0; k < fz; ++k) {
				grid[i][j][k].clear();
			}
		}
	}
	if(fx==1 && fy==1 && fz==1){
		for(int i=0;i<H2O+add;i++){
			grid[0][0][0].push_back(i);	
		}
		return;
	}
	double dx=box.x/fx;
        double dy=box.y/fy;
        double dz=box.z/fz;
	for(int i=0;i<H2O;i++){
		double x = H2O_molecule[i][0].x;
		double y = H2O_molecule[i][0].y;
		double z = H2O_molecule[i][0].z;
		if (x < 0) x += box.x;if (y < 0) y += box.y;if (z < 0) z += box.z;
		if (x >= box.x) x -= box.x;if (y >= box.y) y -= box.y;if (z >= box.z) z -= box.z;
		int xi = int(x / dx);int yi = int(y / dy);int zi = int(z / dz);
		if (xi == fx) xi = fx-1;if (yi == fy) yi = fy-1;if (zi == fz) zi = fz-1;
		double left_x_boundary = xi * dx;
		double left_y_boundary = yi * dy;
		double left_z_boundary = zi * dz;
		bool x_bound = (x <= left_x_boundary + rcut && xi > 0);
		bool y_bound = (y <= left_y_boundary + rcut && yi > 0);
		bool z_bound = (z <= left_z_boundary + rcut && zi > 0);
		bool x_pbc = (x <= rcut);
		bool y_pbc = (y <= rcut);
		bool z_pbc = (z <= rcut);
		for(int offset_x = 0; offset_x >= (x_bound ? -1 : 0); offset_x--) {
 	   		for(int offset_y = 0; offset_y >= (y_bound ? -1 : 0); offset_y--) {
        			for(int offset_z = 0; offset_z >= (z_bound ? -1 : 0); offset_z--) {
            				int nx = xi + offset_x;
            				int ny = yi + offset_y;
            				int nz = zi + offset_z;
           				grid[nx][ny][nz].push_back(i);
					if(y_pbc && z_pbc)
                                                grid[nx][fy-1][fz-1].push_back(i);
                                        if(x_pbc && z_pbc)
                                                grid[fx-1][ny][fz-1].push_back(i);
                                        if(x_pbc && y_pbc)
                                                grid[fx-1][fy-1][nz].push_back(i);
					if (x_pbc)
                                                grid[fx - 1][ny][nz].push_back(i);
                                        if (y_pbc)
                                                grid[nx][fy - 1][nz].push_back(i);
                                        if (z_pbc)
                                                grid[nx][ny][fz - 1].push_back(i);
					
        			}
    			}
		}
		if(x_pbc && y_pbc && z_pbc)
			grid[fx - 1][fy-1][fz-1].push_back(i);
	}
	//clock_gettime(CLOCK_MONOTONIC, &end);
        //cout << "Build Grid: " << get_elapsed_time(start, end) << " sec" << endl;
	if(!cal_additive)
		return;
	int *used = new int[fx*fy*fz];
	fill(used, used + fx*fy*fz, 0);
	int current=0;
	for(int i=0;i<add;i++){
		current++;
		for(int j=0;j<add_atoms.size();j++){
			int a=add_atoms[j];
			double x = add_molecule[i][a].x;
                	double y = add_molecule[i][a].y;
                	double z = add_molecule[i][a].z;
			if (x < 0) x += box.x;if (y < 0) y += box.y;if (z < 0) z += box.z;
                	if (x >= box.x) x -= box.x;if (y >= box.y) y -= box.y;if (z >= box.z) z -= box.z;
                	int xi = int(x / dx);int yi = int(y / dy);int zi = int(z / dz);
                	if (xi == fx) xi = fx-1;if (yi == fy) yi = fy-1;if (zi == fz) zi = fz-1;
			double left_x_boundary = xi * dx;
                	double left_y_boundary = yi * dy;
                	double left_z_boundary = zi * dz;
                	bool x_bound = (x <= left_x_boundary + rcut && xi > 0);
                	bool y_bound = (y <= left_y_boundary + rcut && yi > 0);
                	bool z_bound = (z <= left_z_boundary + rcut && zi > 0);
                	bool x_pbc = (x <= rcut);
                	bool y_pbc = (y <= rcut);
                	bool z_pbc = (z <= rcut);
			int index;
                	for(int offset_x = 0; offset_x >= (x_bound ? -1 : 0); offset_x--) {
                        	for(int offset_y = 0; offset_y >= (y_bound ? -1 : 0); offset_y--) {
                                	for(int offset_z = 0; offset_z >= (z_bound ? -1 : 0); offset_z--) {
						int nx = xi + offset_x;
						int ny = yi + offset_y;
                                        	int nz = zi + offset_z;
						index = nx*(fy*fz)+ny*fz+nz;
						if(used[index]!=current){
							used[index]=current;
                                        		grid[nx][ny][nz].push_back(i+H2O);
						}
						index = nx*(fy*fz)+(fy-1)*fz+(fz-1);
                                        	if(y_pbc && z_pbc && used[index]!=current){
							used[index]=current;
							grid[nx][fy-1][fz-1].push_back(i+H2O);
						}
						index = (fx-1)*(fy*fz)+ny*fz+(fz-1);
						if(x_pbc && z_pbc && used[index]!=current){
							used[index]=current;
                                                	grid[fx-1][ny][fz-1].push_back(i+H2O);
						}
						index = (fx-1)*(fy*fz)+(fy-1)*fz+nz;
                                        	if(x_pbc && y_pbc && used[index]!=current){
							used[index]=current;
                                                	grid[fx-1][fy-1][nz].push_back(i+H2O);
						}
						index = (fx-1)*(fy*fz)+ny*fz+nz;
                                        	if (x_pbc && used[index]!=current){
							used[index]=current;
                                                	grid[fx - 1][ny][nz].push_back(i+H2O);
						}
						index = nx*(fy*fz)+(fy-1)*fz+nz;
                                        	if (y_pbc && used[index]!=current){
							used[index]=current;
                                                	grid[nx][fy - 1][nz].push_back(i+H2O);
						}
						index = nx*(fy*fz)+ny*fz+(fz-1);
                                        	if (z_pbc && used[index]!=current){
							used[index]=current;
                                                	grid[nx][ny][fz - 1].push_back(i+H2O);
						}
                                	}
                        	}
                	}
			index = (fx-1)*(fy*fz)+(fy-1)*fz+(fz-1);
                	if(x_pbc && y_pbc && z_pbc && used[index]!=current){
				used[index]=current;
				grid[fx - 1][fy-1][fz-1].push_back(i+H2O);	
			}
		}
	}
	delete[] used;
	/*for (int xi = 0; xi < fx; ++xi) {
		for (int yi = 0; yi < fy; ++yi) {
			for (int zi = 0; zi < fz; ++zi) {
				if (!grid[xi][yi][zi].empty()) {
					cout << "Grid[" << xi << "][" << yi << "][" << zi << "]: ";
					for (size_t k = 0; k < grid[xi][yi][zi].size(); ++k) {
						cout << grid[xi][yi][zi][k] << " ";
					}
					cout << endl;
				}
			}
		}
	}*/
}
point cage::cross_product(point v1, point v2) {
	point n;
	n.x = v1.y * v2.z - v1.z * v2.y;
	n.y = v1.z * v2.x - v1.x * v2.z;
	n.z = v1.x * v2.y - v1.y * v2.x;
	return n;
}
double cage::dihedral(point n1, point n2) {
	double mag_n1 = sqrt(n1.x * n1.x + n1.y * n1.y + n1.z * n1.z);
	double mag_n2 = sqrt(n2.x * n2.x + n2.y * n2.y + n2.z * n2.z);
	double dot = n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
	double cos_phi = dot / (mag_n1 * mag_n2);
	double di_angle=acos(cos_phi)*180/M_PI;
	if(di_angle>90)
                di_angle=180-di_angle;
	return di_angle;
}
point cage::PBC_vector(point a, point b) {
	double dx, dy, dz;
	dx = box_2.x  - a.x;
	dy = box_2.y  - a.y;
	dz = box_2.z  - a.z;
	a.x = box_2.x ;
	a.y = box_2.y ;
	a.z = box_2.z ;
	b.x = b.x + dx;
	b.y = b.y + dy;
	b.z = b.z + dz;
	if (b.x > box.x)
		b.x -= box.x;
	else if (b.x < 0)
		b.x += box.x;
	if (b.y > box.y)
		b.y -= box.y;
	else if (b.y < 0)
		b.y += box.y;
	if (b.z > box.z)
		b.z -= box.z;
	else if (b.z < 0)
		b.z += box.z;
	point v;
	v.x = b.x - a.x;
	v.y = b.y - a.y;
	v.z = b.z - a.z;
	return v;
}
point cage::PBC_average(const vector<point>& p) {
	if (p.empty()) return point();
	point ref = p[0];
	point sum;
	sum.x=0;sum.y=0;sum.z=0;
	for (size_t i = 1; i < p.size(); ++i) {
		point delta;
		delta.x = p[i].x - ref.x;
		if (delta.x >  box_2.x) delta.x -= box.x;
		if (delta.x < -box_2.x) delta.x += box.x;
		delta.y = p[i].y - ref.y;
		if (delta.y >  box_2.y) delta.y -= box.y;
		if (delta.y < -box_2.y) delta.y += box.y;
		delta.z = p[i].z - ref.z;
		if (delta.z >  box_2.z) delta.z -= box.z;
		if (delta.z < -box_2.z) delta.z += box.z;
		sum.x += delta.x;
		sum.y += delta.y;
		sum.z += delta.z;
	}
	point avg;
	avg.x = ref.x + sum.x / p.size();
	avg.y = ref.y + sum.y / p.size();
	avg.z = ref.z + sum.z / p.size();
	if (avg.x > box.x) avg.x -= box.x;
	else if (avg.x < 0) avg.x += box.x;
	if (avg.y > box.y) avg.y -= box.y;
	else if (avg.y < 0) avg.y += box.y;
	if (avg.z > box.z) avg.z -= box.z;
	else if (avg.z < 0) avg.z += box.z;
	return avg;
}
point cage::PBC_average_P(const polyhedron& P) {
	const vector<int>& v=P.vertex;
	const vector< pair<int, int> >& e = P.edge;
	point ref,sum;
	int n=v.size();
	for (int i = 0; i < n; i++) {
		point p;
		if(v[i]>=H2O){
			vector<point> p_vec;
			int esize=e.size();
			set<int> hbond_atoms;
			for(int j=0;j<esize;j++){
				const int& e_first = e[j].first;
				const int& e_second = e[j].second;
				if(e_first < v[i] && e_second < v[i])
					break;
				const int& map1=v[i];
				int map2;
				if(e_first == v[i])
					map2=e_second;
				else if(e_second == v[i])
					map2=e_first;
				else
					continue;;
				hbond_atoms.insert(adj_map[make_pair(map1,map2)]-1);
			}
			for (set<int>::const_iterator it = hbond_atoms.begin(); it != hbond_atoms.end(); ++it) {
				int idx = *it;
				p_vec.push_back(add_molecule[v[i]-H2O][idx]);
			}
			p=PBC_average(p_vec);
		}
		else
			p = H2O_molecule[v[i]][0];
		if(i==0){
			ref=p;
			sum.x = ref.x;sum.y = ref.y;sum.z = ref.z;
			continue;
		}
		point delta;
        	delta.x = p.x - ref.x;
        	if (delta.x >  box_2.x ) delta.x -= box.x;
        	if (delta.x < -box_2.x ) delta.x += box.x;
        	delta.y = p.y - ref.y;
        	if (delta.y >  box_2.y ) delta.y -= box.y;
        	if (delta.y < -box_2.y ) delta.y += box.y;
        	delta.z = p.z - ref.z;
        	if (delta.z >  box_2.z ) delta.z -= box.z;
        	if (delta.z < -box_2.z ) delta.z += box.z;
        	point image;
        	image.x = ref.x + delta.x;
        	image.y = ref.y + delta.y;
        	image.z = ref.z + delta.z;
        	sum.x += image.x;
        	sum.y += image.y;
        	sum.z += image.z;
	}
    	point avg;
    	avg.x = sum.x / n;
    	avg.y = sum.y / n;
    	avg.z = sum.z / n;
	if (avg.x > box.x) avg.x -= box.x;
        else if (avg.x <  0)     avg.x += box.x;
        if (avg.y > box.y) avg.y -= box.y;
        else if (avg.y <  0)     avg.y += box.y;
        if (avg.z > box.z) avg.z -= box.z;
        else if (avg.z <  0)     avg.z += box.z;
    	return avg;	
}
//find ring==>
//only consider 4 5 6 7 8ring
pair<bool,double> cage::check_planarity(int& level,ring& ring,point n1){
	point a,b,c;
	int a_idx,b_idx,c_idx;
        a_idx=ring.vertex[level-2];
        b_idx=ring.vertex[level-1];
        c_idx=ring.vertex[level-0];
        if(a_idx<H2O)
                a=H2O_molecule[a_idx][0];
        else
		a = add_molecule[a_idx - H2O][adj_map[make_pair(a_idx, b_idx)]-1];
        if(b_idx<H2O)
                b=H2O_molecule[b_idx][0];
        else{
		point b1=add_molecule[b_idx-H2O][adj_map[make_pair(b_idx, a_idx)]-1];
                point b2=add_molecule[b_idx-H2O][adj_map[make_pair(b_idx, c_idx)]-1];
		vector<point> b_vec;
		b_vec.push_back(b1);
		b_vec.push_back(b2);
		b=PBC_average(b_vec);
        }
        if(c_idx<H2O)
                c=H2O_molecule[c_idx][0];
        else
		c=add_molecule[c_idx-H2O][adj_map[make_pair(c_idx, b_idx)]-1];
        point v1=PBC_vector(a,b);
        point v2=PBC_vector(b,c);
        point n2=cross_product(v1,v2);
	double di_angle=dihedral(n1,n2);
	if(di_angle>DA_tolerance){
		pair<bool,double> result=make_pair(true,di_angle);
		return result;
	}
	else{
		pair<bool,double> result=make_pair(false,di_angle);
                return result;
	}
}
bool cage::check_interior_angle(int& level,ring& ring){
	int sides=level+2;
        double total_IA = (sides - 2) * 180.0;
        double IA = total_IA / sides;
	size_t ringIAsize=ring.IA.size();
	for(int i=0;i<ringIAsize;i++){
		if(fabs(ring.IA[i]-IA)>IA_tolerance)
			return true;
	}
	return false;
}
bool cage::check_isring(int& level,ring& ring){
	//cal distort
	if(level<6)
		return false;
	vector<int> count(ring.vertex.size(), 0);
        for(int i=0;i<ring.vertex.size();i++){
                for(int j=i+1;j<ring.vertex.size();j++){
                        bool bond=false;
			if(adj_map.find(make_pair(ring.vertex[i], ring.vertex[j])) != adj_map.end())
				bond=true;
                        if(bond){
                                count[i]++;
                                count[j]++;
                                if(count[i]>2)
                                        return true;
                        }
                }
        }
	size_t ring_vertex_size=ring.vertex.size();
	double sum_x=0,sum_y=0,sum_z=0;
	for (int k = 0; k < level+1; ++k) {
		point a,b;
	        int a_idx=ring.vertex[k];
		int b_idx=ring.vertex[(k+1)%ring_vertex_size];
		if(a_idx<H2O)
			a=H2O_molecule[a_idx][0];
		else
			a=add_molecule[a_idx-H2O][0];
		if(b_idx<H2O)
			b=H2O_molecule[b_idx][0];
		else
			b=add_molecule[b_idx-H2O][0];
		point v;
		v=PBC_vector(a,b);
		sum_x+=v.x;sum_y+=v.y;sum_z+=v.z;
                
	}
	if(fabs(sum_x)>0.1)
		return true;
	if(fabs(sum_y)>0.1)
		return true;
	if(fabs(sum_z)>0.1)
		return true;
	//cout<<sum_x<<" "<<sum_y<<" "<<sum_z<<endl;	
	return false;
}
void cage::ring_dfs(int istart,int& ispath,unordered_set<vector<int>, VectorHash>& path,ring& local_ring,vector<int>& visited,point n1,int level,int max_level,vector<ring>& local_rings){
	//i is current j is neighbor
	visited[istart] = 1;
        local_ring.vertex.push_back(istart);
	//print1D(local_ring.vertex);	
	//Cal The Interior Angle
	point a,b,c;
	double angle = 0;
	if(level>=2){
		int a_idx,b_idx,c_idx;
                a_idx=local_ring.vertex[level-2];
                b_idx=local_ring.vertex[level-1];
                c_idx=local_ring.vertex[level];
                if(a_idx<H2O)
                        a=H2O_molecule[a_idx][0];
                else
			a=add_molecule[a_idx - H2O][adj_map[make_pair(a_idx, b_idx)]-1];
                if(b_idx<H2O)
                        b=H2O_molecule[b_idx][0];
                else{
			point b1=add_molecule[b_idx-H2O][adj_map[make_pair(b_idx, a_idx)]-1];
                	point b2=add_molecule[b_idx-H2O][adj_map[make_pair(b_idx, c_idx)]-1];
			vector<point> b_vec;
			b_vec.push_back(b1);
                	b_vec.push_back(b2);
                	b=PBC_average(b_vec);
                }
                if(c_idx<H2O)
                        c=H2O_molecule[c_idx][0];
                else
			c=add_molecule[c_idx-H2O][adj_map[make_pair(c_idx, b_idx)]-1];
		angle=PBC_angle(b,a,c);
		//cout<<angle<<endl;
	}
	pair<bool,double> distort;
	if(level==2){
		point v1=PBC_vector(a,b);
                point v2=PBC_vector(b,c);
                n1=cross_product(v1,v2);
		local_ring.IA.push_back(angle);
	}
	else if (level > 2) {;
		//cout<<distort.second<<endl;
		bool bond=false;
		distort = check_planarity(level, local_ring, n1);
		if (distort.first) {
                        local_ring.vertex.pop_back();
                        visited[istart] = 0;
                        return;
                }
		if (adj_map.find(make_pair(local_ring.vertex[0], istart)) != adj_map.end())
			bond = true;
		if (bond) {
			if (!check_isring(level, local_ring)){
				vector<int> tmp_path=local_ring.vertex;
				sort(tmp_path.begin(),tmp_path.end());
				if(path.insert(tmp_path).second){
					local_rings.push_back(local_ring);
				}
				ispath=1;
			}
			visited[istart] = 0;
			local_ring.vertex.pop_back();
			return;
		}
		else if (level>4 && check_interior_angle(level, local_ring)) {
			local_ring.vertex.pop_back();
			visited[istart] = 0;
			return;
		}
		local_ring.IA.push_back(angle);
	}
	if (level == max_level) {
		local_ring.vertex.pop_back();
		local_ring.IA.pop_back();
		visited[istart] = 0;
		return;
	}
	size_t H2O_adjlist_size=H2O_adjlist[istart].size();
	for(int j=0;j<H2O_adjlist_size;j++){;
                if (!visited[H2O_adjlist[istart][j]] && H2O_adjlist[istart][j] > local_ring.vertex[0]){
			if(!ispath)
				ring_dfs(H2O_adjlist[istart][j],ispath,path,local_ring,visited,n1,level+1,max_level,local_rings);
			if(ispath && level>1 ){
                		if(level==2)
                        		ispath=0;
                        	local_ring.IA.pop_back();
                		local_ring.vertex.pop_back();
				visited[istart] = 0;
                		return;
        		}
		}
        }
	local_ring.vertex.pop_back();
    	visited[istart] = 0;
	if(level>=2)
		local_ring.IA.pop_back();
}
void cage::find_ring(){
	#pragma omp parallel
	{	
		vector<int> visited(H2O+add,0);
        	vector<ring> local_rings;
		ring local_ring; 
		point n1;
		#pragma omp for schedule(dynamic, 1)
        	for(int i=0;i<H2O+add;i++){
			int ispath=0;
			unordered_set<vector<int>, VectorHash> path; 
                	ring_dfs(i,ispath,path,local_ring,visited,n1,0,max_ring,local_rings);
		}
		#pragma omp critical
			allring.insert(allring.end(), local_rings.begin(), local_rings.end());
	}
}
void cage::find_ring_group(const vector<ring >& ring_groups,unordered_map<pair<int, int>, vector<ring* >, pair_hash >& grouped){
        //grouped save the data of rings share the same edge
        for (size_t i = 0; i < ring_groups.size(); ++i) {
                size_t vertex_count = ring_groups[i].vertex.size();
                const ring& current_ring = ring_groups[i];
                //only n nmber of edge ex: i j k l m==>ij jk kl lm mi
                for (size_t j = 0; j < vertex_count; ++j) {
                        int a = current_ring.vertex[j];//Vertex a
                        int b = current_ring.vertex[(j + 1) % vertex_count];//Vertex b
                        grouped[make_pair(a < b ? a : b, a < b ? b : a)].push_back(const_cast<ring*>(&current_ring));
                }
        }
}
void cage::ring_edge(ring &ring){
	//print1D(ring.vertex);
	size_t n = ring.vertex.size();
	ring.edge.resize(n);
	ring.sort_edge.resize(n);
	for (size_t j = 0; j < ring.vertex.size(); ++j) {
		int a = ring.vertex[j];//Vertex a
		int b = ring.vertex[(j + 1) % ring.vertex.size()];//Vertex b
		int small = a < b ? a : b;
        	int large = a < b ? b : a;
		ring.edge[j] = make_pair(small, large);
		ring.sort_edge[j] = make_pair(small, large);
		//cout<<"( "<<edge.first+1<<" "<<edge.second+1<<" )";
	}
	sort(ring.sort_edge.begin(), ring.sort_edge.end());	
	//cout<<endl;
}
//find cup==>
void cage::find_cup_group(const vector<polyhedron*>& cup_groups,unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>& grouped){
	//cout<<"cup grouped"<<endl;
	for (size_t i = 0; i < cup_groups.size(); ++i) {
		polyhedron* cup = cup_groups[i];
		const vector<int>& key = cup->face[0];
		unordered_map<vector<int>, vector<polyhedron*>*, VectorHash>::iterator it = grouped.find(key);
		if (it == grouped.end())
			it = grouped.emplace(key, new vector<polyhedron*>()).first;
		it->second->push_back(cup);
	}
}
void cage::lateral_ring_combine(
	vector<vector<ring*> >& lateral_ring_grouped_vec,
	vector<ring*>& current_combination,
	vector<vector<ring*> >& all_combinations,
	size_t idx,bool& iscup) {
	if (idx == lateral_ring_grouped_vec.size()) {
		const vector<pair<int,int>>& pre_v = current_combination[0]->sort_edge;
        	const vector<pair<int,int>>& vn = current_combination[idx - 1]->sort_edge;
		if(share_edge(pre_v, vn)){
			all_combinations.push_back(current_combination);
			iscup=true;	
		}
		return;
	}
	vector<ring*>& vec = lateral_ring_grouped_vec[idx];
	for (size_t i = 0; i < vec.size(); ++i) {
		//print1D(vec[i]->vertex);
		if(idx!=0){
                        const vector<pair<int,int>>& pre_v = current_combination[idx - 1]->sort_edge;
                        if (!share_edge(pre_v, vec[i]->sort_edge))
                                continue;
                }
		if(!iscup){
                        current_combination[idx] = vec[i];
                        lateral_ring_combine(lateral_ring_grouped_vec, current_combination, all_combinations, idx + 1, iscup);
                }
		if(iscup && idx>0){
			if(idx==1)
				iscup=false;
			return;
		}
	}
}
//Incorrect cups are not excluded from our calculation.
void cage::find_cup(const vector<ring >& ring_groups,
                    const unordered_map<pair<int, int>, vector<ring* >, pair_hash >& grouped,
                    vector<polyhedron*>& cup_group) {
	#pragma omp parallel
	{
		vector<polyhedron*> local_cup_group;
		//timespec start, end;
		size_t ring_groups_size=ring_groups.size();
		#pragma omp for schedule(dynamic,1)
		for (int i = 0; i < ring_groups_size; ++i) {
			//clock_gettime(CLOCK_MONOTONIC, &start);
			const ring& iring = ring_groups[i];
			size_t iring_vertex_size=iring.vertex.size();
			unordered_map<pair<int, int>, vector<ring* >, pair_hash > lateral_ring_grouped;
			for (size_t j = 0; j < iring_vertex_size; ++j) {
				unordered_map<pair<int, int>, vector<ring* >, pair_hash >::const_iterator it = grouped.find(iring.edge[j]);
				if (it != grouped.end()) {
					const vector<ring*>& shared_rings = it->second;
					for (size_t k = 0; k < shared_rings.size(); ++k) {
						ring* shared_ring = shared_rings[k];
						size_t shared_ring_size=shared_ring->vertex.size();
						if (shared_ring->vertex != iring.vertex) {
							size_t shared_count = 0;
							for (size_t l = 0; l < iring_vertex_size; ++l){
								for (size_t m = 0; m < shared_ring_size; ++m) {
									if (iring.vertex[l] == shared_ring->vertex[m]) {
										++shared_count;
										if (shared_count > 2) break;
									}
								}
								if (shared_count > 2) break;
							}
							if (shared_count == 2) 
								lateral_ring_grouped[iring.edge[j]].push_back(shared_ring);

						}
					}
				} 
				else 
					break;
			}
			//clock_gettime(CLOCK_MONOTONIC, &end);
        		//cout << "Built lateral grouped: " << get_elapsed_time(start, end) << " sec" << endl;
			//clock_gettime(CLOCK_MONOTONIC, &start);
			if (lateral_ring_grouped.size() == iring_vertex_size) {
				//cout<<"iring: ";
				//print1D(iring.vertex);
				bool iscup=false;
				vector<ring*> current_combination(iring_vertex_size);
				vector<vector<ring*> > all_combinations;
				vector<vector<ring*> > lateral_ring_grouped_vec;
				for (size_t j = 0; j < iring_vertex_size; ++j) {
					pair<int, int> key = iring.edge[j];
    					const vector<ring*>& rings = lateral_ring_grouped[key];
    					vector<ring*> ptr_vec;
    					ptr_vec.reserve(rings.size());
					for (ring* rp : rings) 
						 ptr_vec.push_back(rp);
					lateral_ring_grouped_vec.push_back(move(ptr_vec));
				}
				lateral_ring_combine(lateral_ring_grouped_vec, current_combination, all_combinations, 0, iscup);
				int cup_count = 0;
				//cout<<all_combinations.size()<<endl;
				for (size_t c = 0; c < all_combinations.size(); ++c) {
					polyhedron* P = new polyhedron();
					P->face.push_back(iring.vertex);
					unordered_set<int> uni_vertice;
					size_t j_size=all_combinations[c].size();
					for (int j = 0; j < j_size; ++j) {
						size_t k_size=all_combinations[c][j]->vertex.size();
						P->face.push_back(all_combinations[c][j]->vertex);
						for (size_t k = 0; k < k_size; ++k)
							uni_vertice.insert(all_combinations[c][j]->vertex[k]);
					}
					P->vertex.assign(uni_vertice.begin(), uni_vertice.end());
					//print1D(P->vertex);
					local_cup_group.push_back(P);
					cup_count++;
				}
			}
			//clock_gettime(CLOCK_MONOTONIC, &end);
                        //cout << "All combination end: " << get_elapsed_time(start, end) << " sec" << endl;
		}
		#pragma omp critical
			cup_group.insert(cup_group.end(), local_cup_group.begin(), local_cup_group.end());
	}
}
bool cage::share_edge(const vector<pair<int,int> >& v1, const vector<pair<int,int> >& v2){
	size_t i = 0, j = 0;
	size_t n1 = v1.size(), n2 = v2.size();
	while (i < n1 && j < n2) {
    		if (v1[i] == v2[j])
			return true;
		else if (v1[i] < v2[j]) 
        		++i;	
		else 
        		++j;
	}
	return false;	
}
//find cage==>
pair<bool,vector<polyhedron*> > cage::containecup(const vector<int>& face1,const vector<int>& face2,unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>& grouped){
	vector<polyhedron*> matching_cups;
	unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>::const_iterator it = grouped.find(face1);
        if (it != grouped.end() && it->second != NULL && !it->second->empty() ) {
		for (size_t i = 0; i < it->second->size(); ++i) {
			polyhedron* cup = it->second->at(i);
			if(contain_face(face2,cup->face))
				matching_cups.push_back(cup);
		}
	}
	if (!matching_cups.empty()){
		//cout<<"match "<<matching_cups.size()<<endl;
        	return make_pair(true, matching_cups); 
	}
	else
        	return make_pair(false,vector<polyhedron*>());	
}
bool cage::contain_face(const vector<int>& vec, const vector<vector<int> >& vecs) {
	for (vector<vector<int> >::const_iterator it = vecs.begin(); it != vecs.end(); ++it) {
        	if (vec.size() == it->size() && vec[0] == (*it)[0] && equal(vec.begin() + 1, vec.end(), it->begin() + 1)) 
			return true;
    	}
    	return false;
}
string cage::intToString(int number) {
	ostringstream oss;
	oss << number;
	return oss.str();
}
string cage::mapToString(const map<int, int>& face_count) {
    string result;
    for (map<int, int>::const_iterator it = face_count.begin(); it != face_count.end();it++) 
        result += intToString(it->first) + "(" + intToString(it->second) + ")";
    //cout<<result<<" "<<result.length()<<endl;
    return result;
}
void cage::lateral_cup_combine(const unordered_map<vector<int>, vector<polyhedron*> , VectorHash>& lateral_cup_grouped,
     vector<polyhedron*>& current_combination,
     vector<vector<polyhedron*> >& all_combinations,
     unordered_map<vector<int>, vector<polyhedron*> , VectorHash>::const_iterator it,
     unordered_map<vector<int>, vector<polyhedron*> , VectorHash>::const_iterator end,int& type){
	if (it == end) {
		all_combinations.push_back(current_combination);
		return;
	}
	const vector<polyhedron*>& polyhedra = it->second;
	size_t psize=polyhedra.size();
	for (size_t i = 0; i < psize; ++i) {
		if(polyhedra[i]->found > 0)
			continue;
		current_combination.push_back(polyhedra[i]);
		unordered_map<vector<int>, vector<polyhedron*> , VectorHash>::const_iterator next_it = it;
		++next_it;
		lateral_cup_combine(lateral_cup_grouped, current_combination, all_combinations, next_it, end, type);
		current_combination.pop_back();
	}
}
void cage::find_cage_ICO(int layer, int max_layer,
    	unordered_set<vector<int>, VectorHash>& uni_face, unordered_set<int>& uni_vertice,
    	vector<vector<int> >& remaining_face,
    	vector<vector<int> >& pre_layer_face,
    	unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>& grouped, int& iscage,
	vector<polyhedron>& cage,int& type) {
	if (layer >= 3) {
		int total_vertex = 0;
		for (unordered_set<vector<int> >::const_iterator it = uni_face.begin(); it != uni_face.end(); ++it)
			total_vertex += it->size();
		if (uni_vertice.size() * 3 <= total_vertex + type ) {
			polyhedron P;
			P.face.assign(uni_face.begin(), uni_face.end());
			for (int f = 0; f < P.face.size(); f++) {
				P.total += P.face[f].size();
				for (int v = 0; v < P.face[f].size(); v++) {
					int a = P.face[f][v];
					int b = P.face[f][(v + 1) % P.face[f].size()];
					pair<int, int> edge = make_pair(min(a, b), max(a, b));
					P.vertex.push_back(a);
					P.edge.push_back(edge);
				}
			}
			//cout<<layer<<" Find cage!"<<endl;
			P.perfect=SEC_cage_identify(P);
			if(P.perfect>0){
                                iscage=1;
				cage.push_back(P);
			}
                        else if(P.perfect==0){
				iscage=2;
				if(type>0)
					cage.push_back(P);	
			}
			else
				iscage=2;
                        return;
		}
	}
	if (layer >= max_layer || remaining_face.empty()) {
		//cout << "reach max or empty! " << layer << " " << remaining_face.size() << endl;
		return;
	}
	if (remaining_face.size() > 12) {
		//cout << "too many remaining_face!" << endl;
		return;
	}
	unordered_map<vector<int>, vector<polyhedron*> , VectorHash> lateral_cup_grouped;
	// each face of pre_layer must have a cup to form a cage
	for (int j = 0; j < remaining_face.size(); j++) {
		pair<bool, vector<polyhedron*> > result = containecup(remaining_face[j], pre_layer_face[j], grouped);
		if (result.first)
			lateral_cup_grouped[remaining_face[j]] = result.second;
	}
	if (layer == 2 && lateral_cup_grouped.size() < 4) return;
	// find all possible combinations to form new layer
	int total_product = 1;
	for (unordered_map<vector<int>, vector<polyhedron*> , VectorHash>::iterator it = lateral_cup_grouped.begin(); it != lateral_cup_grouped.end(); ++it) {
		const vector<polyhedron*>& value = it->second;
		total_product *= value.size();
		if (total_product > 10) {
			//cout <<total_product <<"too many combinations!" << endl;
			return;
		}
	}
	//cout<<layer<<" "<<total_product<<endl;
	vector<polyhedron*> current_combination;
	vector<vector<polyhedron*> > all_combinations;
	// each face of pre_layer must have a cup to form map
	lateral_cup_combine(lateral_cup_grouped, current_combination, all_combinations, lateral_cup_grouped.begin(), lateral_cup_grouped.end(), type);
	unordered_map<vector<int>, vector<polyhedron*> , VectorHash>::const_iterator it = lateral_cup_grouped.begin();
	// the next layer pre_layer_face is the key of map
	vector<vector<int> > key_vec;
	for (it = lateral_cup_grouped.begin(); it != lateral_cup_grouped.end(); ++it)
		key_vec.push_back(it->first);
	//cout << "all_combinations: " << all_combinations.size() << endl;
	// [[cup1 cup2 cup3 ...] [cup1 cup2 cup3 ...] ...]
	size_t csize=all_combinations.size();
	for (size_t c = 0; c < csize; ++c) {
		//cout << "layer" << layer << " combination" << c << endl;
		vector<vector<int> > tmp_remaining_face, tmp_pre_layer_face;
		unordered_set<vector<int>, VectorHash> tmp_uni_face=uni_face;
		unordered_set<int> tmp_uni_vertice=uni_vertice;
		// [cup1 cup2 cup3 ...]
		//cout<<layer<<" "<<"combination: "<<c<<endl;
		size_t ksize=all_combinations[c].size();
		for (int k = 0; k < ksize; ++k) {
			vector<int> key = key_vec[k];
			size_t lsize=all_combinations[c][k]->face.size();
			//print1D(all_combinations[c][k]->vertex);
			//print2D(all_combinations[c][k]->face);
			for (int l = 0; l < lsize; l++) {
				if (uni_face.insert(all_combinations[c][k]->face[l]).second) {
					uni_vertice.insert(all_combinations[c][k]->vertex.begin(), all_combinations[c][k]->vertex.end());
					tmp_remaining_face.push_back(all_combinations[c][k]->face[l]);
					tmp_pre_layer_face.push_back(key);
				}
			}
		}
		find_cage_ICO(layer + 1, max_layer, uni_face, uni_vertice, tmp_remaining_face, tmp_pre_layer_face, grouped, iscage, cage, type);
		if (iscage) {
			if(iscage==1 || type>0){
				for (int k = 0; k < ksize; ++k)
					all_combinations[c][k]->found=1;
			}
			//cout << "layer" << layer << " find cage" << endl;
			return;
		}
		uni_face = tmp_uni_face;
		uni_vertice = tmp_uni_vertice;
	}
}
void cage::find_cage(vector<polyhedron*>& cup_group, vector<polyhedron>& cage, unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>& grouped,int& type) {
	//#pragma omp parallel 
	{
		//vector<polyhedron> local_cage;
		size_t cup_group_size=cup_group.size();
		//#pragma omp for schedule(static,1)
		//for (int i = 0; i < 1; i++) {
		for (int i = 0; i < cup_group_size; i++) {
			if (cup_group[i]->found>0) continue;
			vector<vector<int> > remaining_face, pre_layer_face;
			unordered_set<int> uni_vertice;
			unordered_set<vector<int>, VectorHash> uni_face;
			for (int j = 1; j < cup_group[i]->face.size(); j++) {
				remaining_face.push_back(cup_group[i]->face[j]);
				pre_layer_face.push_back(cup_group[i]->face[0]);
			}
			uni_vertice.insert(cup_group[i]->vertex.begin(), cup_group[i]->vertex.end());
			uni_face.insert(cup_group[i]->face.begin(), cup_group[i]->face.end());
			int layer = 2, max_layer = 4;
			int iscage = 0;
			//cout<<"Base cup"<<endl;
			//print1D(cup_group[i]->vertex);
			find_cage_ICO(layer, max_layer, uni_face, uni_vertice, remaining_face, pre_layer_face, grouped, iscage, cage,type);
		}
	}
}
double cage::PBC_dist(point a,point b){
	double dx,dy,dz;
	dx=fabs(a.x-b.x);
	dy=fabs(a.y-b.y);
	dz=fabs(a.z-b.z);
	if(dx>box_2.x)
		dx=box.x-dx;
	if(dy>box_2.y)
                dy=box.y-dy;
	if(dz>box_2.z)
                dz=box.z-dz;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
//other==>
double cage::PBC_angle(point a,point b,point c){
	//a O1(donor) b O2(aceptor) c H1
        double dx,dy,dz;
        dx=box_2.x-a.x;dy=box_2.y-a.y;dz=box_2.z-a.z;    
        a.x=box_2.x;a.y=box_2.y;a.z=box_2.z;
        b.x=b.x+dx;b.y=b.y+dy;b.z=b.z+dz;
        if(b.x>box.x)
                b.x-=box.x;
        else if(b.x<0)
                b.x+=box.x;
        if(b.y>box.y)
                b.y-=box.y;
        else if(b.y<0)
                b.y+=box.y;
        if(b.z>box.z)
                b.z-=box.z;
        else if(b.z<0)
                b.z+=box.z;
        c.x=c.x+dx;c.y=c.y+dy;c.z=c.z+dz;
        if(c.x>box.x)
                c.x-=box.x;
        else if(c.x<0)
                c.x+=box.x;
        if(c.y>box.y)
                c.y-=box.y;
        else if(c.y<0)
                c.y+=box.y;
        if(c.z>box.z)
                c.z-=box.z;
        else if(c.z<0)
                c.z+=box.z;
	//a O1(donor) b O2(aceptor) c H1
        point v1,v2;//vector v1 O1-O2, v2 O1-H1, 
        v2.x=c.x-a.x;v2.y=c.y-a.y;v2.z=c.z-a.z;
        v1.x=b.x-a.x;v1.y=b.y-a.y;v1.z=b.z-a.z;
	return acos((v1.x*v2.x + v1.y*v2.y + v1.z*v2.z) / sqrt( (v1.x*v1.x + v1.y*v1.y + v1.z*v1.z) * (v2.x*v2.x + v2.y*v2.y + v2.z*v2.z) ) )*(180.0 / M_PI);
}
//=============================other function=====================================================================
void cage::built_hbond_map() {
	double distance;
	bool Hbond;
        hbond hb;
	int id1,id2;
	for (size_t i = 0; i < H2O_adjlist.size(); ++i) 
                H2O_adjlist[i].clear();
	//#pragma omp parallel for schedule(static) private(distance, Hbond, id1, id2)
	for(int xi = 0; xi < fx; xi++) {
		for(int yi = 0; yi < fy; yi++) {
			for(int zi = 0; zi < fz; zi++) {
				vector<int> &cell = grid[xi][yi][zi];
				size_t cell_size=cell.size();
				//sort(grid[xi][yi][zi].begin(),grid[xi][yi][zi].end());
				for(int j = 0; j < cell_size; j++) {
					for(int k = j + 1; k < cell_size; k++) {
						id1 = cell[j];id2 = cell[k];
						if(id1<id2) {
							if(id1<H2O && id2<H2O){
								//water-water
								double distance = PBC_dist(H2O_molecule[id1][0], H2O_molecule[id2][0]);
								Hbond = false;
								if(distance <= rcut) {
									if(theta<0)
										Hbond = true;
									else	
										Hbond = PBC_angle(H2O_molecule[id1][0], H2O_molecule[id2][0], H2O_molecule[id1][1]) < theta ||
				 							PBC_angle(H2O_molecule[id1][0], H2O_molecule[id2][0], H2O_molecule[id1][2]) < theta ||
				 							PBC_angle(H2O_molecule[id2][0], H2O_molecule[id1][0], H2O_molecule[id2][1]) < theta ||
				 							PBC_angle(H2O_molecule[id2][0], H2O_molecule[id1][0], H2O_molecule[id2][2]) < theta;	
								}
								if (Hbond) {
									adj_map[make_pair(id1, id2)] = 1;
									adj_map[make_pair(id2, id1)] = 1;
									H2O_adjlist[id1].push_back(id2);
									H2O_adjlist[id2].push_back(id1);
								}
							}
							else if(id1<H2O && id2>=H2O){
								//water-additive
								for (int i = 0; i < hbond_H2O_add.size(); i++){
									hbond hb=hbond_H2O_add[i];
                                					if(hb.dtype=="w"){//donor w acceptor a
                                        					distance=PBC_dist(H2O_molecule[id1][hb.d],add_molecule[id2-H2O][hb.a]);
                                        					if(distance<=rcut){
                                                					if(theta<0 || PBC_angle(H2O_molecule[id1][hb.d],add_molecule[id2-H2O][hb.a],H2O_molecule[id1][hb.dh])<hb.th){
                                                        					adj_map[make_pair(id1, id2)] = 1;
                                                        					adj_map[make_pair(id2, id1)] = hb.a+1;
                                                        					H2O_adjlist[id1].push_back(id2);
                                                        					H2O_adjlist[id2].push_back(id1);
                                                        					break;
                                                					}
                                        					}
                                 					}
									else{//donor a acceptor w
										distance=PBC_dist(H2O_molecule[id1][hb.a],add_molecule[id2-H2O][hb.d]);
										if(distance<=rcut){
											if(theta<0 || PBC_angle(add_molecule[id2-H2O][hb.d],H2O_molecule[id1][hb.a],add_molecule[id2-H2O][hb.dh])<hb.th){
												adj_map[make_pair(id1, id2)] = 1;
												adj_map[make_pair(id2, id1)] = hb.d+1;
                                                        					H2O_adjlist[id1].push_back(id2);
                                                        					H2O_adjlist[id2].push_back(id1);
                                                        					break;
                                                					 }
                                       						}
									}
								}
							}
							else if(id1>=H2O && id2>=H2O){
								//additive-additive
								Hbond=false;	
								//donor ja acceptor ka
								for (int i = 0; i < hbond_add_add.size(); i++) {
									hb=hbond_add_add[i];
									distance=PBC_dist(add_molecule[id1-H2O][hb.d],add_molecule[id2-H2O][hb.a]);
									if(distance<=rcut){
										if(theta<0 || PBC_angle(add_molecule[id1-H2O][hb.d],add_molecule[id2-H2O][hb.a],add_molecule[id1-H2O][hb.dh])<hb.th){
											adj_map[make_pair(id1, id2)] = 1+hb.d;
                                                					adj_map[make_pair(id2, id1)] = 1+hb.a;
                                                					Hbond=true;
                                                					H2O_adjlist[id1].push_back(id2);
                                                					H2O_adjlist[id2].push_back(id1);
                                                					break;
										}
									}
								}
								if(Hbond)continue;
								//donor ka acceptor ja
								for (int i = 0; i < hbond_add_add.size(); i++) {
									hb=hbond_add_add[i];
                                					distance=PBC_dist(add_molecule[id2-H2O][hb.d],add_molecule[id1-H2O][hb.a]);
                                					if(distance<=rcut){
										if(theta<0 || PBC_angle(add_molecule[id2-H2O][hb.d],add_molecule[id1-H2O][hb.a],add_molecule[id2-H2O][hb.dh])<hb.th){
											adj_map[make_pair(id1, id2)] = 1+hb.a;
                                                					adj_map[make_pair(id2, id1)] = 1+hb.d;
                                                					Hbond=true;
                                                					H2O_adjlist[id1].push_back(id2);
                                               				 		H2O_adjlist[id2].push_back(id1);
                                                					break;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	for (size_t i = 0; i < H2O_adjlist.size(); i++) {
                sort(H2O_adjlist[i].begin(), H2O_adjlist[i].end());
                H2O_adjlist[i].erase(unique(H2O_adjlist[i].begin(), H2O_adjlist[i].end()), H2O_adjlist[i].end());
        }
	/*for (int i = 0; i < H2O_adjlist.size(); ++i) {
                cout<<i+1<<" [ ";
                for (int j = 0; j < H2O_adjlist[i].size(); ++j) {
                        cout << H2O_adjlist[i][j]+1 << ' ';
                }
                cout <<" ]" <<endl;
        }*/
}
pair<bool,int> cage::define_molecule(string *type, int *atom_number,string filename,int &molecule_num, bool additive) {
	int num;
	*atom_number = 1;
        ifs.open(filename.c_str());
	int file_format;
	pair<bool,int> result;
	//determine the format
	string line,element;
	getline(ifs, line);getline(ifs, line);getline(ifs, line);
	vector<string> elements;
	istringstream iss(line);
	while (iss >> element) {
		elements.push_back(element);
	}
	file_format=elements.size();
	ifs.seekg(0, ifstream::beg);
	if(file_format==9)
		cout<<".gro has velocity information."<<endl;
	else if(file_format==6)
                cout<<".gro has no velocity information."<<endl;
	else{
		cerr<<"Wrong .gro file format!"<<endl;
		result=make_pair(false,file_format);
		return result;
	}
	//define molecule type
	string current_molecule, next_molecule, junk,atom;
	getline(ifs, line);getline(ifs, line);//first two line in gro file
	ifs >> num >> *type; 
	ifs.seekg(0, ifstream::beg);//back to file head
	//define molecule atom number
	getline(ifs, line);getline(ifs, line);
	iss.clear(); 
	iss.str(line); 
	iss>>molecule_num;
	ifs >> current_molecule;
	ifs >> atom >> junk >> junk >> junk >> junk;
	if(file_format==9)
		ifs >> junk >> junk >> junk;//first molecule atom in gro file
	cout<<atom<<" ";
	if(additive)add_name.push_back(atom);
	while (1) {
		ifs >> next_molecule;
		if (next_molecule != current_molecule)
			break;
		ifs >> atom >> junk >> junk >> junk >> junk; 
		if(file_format==9)
			ifs>>junk >> junk >> junk;
		*atom_number= *atom_number+1;
		cout<<atom<<" ";
		if(additive)add_name.push_back(atom);
	}
	cout<<endl<<*type<<" "<<*atom_number<<" atom(s)"<<endl;
	ifs.seekg(0, ifstream::beg);
	//define molecule number
	molecule_num/=*atom_number;
        cout<<"Number of molecules: "<<molecule_num<<endl<<endl;
        ifs.seekg(0, ifstream::beg);
	ifs.close();
	result=make_pair(true,file_format);
	return result;
}
//=============================function of output and input====================================================
void cage::input_guest(int atom_number,point *molecule,int molecule_num) {
	string line, junk;
	//cout<<"Wait for reading guest.gro ..."<<endl;
	getline(ifs_guest, line);getline(ifs_guest, line);
	for (int j = 0; j < molecule_num; j++) {
		for (int k = 0; k < atom_number; k++) {
			if(k<1){
				//skip atom name and resid num
				ifs_guest >> junk;ifs_guest.seekg(12, ios::cur); 
				ifs_guest >> molecule[j].x >> molecule[j].y >> molecule[j].z;
				if(file_format2==9)				
					ifs_guest >> junk>>junk>>junk;
			}
			else{
				ifs_guest>>junk;ifs_guest.seekg(12, ios::cur);
				ifs_guest>>junk>>junk>>junk;
				if(file_format2==9)
					ifs_guest>>junk>>junk>>junk;
			}
		}
	}
	if (!(ifs_guest >> junk >> junk >> junk)) {
                cerr << "\033[31mError: Failed to read box vector. Possible malformed .gro(guest) file at end.\033[0m" << endl;
        }
	ifs_guest.ignore(numeric_limits<streamsize>::max(), '\n');
}
bool cage::input_H2O(int atom_number,point **molecule,int molecule_num) {
        string line, junk;
        getline(ifs_water, line);
	if (ifs_water.eof()) return true;
	getline(ifs_water, line);
	for (int j = 0; j < molecule_num; j++) {
		for (int k = 0; k < atom_number; k++) {
			if(k<3){
				ifs_water>>junk;ifs_water.seekg(12, ios::cur);
				ifs_water >> molecule[j][k].x>>molecule[j][k].y>>molecule[j][k].z;
				if(file_format1==9)
					ifs_water >> junk>>junk>>junk;
			}
			else{
				ifs_water>>junk;ifs_water.seekg(12, ios::cur);
				ifs_water>>junk>>junk>>junk;
				if(file_format1==9)
					ifs_water>>junk>>junk>>junk;
			}
		}		
	}
	if (!(ifs_water >> box.x >> box.y >> box.z)){
		cerr << "\033[31mError: Failed to read box vector. Possible malformed .gro(H2O) file at end.\033[0m" << endl;
		return true;
	}
	box_2.x=box.x/2;box_2.y=box.y/2;box_2.z=box.z/2;
	ifs_water.ignore(numeric_limits<streamsize>::max(), '\n');
	return false;
}
void cage::input_add(int atom_number,point **molecule,int molecule_num) {
        string line, junk;
	getline(ifs_add, line);getline(ifs_add, line);
	for (int j = 0; j < molecule_num; j++) {
		for (int k = 0; k < atom_number; k++) {
			ifs_add>>junk;ifs_add.seekg(12, ios::cur);
			ifs_add >> molecule[j][k].x>>molecule[j][k].y>>molecule[j][k].z;
			if(file_format3==9)
				ifs_add >> junk>>junk>>junk;
			//cout<<molecule[j][k].x<<" "<<molecule[j][k].y<<" "<<molecule[j][k].z<<endl;
		}
	}
	if (!(ifs_add >> junk >> junk >> junk)) {
		cerr << "\033[31mError: Failed to read box vector. Possible malformed .gro(additive) file at end.\033[0m" << endl;
	}
	ifs_add.ignore(numeric_limits<streamsize>::max(), '\n');
}
bool cage::input_hbond() {
	string line;
	while (getline(ifs_hbond, line)) {
		size_t pos = line.find_first_not_of(" \t");
		if (pos == string::npos || line[pos] != '!') continue;
		istringstream iss(line.substr(pos + 1));
		hbond hb;
		string d_str, a_str, dh_str;
		if (!(iss >> hb.dtype >> hb.atype >> d_str >> a_str >> dh_str >> hb.th)) {
			cerr << "Format error in line! " << line << "\n";
			cerr << "Example:! w a 0 1 1 35 "<<endl;
			return false;
		}
		//cout<<hb.dtype<<" "<<hb.atype<<" "<<d_str<<" "<<a_str<<" "<<dh_str<<" "<<hb.r<<" "<<hb.th<<endl;
		for (size_t i = 0; i < d_str.size(); ++i)
			if (!isdigit(d_str[i])) {
				cerr << "Invalid integer d in line: " << line << "\n";
				return false;
			}
		for (size_t i = 0; i < a_str.size(); ++i)
			if (!isdigit(a_str[i])) {
				cerr << "Invalid integer a in line: " << line << "\n";
				return false;
			}
		for (size_t i = 0; i < dh_str.size(); ++i)
			if (!isdigit(dh_str[i])) {
				cerr << "Invalid integer dh in line: " << line << "\n";
				return false;
			}
		hb.d = atoi(d_str.c_str());
		hb.a = atoi(a_str.c_str());
		hb.dh = atoi(dh_str.c_str());
		if ((hb.dtype != "w" && hb.dtype != "a") || (hb.atype != "w" && hb.atype != "a") || (hb.dtype == "w" && hb.atype == "w")) {
    			cerr << "Invalid dtype/atype in line! (a or w, but not both w) " << line << "\n";
    			return false;
		}
		/*if (hb.r <= 0 || hb.th <= 0) {
			cerr << "Invalid r or theta in line! (r and theta should > 0) " << line << "\n";
			return false;
		}*/
		//water atom index
		if(hb.dtype == "w"){
			if(hb.dh>2){
				cerr << "The H2O hydrogen index is out of range! (1-2)" << line << "\n";
				return false;
			}
			if(hb.d!=0){
				cerr << "The H2O oxygen index must be 1!" << line << "\n";
				return false;
			}
		}
		if(hb.atype == "w" && hb.a!=0){
			if(hb.a!=0){
                                cerr << "The H2O oxygen index must be 1!" << line << "\n";
                                return false;
                        }
		}
		//additive atom index
		if(hb.dtype == "a"){
			if(hb.d > (add_atom-1) || hb.dh > (add_atom-1)){
				cerr << "The additive atom index is out of range! (<"<<add_atom-1<<")"<< line << "\n";
                                return false;
			}
		}
		if(hb.atype == "a"){
			if(hb.a > (add_atom-1)){
                                cerr << "The additive atom index is out of range! (<"<<add_atom-1<<")"<< line << "\n";
                                return false;
                        }
		}
		if(hb.dtype=="a" && hb.atype=="a")
			hbond_add_add.push_back(hb);
		else
			hbond_H2O_add.push_back(hb);
	}
	cout<<"Water-additive Hbonds:"<<endl;
	for (int i = 0; i < hbond_H2O_add.size(); i++) {
		cout << hbond_H2O_add[i].dtype << " "
		 << hbond_H2O_add[i].atype << " "
		 << hbond_H2O_add[i].d << " "
		 << hbond_H2O_add[i].a << " "
		 << hbond_H2O_add[i].dh << " "
		 //<< hbond_H2O_add[i].r << " "
		 << hbond_H2O_add[i].th << "\n";
	}
	cout<<endl;
	cout<<"Additive-Additive Hbonds:"<<endl;
	for (int i = 0; i < hbond_add_add.size(); i++) {
		cout << hbond_add_add[i].dtype << " "
		 << hbond_add_add[i].atype << " "
		 << hbond_add_add[i].d << " "
		 << hbond_add_add[i].a << " "
		 << hbond_add_add[i].dh << " "
		 //<< hbond_add_add[i].r << " "
		 << hbond_add_add[i].th << "\n";
	}
	cout<<endl;
	cout<<"# of possible w-a interaction: "<<hbond_H2O_add.size()<<endl<<"# of possible a-a interaction: "<<hbond_add_add.size()<<endl;
        cout<<endl;
	hbond hb;
        for(int i=0;i<hbond_H2O_add.size();i++){
                hb=hbond_H2O_add[i];
                if(hb.dtype=="a")
                        add_atoms.push_back(hb.d);
                else
                        add_atoms.push_back(hb.a);
        }
        for(int i=0;i<hbond_add_add.size();i++){
                hb=hbond_add_add[i];
                add_atoms.push_back(hb.d);
                add_atoms.push_back(hb.a);
        }
	sort(add_atoms.begin(), add_atoms.end());
        vector<int>::iterator new_end = unique(add_atoms.begin(), add_atoms.end());
	add_atoms.erase(new_end, add_atoms.end());
	return true;
}
bool comparePolygons(const polyhedron& a, const polyhedron& b) {
    	if (a.face.size() != b.face.size())
                return a.face.size() < b.face.size();
	if (a.type != b.type)
        	return a.type < b.type;
	if (a.type.length() != b.type.length())
                return a.type.length() < b.type.length();
	if (a.vertex[0] != b.vertex[0])
                return a.vertex[0] < b.vertex[0];
	return false;
}
struct BySecondDesc {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) const {
		return a.second > b.second;
	}
};
int cage::SEC_cage_identify(polyhedron &P){
	vector<int>& vertices = P.vertex; 
	vector<pair<int, int> >& edges = P.edge; 
	unordered_map<int, int> vertex_count;
	unordered_map<pair<int, int>, int, pair_hash> edge_count;
	size_t vsize=vertices.size(),esize=edges.size();
	for (int i = 0; i < vsize; ++i) 
		vertex_count[vertices[i]]++;
	for (int i = 0; i < esize; ++i) 
	   	edge_count[edges[i]]++;
	sort(P.vertex.begin(), P.vertex.end());
	P.vertex.erase(unique(P.vertex.begin(), P.vertex.end()), P.vertex.end());
	sort(P.edge.begin(), P.edge.end());
	P.edge.erase(unique(P.edge.begin(), P.edge.end()), P.edge.end());
	sort(P.edge.begin(), P.edge.end(), BySecondDesc());	
	//euler law F-E+V=2
	//SEC return 2 NSEC return 1 IC return 0 not cage return -1
	if(P.face.size()-P.edge.size()+P.vertex.size()!=2)
		return -1;
	bool has_one_vertex = false;
	bool has_low_vertex = false;//vertex<3
	bool has_high_vertex = false;//vertex>3
	bool has_low_edge = false;//edge<2
	for (unordered_map<int, int>::iterator it = vertex_count.begin(); it != vertex_count.end(); ++it) {
		if (it->second == 1)
			has_one_vertex = true;
		if (it->second < 3)
			has_low_vertex = true;
		if (it->second > 3)
			has_high_vertex = true;
	}
	for (unordered_map<pair<int, int>, int, pair_hash>::iterator it = edge_count.begin(); it != edge_count.end(); ++it) {
		if (it->second < 2)
			has_low_edge = true;
	}
	if ( has_one_vertex || has_low_edge)
		return -1;
	if (has_low_vertex)
		return 0;
	if (has_high_vertex)
		return 1;
	return 2;
}
//find cage cluster
void cage::built_cage_matrix(vector<polyhedron>& ALLcage){
        size_t size = ALLcage.size();
	unordered_map< vector<int>, vector<int>, VectorHash > face_to_cage;
	for (int i = 0; i < (int)size; ++i) {
		for (size_t j = 0; j < ALLcage[i].face.size(); ++j) {
			vector<int> f = ALLcage[i].face[j];
			face_to_cage[f].push_back(i);
		}
	}
	unordered_map< vector<int>, vector<int>, VectorHash >::iterator it;
	for (it = face_to_cage.begin(); it != face_to_cage.end(); ++it) {
		vector<int> &cages = it->second;
		for (size_t i = 0; i < cages.size(); ++i) {
			for (size_t j = i + 1; j < cages.size(); ++j) {
				int a = cages[i], b = cages[j];
				cage_adjlist[a].push_back(b);
				cage_adjlist[b].push_back(a);
			}
		}
	}
}
void cage::cage_dfs(int i, vector<int>& cluster, vector<int>& visited) {
	vector<int> stk;
	stk.push_back(i);
	while (!stk.empty()) {
		int node = stk.back();
		stk.pop_back();
		if (visited[node])
			continue;
		visited[node] = 1;
		cluster.push_back(node);
		for (int j = 0; j < cage_adjlist[node].size(); ++j) {
			if (!visited[cage_adjlist[node][j]])
				stk.push_back(cage_adjlist[node][j]);
		}
	}
}
void cage::find_cluster(vector<vector<int> >& groups,const int& n) {
	vector<int> visited(n, 0);
	vector<int> cluster;
	for (int i = 0; i < n; ++i) {
		if (!visited[i]) {
			cluster.clear();
			cage_dfs(i, cluster, visited);
			groups.push_back(cluster);
		}
	}
}
struct RingALLCompare {
	bool operator()(const ring& a, const ring& b) const {
		if (a.vertex.size() != b.vertex.size())
			return a.vertex.size() < b.vertex.size();
		return a.vertex < b.vertex;
	}
};
struct PolyhedronAllCompare {
	bool operator()(const polyhedron* a, const polyhedron* b) const {
		if (a->vertex.size() != b->vertex.size())
			return a->vertex.size() < b->vertex.size();
		return a->vertex < b->vertex; 
	}
};
//output function
void cage::cal_and_output(int i){
	timespec start, end;
	cout<<"Current frame "<<i<<" ..."<<endl;	
	ofs_detail<<"Frame: "<<i<<endl;
	ofs_ring_detail<<"Frame: "<<i<<endl;
	//reset parameter
	cage512 = 0, cage51262 = 0, cage51263 = 0, cage51264 = 0, cage51265 = 0, cage51266=0, cage51268=0;
	cage4151062 = 0, cage4151063 = 0, cage4151064 = 0, cage4151065 = 0, cage425864 = 0, cage425863 = 0;
 	cage425862 = 0, cage425861 = 0, cage435663 = 0, cage435664 = 0, other=0;
	adj_map.clear();
	if(cal_guest)
                fill(isguest.begin(), isguest.end(), 0);
	fill(H2O_color.begin(), H2O_color.end() ,10);
	//==============================built H-bond map===============================
	built_hbond_map();
	//==============================find ring======================================
	find_ring();
	__gnu_parallel::sort(allring.begin(), allring.end(), RingALLCompare());
	#pragma omp parallel for schedule(static)
	for(int j=0;j<allring.size();j++)
		ring_edge(allring[j]);
	unordered_map<pair<int, int>, vector<ring* >, pair_hash > ring_grouped;
	find_ring_group(allring, ring_grouped);
	write_ring_summary(i);
	//==============================find cup=======================================
	vector<polyhedron*> cup_group;
	find_cup(allring,ring_grouped,cup_group);
	allring.clear();
	ring_grouped.clear();
	__gnu_parallel::sort(cup_group.begin(), cup_group.end(), PolyhedronAllCompare());
	unordered_map<vector<int>, vector<polyhedron*> *, VectorHash> cup_grouped;
	find_cup_group(cup_group,cup_grouped);
	ofs_detail<<"cup: "<<cup_group.size()<<endl;
	//==============================find cage======================================
	vector<int> cage_form(H2O + add, 0);
	int fill=0; 
	vector<polyhedron> cage,SEC_cage,IC_cage,non_SEC_cage;
	int type=0;//type=0 find complete cage, type=2 find IC (distort tolerance T+2)
	find_cage(cup_group,cage,cup_grouped,type);//Find SECs and non-SECs
	type=2;
	find_cage(cup_group,cage,cup_grouped,type);//Find ICs
	classify_and_process_cages(i, fill, cage, SEC_cage, non_SEC_cage, IC_cage, cage_form);		
	process_and_write_clusters(i,SEC_cage,non_SEC_cage,IC_cage);
	output_visual(i,SEC_cage,non_SEC_cage,IC_cage);
	int sec_size = SEC_cage.size(),non_sec_size = non_SEC_cage.size(),ic_size = IC_cage.size();
	write_cage_details(i,cage_form, SEC_cage, non_SEC_cage, IC_cage);
	write_cage_summary(i,fill, sec_size, non_sec_size, ic_size);
	//=================================release memory=============================
	for (size_t j = 0; j < cup_group.size(); ++j)
		delete cup_group[j];
	cup_group.clear();	
	unordered_map<vector<int>, vector<polyhedron*> *, VectorHash>::iterator it = cup_grouped.begin();
	for (it; it != cup_grouped.end(); ++it)
		delete it->second;
	//-*/
}
void cage::print1D(vector<int> v){
	cout<<"[ ";
        for(int i=0;i<v.size();i++)
                cout<<v[i]+1<<" ";
        cout<<"] "<<endl;
}
void cage::printP(polyhedron p,ofstream& ofs_summary){
	
	ofs_summary<<p.cluster<<" Type: "<<p.type<<" t: "<<p.total<<" F: "<<p.face.size()<<" E: "<<p.edge.size()<<" V:"<<p.vertex.size()<<endl;
	ofs_summary << fixed << setprecision(3);
	ofs_summary<<" CP [ "<<p.cp.x<<" "<<p.cp.y<<" "<<p.cp.z<<" ]"<<endl;
	ofs_summary<<" V [ ";
	for(int i=0;i<p.vertex.size();i++)
		ofs_summary<<p.vertex[i]+1<<" ";
	ofs_summary<<"] "<<endl;

	if(p.guest_idx.size()>0){
                ofs_summary<<" g [ ";
                for(int i=0;i<p.guest_idx.size();i++)
                        ofs_summary<<p.guest_idx[i]<<" ";
                ofs_summary<<"] "<<endl;
        }

	if(p.add_list.size()>0){
		sort(p.add_list.begin(), p.add_list.end(), [](const pair<int, double>& a, const pair<int, double>& b) {return a.second < b.second;});
		int limit= p.add_list.size() < 3 ? p.add_list.size() : 3;
                ofs_summary<<" a [ ";
                for (int i = 0; i < limit; i++)
                        ofs_summary << p.add_list[i].first<<": " <<fixed << setprecision(3) << p.add_list[i].second << " ";
                ofs_summary << "]" << endl;
        }
	//if(cal_guest || cal_additive)
	//	ofs_summary << cp << " [ " << fixed << setprecision(3) << P.cp.x << " " << P.cp.y << " " << P.cp.z << " ]" << endl; 
	
	/*ofs_summary<<"edge: [";
        for(int i=0;i<p.edge.size();i++)
                ofs_summary<<"("<<p.edge[i].first+1<<" "<<p.edge[i].second+1<<")";
        ofs_summary<<" ]"<<endl;*/
}
void cage::printP2(polyhedron P){
	cout<<"Polyhedron: ";
	cout<<"Type: "<<P.type<<" total: "<<P.total<<" Face: "<<P.face.size()<<" Edge: "<<P.edge.size()<<" Vertex:"<<P.vertex.size()<<endl;
        //sort(P.vertex.begin(),P.vertex.end());
        //sort(P.edge.begin(),P.edge.end());
        cout<<"face: "<<endl;
	print2D(P.face);
	cout<<"edge: "<<endl;
        cout<<"[ ";
        for(int i=0;i<P.edge.size();i++)
                cout<<"("<<P.edge[i].first+1<<" "<<P.edge[i].second+1<<")";
        cout<<" ]"<<endl;
	cout<<"vertex: "<<endl;
        print1D(P.vertex);
	cout<<endl;
}
void cage::print2D(vector<vector<int> > v){
	cout<<"[ ";
        for(int i=0;i<v.size();i++){
		cout<<"[ ";	
                for(int j=0;j<v[i].size();j++)
                        cout<<v[i][j]+1<<" ";
                cout<<"] "<<endl;
        }
	cout<<" ]"<<endl;
}
void cage::define_cage(polyhedron &P,int frame,int perfect,int &fill){
	//define the cage type
	map<int, int> face_counts;
	double d_from_cp;
	int face_size=P.face.size();
	for(int l=0;l<face_size;l++)
		face_counts[P.face[l].size()]++;
	P.type=mapToString(face_counts);
	//built the neighbor additive list for every cage
	P.cp=PBC_average_P(P);
	if(cal_guest){
		if (face_size > 16)
			d_from_cp = 0.35;
		else if (face_size > 10)
			d_from_cp = 0.2;
		else if (face_size > 8)
			d_from_cp = 0.1;
		else
			d_from_cp = 0.05;	
	}
	if(cal_additive){
		for(int i=0;i<add;i++){
			double min = numeric_limits<double>::max();
			for(int k=0;k<add_atom;k++){
				double d=PBC_dist(P.cp,add_molecule[i][k]);
				if(d<min)
					min=d;
			}
			P.add_list.push_back(make_pair(i+H2O+1,min));
		}
	}
        //Calculate the filled cages.
	if(perfect==2){
		double d;
		//set<int> unique_vertex(P.vertex.begin(), P.vertex.end());
		//=========Count of small cage==============
		bool filled=false;
		if(cal_guest){
			//cout<<"#type "<<P.type<<" "<<P.vertex.size()<<endl;
			//print1D(P.vertex);
			for(int l=0;l<guest;l++){
				if(isguest[l])continue;
				if(PBC_dist(P.cp,guest_molecule[l])<d_from_cp){
					filled=true;
					P.guest_idx.push_back(l+1);//guest resid
					isguest[l]=1;
					if(face_size<=16)
						break;
				}
			}
			if(filled)
				fill++;
		}		   
		if (P.type == "5(12)") {
			cage512++;
			P.color_tier=4;
		} else if (P.type == "5(12)6(2)") {
			cage51262++;
			P.color_tier=5;
		} else if (P.type == "4(1)5(10)6(2)") {
                        cage4151062++;
			P.color_tier=6;
                } else if (P.type == "5(12)6(3)") {
			cage51263++;
			P.color_tier=1;
		} else if (P.type == "5(12)6(4)") {
			cage51264++;
			P.color_tier=3;
		} else if (P.type == "5(12)6(5)") {
			cage51265++;
			P.color_tier=7;
		} else if (P.type == "5(12)6(6)") {
                        cage51266++;
			P.color_tier=7;
		} else if (P.type == "5(12)6(8)"){
			cage51268++;
			P.color_tier=2;
                } else if (P.type == "4(1)5(10)6(3)") {
			cage4151063++;
			P.color_tier=7;
		} else if (P.type == "4(1)5(10)6(4)") {
			cage4151064++;
			P.color_tier=7;
		} else if (P.type == "4(1)5(10)6(5)") {
			cage4151065++;
			P.color_tier=7;
		} else if (P.type == "4(2)5(8)6(1)") {
			cage425864++;
			P.color_tier=7;
		} else if (P.type == "4(2)5(8)6(2)") {
			cage425863++;
			P.color_tier=7;
		} else if (P.type == "4(2)5(8)6(3)") {
			cage425862++;
			P.color_tier=7;
		} else if (P.type == "4(2)5(8)6(4)") {
			cage425861++;
			P.color_tier=7;
		} else if (P.type == "4(3)5(6)6(3)") {
			cage435663++;
			P.color_tier=7;
		} else if (P.type == "4(3)5(6)6(4)") {
			cage435664++;
			P.color_tier=7;
		} else {
			other++;
			P.color_tier=7;
		}
	}
	else {
		if(P.perfect==1)
			P.color_tier=0;
		else
			P.color_tier=8;
                if(cal_guest){
                        for(int l=0;l<guest;l++){
				if(isguest[l])continue;
                                if(PBC_dist(P.cp,guest_molecule[l])<d_from_cp){
					fill++;
					P.guest_idx.push_back(l+1);
					isguest[l]=1;
					if(face_size<=16)
                                                break;
                                }
                        }
                }
        }
	/*cout<<P.type<<" F: "<<P.face.size()<<" V: "<<P.vertex.size()<<endl;	
	for(int v=0;v<P.vertex.size();v++){
		if(P.vertex.size()<H2O)
			cout<<PBC_dist(H2O_molecule[P.vertex[v]][0],P.cp)<<" ";
	}
	cout<<endl;*/
}
//some output function
void cage::write_ring_summary(int i){
	vector<int> ring4_form(H2O + add, 0);
	vector<int> ring5_form(H2O + add, 0);
	vector<int> ring6_form(H2O + add, 0);
	vector<int> ring_count(13, 0);
	vector<int> aring_count(13, 0);
	size_t allring_size = allring.size();
	for (int j = 0; j < allring_size; ++j) {
		
		ofs_ring_detail<<"[ ";
                for(int k=0; k < allring[j].vertex.size(); ++k)
                        ofs_ring_detail<<allring[j].vertex[k]+1<<" ";
                ofs_ring_detail<<"]"<<endl;

		int vsize = allring[j].vertex.size();
		ring_count[vsize]++;
		if (vsize == 4) {
			for (int k = 0; k < vsize; ++k)
				ring4_form[allring[j].vertex[k]]++;
		} else if (vsize == 5) {
			for (int k = 0; k < vsize; ++k)
				ring5_form[allring[j].vertex[k]]++;
		} else if (vsize == 6) {
			for (int k = 0; k < vsize; ++k)
				ring6_form[allring[j].vertex[k]]++;
		}
	}
	/*if (cal_additive) {
                for (size_t j = 0; j < allring.size(); ++j) {
                        ring tmpring = allring[j];
			int tmpring_size=tmpring.vertex.size();
			for(int k=0;k<tmpring_size;k++){
				if(tmpring.vertex[k]>=H2O){
                                	aring_count[tmpring.vertex.size()]++;
					break;
				}
			}
                }
        }*/
	ofs_ring_count << i << " ";
	for (int j = 0; j < H2O + add; ++j)
		ofs_ring_count << ring4_form[j] << "," << ring5_form[j] << "," << ring6_form[j] << " ";
	for (int j = 4; j <= max_ring + 1; ++j)
		ofs_detail << j << "r: " << ring_count[j] << endl;
		//ofs_detail << j << "r: " << ring_count[j] - aring_count[j] << endl;
	/*if (cal_additive) {
		for (int j = 4; j <= max_ring + 1; ++j)
                        ofs_detail << j << "ar: " << aring_count[j] << endl;
	}*/
}
void cage::classify_and_process_cages(int i,int& fill,vector<polyhedron> &cage,vector<polyhedron> &SEC_cage,vector<polyhedron> &non_SEC_cage,vector<polyhedron> &IC_cage,vector<int> &cage_form) {
	size_t cagesize = cage.size();
	for (size_t j = 0; j < cagesize; ++j) {
		if (cage[j].perfect == 2)
			SEC_cage.push_back(cage[j]);
		else if (cage[j].perfect == 0)
			IC_cage.push_back(cage[j]);
		else if (cage[j].perfect == 1)
			non_SEC_cage.push_back(cage[j]);
	}
	cage.clear();
	ofs_detail << "SEC: " << SEC_cage.size() << endl;
	ofs_detail << "NSEC: " << non_SEC_cage.size() << endl;
	ofs_detail << "IC: " << IC_cage.size() << endl;
	for (size_t j = 0; j < SEC_cage.size(); ++j) {
		define_cage(SEC_cage[j], i, SEC_cage[j].perfect, fill);
		for (size_t k = 0; k < SEC_cage[j].vertex.size(); ++k)
			cage_form[SEC_cage[j].vertex[k]]++;
		// printP2(SEC_cage[j]);
	}
	for (size_t j = 0; j < non_SEC_cage.size(); ++j) {
		define_cage(non_SEC_cage[j], i, non_SEC_cage[j].perfect, fill);
		for (size_t k = 0; k < non_SEC_cage[j].vertex.size(); ++k)
			cage_form[non_SEC_cage[j].vertex[k]]++;
		// printP2(non_SEC_cage[j]);
	}
	for (size_t j = 0; j < IC_cage.size(); ++j) {
		define_cage(IC_cage[j], i, IC_cage[j].perfect, fill);
		for (size_t k = 0; k < IC_cage[j].vertex.size(); ++k)
			cage_form[IC_cage[j].vertex[k]]++;
		// printP2(IC_cage[j]);
	}
	__gnu_parallel::sort(SEC_cage.begin(), SEC_cage.end(), comparePolygons);
	__gnu_parallel::sort(non_SEC_cage.begin(), non_SEC_cage.end(), comparePolygons);
	__gnu_parallel::sort(IC_cage.begin(), IC_cage.end(), comparePolygons);
}
void cage::process_and_write_clusters(int i,const vector<polyhedron> &SEC_cage,const vector<polyhedron> &non_SEC_cage,const vector<polyhedron> &IC_cage) {
	vector<polyhedron> ALLcage;
	ALLcage.insert(ALLcage.end(), SEC_cage.begin(), SEC_cage.end());
	ALLcage.insert(ALLcage.end(), non_SEC_cage.begin(), non_SEC_cage.end());
	if (cl == "yes")
		ALLcage.insert(ALLcage.end(), IC_cage.begin(), IC_cage.end());
	cage_adjlist.resize(ALLcage.size());
	built_cage_matrix(ALLcage);
	vector<vector<int>> groups;
	find_cluster(groups, ALLcage.size());
	ofs_cluster << i << " ";
	if (groups.size() == 0)
		ofs_cluster << 0;
	sort(groups.begin(), groups.end(),[](const vector<int> &a, const vector<int> &b) {return a.size() > b.size();});
	for (size_t j = 0; j < groups.size(); ++j) {
		ofs_cluster << groups[j].size() << " ";
		for (size_t k = 0; k < groups[j].size(); ++k) {
			if (SEC_cage.size() > 0 && groups[j][k] <= SEC_cage.size() - 1)
				const_cast<polyhedron&>(SEC_cage[groups[j][k]]).cluster = j + 1;
			else if (non_SEC_cage.size() > 0 && groups[j][k] <= SEC_cage.size() + non_SEC_cage.size() - 1)
				const_cast<polyhedron&>(non_SEC_cage[groups[j][k] - SEC_cage.size()]).cluster = j + 1;
			else
				const_cast<polyhedron&>(IC_cage[groups[j][k] - SEC_cage.size() - non_SEC_cage.size()]).cluster = j + 1;
		}
	}
	ofs_cluster << endl;
	for (size_t j = 0; j < cage_adjlist.size(); ++j)
		cage_adjlist[j].clear();
}
void cage::write_cage_details(const int& i,const vector<int> &cage_form,const vector<polyhedron> &SEC_cage,const vector<polyhedron> &non_SEC_cage,const vector<polyhedron> &IC_cage){
	if (SEC_cage.size() > 0)
		ofs_detail << "=========================SECs=========================" << endl;
	for (size_t j = 0; j < SEC_cage.size(); ++j) {
		if (SEC_cage[j].vertex.back() < H2O)
			ofs_detail << "#cage" << j + 1 << " ";
		else
			ofs_detail << "#a-cage" << j + 1 << " ";
		printP(SEC_cage[j], ofs_detail);
	}
	if (non_SEC_cage.size() > 0)
		ofs_detail << "=========================non-SECs=========================" << endl;
	for (size_t j = 0; j < non_SEC_cage.size(); ++j) {
		if (non_SEC_cage[j].vertex.back() < H2O)
			ofs_detail << "@cage" << j + 1 << " ";
		else
			ofs_detail << "@a-cage" << j + 1 << " ";
		printP(non_SEC_cage[j], ofs_detail);
	}
	if (IC_cage.size() > 0)
		ofs_detail << "=========================ICs=========================" << endl;
	for (size_t j = 0; j < IC_cage.size(); ++j) {
		if (IC_cage[j].vertex.back() < H2O)
			ofs_detail << "!cage" << j + 1 << " ";
		else
			ofs_detail << "!a-cage" << j + 1 << " ";
		printP(IC_cage[j], ofs_detail);
	}
	ofs_crystallinity << i << " ";
	int cage_form_sum = 0;
	for (size_t j = 0; j < cage_form.size(); ++j) {
		ofs_crystallinity << cage_form[j] << " ";
		cage_form_sum += cage_form[j];
	}
	ofs_crystallinity << fixed << setprecision(3) << static_cast<double>(cage_form_sum) / (H2O + add) << endl;
}
void cage::write_cage_summary(int& i,int& fill, int& SEC_cage_size, int& non_SEC_cage_size, int& IC_cage_size){
	char buffer[256];
	sprintf(buffer, "%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d",
		i,
		SEC_cage_size,
		non_SEC_cage_size,
		IC_cage_size,
		other,
		cage512, cage51262, cage51263, cage51264, cage51265, cage51266, cage51268,
		cage4151062, cage4151063, cage4151064, cage4151065,
		cage425861, cage425862, cage425863, cage425864,
		cage435663, cage435664);

	ofs_summary << buffer << endl;
	if (cal_guest) {
		int cage_size = SEC_cage_size + non_SEC_cage_size + IC_cage_size;
		if (cage_size > 0) {
			double vac = static_cast<double>(fill) / cage_size;
			sprintf(buffer, "%-8d%-8.3f%-8d%-8d",i, 1 - vac, cage_size, fill);
			ofs_occupancy << buffer << endl;
		}
		else {
			sprintf(buffer, "%-8d%-8c%-8d%-8d",i, 'x', cage_size, fill);
			ofs_occupancy << buffer << endl;
		}
	}
}
void cage::output_visual(int i,const vector<polyhedron> &SEC_cage,const vector<polyhedron> &non_SEC_cage,const vector<polyhedron> &IC_cage){
	for (size_t j = 0; j < SEC_cage.size(); ++j) {
		for(int k=0; k< SEC_cage[j].vertex.size();k++){
			if(SEC_cage[j].vertex[k]>=H2O)
				break;
			if(H2O_color[SEC_cage[j].vertex[k]]>SEC_cage[j].color_tier){
				H2O_color[SEC_cage[j].vertex[k]]=SEC_cage[j].color_tier;
			}
		}
        }
        for (size_t j = 0; j < non_SEC_cage.size(); ++j) {
		for(int k=0; k< non_SEC_cage[j].vertex.size();k++){
			if(non_SEC_cage[j].vertex[k]>=H2O)
                                break;
                        if(H2O_color[non_SEC_cage[j].vertex[k]]>non_SEC_cage[j].color_tier)
                                H2O_color[non_SEC_cage[j].vertex[k]]=non_SEC_cage[j].color_tier;
                }
        }
        for (size_t j = 0; j < IC_cage.size(); ++j) {
		for(int k=0; k< IC_cage[j].vertex.size();k++){
                        if(IC_cage[j].vertex[k]>=H2O)
                                break;
			if(H2O_color[IC_cage[j].vertex[k]]>IC_cage[j].color_tier)
                                H2O_color[IC_cage[j].vertex[k]]=IC_cage[j].color_tier;
                }
        }
	char text[100];
	sprintf(text,"Generated by smooth trjconv : ");
        ofs_visual_gro<<text<<" t= "<<i<<endl;
	int total_molecule=H2O;
	if(cal_guest)
		total_molecule+=guest;
	if(cal_additive)
		total_molecule+=add*add_atom;
	ofs_visual_gro<<" "<<total_molecule<<endl;
	string type,atom_name;
	//H2O
	atom_name="OICE";
	type="H2O";	
	for (int j = 0; j < H2O; j++) {
		int id1 = j + 1;
		if (id1 > 99999) id1 %= 100000;
		sprintf(text, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%",id1, type.c_str() , atom_name.c_str(), id1,H2O_molecule[j][0].x ,H2O_molecule[j][0].y ,H2O_molecule[j][0].z);
		ofs_visual_gro<<text<<endl;
		ofs_visual_index<<H2O_color[j]<<" ";
	}
	//guest
	if(cal_guest){
		atom_name="C";
        	type="guest";
		for (int j = 0; j < guest; j++) {
			int id1 = j + 1;
			if (id1 > 99999) id1 %= 100000;
			sprintf(text, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%",id1, type.c_str() , atom_name.c_str(),id1,guest_molecule[j].x ,guest_molecule[j].y,guest_molecule[j].z);
			ofs_visual_gro<<text<<endl;
			ofs_visual_index<<9<<" ";
        	}
	}	
	//additive
	if(cal_additive){
		atom_name="a";
                type="add";
                for (int j = 0; j < add; j++) {
			for(int k=0;k<add_atom;k++){
				int id1 = j + 1;
                                int id2 = j * add_atom + k + 1;
                                if (id1 > 99999) id1 %= 100000;
                                if (id2 > 99999) id2 %= 100000;
                                sprintf(text, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%",id1, type.c_str() , add_name[k].c_str(), id2,add_molecule[j][k].x ,add_molecule[j][k].y ,add_molecule[j][k].z);
                        	ofs_visual_gro<<text<<endl;
                        	ofs_visual_index<<10<<" ";
			}
                }

	}
	sprintf(text, "%10.5f%10.5f%10.5f", box.x, box.y, box.z);
        ofs_visual_gro<<text<<endl;	
	ofs_visual_index<<endl;
}
