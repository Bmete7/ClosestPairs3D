/*
@Author: Burak Mete
*/

#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
using namespace std;

class Point {
	int x;
	int y;
	int z;
public:
	Point(int x = 0, int y = 0, int z = 0) : x(x), y(y), z(z) { ; };
	int getX()const { return x; };
	int getY()const { return y; };
	int getZ()const { return z; };
};

class Space {
	int counter;
public:
	vector<Point> points;
	Space(int s) { points.resize(s); counter = 0; };
	void addPoints(const Point &p) { points[counter] = p; ++counter; };
	void print();
	int getLength() { return counter; };

};

int partitionX(vector<Point>&p, int first, int last) { // quicksort conquer step for sorting X values
	int pivotVal = p[last].getX();
	int i = first - 1;
	for (int j = first; j < last; ++j) {
		if (p[j].getX() < pivotVal) {
			++i;
			Point p1 = p[i];
			p[i] = p[j];
			p[j] = p1;
		}
	}
	Point p1 = p[last];
	p[last] = p[++i];
	p[i] = p1;
	return i;
}

int partitionY(vector<Point>&p, int first, int last) { // quicksort conquer step for sorting Y values
	int pivotVal = p[last].getY();
	int i = first - 1;
	for (int j = first; j < last; ++j) {
		if (p[j].getY() < pivotVal) {
			++i;
			Point p1 = p[i];
			p[i] = p[j];
			p[j] = p1;
		}
	}
	Point p1 = p[last];
	p[last] = p[++i];
	p[i] = p1;
	return i;
}

void quickSort(string coor, vector<Point>&p,int first,int last) { // Quicksort divide step for both X and Y values, takes string arguement for sorting axis
	if (last > first) {
		int piv = 0;
		if(coor == "x")
			piv = partitionX(p, first, last);
		else if (coor == "y")
			piv = partitionY(p, first, last);
		
		quickSort(coor,p, first, piv - 1);
		quickSort(coor,p, piv+1, last);
	}
}


double find_distance_3d(const Point&p1, const Point &p2) {
	return sqrt(pow(p1.getX() - p2.getX(), 2) + pow(p1.getY() - p2.getY(), 2) + pow(p1.getZ() - p2.getZ(), 2));
}

double find_distance_2d(const Point&p1, const Point &p2) {
	return sqrt(pow(p1.getX() - p2.getX(), 2) + pow(p1.getY() - p2.getY(), 2));
}

double closest_pair_3D(const vector <Point>& p, const vector <Point>& py, int &calc) {  

	if (p.size() == 1) return 99999; // might as well +inft, but i choose that arbitrary max value

	else if (p.size() == 2) { // 
		calc++;
		return find_distance_3d(p[0], p[1]);
	}
	else {
		int mid = p.size() / 2;

		vector<Point>p1(p.begin(), p.begin() + mid);
		vector<Point>p2(p.begin() + mid, p.end());

		vector<Point>py1 = py;
		vector<Point>py2 = py;
		int left_size = 0;
		int right_size = 0;
		for (int i = 0; i < py.size(); ++i) {
			if (py[i].getX() <= p[mid].getX()) {
				py1[left_size++] = py[i];
			}
			else {
				py2[right_size++] = py[i];
			}
		}
		py1.resize(left_size);
		py2.resize(right_size);

		double delta_1 = closest_pair_3D(p1, py1,calc); // t(n/2)
		double delta_2 = closest_pair_3D(p2, py2,calc); // t(n/2)
		double delta = min(delta_1, delta_2);
		int counter = 0;
		vector<Point>strip = py;
		for (int i = 0; i < py.size(); ++i) {
			if (py[i].getX() >(p[mid].getX() - delta) && py[i].getX() < (p[mid].getX() + delta)) {
				strip[counter++] = py[i];
			}
		}  
		strip.resize(counter); 
		for (int i = 0; i < strip.size(); ++i) {
			for (int j = i + 1; j < strip.size() && abs(strip[j].getY() - strip[i].getY()) < delta; ++j) {
				delta = min(find_distance_3d(strip[i], strip[j]), delta);
				calc++;
			}
		}
		return delta;
	}
}

void Space::print() {
	for (int i = 0; i < counter; ++i) {
		cout << points[i].getX() << " " << " " << points[i].getY() << " " << points[i].getZ() << endl;
	}
}

int main(int argc, char * argv[]) {
	ifstream ii;
	const char* file_name= argv[1];
	ii.open(file_name);
	
	string input_s;
	int input_size;
	getline(ii, input_s);
	istringstream(input_s) >> input_size;
	int x, y, z;
	Space s1(input_size);

	for (int i = 0; i < input_size; ++i) { 
		string xx, yy, zz;
		getline(ii, xx, ' ');
		getline(ii, yy, ' ');
		getline(ii, zz, '\n');
		istringstream(xx) >> x;
		istringstream(yy) >> y;
		istringstream(zz) >> z;
		Point p1(x, y, z);
		s1.addPoints(p1);
	}
	
	quickSort("x", s1.points, 0, s1.points.size() - 1); // firstly, sort the array by x coordinates
	vector<Point> Py = s1.points;
	quickSort("y",Py, 0, s1.points.size() - 1); // sort the array elsewhere by y coordinates, (reduces lgn)
	int dist_calc = 0;
	clock_t beginClock = clock();
	cout << "The distance is " << closest_pair_3D(s1.points, Py,dist_calc) << endl;
	clock_t endClock = clock();

	cout << "Time elapsed:\t" << (double(endClock - beginClock) / CLOCKS_PER_SEC)* 1000 << " ms" << endl; 
	cout << "Number of total distance calculations is " << dist_calc << endl;
	return 0;
}
