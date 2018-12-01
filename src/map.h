#ifndef MAP_H_
#define MAP_H_

class Map {
public:
	
	struct single_landmark_s{

		int id_i ;
		float x_f;
		float y_f;
	};

	std::vector<single_landmark_s> landmark_list ; // List of landmarks in the map

};



#endif /* MAP_H_ */
