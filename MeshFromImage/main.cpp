#include <cvheads.h>
#include <cvtool.h>
#include <stlhandler.h>
#include <eitsheads.h>
//#include <fem.h>

using namespace EITS;

//DEFALT values
#define SIZE_ELECTRODE_DEFAULT 1.8
#define NUM_ELECTRODE_DEFAULT 16
#define SIZE_SHEET_DEFAULT cv::Size(300, 300)
#define SIZE_SHEET_MONITOR_DEFAULT cv::Size(800, 800)
#define NUM_SPLIT_DEFAULT 20

//Window name
cv::String wname_video = "video";
cv::String wname_input = "input";
cv::String wname_trim = "trim";
cv::String wname_mesh = "mesh";
cv::String wname_electrode = "electrode";

//Point stack for homography transformation
int index_corner_cur = 0;
int index_electrode_cur = 0;
int num_click = 0;
std::vector<cv::Point2f> pos_courner;
std::vector<cv::Point2f> pos_object;
cv::Mat H;

//Image data
cv::Mat img;
cv::Mat img_org;
cv::Mat img_resize;
cv::Mat img_trim;
cv::Mat img_trim_org;
cv::Mat img_trim_gray;
cv::Mat img_trim_binary;
cv::Mat img_trim_binary_rgb;
cv::Mat img_trim_dst;

//Electrode information
int num_electrode = NUM_ELECTRODE_DEFAULT;
double size_electrode = SIZE_ELECTRODE_DEFAULT;
//Objecct size [mm]
cv::Size size_sheet = SIZE_SHEET_DEFAULT;
cv::Size size_sheet_monitor = SIZE_SHEET_MONITOR_DEFAULT;
double scale_r2m;
int num_split = NUM_SPLIT_DEFAULT;

//Clicked point stack for electrode
std::vector<cv::Vec2d> pos_electrode;
//Electrode label
int label = SENSING;

//Mesh node
std::vector<cv::Point2f> points;
//Detected point stack for electrode
std::vector<cv::Point2f> points_electrode;

//Threshold range
int range[6] = {255, 0, 255, 0, 255, 0};

//File name
char filename_input[256];
char filename_fem2d[256];
char filename_e2es[256];

bool is_selected_courners = false;
bool is_selected_electrodes = false;
//Image mode (IMAGE or CAMERA)
int mode = IMAGE;


//Set threshold from selected region
void setRange(int *_range, cv::Rect _roi, cv::Mat _img)
{
	range[0] = range[2] = range[4] = 255;
	range[1] = range[3] = range[5] = 0;
	cv::Mat img_tmp = _img(_roi);
	for(int i = 0; i<img_tmp.size().area(); i++){
		if(img_tmp.at<cv::Vec3b>(i)[0] < range[0])range[0] = img_tmp.at<cv::Vec3b>(i)[0]; 
		if(img_tmp.at<cv::Vec3b>(i)[0] > range[1])range[1] = img_tmp.at<cv::Vec3b>(i)[0]; 
		if(img_tmp.at<cv::Vec3b>(i)[1] < range[2])range[2] = img_tmp.at<cv::Vec3b>(i)[1]; 
		if(img_tmp.at<cv::Vec3b>(i)[1] > range[3])range[3] = img_tmp.at<cv::Vec3b>(i)[1]; 
		if(img_tmp.at<cv::Vec3b>(i)[2] < range[4])range[4] = img_tmp.at<cv::Vec3b>(i)[2]; 
		if(img_tmp.at<cv::Vec3b>(i)[2] > range[5])range[5] = img_tmp.at<cv::Vec3b>(i)[2]; 
	}
#ifdef _DEBUG
	std::cout<<range[0]<<","<<range[1]<<","<<range[2]<<","<<range[3]<<","<<range[4]<<","<<range[5]<<std::endl;
#endif
}

//Generate binarized image
void genChromaKey(int *_range, cv::Mat *_img, cv::Mat *_img_binary)
{
	int alpha = 15;
	for(int i = 0; i<_img->size().area(); i++){
		if(_img->at<cv::Vec3b>(i)[0] < range[0] - alpha || _img->at<cv::Vec3b>(i)[0] > range[1] + alpha
		||_img->at<cv::Vec3b>(i)[1] < range[2] - alpha || _img->at<cv::Vec3b>(i)[1] > range[3] + alpha 
		||_img->at<cv::Vec3b>(i)[2] < range[4] - alpha || _img->at<cv::Vec3b>(i)[2] > range[5] + alpha ){
			_img_binary->at<uchar>(i) = 0;
		}else{
			_img_binary->at<uchar>(i) = 255;
		}
	}
}

//STEP1: Click four object's corners
void onClickCorner(int _event, int _x, int _y, int _flags, void* _data)
{
	char text[256];

	if (_flags == cv::EVENT_FLAG_LBUTTON) {
		img_org.copyTo(img_resize);
		pos_courner[index_corner_cur] =cv::Point2i(_x, _y);
	}
	if (_event == cv::EVENT_LBUTTONUP) {
		index_corner_cur = (index_corner_cur + 1) % 4;

		if(index_corner_cur == 0){
			H = cv::findHomography(pos_object, pos_courner);
			img_trim = CVTOOL::trimImageWithHomography(size_sheet_monitor, img_org, H);
			cv::cvtColor(img_trim, img_trim_gray, cv::COLOR_BGR2GRAY);	
			is_selected_courners = true;
			cv::imshow(wname_trim, img_trim_gray);
		}
    }

	for(int i = 0 ; i < 4; i++){
		sprintf_s(text,"%d", i);
		cv::putText(img_resize, text, pos_courner[i], cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::SCALAR_WHITE);
		cv::circle(img_resize, pos_courner[i], 5, cv::SCALAR_RED);
	}
	cv::imshow(wname_input, img_resize);
}

//STEP2: UI for selecting electrode region to binarize
void onSelectRegion(int _event, int _x, int _y, int _flags, void* _data)
{
	cv::Rect *roi_selected = (cv::Rect *)_data;
	if (_event == cv::EVENT_LBUTTONDOWN) {
		roi_selected->x = _x;
		roi_selected->y = _y;
	}
	if (_flags == cv::EVENT_FLAG_LBUTTON) {
		roi_selected->width = _x - roi_selected->x;
		roi_selected->height = _y - roi_selected->y;
		if (!img_trim.empty()) {
			cv::cvtColor(img_trim, img_trim_gray, cv::COLOR_BGR2GRAY);
			cv::rectangle(img_trim_gray, *roi_selected, cv::Scalar(255, 255, 255));
			cv::imshow(wname_trim, img_trim_gray);
			is_selected_electrodes = true;
		}
	}
	if (_event == cv::EVENT_LBUTTONUP) {
		if (roi_selected->area()>0) {
			points.clear();
			img_trim_gray.copyTo(img_trim_binary);

			setRange(range, *roi_selected, img_trim);
			genChromaKey(range, &img_trim, &img_trim_binary);

			//Dilate and erode to reduce noise
			cv::dilate(img_trim_binary, img_trim_binary, cv::Mat(), cv::Point(-1, -1), 1);
			cv::erode(img_trim_binary, img_trim_binary, cv::Mat(), cv::Point(-1, -1), 1);
			cv::cvtColor(img_trim_binary, img_trim_binary_rgb, CV_GRAY2BGR);

			//Draw contours
			std::vector< std::vector<cv::Point >> contours;
			cv::findContours(img_trim_binary, contours, CV_RETR_LIST, CV_CHAIN_APPROX_TC89_KCOS);
			cv::Rect roi_sheet = cv::Rect(0, 0, size_sheet.width + 1, size_sheet.height + 1);
			bool is_close = false;

			for (int i = 0; i < contours.size(); i++) {
				if (contours[i].size()> 5) {
					for (int j = 0; j < contours[i].size(); j++) {
						is_close = false;
						for (int k = 0; k < j; k++) {
							if (k != j) {
								cv::Vec2d dist = cv::Vec2d(contours[i][j].x - contours[i][k].x, contours[i][j].y - contours[i][k].y);
								if (dist.ddot(dist)<4) {
									is_close = true;
								}
							}
						}
						if (!is_close) {
							if (roi_sheet.contains(cv::Point2f(contours[i][j].x / scale_r2m, contours[i][j].y / scale_r2m))) {
								points.push_back(cv::Point2f(contours[i][j].x / scale_r2m, contours[i][j].y / scale_r2m));
								cv::circle(img_trim_binary_rgb, contours[i][j], 2, cv::Scalar(0, 0, 255), -1);
							}
						}
					}
				}
			}

			cv::imshow(wname_mesh, img_trim_binary_rgb);
			std::cout << "   : Press any key if electrodes are detected correctly." << std::endl;

		}
	}
}

//STEP 3: Click electrod for setting electrode label
void onClickElectrode(int _event, int _x, int _y, int _flags, void* _data)
{
	char text[256];
	img_trim_binary_rgb.copyTo(img_trim_org);
	img_trim_binary_rgb.copyTo(img_trim_dst);
	if (_flags == cv::EVENT_FLAG_LBUTTON) {
		pos_electrode[index_electrode_cur] = cv::Vec2d(_x, _y);
	}
	if (_event == cv::EVENT_LBUTTONUP) {
		index_electrode_cur = (index_electrode_cur + 1) % num_electrode;
		num_click++;
	}
	for(int i = 0 ; i < num_electrode; i++){
		sprintf_s(text,"%d", i);
		cv::putText(img_trim_org, text, cv::Point2d(pos_electrode[i].val[0],pos_electrode[i].val[1]), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::SCALAR_WHITE);
		cv::circle(img_trim_org,cv::Point2d(pos_electrode[i].val[0],pos_electrode[i].val[1]), 5, cv::SCALAR_BLUE);
	}
	cv::imshow(wname_mesh, img_trim_org);
}

bool saveElectrodeSetting(const char* _filename, StlMesh *_stl)
{
	std::cout << "Saving electrode setting...";
	int pre = 0;
	std::ofstream file;
	std::vector<cv::Vec2d> pos_ele;
	std::vector<int> label_ele;
	std::vector<int> index_ele;
	file.open(_filename, std::ios::out|std::ios::binary);
	if (!_stl->getIsLoaded() || !file.is_open()){
		std::cout << "[FAIL]" << std::endl;
		return false;
	}

	//Set index for electrode nodes
	std::vector<int> index_edge;
	std::vector<int> is_matching;
	std::vector<int> is_used;
	is_matching.resize(points_electrode.size(), false);
	is_used.resize(_stl->getNumNode(), false);
	index_edge.resize(points_electrode.size());
	for(int i=0;i<_stl->getNumNode();i++){
		for(int j=0;j<points_electrode.size();j++){
			if(!is_matching[j]){
				if(_stl->getVertex(i).x == points_electrode[j].x - size_sheet.width / 2 &&
					_stl->getVertex(i).y == points_electrode[j].y - size_sheet.height / 2 ){
					index_edge[j] = i;
						is_matching[j] = true;
						is_used[i] = true;
				}
			}
		}
	}

	double dir_min;
	int index_min;
	cv::Point t_pos_electrode;
	ColorMap color;
	color.setParam(0, num_electrode-1);
	for(int j=0;j<points_electrode.size();j++){
		dir_min = 1000;
		for(int e = 0;e<pos_electrode.size();e++){
			cv::Vec2d dir = cv::Vec2d(points_electrode[j].x, points_electrode[j].y) - pos_electrode[e]/scale_r2m;
			if(sqrt(dir.dot(dir)) < dir_min){
				dir_min = sqrt(dir.dot(dir));
				index_min = e;
			}
		}
		t_pos_electrode = points_electrode[j] * scale_r2m;
		Vector3d t_color = color.getColorMapForCV(index_min, MAP_JET);
		cv::circle(img_trim_dst, t_pos_electrode, 4, cv::Scalar(t_color.x, t_color.y, t_color.z), -1);
		label_ele.push_back(index_min);
		index_ele.push_back(index_edge[j]);
	}

	for(int i=0;i<_stl->getNumNode();i++){
		if(!is_used[i]){
			pos_ele.push_back(cv::Vec2d(_stl->getVertex(i).x, _stl->getVertex(i).y));			
			for(pre = 0;pre<pos_electrode.size();pre++){
				cv::Vec2d dir = pos_ele[pos_ele.size()-1] + cv::Vec2d(size_sheet.width/2,size_sheet.height/2) - pos_electrode[pre]/scale_r2m;
				if(sqrt(dir.dot(dir)) < size_electrode){
					t_pos_electrode = cv::Point((_stl->getVertex(i).x + size_sheet.width / 2)*scale_r2m, (_stl->getVertex(i).y + size_sheet.height / 2) * scale_r2m);
					Vector3d t_color = color.getColorMapForCV(pre, MAP_JET);
					cv::circle(img_trim_dst, t_pos_electrode, 4, cv::Scalar(t_color.x, t_color.y, t_color.z), -1);
					label_ele.push_back(pre);
					index_ele.push_back(i);
					break;
				}
			}
		}
	}
	cv::imshow(wname_electrode, img_trim_dst);

	file<<pos_electrode.size()<<","<<label<<","<<index_ele.size()<<std::endl;
	//Node indexCStateCElectrode index
	for(int i = 0;i<index_ele.size();i++){
		file<<index_ele[i]<<","<<label<<","<<label_ele[i];
		if(i!=index_ele.size()-1)file<<std::endl;
	}
	file.clear();
	file.close();

	std::cout << "[OK]" << std::endl;
	return true;
}

bool loadProject(char *_filename)
{
	std::cout<<"Loading project...";
	char buf[256];
	std::ifstream file;
	file.open(_filename, std::ios::in);
	if(!file.is_open()){
		std::cout<<"[FAIL]"<<std::endl;
		return false;
	}
	while(!file.eof()){
		file>>buf;
		if(strstr(buf,"input:"))
			sscanf_s(buf,"input:%s", filename_input, (unsigned)sizeof(filename_input));
		else if(strstr(buf,"fem2d:"))
			sscanf_s(buf,"fem2d:%s", filename_fem2d, (unsigned)sizeof(filename_fem2d));
		else if(strstr(buf,"e2es:"))
			sscanf_s(buf,"e2es:%s", filename_e2es, (unsigned)sizeof(filename_e2es));
		else if(strstr(buf,"width:"))
			sscanf_s(buf,"width:%d", &size_sheet.width);
		else if(strstr(buf,"height:"))
			sscanf_s(buf,"height:%d", &size_sheet.height);
		else if(strstr(buf,"nelectrode:"))
			sscanf_s(buf,"nelectrode:%d", &num_electrode);
		else if(strstr(buf,"selectrode:"))
			sscanf_s(buf,"selectrode:%lf", &size_electrode);
		else if(strstr(buf,"nsplit:"))
			sscanf_s(buf,"nsplit:%d", &num_split);
		else if(strstr(buf,"mode:"))
			sscanf_s(buf,"mode:%d", &mode);
	}
	file.close();
	std::cout << "[OK]" << std::endl;
	return true;
}

int main(int _argc, char *_argv[])
{
	std::cout<<"HELLO:)"<<std::endl;
	std::cout<<"********************* MeshFromImage v1.0.0 *********************"<<std::endl;

	//Load project file
	char filename[256];
	if(_argc > 2)
		sscanf_s(_argv[1], "%s", filename, (unsigned)sizeof(filename));
	else 	sprintf_s(filename, "../data/sample.meshgen");
	if(!loadProject(filename)){
		std::cout<<"DEFAULT SETTINGS"<<std::endl;
		sprintf_s(filename_input, "../data/sample.jpg");
		sprintf_s(filename_fem2d, "../data/sample.fem2d");
		sprintf_s(filename_e2es, "../data/sample.e2es");
	}
	std::cout<<"Mode (0:image, 1:camera) : "<<mode<<std::endl;
	std::cout<<"Input image name : "<<filename_input<<std::endl;
	std::cout<<"Output mesh name : "<<filename_fem2d<<std::endl;
	std::cout<<"Output electrode setting name : "<<filename_e2es<<std::endl;
	std::cout<<"Sheet size W*H [mm] : "<<size_sheet<<std::endl;
	std::cout<<"Number of electrode : "<<num_electrode<<std::endl;
	std::cout<<"Size of electrode [mm] : "<<size_electrode<<std::endl;
	std::cout<<"Number of split (Hrizontal) : "<<num_split<<std::endl;
	std::cout<<"*******************************************************"<<std::endl;

	//Load or capture an image
	if(mode == IMAGE){
		std::cout << "Loading image ...";
		img = cv::imread(filename_input);
		if (img.data == NULL) {
			std::cout << "[FAIL]" << std::endl;
			Sleep(2000);
			return -1;
		}
		std::cout << "[OK]" << std::endl;
	}else if(mode == CAMERA){
		std::cout << "Please type any key to capture an image ..." << std::endl;
		cv::VideoCapture cam(0);
		if(cam.isOpened()){
			while(true){
				cam >> img;
				cv::imshow(wname_video, img);
				if(cv::waitKey(30)>=0)break;
			}
			cv::destroyWindow(wname_video);
			char filename[256];
			sprintf_s(filename, "%s.capture.jpg", filename_input);
			if(cv::imwrite(filename, img))
				std::cout << "[OK]" << std::endl;
			else {
				std::cout << "[FAIL]" << std::endl;
				Sleep(2000);
				return -1;
			}
		}else{
			std::cout<<"[ERROR] No camera connection :("<<std::endl;
			Sleep(2000);
			return -1;
		}
	}
	if(img.empty()){
		std::cout<<"[ERROR] Bad image :("<<std::endl;
		Sleep(2000);
		return -1;
	}

	//Scaling for visualization
	scale_r2m = (double)(size_sheet_monitor.width-1) / (size_sheet.width);
	size_sheet_monitor.height = scale_r2m * (size_sheet.height)+1;

	pos_object.push_back(cv::Point(0,0));
	pos_object.push_back(cv::Point(size_sheet_monitor.width,0));
	pos_object.push_back(cv::Point(size_sheet_monitor.width,size_sheet_monitor.height));
	pos_object.push_back(cv::Point(0,size_sheet_monitor.height));

	pos_courner.resize(4);
	pos_electrode.resize(num_electrode);

	//Resize input image
	cv::Size size_new;
	cv::Rect roi_selected;
	size_new.width = cv::SIZE_SUPERVGA.width;
	size_new.height = img.size().height*cv::SIZE_SUPERVGA.width/img.size().width;
	cv::resize(img, img_resize, size_new);
	cv::resize(img, img_org, size_new);
	roi_selected = cv::Rect(0, 0, size_sheet_monitor.width, size_sheet_monitor.height);

	std::cout << "STEP 1 : Please click four object's courners on window " << wname_input << std::endl;
	cv::imshow(wname_input, img_resize);
	cv::setMouseCallback(wname_input, onClickCorner);
	while (!is_selected_courners) {
		cv::waitKey(10);
	}

	std::cout << "STEP 2 : Please drag electrode regoin on the window " << wname_trim << std::endl;
	cv::namedWindow(wname_trim, CV_WINDOW_AUTOSIZE);
	cv::setMouseCallback(wname_trim, onSelectRegion, (void*)&roi_selected);
	while (!is_selected_electrodes) {
		cv::waitKey(0);
	}

	//Electrode points
	points_electrode = points;

	//Grid points
	double dl = (double)(size_sheet.width) / (num_split-1);
	for(int i = 0; i < num_split; i++){
		for(int j = 0; j < num_split; j++){
			points.push_back(cv::Point2f(i*dl, j*dl));
		}
	}

	//Delaunay triangulation
	cv::Subdiv2D subdiv;
	subdiv.initDelaunay(cv::Rect(0, 0, size_sheet.width+1, size_sheet.height+1));
	subdiv.insert(points);
	std::vector<int> idx;
	std::vector<std::vector<cv::Point2f>> facet_lists;
	std::vector<cv::Point2f> facet_centers;
	subdiv.getVoronoiFacetList(idx, facet_lists, facet_centers);
	std::vector<cv::Vec4f> edge_list;
	subdiv.getEdgeList(edge_list);
	std::vector<cv::Vec6f> triangles;
	subdiv.getTriangleList(triangles);

	//Draw mesh
	for(auto it = triangles.begin(); it != triangles.end(); it++)
	{
		cv::Vec6f &vec = *it;
			if(vec[0]>=0&&vec[1]>=0&&vec[2]>=0&&vec[3]>=0&&vec[4]>=0&&vec[5]>=0
				&&vec[0]<=size_sheet.width&&vec[1]<=size_sheet.height
				&&vec[2]<=size_sheet.width&&vec[3]<=size_sheet.height
				&&vec[4]<=size_sheet.width&&vec[5]<=size_sheet.height){
					cv::Point p1(scale_r2m * vec[0], scale_r2m * vec[1]);
					cv::Point p2(scale_r2m * vec[2], scale_r2m * vec[3]);
					cv::Point p3(scale_r2m * vec[4], scale_r2m * vec[5]);
					cv::line(img_trim_binary_rgb, p1, p2, cv::Scalar(64,255,128));
					cv::line(img_trim_binary_rgb, p2, p3, cv::Scalar(64,255,128));
					cv::line(img_trim_binary_rgb, p3, p1, cv::Scalar(64,255,128));
			}
	}

	std::cout << "STEP 3 : Please assign labels to all electrodes by click >> " << wname_mesh << std::endl;
	cv::imshow(wname_mesh, img_trim_binary_rgb);
	cv::setMouseCallback(wname_mesh, onClickElectrode);
	//Selected all electrodes?
	while(num_click < num_electrode){
		cv::waitKey(0);
		std::cout<<num_click<<"/"<<num_electrode<<std::endl;
	}

	//Save mesh (fem2d same as stl format)
	StlMesh *stl = new StlMesh;
	stl->deleteMaterial();
	stl->deleteMesh();
	stl->setNumNode(points.size());
	stl->setNumFacet(triangles.size());
	stl->setNumNormal(triangles.size());
	stl->setNumMaterial(triangles.size());
	stl->setNumLine(3 * triangles.size());
	stl->newMesh();
	stl->newMaterial();

	//Set mesh contents (CoG is origin)
	int i = 0;
	for(auto it = triangles.begin(); it != triangles.end(); it++)
	{
		cv::Vec6f &vec = *it;
		Facet t_facet;
		if(vec[0]>=0&&vec[1]>=0&&vec[2]>=0&&vec[3]>=0&&vec[4]>=0&&vec[5]>=0
			&&vec[0]<=size_sheet.width&&vec[1]<=size_sheet.height
			&&vec[2]<=size_sheet.width&&vec[3]<=size_sheet.height
			&&vec[4]<=size_sheet.width&&vec[5]<=size_sheet.height){
			t_facet.vertex[0] = Vector3d(vec[0]-size_sheet.width/2.0, vec[1]-size_sheet.height/2.0, 0);
			t_facet.vertex[1] = Vector3d(vec[2]-size_sheet.width/2.0, vec[3]-size_sheet.height/2.0, 0);
			t_facet.vertex[2] = Vector3d(vec[4]-size_sheet.width/2.0, vec[5]-size_sheet.height/2.0, 0);
			stl->setFacet(i, t_facet);
			i++;
		}
	}
	stl->setNumFacet(i);
	stl->setup();

	if (stl->save(filename_fem2d))
		saveElectrodeSetting(filename_e2es, stl);

	std::cout << "Press any key to finish" << std::endl;
	cv::waitKey(0);
	std::cout<<"GOOD BYE:("<<std::endl;
	cv::destroyAllWindows();

	Sleep(2000);

	return 0;
}