
class MeshWriter{

public:

        MeshWriter() = default;                               //default constructor
        MeshWriter(const MeshWriter& v) = delete;           //copy constructor
        MeshWriter(MeshWriter&& v) = delete;                //move constructor
        MeshWriter& operator=(const MeshWriter& v) = delete;//copy assignment operator
        MeshWriter& operator=(MeshWriter&& v) = delete;     //move assignment operator 


        static void writeVTK(const std::string &filename,
			 const int mesh[3],
                         const std::vector<std::shared_ptr<LBM_node>> &nodes);
	static void writeflagVTK(const std::string &filename,
                         const int mesh[3],
                         const std::vector<std::shared_ptr<IB_node>> &nodes);

};



