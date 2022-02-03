/*********** XGC *************/

#include <string>

void heatload();
void init();    // initialization
void load_data();  //receving data from XGC1
void heatload_calc(); // calculate heatload
void output();    // output graphs or data for graphs

//extern "C" void set_test_type(int test_type);

int main(int argc, char* argv[])
{
    // Parse command line arguments
    //set_test_type(0);
    if (argc>2){
        printf("ERROR: Too many command line arguments. Available options are '--test', '--update-test', or neither.\n");
        exit(1);
    }
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--test"){
            //set_test_type(1);
        } else if (std::string(argv[i]) == "--update-test") {
            //set_test_type(2);
        } else {
            printf("ERROR: Unknown command line argument. Available options are '--test', '--update-test', or neither.\n");
            exit(1);
        }
    }

    // run actual routine
    heatload();
}

void heatload()
{
    init();

    // put loop?
    load_data();

    heatload_calc();

    output();

}