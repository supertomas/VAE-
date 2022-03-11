#include"bits/stdc++.h"
#include <sstream>
#include <fstream>
using namespace::std;
#define MAX_NUM_NODES 200000
int node_count = 0;
long double t = 0.15;

#pragma warning(disable : 4996)
class Node {
public:
    std::string name;
    std::string sequence = "";

    Node* left;
    Node* right;
    long double left_length = 0;
    long double right_length = 0;
};

std::array<class Node, MAX_NUM_NODES> NODES; // （仮）
int num_used_node = 0;
int num_left_parenthesis = 0;
int num_right_parenthesis = 0;
std::string sequence_name;
std::map <std::string, std::string>name_map;
std::map <std::string, std::string>::iterator pp;




int find_central_comma(std::string Phy_str) {
    int comma_pos = 0;
    //return comma position, which is outside of parenthesis
    //(or:return strings of both sides splitted with the comma)
    //input:
    //      std::string Phy_str:strings of phylogenetic tree
    //return:
    //      int:position of the comma

    int Phy_size = Phy_str.length();
    for (int char_pos = 0; char_pos < Phy_size; ++char_pos) {
        if (Phy_str[char_pos] == '(') num_left_parenthesis++;
        else if (Phy_str[char_pos] == ')') num_right_parenthesis++;

        if (num_left_parenthesis == num_right_parenthesis && Phy_str[char_pos] == ',') {
            comma_pos = char_pos;
            return comma_pos;
        }
        //this code should not be implemented in this function
        /*if(num_left_parenthesis == num_right_parenthesis && Phy_str[char_pos] == ':'){
            colon_pos = char_pos;
            return colon_pos;
        }*/
    }
    return -1;
}

std::string exclude_outer_parenthesis(std::string Phy_str) {
    if (Phy_str[0] == '(' && Phy_str[Phy_str.length() - 1] == ')') {
        return Phy_str.substr(1, Phy_str.length() - 2);
    }
    return Phy_str;
}


int sequence_part() {
    std::ifstream file_open("BLAT_data_upper1000_sequence_remove_gap_10percent_211124.a2m");
    std::string str;
    sequence_name = "";
    if (file_open) {
        while (true) {
            getline(file_open, str);
            if (str[0] == '>') {
                sequence_name = str.substr(1, str.size() - 1);
                name_map.insert(std::pair<std::string, std::string>(sequence_name, ""));
            }
            else name_map[sequence_name] += str;
            if (file_open.eof())
                break;
        }
    }
    return 0;
}
Node* new_node() {
    num_used_node++;
    return &NODES[num_used_node];
}
int sequence_length;

Node* separate(std::string Phy_str) {

    Node* node = new_node();

    Phy_str = exclude_outer_parenthesis(Phy_str);

    //when input Phy_str is "name:length"
    int comma_pos = find_central_comma(Phy_str);
    if (comma_pos == -1) {
        int colon_pos = Phy_str.find(':');
        std::string name_string = Phy_str.substr(0, colon_pos);
        std::string name_of_sequence;
        if (name_string[0] == '\'')
        {
            name_of_sequence = name_string.substr(1, name_string.length() - 2);
        }
        else
        {
            name_of_sequence = name_string;
        }

        pp = name_map.find(name_of_sequence);
        if (pp != name_map.end()) {
            node->name = name_of_sequence;
            node->sequence = pp->second;
            sequence_length = node->sequence.size();
        }
        return node;
    }
    else {
        std::string left_str = Phy_str.substr(0, comma_pos);

        Node* leftnode = separate(left_str.substr(0, left_str.rfind(':')));
        node->left = leftnode;
        node->left_length = std::stod(left_str.substr(left_str.rfind(':') + 1));

        std::string right_str = Phy_str.substr(comma_pos + 1);

        Node* rightnode = separate(right_str.substr(0, right_str.rfind(':')));
        node->right = rightnode;
        node->right_length = std::stod(right_str.substr(right_str.rfind(':') + 1));

        return node;
    }

}

Node* call_node(Node* node, std::string commands) {
    Node* current_node = node;
    for (char& command : commands) {
        if (command == 'L') current_node = current_node->left;
        else if (command == 'R') current_node = current_node->right;
    }
    return current_node;
}

Node* set_cp_node(Node* node, std::string commands, long double t, int cp_set_mode) {
    /*
    set check point node & return setted check point node pointer
    cp_set_mode 1 -> put check point node at left branch
    cp_set_mode 2 -> put check point node at right branch
    */
    Node* current_node = call_node(node, commands);
    Node* cp_node = new_node();
    switch (cp_set_mode) {
    case 1://left branch case
        cp_node->left = current_node->left;
        current_node->left = cp_node;
        cp_node->left_length = current_node->left_length - t;
        cp_node->right_length = 0;
        current_node->left_length = t;
        return cp_node;

    case 2://right branch case
        cp_node->left = current_node->right;
        current_node->right = cp_node;
        cp_node->left_length = current_node->right_length - t;
        cp_node->right_length = 0;
        current_node->right_length = t;
        return cp_node;
    }
    return nullptr;
}

const int AA_variety = 21;
long double AAmatrix_20[AA_variety][AA_variety]{};
long double AAmatrix_19[AA_variety][AA_variety]{};

int make_matrix() {
    int a = 0;
    int c = 0;
    int d = 0;
    int e = 0;
    int f = 0;
    int g = 0;
    int h = 0;
    int i = 0;
    int k = 0;
    int l = 0;
    int m = 0;
    int n = 0;
    int p = 0;
    int q = 0;
    int r = 0;
    int s = 0;
    int t = 0;
    int v = 0;
    int w = 0;
    int y = 0;
    int gap = 0;

    std::ifstream file_open("BLAT_data_upper1000_sequence_remove_gap_10percent_211124.a2m");
    std::string str;
    std::string Phy_str;

    if (file_open) {
        while (true) {
            getline(file_open, str);
            if (str[0] == '>') {
                Phy_str += "";
            }
            else {
                Phy_str += str;
            }
            if (file_open.eof())
                break;
        }
    }
    for (char j : Phy_str) {
        if (j == 'A')a++;
        else if (j == 'C')c++;
        else if (j == 'D')d++;
        else if (j == 'E')e++;
        else if (j == 'F')f++;
        else if (j == 'G')g++;
        else if (j == 'H')h++;
        else if (j == 'I')i++;
        else if (j == 'K')k++;
        else if (j == 'L')l++;
        else if (j == 'M')m++;
        else if (j == 'N')n++;
        else if (j == 'P')p++;
        else if (j == 'Q')q++;
        else if (j == 'R')r++;
        else if (j == 'S')s++;
        else if (j == 'T')t++;
        else if (j == 'V')v++;
        else if (j == 'W')w++;
        else if (j == 'Y')y++;
        else gap++;
    }
    std::cout << Phy_str.size() << std::endl;

    FILE* matrix_file_open;
    matrix_file_open = fopen("dayhoff_matrix.txt", "r");
    double aamat[20 * 20];
    int num_aa_codes = 20;

    //load matrix value
    int i1, j1;
    for (i1 = 0; i1 < 20; ++i1) {
        for (j1 = 0, aamat[i1 * num_aa_codes + i1] = 0; j1 < i1; ++j1) {
            fscanf(matrix_file_open, "%lf", &aamat[i1 * num_aa_codes + j1]);
            aamat[j1 * num_aa_codes + i1] = aamat[i1 * num_aa_codes + j1];
            std::cout << aamat[i1 * num_aa_codes + j1] << "  ";
        }
        std::cout << std::endl;
    }
    for (int i1 = 0; i1 < AA_variety; i1++) {
        for (int j = 0; j < AA_variety; j++) {
            AAmatrix_20[i1][j] = aamat[i1 * 20 + j];
        }
    }

    for (int i = 0; i < AA_variety; i++) {
        AAmatrix_20[i][20] = 1e-5;
        AAmatrix_20[20][i] = AAmatrix_20[i][20];
    }

    for (int i = 0; i < AA_variety; ++i) {
        for (int j = 0; j < AA_variety; ++j) {
            if (j != i) {
                AAmatrix_20[i][i] -= AAmatrix_20[i][j];
            }
            else if (j == i)
                AAmatrix_20[i][i] += 0;
        }
    }
    AAmatrix_20[20][20] = -20e-5;
    return 0;
}
int make_matrix_lessAA() {

    for (int i = 0; i < AA_variety; ++i) {
        for (int j = 0; j < AA_variety; ++j) {
            AAmatrix_19[i][j] = AAmatrix_20[i][j];
        }
    }
    int exclude_AA;
    exclude_AA = 17;
    for (int i = 0; i < AA_variety; ++i) {
        AAmatrix_19[i][i] += AAmatrix_19[exclude_AA][i];
        AAmatrix_19[exclude_AA][i] = 0;
    }
    for (auto& j : AAmatrix_19) {
        j[exclude_AA] = 0;
    }
    for (auto& k : AAmatrix_19) {
        for (long double i : k) {
            std::cout << i << "  ";
        }
        std::cout << std::endl;
    }
    for (auto& l : AAmatrix_20) {
        for (long double i : l) {
            std::cout << i << "  ";
        }
        std::cout << std::endl;
    }
    return 0;
}
long double identity_matrix[AA_variety][AA_variety]{};
long double make_calculation_matrix20[AA_variety][AA_variety]{};
long double make_calculation_matrix19[AA_variety][AA_variety]{};

void make_calculationmatrix(long double t, int mode) {
    int n = 32;

    for (int i = 0; i < AA_variety; i++) {
        for (int j = 0; j < AA_variety; j++) {
            if (i == j)identity_matrix[i][j] = 1;
            else if (i != j)identity_matrix[i][j] = 0;
        }
    }

    switch (mode) {
    case 20:
        for (int i = 0; i < AA_variety; i++) {
            for (int j = 0; j < AA_variety; j++) {
                make_calculation_matrix20[i][j] = identity_matrix[i][j] + ((t * AAmatrix_20[i][j]) / (n * AAmatrix_20[i][i] * (-1)));
            }
        }
        break;
    case 19:
        for (int i = 0; i < AA_variety; i++) {
            for (int j = 0; j < AA_variety; j++) {
                make_calculation_matrix19[i][j] = identity_matrix[i][j] + ((t * AAmatrix_19[i][j]) / (n * AAmatrix_19[i][i] * (-1)));
            }
        }
        break;
    }
}

int matrix_calculation(long double t, long double substitution_rate_matrix[AA_variety][AA_variety], int mode) {
    int i, j, k;
    int m;
    int n = 32;
    long double copy_calculation_matrix20[AA_variety][AA_variety]{};
    long double copy_calculation_matrix19[AA_variety][AA_variety]{};

    switch (mode) {
    case 20:

        for (i = 0; i < AA_variety; i++) {
            for (int j = 0; j < AA_variety; j++) {
                copy_calculation_matrix20[i][j] = make_calculation_matrix20[i][j];
            }
        }

        for (m = 2; m < n + 1; m++) {
            for (i = 0; i < AA_variety; i++) {
                for (j = 0; j < AA_variety; j++) {
                    substitution_rate_matrix[i][j] = 0;
                    for (k = 0; k < AA_variety; k++) {
                        substitution_rate_matrix[i][j] = substitution_rate_matrix[i][j] + make_calculation_matrix20[i][k] * copy_calculation_matrix20[k][j];
                    }
                }
            }

            for (i = 0; i < AA_variety; i++) {
                for (j = 0; j < AA_variety; j++) {
                    copy_calculation_matrix20[i][j] = substitution_rate_matrix[i][j];
                }
            }
        }

        break;
    case 19:

        for (i = 0; i < AA_variety; i++) {
            for (int j = 0; j < AA_variety; j++) {
                copy_calculation_matrix19[i][j] = make_calculation_matrix19[i][j];
            }
        }

        for (m = 2; m < n + 1; m++) {
            for (i = 0; i < AA_variety; i++) {
                for (j = 0; j < AA_variety; j++) {
                    substitution_rate_matrix[i][j] = 0;
                    for (k = 0; k < AA_variety; k++) {
                        substitution_rate_matrix[i][j] = substitution_rate_matrix[i][j] + make_calculation_matrix19[i][k] * copy_calculation_matrix19[k][j];
                    }
                }
            }

            for (i = 0; i < AA_variety; i++) {
                for (j = 0; j < AA_variety; j++) {
                    copy_calculation_matrix19[i][j] = substitution_rate_matrix[i][j];
                }
            }
        }
        break;
    }
    return 0;
}


bool is_external_node(Node* node) {
    return node->left_length == 0 && node->right_length == 0;
}

bool is_checkpoint_node(Node* node) {
    return node->left_length != 0 && node->right_length == 0;
}

void prior_AA_ancestral_uniform(long double prior_AA_array[], int site) {
    for (int i = 0; i < AA_variety; ++i) {
        prior_AA_array[i] = 1.0 / AA_variety;
        //if (i != 17)prior_AA_array[i] = 1.0 / 20;
        //else if (i == 17)prior_AA_array[i] = 0;
        std::cout << prior_AA_array[i] << std::endl;
    }
}

void normalize_prob(long double probabilities[]) {
    long double sum = 0;
    for (int i = 0; i < AA_variety; ++i) {
        sum += probabilities[i];
    }
    for (int i = 0; i < AA_variety; ++i) {
        probabilities[i] /= sum;
    }
}


long double* Likelihood(Node* node, long double likelihood_ret[], int site, int mode = 20) {
    std::map<char, int> mp;

    if (is_external_node(node)) {
            mp.insert({ 'A', 0 });
            mp.insert({ 'R', 1 });
            mp.insert({ 'N', 2 });
            mp.insert({ 'D', 3 });
            mp.insert({ 'C', 4 });
            mp.insert({ 'Q', 5 });
            mp.insert({ 'E', 6 });
            mp.insert({ 'G', 7 });
            mp.insert({ 'H', 8 });
            mp.insert({ 'I', 9 });
            mp.insert({ 'L', 10 });
            mp.insert({ 'K', 11 });
            mp.insert({ 'M', 12 });
            mp.insert({ 'F', 13 });
            mp.insert({ 'P', 14 });
            mp.insert({ 'S', 15 });
            mp.insert({ 'T', 16 });
            mp.insert({ 'W', 17 });
            mp.insert({ 'Y', 18 });
            mp.insert({ 'V', 19 });
            mp.insert({ '-', 20 });


            auto itr = mp.find(node->sequence[site]);
            for (int i = 0; i < AA_variety; ++i) likelihood_ret[i] = 0;
            if (itr != mp.end()) likelihood_ret[itr->second] = 1;
            node_count++;

            //std::cout << " : " << node_count << std::endl;

            return likelihood_ret;

        }

        else if (is_checkpoint_node(node)) {
            long double leftL[AA_variety]{};
            Likelihood(node->left, leftL, site, 20);
            long double substitution_rate_matrix_left[AA_variety][AA_variety]{};
            make_calculationmatrix(node->left_length, mode);
            matrix_calculation(node->left_length, substitution_rate_matrix_left, 20);
            long double likelihood_ret_left;
            long double likelihood_ret_right;
            for (int i = 0; i < AA_variety; ++i) {
                likelihood_ret_left = 0;
                likelihood_ret_right = 0;
                for (int j = 0; j < AA_variety; ++j) {
                    likelihood_ret_left += substitution_rate_matrix_left[i][j] * leftL[j];
                    likelihood_ret_right = 1;
                }
                likelihood_ret[i] = likelihood_ret_left * likelihood_ret_right;
            }

            return likelihood_ret;
        }
        else {
            long double leftL[AA_variety]{};
            long double rightL[AA_variety]{};
            Likelihood(node->left, leftL, site, mode);
            Likelihood(node->right, rightL, site, mode);
            long double substitution_rate_matrix_left[AA_variety][AA_variety]{};
            long double substitution_rate_matrix_right[AA_variety][AA_variety]{};
            make_calculationmatrix(node->left_length, mode);
            make_calculationmatrix(node->right_length, mode);
            matrix_calculation(node->left_length, substitution_rate_matrix_left, mode);
            matrix_calculation(node->right_length, substitution_rate_matrix_right, mode);

            long double likelihood_ret_left;
            long double likelihood_ret_right;
            for (int i = 0; i < AA_variety; ++i) {
                likelihood_ret_left = 0;
                likelihood_ret_right = 0;
                for (int j = 0; j < AA_variety; ++j) {
                    likelihood_ret_left += substitution_rate_matrix_left[i][j] * leftL[j];
                    likelihood_ret_right += substitution_rate_matrix_right[i][j] * rightL[j];
                }
                if (likelihood_ret[i] < 0) {
                    std::cout << "error" << std::endl;
                }
                likelihood_ret[i] = likelihood_ret_left * likelihood_ret_right;
            }
            return likelihood_ret;
        }
}


int main() {
    const char* fileName = "BLAT_20111_upper1001_seq.csv";
    std::ofstream ofs(fileName);
    const char* fileName2 = "BLAT_upper1001_seq_220111.txt";
    ofstream ofs2(fileName2);
    if (!ofs)
    {
        std::cout << "ファイルが開けませんでした。" << std::endl;
        std::cin.get();
        return 0;
    }
     if (!ofs2)
    {
        std::cout << "ファイルが開けませんでした。" << std::endl;
        std::cin.get();
        return 0;
    }
    std::ifstream file_open("BLAT_data_upper1000_sequence_remove_gap_10percent_include_root_211124.a2m.treefile");
    std::string str;
    std::string Phy_str;
    if (file_open) {
        while (true) {
            getline(file_open, str);
            Phy_str += str;
            if (file_open.eof())
                break;
        }
    }
    if (Phy_str[Phy_str.length() - 1] == ';') {
        Phy_str = Phy_str.substr(0, Phy_str.length()-1);
    }
    sequence_part();
    Node* root = separate(Phy_str);
    int count = 0;
    //Node* cp_node_L = set_cp_node(root, "", t, 1);
    //Node* cp_node_R = set_cp_node(root, "", t, 2);
    make_matrix();
    make_matrix_lessAA();
    long double likelihood_ret[AA_variety];
    int mode = 20;
    long double prior_AA_array[AA_variety];
    prior_AA_ancestral_uniform(prior_AA_array, 1);
    long double posterior_AA_array[AA_variety]{};
    std::string estimated_ancestral_sequence;

    for (int j = 0; j < sequence_length; ++j) {
        Likelihood(root, likelihood_ret, j, mode);
        std::cout << "=====LOCUS[" << j + 1 << "]=====" << std::endl;

        char st[AA_variety] = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
                       'V', '-' };
        char st_exclude_gap[20] = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T',
        'W', 'Y', 'V' };
        std::vector<long double> vec(20);
        //calculate posterior probabilty
        for (int i = 0; i < AA_variety; ++i) {
            posterior_AA_array[i] = prior_AA_array[i] * likelihood_ret[i];
        }
        normalize_prob(posterior_AA_array);

        for (int i = 0; i < vec.size(); ++i) {
            vec[i] = posterior_AA_array[i];
        }
        auto iter = std::max_element(vec.begin(), vec.end());
        int index = std::distance(vec.begin(), iter);
        std::cout << st_exclude_gap[index] << " = " << vec[index] << std::endl;
        ofs << j + 1 << endl;
        ofs2 << j+1 << "   " << j+1 << "   ";
        for (int i = 0; i < AA_variety; ++i) {
            std::cout << st[i] << " = ";
            std::cout << std::fixed;
            std::cout << std::setprecision(6) << posterior_AA_array[i] << "   ";
            ofs << st[i] << "," << posterior_AA_array[i] << ",";
            ofs2 << st[i] << "(" << posterior_AA_array[i] << ")  ";
            if (i % 7 == 6) {
                ofs << endl;
                std::cout << std::endl;
            }


        }
        estimated_ancestral_sequence += st_exclude_gap[index];
        std::cout << j + 1 << " : " << node_count << std::endl;
        node_count = 0;
        ofs2 << endl;
    }
    return 0;
}

